

use super::op::CompressMark;
use std::convert::AsRef;
use std::marker::PhantomData;
use std::ops::Range;
use std::mem::size_of;
use std::slice::{ Iter, from_raw_parts, from_raw_parts_mut };


/* Precursor

keeps offsets in Vec, then convert it to slice in Vec.
(Is there better (safer) way? This implementation can't force user to use the same
vector on creating and transforming it. It's better if we can put lifetime paramater
to match it to that of vec.)
*/
#[derive(Copy, Clone, Debug, Default)]
pub(super) struct SlicePrecursor<T> {
	/* just equivalent to Range<usize> */
	ofs: usize,
	len: usize,
	_marker: PhantomData<T>			/* not effectively used for now. what to do with this? */
}

#[allow(dead_code)]
impl<'a, T> SlicePrecursor<T> {
	pub(super) fn compose(range: &Range<usize>) -> Self {
		SlicePrecursor::<T> {
			ofs: range.start,
			len: range.end - range.start,
			_marker: PhantomData::<T>
		}
	}

	pub(super) fn finalize_raw(&self, base: *const u8) -> &'a [T] {
		let ptr = base.wrapping_add(self.ofs) as *const T;
		let cnt = self.len / size_of::<T>();
		unsafe { from_raw_parts(ptr, cnt) }
	}

	pub(super) fn finalize(&self, v: &'a Vec<u8>) -> &'a [T] {
		let base = v.as_ptr() as *const u8;
		self.finalize_raw(base)
	}
}


/* SimdAlignedU8

16-byte aligned array for SSE and NEON
*/
#[repr(align(16))]
pub(super) struct SimdAlignedU8 {
	v: [u8; 16]
}

impl SimdAlignedU8 {
	pub(super) const fn new(v: &[u8; 16]) -> Self {
		SimdAlignedU8 { v: *v }
	}
}

impl AsRef<[u8; 16]> for SimdAlignedU8 {
	fn as_ref(&self) -> &[u8; 16] {
		&self.v
	}
}


/* Writer and reserve_to

This trait and functions is for vectorized array storing. The `reserve_to` function
allocates `count` elements in the writer `T`, which is typically Vec<u8> and passes
it to a clsure for use as a destination array.

The return value `usize` of the closure is the actual number of elements stored to
the array. It is used for shrinking the base buffer to the exact length that valid
elements stored.
*/
pub(super) trait Writer<T: Sized + Copy + Default> {
	type T;

	/*
	`func: FnOnce(&mut [U], &[T]) -> usize` supposed to take working buffer of
	size `count` and corresponding base array at the head, do the work, and
	return object count that is actually used in the procedure.
	*/
	fn reserve_to<U, F>(&mut self, count: usize, func: F) -> Range<usize>
		where U: Copy,
		      F: FnOnce(&mut [U], &[T]) -> usize;
}

/* implementation for Vec<u8> */
impl<T: Sized + Copy + Default> Writer<T> for Vec<T> {

	type T = u8;

	fn reserve_to<U, F>(&mut self, count: usize, func: F) -> Range<usize>
		where U: Copy,
		      F: FnOnce(&mut [U], &[T]) -> usize
	{
		let base_offset = (self.len() + size_of::<U>() - 1) / size_of::<U>() * size_of::<U>();
		let request_size = size_of::<U>() * count;

		/* extend array to hold requested count of U */
		self.resize(base_offset + request_size, T::default());			/* can be mem::uninitialized() ? */

		/* split resized buffer into base (immutable) and working area (mutable) */
		let (base, work) = self.split_at_mut(base_offset);

		/* convert &mut [u8] -> &mut [U] */
		let ptr = work.as_mut_ptr();
		let work = unsafe { from_raw_parts_mut(ptr as *mut U, count) };

		/* do the work */
		let consumed_count = func(work, base);

		/* save the length */
		assert!(consumed_count <= count);
		let consumed_size = size_of::<U>() * consumed_count;
		self.resize(base_offset + consumed_size, T::default());

		/* returns base_offset (count in bytes) for composing packed data structure */
		Range::<usize> {
			start: base_offset,
			end:   base_offset + consumed_size
		}
	}
}


/* PeekFold iterator

Similar to `Iterator::try_fold`, but it doesn't consume the last-peeked element.
The iteration can be resumed from the element of the last failure.
*/
pub trait PeekFold<T: Sized> {
	fn peek_fold<A, F>(&mut self, init: A, func: F) -> A
		where Self: Sized,
		      A: Copy,
		      F: FnMut(A, &T) -> Option<A>;
}

/* this implementation uses `as_slice()` for peeking the head element, assuming T being slice::Iter */
impl<'a, T: Sized> PeekFold<T> for Iter<'a, T> {
	fn peek_fold<A, F>(&mut self, init: A, mut func: F) -> A
		where Self: Sized,
		      A: Copy,
		      F: FnMut(A, &T) -> Option<A>
	{
		let mut accum: A = init;
		loop {
			/* peek */
			let x = match self.as_slice().first() {
				None => { return accum; },
				Some(x) => x
			};

			/* updated accumulator is discarded when func returned false */
			let next_accum = match func(accum, x) {
				None => { return accum; },
				Some(x) => x
			};

			/* if continuous flag is true, overwrite (update) accumulator and forward iterator */
			self.next();
			accum = next_accum;
		}
	}
}


/* convert ascii-printed number to u64, overflow is ignored. */
pub(super) fn atoi_unchecked(v: &mut Iter<u8>) -> u64 {
	v.peek_fold(0, |a, x| {
		let m = (*x as u64).wrapping_sub('0' as u64);
		if m >= 10 { return None; }

		return Some(10 * a + m);
	})
}

#[allow(dead_code)]
pub(super) fn isnum(c: u8) -> bool {
	return (c as u64).wrapping_sub('0' as u64) < 10 as u64;
}

/*
transcode { 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't' } -> { 0x0, 0x1, 0x2, 0x3, 0x0, 0x1, 0x2, 0x3 }
and then mark mismatch.
*/
pub(super) fn transcode_base_unchecked(c: u8) -> u32 {
	let c = c as u32;
	let b2 = ((c>>1) - (c>>3)) & 0x03;
	return CompressMark::Mismatch as u32 + b2;
}

#[allow(dead_code)]
pub(super) fn encode_base_unchecked(c: char) -> u8 {
	match c {
		'A' | 'a' => 0x01,
		'C' | 'c' => 0x02,
		'G' | 'g' => 0x04,
		'T' | 't' => 0x08,
		_ => 0x00
	}
}

#[allow(dead_code)]
pub(super) fn decode_base_unchecked(c: u32) -> char {
	/* one-hot encoding */
	match "-AC-G---T-------".bytes().nth(c as usize) {
		None    => '-',
		Some(x) => x as char
	}
}


/* all the remainings are unittests */
#[cfg(test)]
mod test_utils {
	use super::{ PeekFold, atoi_unchecked, isnum };

	macro_rules! test_atoi_unchecked_impl {
		( $str: expr, $( $num: expr ),* ) => ({
			let s = $str;
			let mut it = s.as_bytes().iter();
			$({
				let n = atoi_unchecked(&mut it);
				assert_eq!(n, $num);

				(&mut it).peek_fold(0, |_, &x| { if isnum(x) { None } else { Some(0) } });
			})*
		})
	}

	#[test]
	fn test_atoi_unchecked() {
		test_atoi_unchecked_impl!("0", 0);
		test_atoi_unchecked_impl!("10", 10);
		test_atoi_unchecked_impl!("-10", 0);

		/* the following also work tests for PeekFold iterator */
		test_atoi_unchecked_impl!("10M11", 10, 11);
		test_atoi_unchecked_impl!("X0M1X222222MMMMMMM1234XXXX", 0, 0, 1, 222222, 1234);
	}

	#[test]
	fn test_isnum() {
		/* trivial */
		assert_eq!(isnum('0' as u8), true);
		assert_eq!(isnum('9' as u8), true);
		assert_eq!(isnum('-' as u8), false);
		assert_eq!(isnum(' ' as u8), false);
	}
}

