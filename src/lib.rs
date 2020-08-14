
/**
@file udon/lib.rs
@brief udon implementation

@author Hajime Suzuki
@license MIT

@detail
Udon is a transcoding and indexing data structure / algorithm for BAM CIGAR / MD strings.
It enables fast querying of subalignments with arbitrary span, with optional scaling for visualization.

See SAM / BAM specification for the definition of the CIGAR / MD strings:
https://samtools.github.io/hts-specs/
*/

#[macro_use]
extern crate bitfield;
extern crate static_assertions;

use std::io::Write;
use std::marker::PhantomData;
use std::mem::{ size_of, transmute };
use std::ops::Range;
use std::ptr::copy_nonoverlapping;
use std::slice::{ Iter, IterMut, from_raw_parts, from_raw_parts_mut };
use std::str::from_utf8;


/* architecture-dependent stuffs */
#[cfg(all(target_arch = "x86_64"))]
use core::arch::x86_64::*;

#[cfg(all(target_arch = "aarch64"))]
use core::arch::aarch64::*;				/* requires nightly */


/* logging */
#[macro_use]
extern crate log;

#[cfg(test)]
use std::println as debug;




/* bam CIGAR (raw CIGAR) related structs:
Transcoder API takes an array of ops, defined as the the following:

```
struct {
	uint32_t op : 4;
	uint32_t len : 28;
};
```

the detailed specification is found in SAM/BAM spec pdf.
*/

/* bam CIGAR operation code */
#[repr(u32)]
#[allow(dead_code)]
enum CigarOp {
	Match    = 0,
	Ins      = 0x01,
	Del      = 0x02,
	Unknown  = 0x03,
	SoftClip = 0x04,
	HardClip = 0x05,
	Pad      = 0x06,
	Eq       = 0x07,
	Mismatch = 0x08
}

/* bam CIGAR element definition */
bitfield!{
	#[derive(Copy, Clone, Default)]
	pub struct Cigar(u32);
	pub op, _:  3,  0;
	pub len, _: 31, 4;
}




/* transcoded CIGAR and index

Transcoded CIGAR is encoded in a run-length manner, as in the original CIGAR string.
Each chunk for the run-length compression is represented by an op, which is 8-bit unsigned.

Op consists of two bitfields: total chunk length (columns) on reference sequence (5bit),
and leading event(s) (3bit). The leading event is one of { insertion, deletion, mismatch }
and the remainder till the end of the chunk is filled with matches. That is, any a chunk
has one of the following composition:

1. insertion of arbitrary length, and trailing matches up to 30 bases.
2. deletion up to three bases, and trailing matches up to 30 - #deletions.
3. a mismatch, and trailing matches up to 29 bases.


The 3-bit leading is defined to hold the events above as the following:

* for insertion, place a marker (0b000) that there is inserted bases between the chunk
    and the previous chunk.
* for deletion, encode the number of deleted bases as 0b001 ~ 0b011 (1 to 3).
* for mismatch, encode a base on the query sequence as 0b100 (A) ~ 0b111 (T).


The transcoded CIGAR stream is indexed at reference-side positions at regular intervals.
The interval is defined as `BLOCK_PITCH`, which is 256 by default. Pointer to the op array
and the insertion sequence array for each block is stored in `Block`. Since block boundary
is not always aligned to chunk (op) boundary, an offset within chunk from the block
boundary is also stored in `op_skip` field.
*/


/* block pitch, see above */
const BLOCK_PITCH: usize = 256;

/* block object, see above */
bitfield!{
	#[derive(Copy, Clone, Debug, Default)]
	pub struct Block(u64);
	pub ins_offset, set_ins_offset: 28, 0;
	pub op_offset,  set_op_offset:  58, 29;
	pub op_skip,    set_op_skip:    63, 59;
}

/* op trailing events */
#[repr(u32)]
enum CompressMark {
	Ins      = 0x00,

	/* mismatch base for 'A': 0x04, .., 'T': 0x07 */
	Mismatch = 0x04,
}

/* event for output (expanded) array */
#[repr(u32)]
#[allow(dead_code)]
enum Op {
	MismatchA = CompressMark::Mismatch as u32 | 0x00,
	MismatchC = CompressMark::Mismatch as u32 | 0x01,
	MismatchG = CompressMark::Mismatch as u32 | 0x02,
	MismatchT = CompressMark::Mismatch as u32 | 0x03,
	Del       = 0x08,
	Ins       = 0x10
}


/* Transcoded CIGAR and index object container

everything is laid out in a single flat memory block, thus constructed via
some unsafe operations. see UdonPrecursor and Precursor traits.
*/
#[derive(Copy, Clone, Debug, Default)]
pub struct Udon<'a> {
	size: usize,				/* object size for copying; for Clone */
	ref_span: usize,			/* reference side span */
	op: &'a [u8],
	block: &'a [Block],
	ins: &'a [u8]
}


/* UdonPrecursor

Precursor is type that has some invalid pointers (slices) but the others are all sane.
The invalid pointer has an offset from a certain memory block, which may be reallocated
along with the object construction. The offset is converted to valid pointer, by adding
base pointer to the offset, at the very last of the object construction.
*/
#[repr(C)]
#[derive(Copy, Clone, Debug, Default)]
pub struct UdonPrecursor {
	size: usize,
	ref_span: usize,
	op: SlicePrecursor<u8>,
	block: SlicePrecursor<Block>,
	ins: SlicePrecursor<u8>
}
assert_eq_size!(Udon, UdonPrecursor);


/* Precursor

keeps offsets in Vec, then convert it to slice in Vec.
(Is there better (safer) way? This implementation can't force user to use the same
vector on creating and transforming it. It's better if we can put lifetime paramater
to match it to that of vec.)
*/
#[derive(Copy, Clone, Debug, Default)]
struct SlicePrecursor<T> {
	/* just equivalent to Range<usize> */
	ofs: usize,
	len: usize,
	_marker: PhantomData<T>
}

#[allow(dead_code)]
impl<'a, T> SlicePrecursor<T> {
	fn compose(range: Range<usize>) -> Self {
		SlicePrecursor::<T> {
			ofs: range.start,
			len: range.end - range.start,
			_marker: PhantomData::<T>
		}
	}

	fn finalize_raw(self, base: *const u8) -> &'a [T] {
		let ptr = base.wrapping_add(self.ofs) as *const T;
		let cnt = self.len / size_of::<T>();
		unsafe { from_raw_parts(ptr, cnt) }
	}

	fn finalize(self, v: &'a Vec<u8>) -> &'a [T] {
		let base = v.as_ptr() as *const u8;
		Self::finalize_raw(self, base)
	}
}

/*
trait SlicePrecursor<'a, T> {
	fn compose(range: Range<usize>) -> Self;
	fn finalize(self, base: *const u8) -> Self;
}

impl<'a, T> Precursor<'a, T> for &'a [T] {
	fn compose(range: Range<usize>) -> Self {
		assert!(((range.end - range.start) % size_of::<T>()) == 0);

		let ptr   = range.start as *const T;
		let count = (range.end - range.start) / size_of::<T>();

		unsafe { from_raw_parts(ptr, count) }
	}

	fn finalize(self, base: *const u8) -> &'a [T] {
		let offset = self.as_ptr() as usize;

		let ptr   = base.wrapping_add(offset) as *const T;
		let count = self.len();

		unsafe { from_raw_parts(ptr, count) }
	}
}
*/


/* SimdAlignedU8

16-byte aligned array for SSE and NEON
*/
#[repr(align(16))]
struct SimdAlignedU8 {
	v: [u8; 16]
}


/* Writer and reserve_to

This trait and functions is for vectorized array storing. The `reserve_to` function
allocates `count` elements in the writer `T`, which is typically Vec<u8> and passes
it to a clsure for use as a destination array.

The return value `usize` of the closure is the actual number of elements stored to
the array. It is used for shrinking the base buffer to the exact length that valid
elements stored.
*/
trait Writer<T: Sized + Copy + Default> {
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
		debug!("base_offset({}), request_size({})", base_offset, request_size);

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

Similar to Iterator::`try_fold`, but it doesn't consume the last-peeked element.
The iteration can be resumed from the element of the last failure.
*/
trait PeekFold<T: Sized> {
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
fn atoi_unchecked(v: &mut Iter<u8>) -> u64 {
	v.peek_fold(0, |a, x| {
		let m = (*x as u64).wrapping_sub('0' as u64);
		if m >= 10 { return None; }

		return Some(10 * a + m);
	})
}

#[allow(dead_code)]
fn isnum(c: u8) -> bool {
	return (c as u64).wrapping_sub('0' as u64) < 10 as u64;
}

/*
transcode { 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't' } -> { 0x0, 0x1, 0x2, 0x3, 0x0, 0x1, 0x2, 0x3 }
and then mark mismatch.
*/
fn transcode_base_unchecked(c: u8) -> u32 {
	let c = c as u32;
	let b2 = ((c>>1) - (c>>3)) & 0x03;
	return CompressMark::Mismatch as u32 + b2;
}

#[allow(dead_code)]
fn encode_base_unchecked(c: char) -> u8 {
	match c {
		'A' | 'a' => 0x01,
		'C' | 'c' => 0x02,
		'G' | 'g' => 0x04,
		'T' | 't' => 0x08,
		_ => 0x00
	}
}

#[allow(dead_code)]
fn decode_base_unchecked(c: u32) -> char {
	match "ACGT".bytes().nth(c as usize) {
		None    => '-',
		Some(x) => x as char
	}
}

/*
op -> len conversion
*/
fn op_len(x: u8) -> usize {
	let len = x as usize & 0x1f;
	len - (len == 31) as usize
}

fn op_marker(x: u8) -> u32 {
	(x>>5) as u32
}


#[cfg(test)]
mod test_utils {
	use crate::{ atoi_unchecked, isnum, PeekFold };

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
		assert_eq!(isnum('0' as u8), true);
		assert_eq!(isnum('9' as u8), true);
		assert_eq!(isnum('-' as u8), false);
		assert_eq!(isnum(' ' as u8), false);
	}
}



/* strip_clips

bam CIGAR string sometimes has clipping ops at the head and tail, and they need removed
before conversion to augumented CIGAR. This function does this, and returns cigar slice
without clips along with clipped lengths.

Clipped length is defined as the number of bases to be removed from the matched query
sequence. Thus (head, tail) = (0, 0) is returned for soft-clipped alignment.
*/
#[derive(Copy, Clone, Default)]
struct QueryClip {
	head: usize,
	tail: usize
}

fn strip_clips(cigar: &[Cigar]) -> Option<(&[Cigar], QueryClip)> {

	/* cigar length must be at least one */
	if cigar.len() == 0 { return None; }

	let mut cigar = cigar;
	let mut clip = QueryClip {
		head: 0,
		tail: 0
	};

	/* strip head */ {
		let head = cigar[0];
		if head.op() == CigarOp::SoftClip as u32 {
			/* query string is not clipped so we have to do it */
			clip.head = head.len() as usize;

			let (_, body) = cigar.split_first()?;
			cigar = body;
		} else if head.op() == CigarOp::HardClip as u32 {
			let (_, body) = cigar.split_first()?;
			cigar = body;
		}
	}

	/* strip tail */ {
		let tail = cigar[cigar.len() - 1];
		if tail.op() == CigarOp::SoftClip as u32 {
			clip.tail = tail.len() as usize;

			let (_, body) = cigar.split_last()?;
			cigar = body;
		} else if tail.op() == CigarOp::HardClip as u32 {
			let (_, body) = cigar.split_last()?;
			cigar = body;
		}
	}
	return Some((cigar, clip));
}


/* copy 4bit-encoded bases, for storing insertions */
fn copy_packed_nucl(src: &[u8], dst: &mut [u8], ofs: usize, len: usize) {

	/* very very very very very very very very very naive way */
	for i in 0 .. len {
		let pos = ofs + i;

		/* bam packed nucl is big endian */
		let c = if (pos & 0x01) == 0 {
			src[pos / 2]>>4
		} else {
			src[pos / 2] & 0x0f
		};

		/* store to dst */
		if (i & 0x01) == 0 {
			dst[i / 2]  = c;
		} else {
			dst[i / 2] |= c<<4;
		}
	}
}


/* transcoder state object */
pub struct UdonBuilder<'a> {
	buf: Vec<u8>,
	ins: Vec<u8>,
	cigar: Iter<'a, Cigar>,
	mdstr: Iter<'a, u8>,
	query: &'a [u8],
	qofs: usize,
	rspan: usize
}

impl<'a, 'b> UdonBuilder<'a> {

	/*
	 * input parsers
	 */

	/* imitate iterator on 4-bit packed nucleotide sequence */
	fn next_base(&mut self) -> u32 {
		let ofs = self.qofs;
		self.qofs += 1;

		debug!("qlen({}), ofs({})", self.query.len(), ofs);
		assert!(ofs / 2 < self.query.len(), "qlen({}), ofs({})", self.query.len(), ofs);

		let c = self.query[ofs / 2];
		if (ofs & 0x01) == 0x01 {
			return transcode_base_unchecked(c & 0x0f);
		}
		return transcode_base_unchecked(c>>4);
	}

	/*
	We need cigar iterator peekable, but Iter::Peekable implementation is redundant
	for slice::Iter. What we really need is peeking the head by .as_slice()[0].
	*/
	fn peek_cigar(&self) -> Option<Cigar> {
		if self.cigar.as_slice().len() == 0 {
			return None;
		}

		return Some(self.cigar.as_slice()[0]);
	}

	/* MD string handling: forward pointer along with atoi */
	fn eat_md_eq(&mut self) -> usize {
		return atoi_unchecked(&mut self.mdstr) as usize;
	}

	/* make MD iterator peekable */
	fn is_double_mismatch(&self) -> bool {
		if self.mdstr.as_slice().len() < 2 {
			return false;
		}
		return self.mdstr.as_slice()[1] == '0' as u8;
	}

	/* writers, on top of impl Writer for Vec<u8> */
	fn push_op(&mut self, match_len: usize, marker: u32) {
		assert!(marker    < 0x08);
		assert!(match_len < 0x20);

		debug!("len({}), marker({})", match_len, marker);

		let op = (marker<<5) as u8 | match_len as u8;
		assert!(op != 0);

		self.buf.reserve_to(1, |arr: &mut [u8], _: &[u8]| -> usize {
			arr[0] = op;
			1
		});
	}
	fn push_match(&mut self, match_len: usize, last_op: u32) {
		let mut op = last_op;
		let mut rem = match_len;

		/* divide into 30-column chunks (maximum chunk length for an op is 30) */
		while rem > 30 {
			self.push_op(31, op);		/* 31 is alias for 30-column chunk with "continuous" flag */
			op   = 0;
			rem -= 30;
		}
		self.push_op(rem, op);			/* push the last */
	}

	/* forward both cigar and md strings */
	fn eat_del(&mut self) -> Option<u32> {
		let c = self.cigar.next()?;
		assert!(c.op() == CigarOp::Del as u32);

		let len = c.len() as usize;
		self.mdstr.nth(len - 1)?;		/* error if starved */
		self.mdstr.next();				/* not regarded as error for the last element */

		/* 3 columns at maximum for deletion per chunk */
		let mut rem = len;
		while rem > 3 {
			self.push_op(3, 3);			/* deletion length == 3, chunk length == 3 */
			rem -= 3;
		}

		assert!(rem > 0);
		return Some(rem as u32);		/* remainder is concatenated to the next match into the next chunk */
	}
	fn eat_ins(&mut self, xrem: usize) -> Option<usize> {
		assert!(self.qofs <= std::i32::MAX as usize);

		let c = self.cigar.next()?;
		assert!(c.op() == CigarOp::Ins as u32);

		/* forward query offset */
		let ofs = self.qofs;
		let len = c.len() as usize;
		self.qofs += len;

		self.save_ins(ofs, len);		/* use offset before forwarding */
		return Some(xrem - len);		/* there might be remainder */
	}
	fn eat_match(&mut self, xrem: usize, last_op: u32) -> Option<usize> {
		let c = self.peek_cigar()?;
		debug!("c({}, {})", c.op(), c.len());

		let is_valid = c.op() != CigarOp::Ins as u32;
		if is_valid { self.cigar.next()?; }

		let mut crem = if is_valid { c.len() as usize } else { 0 } + last_op as usize;
		let mut xrem = xrem;			/* might continues from the previous cigar op, possibly insertion */
		let mut op = last_op;			/* insertion, deletion, or mismatch */

		debug!("eat_match, crem({}), xrem({})", crem, xrem);

		while xrem < crem {
			debug!("mismatch?, crem({}), xrem({})", crem, xrem);

			/* xrem < crem indicates this op (cigar span) is interrupted by mismatch(es) at the middle */
			self.push_match(xrem, op);
			crem      -= xrem;
			self.qofs += xrem;

			while self.is_double_mismatch() {
				debug!("{:?}", from_utf8(self.mdstr.as_slice()));

				let c = self.next_base();
				self.push_op(1, c);		/* this chunk contains only a single mismatch */
				self.mdstr.nth(1)?;
				crem -= 1;
			}
			debug!("{:?}", from_utf8(self.mdstr.as_slice()));

			op = self.next_base();		/* we only have a single mismatch remaining, will be combined to succeeding matches */
			self.mdstr.nth(0)?;

			/* next match length */
			xrem = self.eat_md_eq() + 1;/* xrem for the next { match, insertion } region, including the last mismatch */
			debug!("updated xrem, crem({}), xrem({})", crem, xrem);
		}

		self.push_match(crem, op);		/* tail match; length is the remainder of crem */
		xrem      -= crem;
		self.qofs += crem;

		debug!("done, crem({}), xrem({}), qofs({})", crem, xrem, self.qofs);
		return Some(xrem);				/* nonzero if insertion follows */
	}

	/* insertion bin handling */
	fn save_ins(&mut self, qofs: usize, qlen: usize) {

		if qlen == 0 { return; }		/* must be corrupted cigar, but safe to ignore */

		/* copy subsequence; 4-bit packed */
		let packet_size = (qlen + 1) / 2 + 1;
		let query = self.query;
		self.ins.reserve_to(packet_size, |arr: &mut [u8], _: &[u8]| -> usize {
			copy_packed_nucl(&query, arr, qofs, qlen);
			arr[arr.len() - 1] = 0;

			packet_size
		});
	}

	/* parse raw CIGAR into augumented CIGAR, takes clip-stripped raw CIGAR, which must be stored in Self. */
	fn parse_cigar(&mut self) -> Option<Range<usize>> {

		let base_offset = self.buf.len();
		let mut xrem = 0;

		/*
		CIGAR might start with insertion (possible in short-read alignment).cigar
		in such case, we need a dummy ins marker at the head of op array,
		because dummy ins is placed for normal CIGAR (that begins with match) and is
		skipped on decoding. we have to tell the decoder the head insertion must not
		be ignored by adding one more insertion.
		*/
		let c = self.peek_cigar()?;
		debug!("c({}, {})", c.op(), c.len());
		if c.op() == CigarOp::Ins as u32 {
			/* op iterator not forwarded here */
			self.push_op(0, CompressMark::Ins as u32);
			self.ins.write(&[0]).ok()?;	/* dummy insertion marker for this */

		} else if c.op() == CigarOp::Match as u32 {
			xrem = self.eat_md_eq();

			/*
			the first op consists of dummy insertion and match
			(the insertion is real one for a CIGAR that starts with insertion. see above.)
			*/
			xrem = self.eat_match(xrem, CompressMark::Ins as u32)?;
			self.ins.write(&[0]).ok()?;	/* dummy insertion marker at the head */
		}

		'outer: loop {
			macro_rules! peek_or_break {
				( $self: expr ) => ({
					match $self.peek_cigar() {
						None    => { break 'outer; },
						Some(x) => x
					}
				});
			}

			/* deletion-match pair */
			loop {
				let c = peek_or_break!(self);
				debug!("c({}, {})", c.op(), c.len());
				if c.op() != CigarOp::Del as u32 { break; }

				/* the CIGAR ends with deletion; must be treated specially */
				if self.cigar.as_slice().len() < 2 { break; }

				/* push deletion-match pair, then parse next eq length */
				let op = self.eat_del()?;
				xrem = self.eat_md_eq() + op as usize;
				debug!("op({}), xrem({})", op, xrem);
				xrem = self.eat_match(xrem, op)?;
			}

			/* it's insertion-match pair when it appeared not be deletion-match */
			let c = peek_or_break!(self);
			debug!("c({}, {})", c.op(), c.len());
			if c.op() != CigarOp::Ins as u32 {
				return None;			/* if not, we regard it broken */
			}

			debug!("op({}), len({}), remaining cigars({})", c.op(), c.len(), self.cigar.as_slice().len());

			/* the CIGAR ends with insertion; must be treated specially */
			if self.cigar.as_slice().len() < 2 { break; }

			/* push insertion-match pair, update eq length remainder */
			xrem = self.eat_ins(xrem)?;
			xrem = self.eat_match(xrem, CompressMark::Ins as u32)?;
		}

		/* CIGAR ends with isolated insertion or deletion */
		if self.cigar.as_slice().len() > 0 {
			let c = self.peek_cigar()?;
			debug!("c({}, {})", c.op(), c.len());

			if c.op() == CigarOp::Del as u32 {
				let op = self.eat_del()?;
				self.push_op(op as usize, op);
			} else if c.op() == CigarOp::Ins as u32 {
				self.eat_ins(xrem)?;
			}
		}

		/* range in output buffer, for composing precursor */
		return Some(Range::<usize> {
			start: base_offset,
			end:   self.buf.len()
		});
	}

	/* accumulate chunk lengths to compute reference span (for computing index table size) */
	fn calc_reference_span(&self, op_range: &Range<usize>) -> usize {
		let ops = &self.buf[op_range.start..op_range.end];

		ops.iter().fold(0, |a, &x| { a + op_len(x) })
	}

	/* construct index for blocks */
	fn forward_ins(ins: &mut Iter<u8>, op_offset: usize, marker: u32) -> usize {
		if marker != 0 { return 0; }
		if op_offset == 0 { return 0; }		/* ignore first ins marker */

		ins.peek_fold(0, |a, &x| {
			if x == 0 { return None; }
			return Some(a + 1);
		})
	}

	fn push_block(dst: &mut IterMut<Block>, ins_offset: usize, op_offset: usize, op_skip: usize) {
		let bin = match dst.next() {
			None => { return; },
			Some(bin) => bin
		};
		bin.set_ins_offset(ins_offset as u64);
		bin.set_op_offset(op_offset as u64);
		bin.set_op_skip(op_skip as u64);
	}

	fn build_index(&mut self, ref_span: usize, op_range: &Range<usize>, ins_range: &Range<usize>) -> Option<Range<usize>> {

		/* block pitch must be 2^n */
		assert!(BLOCK_PITCH.next_power_of_two() == BLOCK_PITCH);

		/* rip buffer for disjoint ownership */
		let buf = &mut self.buf;
		debug!("{:?}", buf);

		/* index for regular pitch on span */
		let block_count = (ref_span + BLOCK_PITCH - 1) / BLOCK_PITCH;
		let range = buf.reserve_to(block_count, |block: &mut [Block], base: &[u8]| {
			debug!("{:?}", base);

			/* src: both are &[u8] */
			let ops = (&base[op_range.start .. op_range.end]).iter();
			let mut ins = (&base[ins_range.start .. ins_range.end]).iter();

			/* dst */
			let mut dst = block.iter_mut();
			Self::push_block(&mut dst, 0, 0, 0);

			/* I want this loop be more lightweight... */
			let mut rpos: usize = 0;		/* ref_pos */
			let mut iofs: usize = 0;		/* ins_offset */
			for (i, &x) in ops.enumerate() {
				let next_boundary = (rpos | (BLOCK_PITCH - 1)) + 1;

				/* forward insertion array and reference-side offset */
				rpos += op_len(x);
				iofs += Self::forward_ins(&mut ins, i, op_marker(x));

				/* if forwarded rpos doexn't exceed the next boundary, just skip this op */
				if rpos < next_boundary { continue; }

				/* boundary found; save block info */
				Self::push_block(&mut dst, iofs, i, rpos & 0x1f);
			}

			let dst = dst.into_slice();
			assert!(dst.len() == 0, "len({})", dst.len());
			block.len()
		});
		return Some(range);
	}

	fn pack_ins(&mut self) -> Option<Range<usize>> {

		/* rip vectors from self for ownership */
		let buf = &mut self.buf;
		let ins = &self.ins;

		/* just copy */
		let range = buf.reserve_to(ins.len(), |arr: &mut [u8], _: &[u8]| {
			arr.copy_from_slice(ins.as_slice());
			ins.len()
		});
		return Some(range);
	}

	/* convert internal object (UdonBuilder) to UdonPrecursor */
	fn finalize(&self, ref_span: usize, op: &Range<usize>, ins: &Range<usize>, block: &Range<usize>) -> UdonPrecursor {
		assert!(op.end  <= ins.start);
		assert!(ins.end <= block.start);

		assert!(((block.end - block.start) % size_of::<Block>()) == 0);

		/* entire range on buffer for this object */
		let range = Range::<usize> {
			start: op.start,
			end:   ins.end
		};

		/* convert range to slice */
		UdonPrecursor {
			/* just save */
			size:     range.end - range.start,
			ref_span: self.rspan,

			/* compose fake slices (dereference causes SEGV) */
			op:    SlicePrecursor::compose(op),
			block: SlicePrecursor::compose(block),
			ins:   SlicePrecursor::compose(ins)
		}
	}
}

/* build UdonPrecursor using UdonBuilder internally */
impl<'a> UdonPrecursor {
	fn build_core(buf: Vec<u8>, cigar: &'a [Cigar], packed_query: &'a [u8], mdstr: &'a [u8]) -> (Vec<u8>, Option<UdonPrecursor>) {
		/* save initial offset for unwinding */
		let base_offset = buf.len();
		debug!("{:?}", buf);

		/* compose working variables */
		let (cigar, qclip) = match strip_clips(cigar) {
			None => { return (buf, None); },		/* just return buffer (nothing happened) */
			Some((cigar, qclip)) => (cigar, qclip)
		};
		let mut state = UdonBuilder::<'a> {
			buf:   buf,				/* move */
			ins:   Vec::new(),
			cigar: cigar.iter(),	/* iterator for slice is a pair of pointers (ptr, tail) */
			mdstr: mdstr.iter(),
			query: packed_query,
			qofs:  qclip.head,		/* initial offset for query sequence, non-zero for soft clipped alignments */
			rspan: 0				/* calcd in build_index */
		};

		/* if error detected, unwind destination vector and return it (so that the vector won't lost) */
		macro_rules! unwrap_or_unwind { ( $expr:expr ) => {
			match $expr {
				/* equivalent to tail error handling in C? (I'm not sure how rustc treat this...) */
				None => {
					debug!("unwinding");
					let mut buf = state.buf;
					buf.resize(base_offset, 0);		/* unwind buffer for cleaning up errored sequence */
					return (buf, None);
				},
				Some(val) => val
			}
		}}
		debug!("a");

		/* parse input op array */
		let op = unwrap_or_unwind!(state.parse_cigar());
		let ref_span = state.calc_reference_span(&op);
		debug!("op({:?}), ref_span({:?})", op, ref_span);
		debug!("{:?}", state.buf);

		/* pack ins vector */
		let ins = unwrap_or_unwind!(state.pack_ins());
		debug!("ins({:?})", ins);
		debug!("{:?}", state.buf);

		/* build index for op and ins arrays */
		let block = unwrap_or_unwind!(state.build_index(ref_span, &op, &ins));
		debug!("block({:?})", block);

		/* everything done; compose precursor */
		let precursor = state.finalize(ref_span, &op, &ins, &block);
		debug!("done, precursor({:?})", precursor);

		/* return back ownership of the buffer */
		let buf = state.buf;
		(buf, Some(precursor))
	}

	pub unsafe fn build(buf: Vec<u8>, cigar: &'a [u32], packed_query: &'a [u8], mdstr: &'a [u8]) -> (Vec<u8>, Option<UdonPrecursor>) {
		
		/* for compatibility with bam streamed reader */
		let cigar = unsafe { transmute::<&'a [u32], &'a [Cigar]>(cigar) };

		let (buf, precursor) = Self::build_core(buf, cigar, packed_query, mdstr);
		match precursor {
			None    => { debug!("None"); },
			Some(_) => { debug!("Some"); }
		};
		(buf, precursor)
	}
}

/* Udon builder and ribbon slicing APIs

We have two choices for building Udon object: `Udon::build` and `UdonVec::append` in a
safe way.

The first funcion `build` has the simplest interface. It's good to start from this one
since the signature, taking CIGAR, MD, and query sequence slices and returning Box<Udon>>
looks quite natural as a Rust-native library.

The others are for memory-optimized applications. Since boxing udon object and creating
lots of them makes heap fragmented, even if the allocator (GlobalAlloc / jemalloc) does
its best. Thinking that the number of objects sometimes reaches as many as several million
for visualization apps, it's important to prevent such fragmentation for smaller memory
consumption and better performance.

`UdonVec` is for this purpose, creating multiple udon objects in a single monolithic
region of memory. It has two states, mutable for building objects in it and immutable
for using the objects. The mutable vec can be converted to immutable one by `freeze` method.


We also provide unsafe APIs, `UdonPrecursor::build` and `Udon::from_precursor`. This two
methods were originally implemented for use in `UdonVec::append` and `UdonVec::freeze` each.
The first function creates `UdonPrecursor` object, which is almost equivalent to `Udon`
except that pointers in slices are not valid. The invalid pointer in slice actually
retains an offset from the memory region. Being offset, not an absolute pointer, it allows
reallocation (on-the-fly extension) of the memory region without breaking consistency,
which is important for creating multiple objects in a single region one by one. The
offsets are converted to valid slices at once in the `freeze` procedure.

The raw unsafe APIs `UdonPrecursor::build` and `Udon::from_precursor` are useful if we
want to pack udon objects in a more complex way, not a flat vector. It is especially
useful when we design a data structure with additional metadata for visualization, by
keeping `Udon` object as one member of larger struct.

Please make sure the `UdonPrecursor::build` and `Udon::from_precursor` are unsafe APIs.
The unsafeness comes from the requirement that a vector precursor was built and a vector
the precursor is freezed are the same one, which can't be forced by the type system.
If the two vector don't match, the created `Udon` objects are definitely broken and
causes exception.
*/
impl<'a, 'b> Udon<'a> {

	/*
	APIs
	*/
	pub unsafe fn from_precursor(buf: &Vec<u8>, precursor: UdonPrecursor) -> Udon<'a> {
		let base: *const u8 = buf.as_ptr();
		Self::from_precursor_raw(base, &precursor)
	}

	pub fn build(cigar: &'b [u32], packed_query: &'b [u8], mdstr: &'b [u8]) -> Option<Box<Udon<'a>>> {
		/*
		Note:
		This function creates box for `Udon` but it has a flaw that the created box has
		a trailer for the object which is not visible from compiler. The current runtime
		does not cause problem around this because the current allocator (`GlobalAlloc`)
		ignores memory layout information. If the future allocator change this behavior
		to explicitly check the sanity of the layout, it collapses by something like
		"unmatched free size" error.

		(so what to do?)
		*/

		/* size_of::<u32>() == size_of::<Cigar>(); for compatibility with bam streamed reader */
		let cigar = unsafe { transmute::<&'b [u32], &'b [Cigar]>(cigar) };

		Self::build_core(cigar, packed_query, mdstr)
	}

	/*
	API internals
	*/
	unsafe fn from_precursor_raw(base: *const u8, precursor: &UdonPrecursor) -> Udon<'a> {
		Udon::<'a> {
			/* just copy */
			size:     precursor.size,
			ref_span: precursor.ref_span,

			/* adjusting pointer */
			op:    precursor.op.finalize_raw(base),
			block: precursor.block.finalize_raw(base),
			ins:   precursor.ins.finalize_raw(base)
		}
	}

	fn build_core(cigar: &'b [Cigar], packed_query: &'b [u8], mdstr: &'b [u8]) -> Option<Box<Udon<'a>>> {
		let mut buf = Vec::<u8>::new();

		/* header always at the head */
		let range = buf.reserve_to(1, |header: &mut [Udon], _: &[u8]| {
			header[0] = Default::default();
			1
		});
		assert!(range.start == 0);
		assert!(range.end == size_of::<Udon>());

		return match UdonPrecursor::build_core(buf, cigar, packed_query, mdstr) {
			(_, None) => None,
			(buf, Some(precursor)) => Some(Self::compose_box(buf, precursor))
		}
	}

	fn compose_box(mut buf: Vec<u8>, precursor: UdonPrecursor) -> Box<Udon<'a>> {

		/* compose pointer-adjusted header on stack */
		let base: *const u8 = buf.as_ptr();
		let header = unsafe { Self::from_precursor_raw(base, &precursor) };

		/* compose box and copy header into it */
		let udon = unsafe {
			/* convert buffer (on heap) to box */
			let mut udon = Box::<Udon<'a>>::from_raw(
				transmute::<*mut u8, *mut Udon<'a>>(buf.as_mut_ptr())
			);

			/* copy header from stack to heap (box) */
			let src = &header as *const Udon<'a>;
			let dst = &mut udon as &mut Udon<'a> as *mut Udon<'a>;
			copy_nonoverlapping(src, dst, 1);

			udon
		};
		udon
	}


	/*
	decoder implementations
	*/
	pub fn ref_span(&self) -> usize {
		self.ref_span
	}

	pub fn decode_into(&self, dst: &mut Vec<u8>, ref_span: Range<usize>) -> Option<usize> {

		/* check sanity of the span (range) and return None if broken */
		self.check_span(&ref_span)?;

		/* linear polling from block head, then decode until the end of the span */
		let (ops, rem) = self.scan_op_array(ref_span.start);
		let used = Self::decode_core(dst, ops, rem, ref_span.end - ref_span.start)?;
		Some(used)
	}

	pub fn decode(&self, ref_span: Range<usize>) -> Option<Vec<u8>> {

		let size = ref_span.end - ref_span.start;
		let mut buf = Vec::<u8>::with_capacity(size);

		let used = self.decode_into(&mut buf, ref_span)?;
		buf.resize(used, 0);

		return Some(buf);
	}

	#[allow(dead_code, unused_variables)]
	pub fn decode_scaled_into(&self, dst: &mut Vec<u8>, ref_span: Range<usize>, scale: f64) -> Option<usize> {
		unimplemented!();
	}

	#[allow(dead_code, unused_variables)]
	pub fn decode_scaled(&self, ref_span: Range<usize>, scale: f64) -> Option<Vec<u8>> {
		unimplemented!();
	}

	const DEL_MASK: SimdAlignedU8 = {
		let mut x = [0u8; 16];
		x[0] = Op::Del as u8;
		x[1] = Op::Del as u8;
		x[2] = Op::Del as u8;
		SimdAlignedU8 { v: x }
	};
	const SCATTER_MASK: SimdAlignedU8 = {
		let mut x = [0x80u8; 16];
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
		SimdAlignedU8 { v: x }
	};
	const IS_DEL_THRESH: SimdAlignedU8 = {
		let mut x = [0xffu8; 16];
		x[0] = 0x1f;
		x[1] = 0x3f;
		x[2] = 0x5f;
		SimdAlignedU8 { v: x }
	};

	#[cfg(all(target_arch = "x86_64"))]
	unsafe fn decode_core_block(dst: &mut [u8], op: u32, ins: u64) -> (usize, u64) {

		/* load constants; expelled out of the innermost loop when inlined */
		let del_mask      = _mm_load_si128(&Self::DEL_MASK.v as *const [u8; 16] as *const __m128i);
		let scatter_mask  = _mm_load_si128(&Self::SCATTER_MASK.v as *const [u8; 16] as *const __m128i);
		let is_del_thresh = _mm_load_si128(&Self::IS_DEL_THRESH.v as *const [u8; 16] as *const __m128i);

		/* keep op on general-purpose register */
		let rop = op as u64;

		/* compute deletion mask */
		let xop = _mm_cvtsi64_si128(rop as i64);						/* copy op to xmm */
		let xop_scattered = _mm_shuffle_epi8(xop, scatter_mask);		/* [op, op, op, 0, ...] */
		let is_del = _mm_cmpgt_epi8(xop_scattered, is_del_thresh);		/* signed comparison */

		/* compute mismatch / insertion mask */
		let marker = if rop == 0 { ins } else { rop>>5 };
		let ins_mismatch = _mm_cvtsi64_si128(marker as i64);

		/* merge deletion / insertion-mismatch vector */
		let merged = _mm_blendv_epi8(ins_mismatch, del_mask, is_del);
		_mm256_storeu_si256(dst as *mut [u8] as *mut __m256i, _mm256_castsi128_si256(merged));

		/*
		compute forward length; 31 is "continuous marker"
		rop-to-rop critical path length is 6
		*/
		let len = rop as usize & 0x1f;
		let next_ins     = if len == 0x1f { 0 } else { Op::Ins as u64 };	/* next ins will be masked if 0x1f */
		let adjusted_len = if len == 0x1f { len - 1 } else { len };
		(adjusted_len, next_ins)
	}

	#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
	unsafe fn decode_core_block(dst: &mut [u8], op: u32, ins: u64) -> (usize, u64) {
		(0, 0)
	}

	fn decode_core(dst: &mut Vec<u8>, ops: &[u8], offset: usize, len: usize) -> Option<usize> {

		/* working variables */
		let mut buf: [u8; 96] = [0; 96];
		let mut ops = ops.iter();
		let mut ins = Op::Ins as u64;
		let mut ofs = offset;
		let mut rem = len;

		if rem > 32 {
			/* decode head block */
			let op = ops.next()?;
			let (block_len, next_ins) = unsafe {
				Self::decode_core_block(&mut buf, *op as u32, ins)
			};

			let block_len = block_len - ofs;
			dst.write(&buf[ofs .. ofs + block_len]).ok()?;		/* I expect it never fails though... */

			rem -= block_len;
			ins = next_ins;			/* just copy to mutable variable for use in the loop */
			ofs = 0;				/* clear offset for tail */

			/* decode body blocks */
			while rem > 32 {
				let op = ops.next()?;
				let (block_len, next_ins) = unsafe {
					Self::decode_core_block(&mut buf, *op as u32, ins)
				};

				dst.write(&buf[0..32]).ok()?;
				rem -= block_len;
				ins = next_ins;
			}
		}

		/* tail */ {
			let end = ofs + rem;
			while ofs < end {
				let op = ops.next()?;
				let (block_len, next_ins) = unsafe {
					Self::decode_core_block(&mut buf[ofs ..], *op as u32, ins)
				};

				ofs += block_len;
				ins = next_ins;
			}

			dst.write(&buf[end - rem .. end]).ok()?;
		}
		Some(len)
	}

	/* Check sanity of the span. If queried span (range) is out of the indexed one, return None */
	fn check_span(&self, range: &Range<usize>) -> Option<()> {

		if range.end < range.start { return None; }
		if range.end > self.ref_span { return None; }
		return Some(());
	}

	/* fetch block head for this pos. */
	fn fetch_ops_block(&self, pos: usize) -> (&[u8], usize) {

		let block_index = pos / BLOCK_PITCH;
		let block_rem   = pos & (BLOCK_PITCH - 1);
		let block = &self.block[block_index];

		let op_offset  = block.op_offset() as usize;
		let op_skip    = block.op_skip() as usize;
		// let ins_offset = block.ins_offset() as usize;

		let ops = &self.op[op_offset ..];
		return (ops, block_rem + op_skip);
	}

	fn fetch_ins_block(&self, pos: usize) -> (&[u8], bool) {

		let block_index = pos / BLOCK_PITCH;
		let block = &self.block[block_index];

		let ins_offset = block.ins_offset() as usize;
		let skip = (block_index == 0) as bool;

		return (&self.ins[ins_offset ..], skip);
	}

	fn scan_op_array(&self, pos: usize) -> (&[u8], usize) {

		/* get block head for this pos */
		let (ops, rem) = self.fetch_ops_block(pos);

		/* linear polling */
		let mut ops = ops.iter();
		let ofs = (&mut ops).peek_fold(0, |a, &x| {
			let len = a + op_len(x);

			/* continue if at least one column remaining */
			if len >= rem { return None; }
			return Some(len);
		});
		return (ops.as_slice(), rem - ofs);		/* is this sound? (I'm not sure...) */
	}


	/* ins */
	fn scan_ins_array(&self, pos: usize) -> Option<&[u8]> {

		let (ops, rem)  = self.fetch_ops_block(pos);
		let (ins, skip) = self.fetch_ins_block(pos);

		/* linear polling on op array */
		let mut ops = ops.iter();
		let (len, count) = (&mut ops).peek_fold((0, skip as usize), |(a, c), &x| {
			let len = a + op_len(x);
			if len >= rem { return None; }

			let is_ins = (op_marker(x) == CompressMark::Ins as u32) as usize;
			return Some((len, c + is_ins));
		});

		/* check if the column has insertion */
		if len > rem { return None; }

		/* linear polling on ins array */
		let mut ins = ins.iter();
		ins.peek_fold(count, |a, &x| {
			if a == 0 { return None; }

			return Some(a - (x == 0) as usize);
		});

		let ins = ins.as_slice();
		let len = ins.clone().iter().peek_fold(0, |a, &x| {
			if x == 0 { return None; }
			return Some(a + 1);
		});
		return Some(&ins[.. len]);
	}

	fn get_ins_core(dst: &mut Vec<u8>, ins: &[u8]) -> usize {
		/* expand packed nucleotide to Vec */
		let range = dst.reserve_to(ins.len() * 2, |arr: &mut [u8], _: &[u8]| -> usize {
			for (i, x) in arr.iter_mut().enumerate() {
				*x = if (i & 0x01) == 0 {
					ins[i / 2] & 0x0f
				} else {
					ins[i / 2]>>4
				};
			}

			let remove_tail = ins[ins.len() - 1] == 0;
			arr.len() - remove_tail as usize
		});

		dst.resize(range.end, 0);
		return range.end - range.start;
	}

	pub fn get_ins_into(&self, dst: &mut Vec<u8>, pos: usize) -> Option<usize> {
		self.check_span(&Range::<usize> { start: 0, end: pos })?;

		let ins = self.scan_ins_array(pos)?;
		let len = Self::get_ins_core(dst, ins);

		Some(len)
	}

	pub fn get_ins(&self, pos: usize) -> Option<Vec<u8>> {
		self.check_span(&Range::<usize> { start: 0, end: pos })?;

		/* fetch ins vector */
		let ins = self.scan_ins_array(pos)?;

		/* expand packed nucleotide to Vec */
		let mut v = Vec::with_capacity(ins.len() * 2);
		Self::get_ins_core(&mut v, ins);

		Some(v)
	}
}


/* UdonVec builder, see the comment above */
struct UdonPrecursorVec {
	buf: Vec<u8>,
	precursors: Vec<u8>
}
struct UdonVec {
	buf: Vec<u8>,
	precursors: Vec<u8>
}

impl<'a, 'b> UdonVec<'a, 'b> {
	/* returns reference side span */
	pub fn append(&mut self, cigar: &'b [u32], packed_query: &'b [u8], mdstr: &'b [u8]) -> Option<usize> {

		/* move buf to push precursor content */
		let (buf, precursor) = unsafe { UdonPrecursor::build(self.buf, cigar, packed_query, mdstr) };

		/* put back */
		self.buf = buf;
		let precursor = precursor?;
		let ref_span = precursor.ref_span;
		self.precursors.push(precursor);

		return Some(ref_span);
	}

	unsafe fn from_precursor_vec(buf: Vec<u8>, precursors: Vec<UdonPrecursor>) -> Vec<Udon<'a>> {

		let base: *const u8 = buf.as_ptr();
		for mut precursor in precursors.iter_mut() {		/* consume? */

			/* compose udon header on stack */
			let header = Udon::<'a>::from_precursor_raw(base, &precursor);

			/* copy back udon on stack to vec */
			let src = &header as *const Udon<'a>;
			let dst = &mut precursor as &mut UdonPrecursor as *mut UdonPrecursor;

			let dst = transmute::<*mut UdonPrecursor, *mut Udon<'a>>(dst);
			copy_nonoverlapping(src, dst, 1);
		}

		/* is there simpler way? */
		let ptr = buf.as_mut_ptr();
		let len = buf.len();
		let cap = buf.capacity();

		Vec::<Udon<'a>>::from_raw_parts(ptr as *mut Udon<'a>, len, cap)

	}

	pub fn freeze(self) -> Vec<Udon<'a>> {

		/* get precursor count before merging two vectors */
		let mut precursors = self.precursors;
		let count = precursors.len() / size_of::<UdonPrecursor>();

		/* buf is consumed */
		let buf = self.buf;
		(&mut precursors).reserve_to(buf.len(), |arr: &mut [u8], _: &[u8]| -> usize {
			arr.copy_from_slice(&buf);
		});

		precursors.forget();
		let 

	}

}


/*
impl Transcode for Udon<'_> {
	fn transcode<'a>(cigar: &[Cigar], packed_query: &[u8], md: &[u8]) -> Option<Box<Udon<'a>>> {

		// let mut builder = UdonBuilder::new(cigar, packed_query, md);
		// let _ = builder.parse_cigar()?;
		// let udon = builder.finalize()?;
		return None;
	}
}
*/


#[cfg(test)]
mod test {
	use std::ops::Range;
	use std::str::from_utf8;
	use crate::{ Udon, UdonPrecursor, CigarOp, encode_base_unchecked };

	macro_rules! cigar {
		[ $( ( $op: ident, $len: expr ) ),* ] => ({
			vec![ $( CigarOp::$op as u32 | (( $len as u32 )<<4) ),* ]
		});
	}

	macro_rules! nucl {
		( $st: expr ) => ({
			let s = $st;
			let mut v = Vec::<u8>::new();
			let mut a: u8 = 0;
			for (i, c) in s.bytes().enumerate() {
				if (i & 0x01) == 0x01 {
					v.push(a | encode_base_unchecked(c as char));
				} else {
					a = encode_base_unchecked(c as char)<<4;
				}
			}

			if (s.len() & 0x01) == 0x01 {
				v.push(a);
			}
			v
		});
	}

	macro_rules! encode_flat {
		( $arr: expr ) => ({
			let mut v = Vec::<u8>::new();
			for &x in arr {
				let c = match x {
					0x00 => 'M', 0x04 => 'A', 0x05 => 'C', 0x06 => 'G', 0x07 => 'T', 0x08 => 'D',
					0x10 => 'M', 0x04 => 'A', 0x05 => 'C', 0x06 => 'G', 0x07 => 'T', 0x08 => 'D'
				};
				v.push(c);
			}
			v
		});
	}

	macro_rules! compare {
		( $cigar: expr, $nucl: expr, $mdstr: expr, $flat: expr ) => ({
			// let v = Vec::<u8>::new();
			let u = Udon::build(&( $cigar ), &( $nucl ), &(( $mdstr ).as_bytes())).unwrap();
			let f = u.decode(Range { start: 0, end: u.ref_span() }).unwrap();
			let d = from_utf8(&f).unwrap();
			assert!(d == $flat);
		});
	}

	#[test]
	fn test_udon() {
		compare!(
			cigar![(Match, 4), (Ins, 1), (Match, 4), (Del, 1), (Match, 2), (Del, 7), (Match, 40)],
			nucl!("ACGTACGTACGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
			"9^A2^ACGTACG20A9A0C0G7",
			"MMMMMMMMDMMDDDDDDDMMMMMMMMMMMMMMMMMMMMAMMMMMMMMMACGMMMMMMM"
		);
		


		// assert!(false);
	}


}




