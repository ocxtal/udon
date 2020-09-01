
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

use std::f64::consts::PI;
use std::io::Write;
use std::marker::PhantomData;
use std::mem::{ size_of, transmute, forget };
use std::ops::{ Range, Add, Mul, AddAssign };
use std::ptr::copy_nonoverlapping;
use std::slice::{ Iter, IterMut, from_raw_parts, from_raw_parts_mut };


/* architecture-dependent stuffs */
#[cfg(all(target_arch = "x86_64"))]
use core::arch::x86_64::*;

#[cfg(all(target_arch = "aarch64"))]
use core::arch::aarch64::*;				/* requires nightly */


/* logging */
#[macro_use]
#[allow(unused_imports)]				/* we need identifier `debug`, but it'll be aliased to println */
extern crate log;

/*
macro_rules! debug {
	( $( $expr: expr ),* ) => ({
		println!("[{}]")
	});
}
*/
// #[cfg(test)]
// use std::println as debug;




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
	MisA = CompressMark::Mismatch as u32 | 0x00,
	MisC = CompressMark::Mismatch as u32 | 0x01,
	MisG = CompressMark::Mismatch as u32 | 0x02,
	MisT = CompressMark::Mismatch as u32 | 0x03,
	Del  = 0x08,
	Ins  = 0x10
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


/* Scaled rendering API


*/
#[derive(Copy, Clone, Debug)]
pub struct UdonPalette {
	/* all in (r, g, b, unused) form */
	background:  [u8; 4],
	del: [u8; 4],
	ins: [u8; 4],
	mismatch: [[u8; 4]; 4]
}

impl Default for UdonPalette {
	fn default() -> UdonPalette {
		UdonPalette {
			background: [0x00, 0x00, 0x00, 0x00],
			del: [0xff, 0xff, 0x22, 0x10],
			ins: [0x22, 0xff, 0x22, 0x10],
			mismatch: [
				[0xcc, 0xcc, 0xff, 0x10],
				[0xff, 0xff, 0x00, 0x10],
				[0xff, 0x33, 0x99, 0x10],
				[0x00, 0xff, 0x99, 0x10]
			]
		}
	}
}


/* UdonUtils

Provides utilities on ribbon (&mut [u32]): blending and gamma correction
*/
pub trait UdonUtils {
	fn append_on_basecolor(&mut self, basecolor: [u8; 4]) -> &mut Self;
	fn correct_gamma(&mut self) -> &mut Self;
}

impl UdonUtils for [u32] {

	fn append_on_basecolor(&mut self, basecolor: [u8; 4]) -> &mut Self {
		for x in self.iter_mut() {
			let mut v = x.to_le_bytes();
			for i in 0 .. 4 {
				v[i] = basecolor[i] - v[i].min(basecolor[i]);
			}
			*x = u32::from_le_bytes(v);
		}
		self
	}

	fn correct_gamma(&mut self) -> &mut Self {
		const GAMMA: [u8; 256] = [
			  0,  12,  21,  28,  33,  38,  42,  46,  49,  52,  55,  58,  61,  63,  66,  68,
			 70,  73,  75,  77,  79,  81,  83,  84,  86,  88,  90,  91,  93,  94,  96,  97,
			 99, 100, 102, 103, 105, 106, 107, 109, 110, 111, 113, 114, 115, 116, 118, 119,
			120, 121, 122, 123, 124, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136,
			137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 146, 147, 148, 149, 150, 151,
			152, 153, 153, 154, 155, 156, 157, 158, 159, 159, 160, 161, 162, 163, 163, 164,
			165, 166, 166, 167, 168, 169, 169, 170, 171, 172, 172, 173, 174, 175, 175, 176,
			177, 178, 178, 179, 180, 180, 181, 182, 182, 183, 184, 184, 185, 186, 186, 187,
			188, 188, 189, 190, 190, 191, 192, 192, 193, 194, 194, 195, 195, 196, 197, 197,
			198, 199, 199, 200, 200, 201, 202, 202, 203, 203, 204, 205, 205, 206, 206, 207,
			207, 208, 209, 209, 210, 210, 211, 211, 212, 213, 213, 214, 214, 215, 215, 216,
			216, 217, 218, 218, 219, 219, 220, 220, 221, 221, 222, 222, 223, 223, 224, 224,
			225, 226, 226, 227, 227, 228, 228, 229, 229, 230, 230, 231, 231, 232, 232, 233,
			233, 234, 234, 235, 235, 236, 236, 237, 237, 238, 238, 238, 239, 239, 240, 240,
			241, 241, 242, 242, 243, 243, 244, 244, 245, 245, 246, 246, 246, 247, 247, 248,
			248, 249, 249, 250, 250, 251, 251, 252, 252, 252, 253, 253, 254, 254, 255, 255
		];

		for x in self.iter_mut() {
			let mut v = x.to_le_bytes();
			for i in 0 .. 4 {
				v[i] = GAMMA[v[i] as usize];
			}
			*x = u32::from_le_bytes(v);
		}
		self
	}
}



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
	_marker: PhantomData<T>			/* not effectively used for now. what to do with this? */
}

#[allow(dead_code)]
impl<'a, T> SlicePrecursor<T> {
	fn compose(range: &Range<usize>) -> Self {
		SlicePrecursor::<T> {
			ofs: range.start,
			len: range.end - range.start,
			_marker: PhantomData::<T>
		}
	}

	fn finalize_raw(&self, base: *const u8) -> &'a [T] {
		let ptr = base.wrapping_add(self.ofs) as *const T;
		let cnt = self.len / size_of::<T>();
		unsafe { from_raw_parts(ptr, cnt) }
	}

	fn finalize(&self, v: &'a Vec<u8>) -> &'a [T] {
		let base = v.as_ptr() as *const u8;
		self.finalize_raw(base)
	}
}


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
	/* one-hot encoding */
	match "-AC-G---T-------".bytes().nth(c as usize) {
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

fn op_is_cont(x: u8) -> bool {
	(x & 0x1f) == 0x1f
}


/* Insertion marker tracker

CompressMark::Ins should be ignored when the succeeding op length is 31,
which actually indicates its chunk length is 30. The length 31 appears when
longer matching region is divided into multiple chunks, so it shouldn't be
regarded as insertion even if the marker is Ins.

InsTracker handles this double-meaning Ins marker. It keeps "ignore" state,
which indicates Ins marker of the next op is treated NOP. It also keeps
track of the inserted sequence array, by eating each variable-length sub-
string on the array every time Ins marker apperared on the op stream. The
current substring (inserted sequence for the current op with Ins marker)
is available via `as_slice()` method.
*/
struct InsTracker<'a> {
	ignore_next: bool,
	len: usize,
	ins: Iter<'a, u8>
}

impl<'a> InsTracker<'a> {
	fn new(last_op: u8, ins: &'a [u8]) -> Self {
		InsTracker::<'a> {
			ignore_next: op_is_cont(last_op),
			len: ins.len(),
			ins: ins.iter()
		}
	}

	fn forward(&mut self, op: u8) -> bool {
		/* check if the current op has ins, and update state for the next op */
		let is_ins = !self.ignore_next && op_marker(op) == CompressMark::Ins as u32;
		self.ignore_next = op_is_cont(op);

		/* if current does not have ins, just return it's not ins */
		if !is_ins { return false; }

		/* has ins, forward iterator */
		self.ins.try_fold(1, |a, &x| {			/* consumes at least one byte */
			let a = a - (x == 0) as usize;
			if a == 0 { return None; }
			return Some(a);
		});
		return true;
	}

	fn get_offset(&self) -> usize {
		self.len - self.ins.as_slice().len()	/* is there better way for tracking the offset? */
	}

	fn as_slice(&self) -> &'a [u8] {
		self.ins.as_slice()
	}
}


/* OpsIterator

*/
#[allow(dead_code)]
struct OpsIter<'a> {
	iter: Iter<'a, u8>,
	rofs: usize
}

trait IntoOpsIterator<'a> {
	fn iter_ops(self, base_rofs: usize) -> OpsIter<'a>;
}

impl<'a> IntoOpsIterator<'a> for &'a [u8] {
	fn iter_ops(self, base_rofs: usize) -> OpsIter<'a> {
		OpsIter {
			iter: self.iter(),
			rofs: base_rofs
		}
	}
}

impl<'a> Iterator for OpsIter<'a> {
	type Item = (u8, usize);

	#[allow(dead_code)]
	fn next(&mut self) -> Option<Self::Item> {
		let x    = self.iter.next()?;
		self.rofs += op_len(*x);
		Some((*x, self.rofs))
	}
}


/* UdonScaler


*/
#[derive(Copy, Clone, Debug, Default)]
struct Color {
	v: [i32; 4]
}

impl From<&[u8; 4]> for Color {
	fn from(val: &[u8; 4]) -> Color {
		let mut x = Color::default();
		for i in 0 .. 4 { x.v[i] = val[i] as i32; }
		x
	}
}

impl From<&Color> for [u8; 4] {
	fn from(val: &Color) -> [u8; 4] {
		let mut x: [u8; 4] = Default::default();
		for i in 0 .. 4 {
			x[i] = val.v[i].min(255).max(0) as u8;
		}
		x
	}
}

impl From<&Color> for u32 {
	fn from(val: &Color) -> u32 {
		let mut x: [u8; 4] = Default::default();
		for i in 0 .. 4 {
			x[i] = val.v[i].min(255).max(0) as u8;
		}
		u32::from_le_bytes(x)
	}
}

impl Add<Color> for Color {
	type Output = Color;
	fn add(self, other: Color) -> Color {
		let mut x = self;
		for i in 0 .. 4 { x.v[i] += other.v[i]; }
		x
	}
}
impl Mul<Color> for Color {
	type Output = Color;
	fn mul(self, other: Color) -> Color {
		let mut x = self;
		for i in 0 .. 4 {
			/* upcast, multiply, then downcast. (expect pmuludq) */
			let n = x.v[i] as i64 * other.v[i] as i64;
			x.v[i] = (n>>24) as i32;		/* I expect this won't overflow but no guarantee */
		}
		x
	}
}

impl Add<i32> for Color {
	type Output = Color;
	fn add(self, other: i32) -> Color {
		let mut x = self;
		for i in 0 .. 4 { x.v[i] += other; }
		x
	}
}
impl Mul<i32> for Color {
	type Output = Color;
	fn mul(self, other: i32) -> Color {
		let mut x = self;
		for i in 0 .. 4 {
			let n = x.v[i] as i64 * other as i64;
			x.v[i] = (n>>24) as i32;
		}
		x
	}
}

/* arithmetic assign */
impl AddAssign<Color> for Color {
	fn add_assign(&mut self, other: Color) { *self = self.add(other); }
}
impl AddAssign<i32> for Color {
	fn add_assign(&mut self, other: i32) { *self = self.add(other); }
}


/* UdonScaler second impl */
#[derive(Default)]
pub struct UdonScaler {
	columns_per_pixel: f64,
	window: f64,
	offset: f64,
	color: [Color; 12],
	table: [Vec<i32>; 17]
}

fn sinc(rad: f64) -> f64 {
	if !rad.is_normal() {
		return 1.0;
	}
	rad.sin() / rad
}

fn sincx(x: f64, order: f64) -> f64 {
	sinc(PI * x) * sinc(PI * (x / order).max(-1.0).min(1.0))
}

fn clip(x: f64, window: f64) -> f64 {
	if x > 0.0 {
		return if x <  window { 0.0 } else { x - window };
	} else {
		return if x > -window { 0.0 } else { x + window };
	}
}

impl UdonScaler {

	/* column -> color index */
	const INDEX: [u8; 32] = {
		let mut index = [0; 32];
		index[Op::MisA as usize] = 1;
		index[Op::MisC as usize] = 2;
		index[Op::MisG as usize] = 3;
		index[Op::MisT as usize] = 4;
		index[Op::Del as usize]  = 5;
		index[Op::Ins as usize | Op::MisA as usize] = 6;
		index[Op::Ins as usize | Op::MisC as usize] = 7;
		index[Op::Ins as usize | Op::MisG as usize] = 8;
		index[Op::Ins as usize | Op::MisT as usize] = 9;
		index[Op::Ins as usize | Op::Del as usize]  = 10;
		index
	};
	fn index(column: u8) -> usize {
		assert!(column < 32, "{}", column);

		Self::INDEX[column as usize] as usize
	}

	fn pick_color(&self, column: u8) -> Color {
		self.color[Self::index(column)]
	}

	fn build_color_table(color: &UdonPalette) -> [Color; 12] {

		let mismatch: [Color; 4] = [
			Color::from(&color.mismatch[0]),
			Color::from(&color.mismatch[1]),
			Color::from(&color.mismatch[2]),
			Color::from(&color.mismatch[3])
		];
		let del = Color::from(&color.del);
		let ins = Color::from(&color.ins);
		let bg  = Color::from(&color.background);

		let mut x = [bg; 12];
		x[Self::index(Op::MisA as u8)] = mismatch[0];
		x[Self::index(Op::MisC as u8)] = mismatch[1];
		x[Self::index(Op::MisG as u8)] = mismatch[2];
		x[Self::index(Op::MisT as u8)] = mismatch[3];
		x[Self::index(Op::Del as u8)]  = del;
		x[Self::index(Op::Ins as u8 + Op::MisA as u8)] = ins * 0x800000 + mismatch[0] * 0x800000;
		x[Self::index(Op::Ins as u8 + Op::MisC as u8)] = ins * 0x800000 + mismatch[1] * 0x800000;
		x[Self::index(Op::Ins as u8 + Op::MisG as u8)] = ins * 0x800000 + mismatch[2] * 0x800000;
		x[Self::index(Op::Ins as u8 + Op::MisT as u8)] = ins * 0x800000 + mismatch[3] * 0x800000;
		x[Self::index(Op::Ins as u8 + Op::Del as u8)]  = ins * 0x800000 + del * 0x800000;

		x
	}

	const WINDOW: f64 = 1.0;

	fn build_coef_table(v: &mut Vec<i32>, i: usize, scale: f64, pitch: f64, width: f64) {
		let span = 2 * ((0.5 * scale).ceil() as usize);
		let offset = (i as f64 - 8.0) / 16.0;

		/* FIXME */
		let center = pitch * pitch * offset + span as f64 / 2.0;

		for j in 0 ..= span {
			let dist = center - j as f64;
			let coef = sincx(clip(dist, width) / pitch, Self::WINDOW);

			/* FIXME */
			let coef = coef.max(0.0);
			v.push((0x01000000 as f64 * coef) as i32);
		}
	}

	pub fn new(color: &UdonPalette, columns_per_pixel: f64) -> UdonScaler {
		let scale = columns_per_pixel.max(1.0);
		let pitch = columns_per_pixel / scale;
		let width = 0.5 * (scale - 1.0);

		let mut x = UdonScaler {
			columns_per_pixel: columns_per_pixel,
			window: scale.ceil() + 1.0,
			offset: (scale.ceil() + 1.0) / 2.0,
			color: Self::build_color_table(&color),
			table: Default::default()
		};

		for i in 0 .. 17 {
			Self::build_coef_table(&mut x.table[i], i, scale, pitch, width);
		}
		x
	}

	/*
	fn calc_coef(dist: f64, shoulder: f64, magnifier: f64) -> f64 {
		let width = magnifier - 1.0;

		if dist.abs() < magnifier / 4.0 {
			return 1.0;
		}
		0.5 * (
			  sincx((dist + shoulder) * magnifier, Self::WINDOW)
			+ sincx((dist - shoulder) * magnifier, Self::WINDOW)
		)
	}

	fn build_coef_table(v: &mut Vec<i32>, i: usize, magnifier: f64, shoulder: f64) {

		// println!("columns_per_pixel({:?}), scale({:?}), width({:?}), shoulder({:?})", columns_per_pixel, scale, pitch - 1.0, shoulder);
		let offset = i as f64 / 16.0;
		let center = offset + Self::WINDOW * magnifier;

		// let start = offset;
		let end  = offset + 2.0 * Self::WINDOW * magnifier;
		let span = end.ceil() as usize;

		// let mut x = Vec::new();
		for j in 0 ..= span {

			let dist = center - j as f64;
			let coef = Self::calc_coef(dist, shoulder, magnifier);
			// v.push((j, dist, format!("{:#.3}", coef)));
			v.push((0x01000000 as f64 * coef) as i32);
		}
		// println!("offset({:#.05}), center({:#.05}), r({:#.3}, {:#.3}), r({:#.3}, {:#.3}), {:?}", offset, center, start, end, start.floor() as i64, end.ceil() as i64, x);
	}

	pub fn new(color: &UdonPalette, columns_per_pixel: f64) -> UdonScaler {
		let scale = columns_per_pixel.max(1.0);
		let magnifier = scale / columns_per_pixel;

		let mut x = UdonScaler {
			columns_per_pixel: columns_per_pixel,
			window: (2.0 * Self::WINDOW) * magnifier,
			color: Self::build_color_table(&color),
			table: Default::default()
		};

		for i in 0 .. 17 {
			Self::build_coef_table(&mut x.table[i],
				i, magnifier, 0.25 / columns_per_pixel
			);
		}
		x
	}
	*/

	fn expected_size(&self, span: usize) -> usize {
		(span as f64 / self.columns_per_pixel) as usize + 2
	}

	fn init(&self, offset_in_pixels: f64) -> (f64, usize) {
		let offset = (offset_in_pixels + 0.5) * self.columns_per_pixel;
		let margin = (offset + self.offset) as usize;
		debug!("init, offset({}, {}), margin({})", offset, self.offset, margin);

		(offset, margin)
	}

	fn scale(&self, dst: &mut Vec<u32>, src: &[u8], offset: f64) -> Option<(f64, usize)> {

		debug!("scale, offset({})", offset);
		for i in 0 .. {
			/*
			offset := offset_in_columns
			*/
			let base = offset + (i as f64 * self.columns_per_pixel);
			let range = Range::<usize> {
				start: base as usize,
				end: (base + self.window + 1.0).ceil() as usize
			};
			debug!("base({}), range({:?})", base, range);

			if range.end > src.len() {
				return Some((base.fract(), range.start));
			}

			let table = &self.table[(base.fract() * 16.0) as usize];
			debug!("frac({}), {:?}", base.fract(), table);

			let mut a = Color::default();
			for (&coef, &column) in (&table).iter().zip((&src[range]).iter()) {
				debug!("col({}), color({:?}), coef({})", column, self.pick_color(column), coef);
				a += self.pick_color(column) * coef;
			}
			debug!("acc({:?})", a);
			dst.push(u32::from(&a));
		}

		return None;
	}
}


/*
/* UdonScaler and its impl */
struct UdonScaler {
	accum: Color,
	prev:  Color,
	color: [Color; 12],
	normalizer: u32
}

impl UdonScaler {

	/* column -> color index */
	const INDEX: [u8; 32] = {
		let mut index = [0; 32];
		index[Op::MisA as usize] = 1;
		index[Op::MisC as usize] = 2;
		index[Op::MisG as usize] = 3;
		index[Op::MisT as usize] = 4;
		index[Op::Del as usize]  = 5;
		index[Op::Ins as usize | Op::MisA as usize] = 6;
		index[Op::Ins as usize | Op::MisC as usize] = 7;
		index[Op::Ins as usize | Op::MisG as usize] = 8;
		index[Op::Ins as usize | Op::MisT as usize] = 9;
		index[Op::Ins as usize | Op::Del as usize]  = 10;
		index
	};
	fn index(column: u8) -> usize {
		assert!(column < 32, "{}", column);

		Self::INDEX[column as usize] as usize
	}

	/* fraction in margin -> coefficient */
	const COEF: [u32; 33] = [
		0x01000000,
		0x00ff47f5,
		0x00fd2228,
		0x00f99580,
		0x00f4ad63,
		0x00ee7987,
		0x00e70db5,
		0x00de8179,
		0x00d4efc8,
		0x00ca7697,
		0x00bf3666,
		0x00b351bf,
		0x00a6ecb7,
		0x009a2c61,
		0x008d3640,
		0x00802fba,
		0x00733d90,
		0x00668354,
		0x005a22e9,
		0x004e3c09,
		0x0042ebda,
		0x00384c89,
		0x002e74f8,
		0x00257873,
		0x001d6681,
		0x00164ab8,
		0x00102ca9,
		0x000b0fdf,
		0x0006f3e8,
		0x0003d473,
		0x0001a97f,
		0x00006790,
		0x00000000
	];
	fn coef(frac: f64) -> (u32, u32) {

		/*
		let index = (frac * 32.0) as usize;
		assert!(index < 33, "frac({}), index({})", frac, index);

		(Self::COEF[index], Self::COEF[32 - index])
		*/

		let coef = (0x01000000 as f64 * frac) as u32;
		(coef, 0x01000000 - coef)

	}


	fn new(color: &UdonPalette, pitch: f64) -> Self {
		// let normalizer = (0x01000000 as f64 * (pitch + 1.0).log(2.7) + 1.0) as u32;
		let normalizer = (0x01000000 as f64 * pitch) as u32;
		debug!("pitch({}), normalizer({})", pitch, (pitch + 1.0).log(2.7));

		UdonScaler {
			accum: Color::default(),
			prev:  Color::default(),
			normalizer: normalizer,
			color: Self::build_color_table(&color)
		}
	}

	fn build_color_table(color: &UdonPalette) -> [Color; 12] {

		let mismatch: [Color; 4] = [
			Color::from(&color.mismatch[0]),
			Color::from(&color.mismatch[1]),
			Color::from(&color.mismatch[2]),
			Color::from(&color.mismatch[3])
		];
		let del = Color::from(&color.del);
		let ins = Color::from(&color.ins);
		let bg  = Color::from(&color.background);

		let mut x = [bg; 12];
		x[Self::index(Op::MisA as u8)] = mismatch[0];
		x[Self::index(Op::MisC as u8)] = mismatch[1];
		x[Self::index(Op::MisG as u8)] = mismatch[2];
		x[Self::index(Op::MisT as u8)] = mismatch[3];
		x[Self::index(Op::Del as u8)]  = del;
		x[Self::index(Op::MisA as u8 + Op::Ins as u8)] = ins * 0x800000 + mismatch[0] * 0x800000;
		x[Self::index(Op::MisC as u8 + Op::Ins as u8)] = ins * 0x800000 + mismatch[1] * 0x800000;
		x[Self::index(Op::MisG as u8 + Op::Ins as u8)] = ins * 0x800000 + mismatch[2] * 0x800000;
		x[Self::index(Op::MisT as u8 + Op::Ins as u8)] = ins * 0x800000 + mismatch[3] * 0x800000;
		x[Self::index(Op::Del as u8 + Op::Ins as u8)]  = ins * 0x800000 + del * 0x800000;

		debug!("{:?}", x);
		x
	}

	fn accumulate(&mut self, column: u8) {
		let color = &self.color[Self::index(column)];

		/* let compiler vectorize this! */
		self.accum += *color;

		debug!("accumulate: color({:?}), accum({:?})", color, self.accum);

		/* copy to save */
		self.prev = *color;
	}

	fn interpolate(&mut self, column: u8, frac: f64) {

		let curr     = &self.color[Self::index(column)];
		let (c0, c1) = Self::coef(frac);

		/* let compiler vectorize this!! */
		self.accum += *curr     * c0;
		self.accum += self.prev * c1;

		debug!("interpolate: color({:?}, {:?}), frac({}), c({}, {}), accum({:?})", self.prev, curr, frac, c0, c1, self.accum);

		/* copy to save */
		self.prev = *curr;
	}

	fn flush(&mut self, dst: &mut Vec<u32>, cnt: usize) -> Option<()> {
		assert!(cnt > 0, "{}", cnt);

		/* I hope this conversions are automatically vectorized!!!! */
		let body = u32::from(&self.prev);
		let tail = u32::from(&(self.accum * self.normalizer));

		debug!("flush: body({:?}), tail({:?})", self.prev, self.accum * self.normalizer);

		for _ in 0 .. cnt - 1 { dst.push(body); }	/* broadcast for scale < 1.0 */
		dst.push(tail);

		/* clear all */
		self.accum = Color::default();
		return Some(());
	}

	fn scale(&mut self, dst: &mut Vec<u32>, src: &[u8], offset: f64, pitch: f64, margin: f64) -> Option<f64> {

		/* returns new offset on success, None on failure */

		debug!("offset({}), pitch({}), margin({})", offset, pitch, margin);
		debug!("{:?}", src);

		assert!(!(offset < 0.0));
		assert!(pitch  > 0.0);
		assert!(margin > 0.0);

		let rev_pitch  = 1.0 / pitch;
		let rev_margin = 1.0 / margin;
		let mut last_bin = 0;
		let mut dst = dst;

		for (i, &column) in src.iter().enumerate() {
			let pos = i as f64 + offset;
			let bin = (pos * rev_pitch).floor();
			let thresh = bin * pitch + 1.0;

			debug!("i({}), pos({}), bin({}), frac({}), thresh({}), column({}), index({}), color({:?})", i, pos, bin, (pos - bin * pitch) * rev_margin, thresh, column, Self::index(column), &self.color[Self::index(column)]);

			/* flush if needed */ {
				let bin = bin as usize;
				if bin != last_bin {
					self.flush(&mut dst, bin - last_bin)?;
				}
				last_bin = bin;
			}

			/* do interpolation if on margin */
			if pos < thresh {
				/* relative position in margin: [0.0, 1.0) */
				self.interpolate(column, (pos - bin * pitch) * rev_margin);
				continue;
			}

			/* just accumulate otherwise */
			self.accumulate(column);
		}
		self.flush(&mut dst, 1);

		/* compute new offset */
		let last_pos = src.len() as f64 + offset;
		let last_bin = (last_pos / pitch).floor();
		return Some(last_pos - last_bin * pitch);
	}
}
*/



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
		/* trivial */
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

	// debug!("{}, {}, {:?}", ofs, len, src);

	/* very very very very very very very very very naive way though fine */
	for i in 0 .. len {
		let pos = ofs + i;

		/* bam packed nucl is big endian */
		let c = if (pos & 0x01) == 0 {
			src[pos / 2]>>4
		} else {
			src[pos / 2] & 0x0f
		};

		/* store to dst in little endian */
		if (i & 0x01) == 0 {
			dst[i / 2]  = c;
		} else {
			dst[i / 2] |= c<<4;
		}

		// debug!("push {} at {}, ({}, {})", c, i, pos, src[pos / 2]);
	}

	// debug!("{:?}", &dst[0 .. (len + 1) / 2]);
}


/* transcoder state object */
pub struct UdonBuilder<'a> {
	buf: Vec<u8>,
	ins: Vec<u8>,
	cigar: Iter<'a, Cigar>,
	mdstr: Iter<'a, u8>,
	query: &'a [u8],
	qofs: usize
}

impl<'a, 'b> UdonBuilder<'a> {

	/*
	 * input parsers
	 */

	/* imitate iterator on 4-bit packed nucleotide sequence */
	fn next_base(&mut self) -> u32 {
		let ofs = self.qofs;
		self.qofs += 1;

		// debug!("qlen({}), ofs({})", self.query.len(), ofs);
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

		/* if the last op is mismatch, MD string ends like "A0" */
		if self.mdstr.as_slice().len() < 3 {
			return false;
		}
		return self.mdstr.as_slice()[1] == '0' as u8;
	}

	/* writers, on top of impl Writer for Vec<u8> */
	fn push_op(&mut self, match_len: usize, marker: u32) {
		assert!(marker    < 0x08);
		assert!(match_len < 0x20);

		// debug!("push_op, len({}), marker({})", match_len, marker);

		let op = (marker<<5) as u8 | match_len as u8;
		// assert!(op != 0);

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

		// debug!("eat_del, len({}), rem({})", len, rem);

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
		// debug!("eat_ins, xrem({}), len({}), qofs({})", xrem, len, self.qofs);

		self.save_ins(ofs, len);		/* use offset before forwarding */
		return Some(xrem);				/* there might be remainder */
	}
	fn eat_match(&mut self, xrem: usize, last_op: u32) -> Option<usize> {

		/* consumes only when the op is match. works as an escape route for del after ins or ins after del */
		let c = self.peek_cigar()?;
		let is_valid = c.op() == CigarOp::Match as u32;
		if is_valid { self.cigar.next()?; }

		/* forward qofs before adjusting xrem with the previous op */
		self.qofs += xrem;

		/* adjust crem and xrem; last_op == 0 for ins, > 0 for del */
		let mut crem = if is_valid { c.len() as usize } else { 0 } + last_op as usize;
		let mut xrem = xrem + last_op as usize;		/* might continues from the previous cigar op, possibly insertion */
		let mut op = last_op;						/* insertion, deletion, or mismatch */

		// debug!("eat_match, crem({}), xrem({})", crem, xrem);

		while xrem < crem {
			/* xrem < crem indicates this op (cigar span) is interrupted by mismatch(es) at the middle */
			self.push_match(xrem, op);
			crem -= xrem;
			// debug!("eat_match mismatch?, crem({}), xrem({}), qofs({})", crem, xrem, self.qofs);

			while self.is_double_mismatch() {
				// debug!("eat_match, crem({}), {:?}", crem, from_utf8(self.mdstr.as_slice()));

				let c = self.next_base();
				self.push_op(1, c);		/* this chunk contains only a single mismatch */
				self.mdstr.nth(1)?;
				crem -= 1;
			}
			// debug!("eat_match, crem({}), {:?}", crem, from_utf8(self.mdstr.as_slice()));

			op = self.next_base();		/* we only have a single mismatch remaining, will be combined to succeeding matches */
			self.mdstr.nth(0)?;

			/*
			xrem for the next { match, insertion } region, including +1 for the last mismatch.
			adjustment is already done on crem, by not decrementing it for the last mismatch.
			*/
			xrem = self.eat_md_eq() + 1;
			self.qofs += xrem - 1;
			// debug!("eat_match, updated xrem, crem({}), xrem({})", crem, xrem);
		}

		self.push_match(crem, op);		/* tail match; length is the remainder of crem */
		xrem      -= crem;
		self.qofs -= xrem;

		// debug!("eat_match, done, crem({}), xrem({}), qofs({})", crem, xrem, self.qofs);
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
		// debug!("c({}, {})", c.op(), c.len());
		if c.op() == CigarOp::Ins as u32 {
			xrem = self.eat_md_eq();

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
				if c.op() != CigarOp::Del as u32 { break; }
				// debug!("op({}), len({}), remaining cigars({})", c.op(), c.len(), self.cigar.as_slice().len());

				/* the CIGAR ends with deletion; must be treated specially */
				if self.cigar.as_slice().len() < 2 { break 'outer; }

				/* push deletion-match pair, then parse next eq length */
				let op = self.eat_del()?;
				xrem = self.eat_md_eq();
				// debug!("eat_del done, op({}), xrem({})", op, xrem);
				xrem = self.eat_match(xrem, op)?;
			}

			/* it's insertion-match pair when it appeared not be deletion-match */
			let c = peek_or_break!(self);
			if c.op() != CigarOp::Ins as u32 {
				return None;			/* if not, we regard it broken */
			}
			// debug!("op({}), len({}), remaining cigars({})", c.op(), c.len(), self.cigar.as_slice().len());

			/* the CIGAR ends with insertion; must be treated specially */
			if self.cigar.as_slice().len() < 2 { break 'outer; }

			/* push insertion-match pair, update eq length remainder */
			xrem = self.eat_ins(xrem)?;
			xrem = self.eat_match(xrem, CompressMark::Ins as u32)?;
		}

		/* CIGAR ends with isolated insertion or deletion */
		if self.cigar.as_slice().len() > 0 {
			let c = self.peek_cigar()?;
			// debug!("c({}, {})", c.op(), c.len());

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

		/* index for regular pitch on span */
		let block_count = (ref_span + BLOCK_PITCH - 1) / BLOCK_PITCH;
		let range = buf.reserve_to(block_count, |block: &mut [Block], base: &[u8]| {

			/* src: both are &[u8] and placed within base */
			let ops = (&base[op_range.start .. op_range.end]).iter_ops(0);
			let mut ins = InsTracker::new(0, &base[ins_range.start .. ins_range.end]);

			/* prepare dst. put head boundary info for simplicity */
			let mut dst = block.iter_mut();
			Self::push_block(&mut dst, 0, 0, 0);

			/* I want this loop be more lightweight... */
			let mut rbnd: usize = BLOCK_PITCH;
			for (i, (x, rpos)) in ops.enumerate() {
				let rem  = rbnd - (rpos - op_len(x));

				/* forward insertion array */
				let iofs = ins.get_offset();
				ins.forward(x);

				/* if forwarded rpos doexn't exceed the next boundary, just skip this op */
				if rpos <= rbnd { continue; }
				rbnd += BLOCK_PITCH;

				/* boundary found; save block info */
				Self::push_block(&mut dst, iofs, i, rem);
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
			ref_span: ref_span,

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
			qofs:  qclip.head		/* initial offset for query sequence, non-zero for soft clipped alignments */
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

		/* parse input op array */
		let op = unwrap_or_unwind!(state.parse_cigar());
		let ref_span = state.calc_reference_span(&op);

		/* pack ins vector */
		let ins = unwrap_or_unwind!(state.pack_ins());

		/* build index for op and ins arrays */
		let block = unwrap_or_unwind!(state.build_index(ref_span, &op, &ins));

		/* everything done; compose precursor */
		let precursor = state.finalize(ref_span, &op, &ins, &block);

		/* return back ownership of the buffer */
		let buf = state.buf;
		(buf, Some(precursor))
	}

	pub unsafe fn build(buf: Vec<u8>, cigar: &'a [u32], packed_query: &'a [u8], mdstr: &'a [u8]) -> (Vec<u8>, Option<UdonPrecursor>) {
		
		/* for compatibility with bam streamed reader */
		let cigar = transmute::<&'a [u32], &'a [Cigar]>(cigar);

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

	fn compose_box(buf: Vec<u8>, precursor: UdonPrecursor) -> Box<Udon<'a>> {
		let mut buf = buf;

		/* compose pointer-adjusted header on stack */
		let base: *mut u8 = buf.as_mut_ptr();
		let header = unsafe { Self::from_precursor_raw(base as *const u8, &precursor) };

		/* compose box and copy header into it */
		let udon = unsafe {
			/* convert buffer (on heap) to box */
			let mut udon = Box::<Udon<'a>>::from_raw(
				transmute::<*mut u8, *mut Udon<'a>>(base)
			);

			/* copy header from stack to heap (box) */
			let src = &header as *const Udon<'a>;
			let dst = &mut udon as &mut Udon<'a> as *mut Udon<'a>;
			copy_nonoverlapping(src, dst, 1);

			udon
		};

		/* heap block inside buf was moved to udon so we have to release buf */
		forget(buf);		/* note: this allowed outside unsafe */

		udon
	}


	/*
	decoder implementations
	*/
	pub fn ref_span(&self) -> usize {
		self.ref_span
	}

	pub fn decode_raw_into(&self, dst: &mut Vec<u8>, ref_span: &Range<usize>) -> Option<usize> {
		self.check_span(&ref_span)?;
		self.decode_core(dst, &ref_span)
	}

	pub fn decode_raw(&self, ref_span: &Range<usize>) -> Option<Vec<u8>> {
		self.check_span(&ref_span)?;

		let size = ref_span.end - ref_span.start;
		let mut buf = Vec::<u8>::with_capacity(size);

		let used = self.decode_core(&mut buf, &ref_span)?;
		buf.resize(used, 0);

		return Some(buf);
	}

	pub fn decode_scaled_into(&self, dst: &mut Vec<u32>, ref_span: &Range<usize>, offset_in_pixels: f64, scaler: &UdonScaler) -> Option<usize> {
		self.check_span(&ref_span)?;
		self.decode_scaled_core(dst, &ref_span, offset_in_pixels, &scaler)
	}

	pub fn decode_scaled(&self, ref_span: &Range<usize>, offset_in_pixels: f64, scaler: &UdonScaler) -> Option<Vec<u32>> {
		self.check_span(&ref_span)?;

		let span = ref_span.end - ref_span.start;
		let size = scaler.expected_size(span);
		let mut buf = Vec::<u32>::with_capacity(size);
		// debug!("ref_span({:?}), span({}), size({})", ref_span, span, size);

		let used = self.decode_scaled_core(&mut buf, &ref_span, offset_in_pixels, &scaler)?;
		buf.resize(used, 0);
		// debug!("used({})", used);

		return Some(buf);
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
	unsafe fn decode_core_block(dst: &mut [u8], op: u8, ins: u32) -> (usize, u32) {

		/* load constants; expelled out of the innermost loop when inlined */
		let del_mask      = _mm_load_si128(&Self::DEL_MASK.v as *const [u8; 16] as *const __m128i);
		let scatter_mask  = _mm_load_si128(&Self::SCATTER_MASK.v as *const [u8; 16] as *const __m128i);
		let is_del_thresh = _mm_load_si128(&Self::IS_DEL_THRESH.v as *const [u8; 16] as *const __m128i);

		/* compute deletion mask */
		let xop = _mm_cvtsi64_si128(op as i64);					/* copy op to xmm */
		let xop = _mm_shuffle_epi8(xop, scatter_mask);			/* [op, op, op, 0, ...] */
		let is_del = _mm_cmpgt_epi8(xop, is_del_thresh);		/* signed comparison */

		/* compute mismatch / insertion mask */
		let marker = if op_marker(op) == CompressMark::Ins as u32 { ins } else { op_marker(op) };
		let ins_mismatch_mask = _mm_cvtsi64_si128(marker as i64);

		/* merge deletion / insertion-mismatch vector */
		let merged = _mm_blendv_epi8(ins_mismatch_mask, del_mask, is_del);

		_mm_storeu_si128(&mut dst[0 .. 16] as *mut [u8] as *mut __m128i, merged);
		_mm_storeu_si128(&mut dst[16 .. 32] as *mut [u8] as *mut __m128i, _mm_setzero_si128());

		/*
		compute forward length; 31 is "continuous marker"
		rop-to-rop critical path length is 6
		*/
		let next_ins     = if op_is_cont(op) { 0 } else { Op::Ins as u32 };	/* next ins will be masked if 0x1f */
		let adjusted_len = op_len(op);

		// debug!("{:#x}, {:#x}, {:#x}, {:#x}", op>>5, ins, marker, op_len(op));
		// debug!("{}, {:?}", adjusted_len, transmute::<__m128i, [u8; 16]>(merged));

		(adjusted_len, next_ins)
	}

	#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
	unsafe fn decode_core_block(dst: &mut [u8], op: u32, ins: u64) -> (usize, u64) {
		(0, 0)
	}

	fn decode_core(&self, dst: &mut Vec<u8>, ref_span: &Range<usize>) -> Option<usize> {

		/*
		we suppose the ref_span is sane.
		linear polling from block head, then decode until the end of the span
		*/
		let len = ref_span.end - ref_span.start;
		let (ops, offset) = self.scan_op_array(ref_span.start);

		/* working variables */
		let mut buf: [u8; 96] = [0; 96];
		let mut ops = ops.iter();
		let mut ins = 0; // Op::Ins as u64;
		let mut ofs = offset;
		let mut rem = len;

		if rem > 32 {
			/* decode head block */
			let op = ops.next()?;
			let (block_len, next_ins) = unsafe {
				Self::decode_core_block(&mut buf, *op, ins)
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
					Self::decode_core_block(&mut buf, *op, ins)
				};

				dst.write(&buf[0 .. block_len]).ok()?;
				rem -= block_len;
				ins = next_ins;
			}
		}

		/* tail */ {
			let end = ofs + rem;
			let mut ofs = 0;
			while ofs < end {
				let op = ops.next()?;
				let (block_len, next_ins) = unsafe {
					Self::decode_core_block(&mut buf[ofs ..], *op, ins)
				};

				ofs += block_len;
				ins = next_ins;
			}

			dst.write(&buf[end - rem .. end]).ok()?;
		}

		Some(len)
	}

	fn decode_scaled_core(&self, dst: &mut Vec<u32>, ref_span: &Range<usize>, offset_in_pixels: f64, scaler: &UdonScaler) -> Option<usize> {

		/* states (working variables) */
		let (mut offset, margin) = scaler.init(offset_in_pixels);
		let mut dst = dst;

		/* buffer */
		let bulk_size: usize = 4 * 1024;		/* 4KB */
		let mut buf = Vec::<u8>::with_capacity(bulk_size + margin);
		for _ in 0 .. margin { buf.push(0); }

		for spos in ref_span.clone().step_by(bulk_size) {
			/* decode a block */
			let used = self.decode_core(&mut buf, &Range::<usize> {
				start: spos,
				end:   ref_span.end.min(spos + bulk_size)
			})?;
			debug!("{:?}", &buf[margin ..]);
			if used < bulk_size {
				for _ in 0 .. margin + 1 { buf.push(0); }
			}

			/* rescale to dst array, forward offset */
			let (next_offset, used) = scaler.scale(&mut dst, &buf, offset)?;
			offset = next_offset;

			if used < buf.len() - used { continue; }

			let (base, tail) = buf.split_at_mut(used);
			let (base, _) = base.split_at_mut(tail.len());
			base.copy_from_slice(tail);
		}
		Some(dst.len())
	}

	/* Check sanity of the span. If queried span (range) is out of the indexed one, return None */
	fn check_span(&self, range: &Range<usize>) -> Option<()> {

		// debug!("span({}, {}), ref_span({})", range.start, range.end, self.ref_span);

		if range.end < range.start { return None; }
		if range.end > self.ref_span { return None; }
		return Some(());
	}

	/* fetch block head for this pos. */
	fn fetch_ops_block(&self, pos: usize) -> (u8, &[u8], usize) {

		let block_index = pos / BLOCK_PITCH;
		let block_rem   = pos & (BLOCK_PITCH - 1);
		let block = &self.block[block_index];

		let op_offset  = block.op_offset() as usize;
		let op_skip    = block.op_skip() as usize;
		// let ins_offset = block.ins_offset() as usize;

		let ops = &self.op[op_offset ..];
		let last_op = if op_offset == 0 {
			0
		} else {
			self.op[op_offset - 1]
		};

		// debug!("pos({}), rem({}), skip({})", pos, block_rem, op_skip);
		return (last_op, ops, block_rem + op_skip);
	}

	fn fetch_ins_block(&self, pos: usize) -> &[u8] {

		let block_index = pos / BLOCK_PITCH;
		let block = &self.block[block_index];

		let ins_offset = block.ins_offset() as usize;
		return &self.ins[ins_offset as usize ..];
	}

	fn scan_op_array(&self, pos: usize) -> (&[u8], usize) {

		/* get block head for this pos */
		let (_, ops, rem) = self.fetch_ops_block(pos);
		// debug!("rem({}), ops({:?})", rem, ops);

		/* linear polling */
		let mut ops = ops.iter();
		let ofs = (&mut ops).peek_fold(0, |a, &x| {
			let len = a + op_len(x);

			/* continue if at least one column remaining */
			if len >= rem { return None; }
			return Some(len);
		});

		// debug!("rem({}), ofs({})", rem, ofs);
		return (ops.as_slice(), rem - ofs);		/* is this sound? (I'm not sure...) */
	}


	/* ins */
	fn scan_ins_array(&self, pos: usize) -> Option<&[u8]> {

		let (last_op, ops, rem) = self.fetch_ops_block(pos);
		let ins = self.fetch_ins_block(pos);

		/* linear polling on op array */
		let mut ops = ops.iter();
		let ins = InsTracker::new(last_op, ins);
		let len = (&mut ops).peek_fold(0, |a, &x| {
			let len = a + op_len(x);
			if len > rem { return None; }
			return Some(len);
		});

		/* if the length doesn't match, it indicates the column doesn't have insertion (so the query is wrong) */
		if len < rem { return None; }

		/* insertion found. determine the valid length of the slice */
		let ins = ins.as_slice();
		let len = ins.iter().peek_fold(0, |a, &x| {
			if x == 0 { return None; }
			return Some(a + 1);
		});
		return Some(&ins[.. len]);
	}

	fn get_ins_core(dst: &mut Vec<u8>, ins: &[u8]) -> usize {
		/* expand packed nucleotide to Vec */
		let range = dst.reserve_to(ins.len() * 2, |arr: &mut [u8], _: &[u8]| -> usize {
			for (i, x) in arr.iter_mut().enumerate() {
				/* little endian */
				*x = if (i & 0x01) == 0 {
					ins[i / 2] & 0x0f
				} else {
					ins[i / 2]>>4
				};
			}

			let remove_tail = arr[arr.len() - 1] == 0;
			// debug!("{}, {}, {}, {:?}", ins.len(), arr.len(), remove_tail, arr);

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

/*
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
*/


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
	use arraytools::ArrayTools;
	use std::ops::Range;
	use std::str::from_utf8;

	#[allow(unused_imports)]
	use crate::{ Udon, UdonPrecursor, CigarOp, UdonPalette, BLOCK_PITCH, encode_base_unchecked, decode_base_unchecked };

	macro_rules! cigar {
		[ $( ( $op: ident, $len: expr ) ),* ] => ({
			vec![ $( CigarOp::$op as u32 | (( $len as u32 )<<4) ),* ]
		});
	}

	macro_rules! nucl {
		( $st: expr ) => ({
			let s = $st;
			println!("{:?}", s);

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
			let mut v = Vec::new();
			for &x in $arr {
				let c = match x {
					0x00 => 'M', 0x04 => 'A', 0x05 => 'C', 0x06 => 'G', 0x07 => 'T', 0x08 => 'D',
					0x10 => 'M', 0x14 => 'A', 0x15 => 'C', 0x16 => 'G', 0x17 => 'T', 0x18 => 'D',
					_ => ' '
				};
				v.push(c as u8);
			}
			v
		});
	}

	macro_rules! encode_ins {
		( $arr: expr ) => ({
			let mut v = Vec::new();
			for &x in $arr {
				v.push(if (x & 0x10) == 0x10 { 'I' } else { '-' } as u8);
			}
			v
		});
	}

	#[allow(unused_macros)]
	macro_rules! decode_nucl {
		( $arr: expr ) => ({
			let mut v = Vec::new();
			for &x in $arr {
				v.push(decode_base_unchecked((x as u32)>>4) as u8);
				if (x & 0x0f) == 0 { continue; }		/* should be the last base */

				v.push(decode_base_unchecked((x as u32) & 0x0f) as u8);
			}
			v
		});
	}

	macro_rules! compare {
		( $cigar: expr, $nucl: expr, $mdstr: expr, $range: expr, $flat: expr, $ins: expr ) => ({
			// let v = Vec::<u8>::new();
			let c = $cigar;
			let n = $nucl;
			let m = $mdstr;
			let u = match Udon::build(&c, &n, &m.as_bytes()) {
				None => {
					assert!(false, "failed to build index");
					return;
				},
				Some(u) => u
			};
			let mut r: Range<usize> = $range;
			if r.start == 0 && r.end == 0 {
				r.end = u.ref_span();
			}

			let a = u.decode_raw(&r).unwrap();
			let f = encode_flat!(&a);
			let d = from_utf8(&f).unwrap();
			assert!(d == $flat, "{:?}, {:?}", d, $flat);

			let j = encode_ins!(&a);
			let i = from_utf8(&j).unwrap();
			assert!(i == $ins, "{:?}, {:?}", i, $ins);
		});
	}

	macro_rules! compare_ins {
		( $cigar: expr, $nucl: expr, $mdstr: expr, $pos: expr, $ins_seq: expr ) => ({
			let c = $cigar;
			let n = $nucl;
			let m = $mdstr;
			let u = match Udon::build(&c, &n, &m.as_bytes()) {
				None => {
					assert!(false, "failed to build index");
					return;
				},
				Some(u) => u
			};

			let i = match u.get_ins($pos) {
				None    => vec!['*' as u8],
				Some(v) => v.iter().map(|x| decode_base_unchecked(*x as u32) as u8).collect()
			};
			let i = from_utf8(&i).unwrap();
			assert!(i == $ins_seq, "{:?}, {:?}", i, $ins_seq);
		});
	}

	const BG: [u8; 4]   = [0x00, 0x00, 0x00, 0x00];
	const DEL: [u8; 4]  = [0x00, 0x00, 0xff, 0x00];
	const INS: [u8; 4]  = [0xff, 0x00, 0xff, 0x00];
	const MISA: [u8; 4] = [0x7f, 0x1f, 0x1f, 0x00];
	const MISC: [u8; 4] = [0x00, 0x00, 0xff, 0x00];
	const MISG: [u8; 4] = [0x00, 0xff, 0x00, 0x00];
	const MIST: [u8; 4] = [0xff, 0x00, 0x00, 0x00];

	macro_rules! compare_color {
		( $cigar: expr, $nucl: expr, $mdstr: expr, $range: expr, $offset: expr, $scale: expr, $ribbon: expr ) => ({
			let c = $cigar;
			let n = $nucl;
			let m = $mdstr;
			let u = match Udon::build(&c, &n, &m.as_bytes()) {
				None => {
					assert!(false, "failed to build index");
					return;
				},
				Some(u) => u
			};
			let mut r: Range<usize> = $range;
			if r.start == 0 && r.end == 0 {
				r.end = u.ref_span();
			}

			let c = UdonPalette {
				background: BG,
				del: DEL,
				ins: INS,
				mismatch: [ MISA, MISC, MISG, MIST ]
			};
			let b = match u.decode_scaled(&r, $offset, $scale, &c) {
				None    => Vec::<u32>::new(),
				Some(v) => v
			};
			let n: Vec::<u32> = $ribbon.iter().map(|x| u32::from_le_bytes(*x)).collect();
			assert!(b == n, "{:?}, {:?}", b, n);
		});
	}

	#[test]
	fn test_udon_build_match() {
		compare!(
			cigar![(Match, 4)],
			nucl!("ACGT"),
			"4",
			Range { start: 0, end: 0 },
			"MMMM",
			"----"
		);
		compare!(
			cigar![(Match, 30)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTAC"),
			"30",
			Range { start: 0, end: 0 },
			"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
			"------------------------------"
		);
		compare!(
			cigar![(Match, 31)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACG"),
			"31",
			Range { start: 0, end: 0 },
			"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
			"-------------------------------"
		);
		compare!(
			cigar![(Match, 32)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGT"),
			"32",
			Range { start: 0, end: 0 },
			"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
			"--------------------------------"
		);
		compare!(
			cigar![(Match, 128)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
			"128",
			Range { start: 0, end: 0 },
			"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
			"--------------------------------------------------------------------------------------------------------------------------------"
		);
	}

	#[test]
	fn test_udon_build_del() {
		compare!(
			cigar![(Match, 4), (Del, 1), (Match, 4)],
			nucl!("ACGTACGT"),
			"4^A4",
			Range { start: 0, end: 0 },
			"MMMMDMMMM",
			"---------"
		);
		compare!(
			cigar![(Match, 4), (Del, 3), (Match, 4)],
			nucl!("ACGTACGT"),
			"4^AGG4",
			Range { start: 0, end: 0 },
			"MMMMDDDMMMM",
			"-----------"
		);
		compare!(
			cigar![(Match, 4), (Del, 4), (Match, 4)],
			nucl!("ACGTACGT"),
			"4^AAAC4",
			Range { start: 0, end: 0 },
			"MMMMDDDDMMMM",
			"------------"
		);
		compare!(
			cigar![(Match, 4), (Del, 11), (Match, 4)],
			nucl!("ACGTACGT"),
			"4^GATAGATAGGG4",
			Range { start: 0, end: 0 },
			"MMMMDDDDDDDDDDDMMMM",
			"-------------------"
		);
	}

	#[test]
	fn test_udon_build_ins() {
		compare!(
			cigar![(Match, 4), (Ins, 1), (Match, 4)],
			nucl!("ACGTACGTA"),
			"8",
			Range { start: 0, end: 0 },
			"MMMMMMMM",
			"----I---"
		);
		compare!(
			cigar![(Match, 4), (Ins, 2), (Match, 4)],
			nucl!("ACGTACGTAC"),
			"8",
			Range { start: 0, end: 0 },
			"MMMMMMMM",
			"----I---"
		);
	}

	#[test]
	fn test_udon_build_mismatch() {
		compare!(
			cigar![(Match, 10)],
			nucl!("ACGTACGTAC"),
			"4T5",
			Range { start: 0, end: 0 },
			"MMMMAMMMMM",
			"----------"
		);
		compare!(
			cigar![(Match, 10)],
			nucl!("ACGTACGTAC"),
			"4T0C4",
			Range { start: 0, end: 0 },
			"MMMMACMMMM",
			"----------"
		);
		compare!(
			cigar![(Match, 10)],
			nucl!("ACGTACGTAC"),
			"4T1A0T2",
			Range { start: 0, end: 0 },
			"MMMMAMGTMM",
			"----------"
		);
	}

	#[test]
	fn test_udon_build_cont_mismatch() {
		/* continuous flag */
		compare!(
			cigar![(Match, 34)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"),
			"30T3",
			Range { start: 0, end: 0 },
			"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMGMMM",
			"----------------------------------"
		);
		compare!(
			cigar![(Match, 64)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
			"60T3",
			Range { start: 0, end: 0 },
			"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMAMMM",
			"----------------------------------------------------------------"
		);
	}

	#[test]
	fn test_udon_build_cont_del() {
		/* continuous flag */
		compare!(
			cigar![(Match, 30), (Del, 4), (Match, 4)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"),
			"30^ACGT4",
			Range { start: 0, end: 0 },
			"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMDDDDMMMM",
			"--------------------------------------"
		);
		compare!(
			cigar![(Match, 60), (Del, 4), (Match, 4)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
			"60^ACGT4",
			Range { start: 0, end: 0 },
			"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMDDDDMMMM",
			"--------------------------------------------------------------------"
		);
	}

	#[test]
	fn test_udon_build_cont_ins() {
		/* continuous flag */
		compare!(
			cigar![(Match, 30), (Ins, 4), (Match, 4)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"),
			"34",
			Range { start: 0, end: 0 },
			"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
			"------------------------------I---"
		);
		compare!(
			cigar![(Match, 60), (Ins, 4), (Match, 4)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
			"64",
			Range { start: 0, end: 0 },
			"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
			"------------------------------------------------------------I---"
		);
	}

	#[test]
	fn test_udon_build_softclip() {
		compare!(
			cigar![(SoftClip, 10), (Match, 10)],
			nucl!("ACGTACGTACGTACGTACGT"),
			"4T5",
			Range { start: 0, end: 0 },
			"MMMMGMMMMM",
			"----------"
		);
		compare!(
			cigar![(Match, 10), (SoftClip, 10)],
			nucl!("ACGTACGTACGTACGTACGT"),
			"4T5",
			Range { start: 0, end: 0 },
			"MMMMAMMMMM",
			"----------"
		);
		compare!(
			cigar![(SoftClip, 10), (Match, 10), (SoftClip, 10)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTAC"),
			"4T5",
			Range { start: 0, end: 0 },
			"MMMMGMMMMM",
			"----------"
		);
	}

	#[test]
	fn test_udon_build_hardclip() {
		compare!(
			cigar![(HardClip, 10), (Match, 10)],
			nucl!("GTACGTACGT"),
			"4T5",
			Range { start: 0, end: 0 },
			"MMMMGMMMMM",
			"----------"
		);
		compare!(
			cigar![(Match, 10), (HardClip, 10)],
			nucl!("ACGTACGTAC"),
			"4T5",
			Range { start: 0, end: 0 },
			"MMMMAMMMMM",
			"----------"
		);
		compare!(
			cigar![(HardClip, 10), (Match, 10), (HardClip, 10)],
			nucl!("GTACGTACGT"),
			"4T5",
			Range { start: 0, end: 0 },
			"MMMMGMMMMM",
			"----------"
		);
	}

	#[test]
	fn test_udon_build_head_del() {
		compare!(
			cigar![(Del, 4), (Match, 4)],
			nucl!("ACGT"),
			"^ACGT4",
			Range { start: 0, end: 0 },
			"DDDDMMMM",
			"--------"
		);
	}

	#[test]
	fn test_udon_build_head_ins() {
		compare!(
			cigar![(Ins, 4), (Match, 4)],
			nucl!("ACGTACGT"),
			"4",
			Range { start: 0, end: 0 },
			"MMMM",
			"I---"
		);
	}

	#[test]
	fn test_udon_build_tail_del() {
		compare!(
			cigar![(Match, 4), (Del, 4)],
			nucl!("ACGT"),
			"4^ACGT",
			Range { start: 0, end: 0 },
			"MMMMDDDD",
			"--------"
		);
	}

	#[test]
	fn test_udon_build_tail_ins() {
		compare!(
			cigar![(Match, 4), (Ins, 4)],
			nucl!("ACGTACGT"),
			"4",
			Range { start: 0, end: 0 },
			"MMMM",
			"----"
		);
	}

	#[test]
	fn test_udon_build_del_ins() {
		/* not natural as CIGAR string but sometimes appear in real data */
		compare!(
			cigar![(Match, 4), (Del, 4), (Ins, 4), (Match, 4)],
			nucl!("ACGTGGGGACGT"),
			"4^CCCC4",
			Range { start: 0, end: 0 },
			"MMMMDDDDMMMM",
			"--------I---"
		);
		compare!(
			cigar![(Match, 4), (Del, 4), (Ins, 4), (Del, 4), (Ins, 4), (Match, 4)],
			nucl!("ACGTGGGGAAAAACGT"),
			"4^CCCC0^TTTT4",		/* is this correct? */
			Range { start: 0, end: 0 },
			"MMMMDDDDDDDDMMMM",
			"------------I---"		/* FIXME: insertion marker lost */
		);
	}

	#[test]
	fn test_udon_build_ins_del() {
		/* also not natural */
		compare!(
			cigar![(Match, 4), (Ins, 4), (Del, 4), (Match, 4)],
			nucl!("ACGTGGGGACGT"),
			"4^CCCC4",
			Range { start: 0, end: 0 },
			"MMMMDDDDMMMM",
			"------------"			/* insertion marker lost */
		);
		compare!(
			cigar![(Match, 4), (Ins, 4), (Del, 4), (Ins, 4), (Del, 4), (Match, 4)],
			nucl!("ACGTGGGGACGT"),
			"4^CCCC0^AAAA4",
			Range { start: 0, end: 0 },
			"MMMMDDDDDDDDMMMM",
			"----------------"		/* again, lost */
		);
	}

	#[test]
	fn test_udon_build_complex() {
		compare!(
			cigar![(SoftClip, 7), (Match, 4), (Ins, 1), (Match, 4), (Del, 1), (Match, 2), (Del, 7), (Match, 40), (HardClip, 15)],
			nucl!("TTTTTTTACGTACGTACGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
			"8^A2^ACGTACG4T9A0C0G23",
			Range { start: 0, end: 0 },
			"MMMMMMMMDMMDDDDDDDMMMMAMMMMMMMMMGTAMMMMMMMMMMMMMMMMMMMMMMM",
			"----I-----------------------------------------------------"
		);
	}

	#[test]
	fn test_udon_decode_match() {
		compare!(
			cigar![(Match, 8)],
			nucl!("ACGTACGT"),
			"8",
			Range { start: 2, end: 6 },
			"MMMM",
			"----"
		);
	}

	#[test]
	fn test_udon_decode_mismatch() {
		compare!(
			cigar![(Match, 8)],
			nucl!("ACGTACGT"),
			"4T3",
			Range { start: 2, end: 6 },
			"MMAM",
			"----"
		);
	}

	#[test]
	fn test_udon_decode_del() {
		compare!(
			cigar![(Match, 4), (Del, 1), (Match, 4)],
			nucl!("ACGTACGT"),
			"4^T4",
			Range { start: 2, end: 7 },
			"MMDMM",
			"-----"
		);
	}

	#[test]
	fn test_udon_decode_ins() {
		compare!(
			cigar![(Match, 4), (Ins, 1), (Match, 4)],
			nucl!("ACGTACGTA"),
			"8",
			Range { start: 2, end: 6 },
			"MMMM",
			"--I-"
		);
	}

	#[test]
	fn test_udon_decode_cont_mismatch() {
		/* mismatch on boundary */
		compare!(
			cigar![(Match, 34)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"),
			"30T3",
			Range { start: 30, end: 34 },
			"GMMM",
			"----"
		);
		compare!(
			cigar![(Match, 64)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
			"60T3",
			Range { start: 60, end: 64 },
			"AMMM",
			"----"
		);

		/* skip one */
		compare!(
			cigar![(Match, 34)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"),
			"30T3",
			Range { start: 31, end: 34 },
			"MMM",
			"---"
		);
	}

	#[test]
	fn test_udon_decode_cont_del() {
		/* deletion on boundary */
		compare!(
			cigar![(Match, 30), (Del, 4), (Match, 4)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"),
			"30^ACGT4",
			Range { start: 30, end: 38 },
			"DDDDMMMM",
			"--------"
		);
		compare!(
			cigar![(Match, 60), (Del, 4), (Match, 4)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
			"60^ACGT4",
			Range { start: 60, end: 68 },
			"DDDDMMMM",
			"--------"
		);

		/* skipping one */
		compare!(
			cigar![(Match, 30), (Del, 4), (Match, 4)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"),
			"30^ACGT4",
			Range { start: 31, end: 38 },
			"DDDMMMM",
			"-------"
		);

		/* leaving one */
		compare!(
			cigar![(Match, 30), (Del, 4), (Match, 4)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"),
			"30^ACGT4",
			Range { start: 29, end: 33 },
			"MDDD",
			"----"
		);
	}

	#[test]
	fn test_udon_decode_cont_ins() {
		/* insertion on boundary */
		compare!(
			cigar![(Match, 30), (Ins, 4), (Match, 4)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"),
			"34",
			Range { start: 30, end: 34 },
			"MMMM",
			"I---"
		);
		compare!(
			cigar![(Match, 60), (Ins, 4), (Match, 4)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
			"64",
			Range { start: 60, end: 64 },
			"MMMM",
			"I---"
		);
		/* skip one */
		compare!(
			cigar![(Match, 30), (Ins, 4), (Match, 4)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"),
			"34",
			Range { start: 31, end: 34 },
			"MMM",
			"---"
		);
	}
	#[test]
	fn test_udon_decode_poll_match() {
		/* test block polling */
		compare!(
			cigar![(Match, BLOCK_PITCH + 8)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}", BLOCK_PITCH + 8),
			Range { start: BLOCK_PITCH - 2, end: BLOCK_PITCH + 6 },
			"MMMMMMMM",
			"--------"
		);
		compare!(
			cigar![(Match, BLOCK_PITCH + 8)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}", BLOCK_PITCH + 8),
			Range { start: BLOCK_PITCH + 2, end: BLOCK_PITCH + 6 },
			"MMMM",
			"----"
		);

		/* longer */
		compare!(
			cigar![(Match, 21 * BLOCK_PITCH + 8)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; 21 * BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}", 21 * BLOCK_PITCH + 8),
			Range { start: 21 * BLOCK_PITCH - 2, end: 21 * BLOCK_PITCH + 6 },
			"MMMMMMMM",
			"--------"
		);
		compare!(
			cigar![(Match, 21 * BLOCK_PITCH + 8)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; 21 * BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}", 21 * BLOCK_PITCH + 8),
			Range { start: 21 * BLOCK_PITCH + 2, end: 21 * BLOCK_PITCH + 6 },
			"MMMM",
			"----"
		);
	}

	#[test]
	fn test_udon_decode_poll_mismatch() {
		compare!(
			cigar![(Match, BLOCK_PITCH + 8)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}T3", BLOCK_PITCH + 4),
			Range { start: BLOCK_PITCH - 2, end: BLOCK_PITCH + 6 },
			"MMMMMMAM",
			"--------"
		);
		compare!(
			cigar![(Match, BLOCK_PITCH + 8)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}T3", BLOCK_PITCH + 4),
			Range { start: BLOCK_PITCH + 2, end: BLOCK_PITCH + 6 },
			"MMAM",
			"----"
		);

		/* mismatch on block boundary */
		compare!(
			cigar![(Match, BLOCK_PITCH + 8)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}T7", BLOCK_PITCH),
			Range { start: BLOCK_PITCH - 2, end: BLOCK_PITCH + 2 },
			"MMAM",
			"----"
		);
		/* mismatch right before block boundary */
		compare!(
			cigar![(Match, BLOCK_PITCH + 8)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}T8", BLOCK_PITCH - 1),
			Range { start: BLOCK_PITCH - 2, end: BLOCK_PITCH + 2 },
			"MAMM",
			"----"
		);
	}

	#[test]
	fn test_udon_decode_poll_mismatch_long() {
		/* much longer */
		compare!(
			cigar![(Match, 321 * BLOCK_PITCH + 8)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; 321 * BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}T3", 321 * BLOCK_PITCH + 4),
			Range { start: 321 * BLOCK_PITCH - 2, end: 321 * BLOCK_PITCH + 6 },
			"MMMMMMAM",
			"--------"
		);
		compare!(
			cigar![(Match, 321 * BLOCK_PITCH + 8)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; 321 * BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}T3", 321 * BLOCK_PITCH + 4),
			Range { start: 321 * BLOCK_PITCH + 2, end: 321 * BLOCK_PITCH + 6 },
			"MMAM",
			"----"
		);
		compare!(
			cigar![(Match, 321 * BLOCK_PITCH + 8)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; 321 * BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}T7", 321 * BLOCK_PITCH),
			Range { start: 321 * BLOCK_PITCH - 2, end: 321 * BLOCK_PITCH + 2 },
			"MMAM",
			"----"
		);
	}

	#[test]
	fn test_udon_decode_poll_del() {
		/* boundary on boundary */
		compare!(
			cigar![(Match, BLOCK_PITCH), (Del, 4), (Match, 4)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}^ACGT4", BLOCK_PITCH),
			Range { start: BLOCK_PITCH - 2, end: BLOCK_PITCH + 6 },
			"MMDDDDMM",
			"--------"
		);
		compare!(
			cigar![(Match, BLOCK_PITCH), (Del, 4), (Match, 4)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}^ACGT4", BLOCK_PITCH),
			Range { start: BLOCK_PITCH, end: BLOCK_PITCH + 6 },
			"DDDDMM",
			"------"
		);
	}

	#[test]
	fn test_udon_decode_poll_del2() {
		/* over boundary */
		compare!(
			cigar![(Match, BLOCK_PITCH - 2), (Del, 4), (Match, 6)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}^ACGT6", BLOCK_PITCH - 2),
			Range { start: BLOCK_PITCH - 2, end: BLOCK_PITCH + 6 },
			"DDDDMMMM",
			"--------"
		);
		compare!(
			cigar![(Match, BLOCK_PITCH - 2), (Del, 4), (Match, 6)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}^ACGT6", BLOCK_PITCH - 2),
			Range { start: BLOCK_PITCH, end: BLOCK_PITCH + 6 },
			"DDMMMM",
			"------"
		);
		compare!(
			cigar![(Match, BLOCK_PITCH - 2), (Del, 4), (Match, 6)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}^ACGT6", BLOCK_PITCH - 2),
			Range { start: BLOCK_PITCH + 2, end: BLOCK_PITCH + 6 },
			"MMMM",
			"----"
		);
	}

	#[test]
	fn test_udon_decode_poll_ins() {
		/* boundary on boundary */
		compare!(
			cigar![(Match, BLOCK_PITCH - 2), (Ins, 4), (Match, 6)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}", BLOCK_PITCH + 4),
			Range { start: BLOCK_PITCH - 2, end: BLOCK_PITCH + 2 },
			"MMMM",
			"I---"
		);
		compare!(
			cigar![(Match, BLOCK_PITCH - 2), (Ins, 4), (Match, 6)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGT")),
			format!("{}", BLOCK_PITCH + 4),
			Range { start: BLOCK_PITCH, end: BLOCK_PITCH + 4 },
			"MMMM",
			"----"
		);
	}

	#[test]
	fn test_udon_decode_query_ins() {
		compare_ins!(
			cigar![(Match, 4), (Ins, 4), (Match, 4)],
			nucl!("ACGTACGTACGT"),
			"8",
			4,
			"ACGT"
		);
		compare_ins!(
			cigar![(Match, 4), (Ins, 4), (Match, 4)],
			nucl!("ACGTACGTACGT"),
			"8",
			3,
			"*"
		);
		compare_ins!(
			cigar![(Match, 4), (Ins, 4), (Match, 4)],
			nucl!("ACGTACGTACGT"),
			"8",
			5,
			"*"
		);
	}

	#[test]
	fn test_udon_decode_query_ins_double() {
		compare_ins!(
			cigar![(Match, 4), (Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!("ACGTACGTACGTGGGGACGT"),
			"12",
			8,
			"GGGG"
		);
		compare_ins!(
			cigar![(Match, 4), (Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!("ACGTACGTACGTGGGGACGT"),
			"12",
			7,
			"*"
		);
		compare_ins!(
			cigar![(Match, 4), (Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!("ACGTACGTACGTGGGGACGT"),
			"12",
			9,
			"*"
		);
	}

	#[test]
	fn test_udon_decode_query_ins_head() {
		compare_ins!(
			cigar![(Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!("CCCCACGTACGTACGT"),
			"8",
			0,
			"CCCC"
		);
		compare_ins!(
			cigar![(Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!("CCCCACGTACGTACGT"),
			"8",
			1,
			"*"
		);
		compare_ins!(
			cigar![(Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!("CCCCACGTGGGGACGT"),
			"8",
			4,
			"GGGG"
		);
	}

	#[test]
	fn test_udon_decode_query_ins_poll() {
		compare_ins!(
			cigar![(Match, BLOCK_PITCH + 4), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGTACGT")),
			format!("{}", BLOCK_PITCH + 8),
			BLOCK_PITCH + 4,
			"ACGT"
		);
		compare_ins!(
			cigar![(Match, BLOCK_PITCH + 4), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGTACGT")),
			format!("{}", BLOCK_PITCH + 8),
			BLOCK_PITCH + 3,
			"*"
		);
		compare_ins!(
			cigar![(Match, BLOCK_PITCH + 4), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGTACGT")),
			format!("{}", BLOCK_PITCH + 8),
			BLOCK_PITCH + 5,
			"*"
		);
	}

	#[test]
	fn test_udon_decode_query_ins_double_poll() {
		compare_ins!(
			cigar![(Match, BLOCK_PITCH + 4), (Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGTACGTGGGGACGT")),
			format!("{}", BLOCK_PITCH + 12),
			BLOCK_PITCH + 8,
			"GGGG"
		);
		compare_ins!(
			cigar![(Match, BLOCK_PITCH + 4), (Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGTACGTGGGGACGT")),
			format!("{}", BLOCK_PITCH + 12),
			BLOCK_PITCH + 7,
			"*"
		);
		compare_ins!(
			cigar![(Match, BLOCK_PITCH + 4), (Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}", from_utf8(&['A' as u8; BLOCK_PITCH]).unwrap(), "ACGTACGTACGTGGGGACGT")),
			format!("{}", BLOCK_PITCH + 12),
			BLOCK_PITCH + 9,
			"*"
		);
	}

	#[test]
	fn test_udon_decode_query_ins_double_poll2() {
		compare_ins!(
			cigar![(Match, 4), (Ins, 4), (Match, BLOCK_PITCH - 8), (Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}{}", "ACGTCCCC", from_utf8(&['A' as u8; BLOCK_PITCH - 8]).unwrap(), "TTTTACGTGGGGACGT")),
			format!("{}", BLOCK_PITCH + 4),
			4,
			"CCCC"
		);
		compare_ins!(
			cigar![(Match, 4), (Ins, 4), (Match, BLOCK_PITCH - 8), (Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}{}", "ACGTCCCC", from_utf8(&['A' as u8; BLOCK_PITCH - 8]).unwrap(), "TTTTACGTGGGGACGT")),
			format!("{}", BLOCK_PITCH + 4),
			BLOCK_PITCH - 4,
			"TTTT"
		);
		compare_ins!(
			cigar![(Match, 4), (Ins, 4), (Match, BLOCK_PITCH - 8), (Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}{}", "ACGTCCCC", from_utf8(&['A' as u8; BLOCK_PITCH - 8]).unwrap(), "TTTTACGTGGGGACGT")),
			format!("{}", BLOCK_PITCH + 4),
			BLOCK_PITCH,
			"GGGG"
		);
	}

	#[test]
	fn test_udon_decode_query_ins_head_double_poll() {
		compare_ins!(
			cigar![(Ins, 4), (Match, 4), (Ins, 4), (Match, BLOCK_PITCH - 8), (Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}{}", "GGGGACGTCCCC", from_utf8(&['A' as u8; BLOCK_PITCH - 8]).unwrap(), "TTTTACGTGGGGACGT")),
			format!("{}", BLOCK_PITCH + 4),
			0,
			"GGGG"
		);
		compare_ins!(
			cigar![(Ins, 4), (Match, 4), (Ins, 4), (Match, BLOCK_PITCH - 8), (Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}{}", "GGGGACGTCCCC", from_utf8(&['A' as u8; BLOCK_PITCH - 8]).unwrap(), "TTTTACGTGGGGACGT")),
			format!("{}", BLOCK_PITCH + 4),
			4,
			"CCCC"
		);
		compare_ins!(
			cigar![(Ins, 4), (Match, 4), (Ins, 4), (Match, BLOCK_PITCH - 8), (Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}{}", "GGGGACGTCCCC", from_utf8(&['A' as u8; BLOCK_PITCH - 8]).unwrap(), "TTTTACGTGGGGACGT")),
			format!("{}", BLOCK_PITCH + 4),
			BLOCK_PITCH - 4,
			"TTTT"
		);
		compare_ins!(
			cigar![(Ins, 4), (Match, 4), (Ins, 4), (Match, BLOCK_PITCH - 8), (Ins, 4), (Match, 4), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}{}", "GGGGACGTCCCC", from_utf8(&['A' as u8; BLOCK_PITCH - 8]).unwrap(), "TTTTACGTGGGGACGT")),
			format!("{}", BLOCK_PITCH + 4),
			BLOCK_PITCH,
			"GGGG"
		);
	}

	#[test]
	fn test_udon_decode_query_ins_head_double_poll2() {
		compare_ins!(
			cigar![(Ins, 4), (Match, 4), (Ins, 4), (Match, BLOCK_PITCH - 8), (Ins, 4), (Match, 8), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}{}", "GGGGACGTCCCC", from_utf8(&['A' as u8; BLOCK_PITCH - 8]).unwrap(), "TTTTACGTGGGGACGTACGT")),
			format!("{}", BLOCK_PITCH + 8),
			BLOCK_PITCH + 4,
			"ACGT"
		);
		compare_ins!(
			cigar![(Ins, 4), (Match, 4), (Ins, 4), (Match, BLOCK_PITCH - 8), (Ins, 4), (Match, 2), (Del, 4), (Match, 2), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}{}", "GGGGACGTCCCC", from_utf8(&['A' as u8; BLOCK_PITCH - 8]).unwrap(), "TTTTACGTGGGGACGT")),
			format!("{}^ACGT8", BLOCK_PITCH - 2),
			BLOCK_PITCH + 4,
			"GGGG"
		);
		compare_ins!(
			cigar![(Ins, 4), (Match, 4), (Ins, 4), (Match, BLOCK_PITCH - 8), (Ins, 4), (Del, 8), (Ins, 4), (Match, 4)],
			nucl!(format!("{}{}{}", "GGGGACGTCCCC", from_utf8(&['A' as u8; BLOCK_PITCH - 8]).unwrap(), "TTTTGGGGACGT")),
			format!("{}^ACGTACGT4", BLOCK_PITCH - 4),
			BLOCK_PITCH + 4,
			"GGGG"
		);
	}

	#[test]
	fn test_udon_decode_scaled() {
		compare_color!(
			cigar![(Match, 4), (Del, 1), (Match, 4)],
			nucl!("ACGTACGT"),
			"4^A4",
			Range { start: 0, end: 0 },
			0.0, 1.0,
			vec![BG, BG, BG, BG, DEL, BG, BG, BG, BG]
		);

		compare_color!(
			cigar![(Match, 4), (Del, 1), (Match, 4)],
			nucl!("ACGTACGT"),
			"4^A4",
			Range { start: 0, end: 0 },
			0.0, 3.0,
			vec![BG, DEL.map(|x| x / 3), BG]
		);

		/* we need more tests but how to do */
		/*
		compare_color!(
			cigar![(Match, 9)],
			nucl!("ACGTACGTA"),
			"G0T0A0C0G0T0A0C0G",
			Range { start: 0, end: 0 },
			1.0 / 3.0, 1.5,
			vec![BG, DEL.map(|x| x / 3), BG]
		);
		*/
	}
}





