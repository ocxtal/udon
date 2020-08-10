
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

// extern crate packed_struct;
// #[macro_use]
// extern crate packed_struct_codegen;
// use packed_struct::prelude::*;

#[macro_use]
extern crate bitfield;

use std::io::Write;
use std::mem::{ size_of, transmute };
use std::ops::Range;
use std::ptr::copy_nonoverlapping;
use std::slice::{ Iter, from_raw_parts, from_raw_parts_mut };


#[cfg(all(target_arch = "x86_64"))]
use core::arch::x86_64::*;

#[cfg(all(target_arch = "aarch64"))]
use core::arch::aarch64::*;				/* requires nightly */




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

/** bam CIGAR operation code */
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

/** bam CIGAR element definition */
bitfield!{
	#[derive(Copy, Clone, Default)]
	pub struct Cigar(u32);
	pub op, _:  4,  0;
	pub len, _: 32, 4;
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
	#[derive(Copy, Clone, Default)]
	pub struct Block(u64);
	pub ins_offset, set_ins_offset: 29, 0;
	pub op_offset,  set_op_offset:  59, 29;
	pub op_skip,    set_op_skip:    64, 59;
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
#[derive(Copy, Clone, Default)]
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
#[derive(Copy, Clone, Default)]
pub struct UdonPrecursor<'a>(Udon<'a>);


/* Precursor

pointer arithmetic for creating precursor, and converting precursor to product
*/
trait Precursor<'a, T> {
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
		let base_offset = (self.len() + size_of::<U>() - 1) / size_of::<U>();
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
		assert!(consumed_count < count);
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

Passes an element fetched from an slice iterator to a closure, and expect
to get a bool variable. The iteration continues when it's true and breaks
otherwise. On breaking the iteration, the last element remains (not fetched)
in the iterator.
*/
trait PeekFold<T: Sized> {
	fn peek_fold<A, F>(&mut self, init: A, func: F) -> A
		where Self: Sized,
		      A: Copy,
		      F: FnMut(A, &T) -> (A, bool);
}

impl<'a, T: Sized> PeekFold<T> for Iter<'a, T> {
	fn peek_fold<A, F>(&mut self, init: A, mut func: F) -> A
		where Self: Sized,
		      A: Copy,
		      F: FnMut(A, &T) -> (A, bool)
	{
		let mut accum: A = init;
		loop {
			let x = match self.as_slice().first() {
				None => { return accum; },
				Some(x) => x
			};

			let (next_accum, cont) = func(accum, x);
			if !cont { return accum; }

			self.next();
			accum = next_accum;
		}
	}
}


/* convert ascii-printed number to u64, overflow is ignored. */
fn atoi_unchecked(v: &mut Iter<u8>) -> u64 {
	v.peek_fold(0, |a, x| {
		let m = (*x as u64) - ('0' as u64);
		let a = 10 * a + m;
		return (a, m < 10);
	})
}

fn isnum(c: u8) -> bool {
	return (c as u64 - '0' as u64) < 10 as u64;
}

/*
transcode { 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't' } -> { 0x0, 0x1, 0x2, 0x3, 0x0, 0x1, 0x2, 0x3 }
and then mark mismatch.
*/
fn encode_base_unchecked(c: u8) -> u32 {
	let c = c as u32;
	let b2 = ((c>>1) - (c>>3)) & 0x03;
	return CompressMark::Mismatch as u32 + b2;
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
fn copy_packed_nucl(src: &[u8], dst: &mut [u8], skip_first: bool) {

	if skip_first == false {
		dst.copy_from_slice(src);
		return;
	}

	/* shift needed */
	for (s, d) in src.windows(2).zip(dst) {
		let mut v: [u8; 2] = [0; 2];
		v.copy_from_slice(s);

		*d = (u16::from_le_bytes(v)>>4) as u8;
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

		let c = self.query[ofs / 2];
		if (ofs & 0x01) == 0x01 {
			return encode_base_unchecked(c & 0x0f);
		}
		return encode_base_unchecked(c>>4);
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
		return !isnum(self.mdstr.as_slice()[1]);
	}

	/* writers, on top of impl Writer for Vec<u8> */
	fn push_op(&mut self, match_len: usize, marker: u32) {
		assert!(marker    < 0x08);
		assert!(match_len < 0x20);

		/*
		 * op == 0x00 indicates 
		 */
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
	fn push_ins(&mut self, qofs: usize, qlen: usize) {

		if qlen == 0 { return; }		/* must be corrupted cigar, but safe to ignore */

		/* copy subsequence; 4-bit packed */
		let copy_len = (qlen + 1) / 2;
		let packet_size = copy_len + 1;

		let query = self.query;
		self.ins.reserve_to(packet_size, |arr: &mut [u8], _: &[u8]| -> usize {
			/* marker: valid bases never generate zero byte thus can be used for boundary marker */
			let (head, mut payload) = match arr.split_first_mut() {
				None => { return 0; },
				Some((head, payload)) => (head, payload)
			};

			/* header; just store zero */
			*head = 0;

			/* payload */
			let skip_head = (qofs & 0x01) == 0x01;
			let len = copy_len + skip_head as usize;
			copy_packed_nucl(&query[qofs .. qofs + len], &mut payload, skip_head);

			packet_size
		});
	}
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

		self.push_ins(ofs, len);		/* use offset before forwarding */
		return Some(xrem - len);		/* there might be remainder */
	}
	fn eat_match(&mut self, xrem: usize, last_op: u32) -> Option<usize> {
		let c = self.cigar.next()?;

		let is_valid = c.op() != CigarOp::Ins as u32;
		if is_valid { self.cigar.next()?; }

		let mut crem = if is_valid { c.len() as usize } else { 0 } + xrem;
		let mut xrem = xrem;			/* might continues from the previous cigar op, possibly insertion */
		let mut op = last_op;			/* insertion, deletion, or mismatch */

		while xrem < crem {
			/* xrem < crem indicates this op (cigar span) is interrupted by mismatch(es) at the middle */
			self.push_match(xrem, op);
			crem      -= xrem;
			self.qofs += xrem;

			while self.is_double_mismatch() {
				let c = self.next_base();
				self.push_op(1, c);		/* this chunk contains only a single mismatch */
				crem -= 1;
			}
			op = self.next_base();		/* we only have a single mismatch remaining, will be combined to succeeding matches */

			/* next match length */
			xrem = self.eat_md_eq();	/* xrem for the next { match, insertion } region */
		}

		self.push_match(crem, op);		/* tail match; length is the remainder of crem */
		xrem      -= crem;
		self.qofs += crem;
		return Some(xrem);				/* nonzero if insertion follows */
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
		if c.op() == CigarOp::Ins as u32 {
			self.push_ins(0, 0);		/* dummy */
		}

		/* normal CIGAR starts with match */
		let c = self.peek_cigar()?;
		if c.op() == CigarOp::Match as u32 {
			xrem = self.eat_md_eq();

			/*
			the first op consists of dummy insertion and match
			(the insertion is real one for a CIGAR that starts with insertion. see above.)
			*/
			xrem = self.eat_match(xrem, CompressMark::Ins as u32)?;
		}

		loop {
			/* deletion-match pair */
			let c = self.peek_cigar()?;
			while c.op() == CigarOp::Del as u32 {

				/* the CIGAR ends with deletion; must be treated specially */
				if self.cigar.as_slice().len() < 2 { break; }

				/* push deletion-match pair, then parse next eq length */
				let op = self.eat_del()?;
				xrem = self.eat_md_eq() + op as usize;
				xrem = self.eat_match(xrem, op)?;
			}

			/* it's insertion-match pair when it appeared not be deletion-match */
			let c = self.peek_cigar()?;
			if c.op() != CigarOp::Ins as u32 {
				return None;
			}

			/* the CIGAR ends with insertion; must be treated specially */
			if self.cigar.as_slice().len() < 2 { break; }

			/* push insertion-match pair, update eq length remainder */
			xrem = self.eat_ins(xrem)?;
			xrem = self.eat_match(xrem, CompressMark::Ins as u32)?;
		}

		/* CIGAR ends with isolated insertion or deletion */
		if self.cigar.as_slice().len() > 0 {
			let c = self.peek_cigar()?;

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

		ops.iter().fold(0, |a, &x| {
			let len = x as usize & 0x1f;
			a + len - (if len == 31 { 1 } else { 0 })
		})
	}

	/* construct index for blocks */
	fn build_index(&mut self, op_range: &Range<usize>, ref_span: usize) -> Option<Range<usize>> {

		/* block pitch must be 2^n */
		assert!(BLOCK_PITCH.next_power_of_two() == BLOCK_PITCH);

		/* rip buffer for disjoint ownership */
		let buf = &mut self.buf;

		/* FIXME: we need to cover entire ref_span */
		let block_count = (ref_span + BLOCK_PITCH - 1) / BLOCK_PITCH;
		let range = buf.reserve_to(block_count, |block: &mut [Block], base: &[u8]| {
			let ops = &base[op_range.start..op_range.end];

			let mut rpos: usize = 0;
			let mut iofs: usize = 0;
			let mut cofs: usize = 0;
			for (i, &op) in ops.iter().enumerate() {
				let next_boundary = (rpos | (BLOCK_PITCH - 1)) + 1;

				/* forward insertion array and reference-side offset */
				let marker = op>>5;
				let len = op as usize & 0x1f;

				iofs += if marker == 0 { 1 } else { 0 };
				rpos += if len == 31 { len - 1 } else { len };

				/* if forwarded rpos doexn't exceed the next boundary, just skip this op */
				if rpos < next_boundary { continue; }

				/* boundary found; save block info */
				block[cofs].set_ins_offset(iofs as u64);
				block[cofs].set_op_offset(i as u64);
				block[cofs].set_op_skip(rpos as u64 & 0x1f);
				cofs += 1;
			}

			assert!(cofs == block_count);
			cofs
		});

		self.rspan = ref_span;
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
	fn finalize(&self, op: Range<usize>, block: Range<usize>, ins: Range<usize>) -> UdonPrecursor<'b> {
		assert!(op.end <= block.start);
		assert!(block.end <= ins.start);

		assert!(((block.end - block.start) % size_of::<Block>()) == 0);

		/* entire range on buffer for this object */
		let range = Range::<usize> {
			start: op.start,
			end:   ins.end
		};

		/* convert range to slice */
		UdonPrecursor::<'b> (
			Udon::<'b> {
				/* just save */
				size:     range.end - range.start,
				ref_span: self.rspan,

				/* compose fake slices (dereference causes SEGV) */
				op:    Precursor::<'b, u8>::compose(op),
				block: Precursor::<'b, Block>::compose(block),
				ins:   Precursor::<'b, u8>::compose(ins)
			}
		)
	}
}

/* build UdonPrecursor using UdonBuilder internally */
impl<'a, 'b> UdonPrecursor<'a> {
	pub fn build(buf: Vec<u8>, cigar: &'a [Cigar], packed_query: &'a [u8], mdstr: &'a [u8]) -> (Vec<u8>, Option<UdonPrecursor<'b>>)
	{
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
			qofs:  qclip.head,		/* initial offset for query sequence, non-zero for soft clipped alignments */
			rspan: 0				/* calcd in build_index */
		};

		/* if error detected, unwind destination vector and return it (so that the vector won't lost) */
		macro_rules! unwrap_or_unwind { ( $expr:expr ) => {
			match $expr {
				/* this block comes tail error handling with goto in C */
				None => {
					let mut buf = state.buf;
					buf.resize(base_offset, 0);		/* unwind buffer for cleaning up errored sequence */
					return (buf, None);
				},
				Some(val) => val
			}
		}}

		/* parse */
		let op = unwrap_or_unwind!(state.parse_cigar());
		let ref_span = state.calc_reference_span(&op);
		let block = unwrap_or_unwind!(state.build_index(&op, ref_span));
		let ins = unwrap_or_unwind!(state.pack_ins());

		/* everything done; compose precursor */
		let precursor = state.finalize(op, block, ins);

		/* return back ownership of the buffer */
		let buf = state.buf;
		return (buf, Some(precursor));
	}
}

impl<'a, 'b> Udon<'a> {

	/*
	 * builders
	 */
	fn from_precursor(buf: &Vec<u8>, precursor: UdonPrecursor<'a>) -> Udon<'a> {
		let base: *const u8 = buf.as_ptr();
		Udon::<'a> {
			/* just copy */
			size:     precursor.0.size,
			ref_span: precursor.0.ref_span,

			/* adjusting pointer */
			op:    precursor.0.op.finalize(base),
			block: precursor.0.block.finalize(base),
			ins:   precursor.0.ins.finalize(base)
		}
	}

	fn compose_box(mut buf: Vec<u8>, precursor: UdonPrecursor<'a>) -> Box<Udon<'a>> {

		/* compose pointer-adjusted header on stack */
		let header = Self::from_precursor(&buf, precursor);

		/* compose box and copy header into it */
		let udon = unsafe {
			/* convert buffer (on heap) to box */
			let mut udon = Box::<Udon<'a>>::from_raw(
				transmute::<*mut u8, *mut Udon<'a>>(buf.as_mut_ptr())
			);

			/* copy header from stack to heap (box) */
			copy_nonoverlapping(
				&header as *const Udon<'a>,
				&mut udon as &mut Udon<'a> as *mut Udon<'a>,
				1
			);
			udon
		};
		udon
	}

	/*
	pub fn from_precursor_vec(buf: &Vec<u8>, precursors: Vec<UdonPrecursor<'a>>) -> Vec<Udon<'a>> {
		unimplemented!();
	}
	*/

	pub fn build(cigar: &'b [Cigar], packed_query: &'b [u8], mdstr: &'b [u8]) -> Option<Box<Udon<'a>>>
	{
		let mut buf = Vec::<u8>::new();

		/* header always at the head */
		let range = buf.reserve_to(1, |header: &mut [Udon], _: &[u8]| {
			header[0] = Default::default();
			1
		});
		assert!(range.start == 0);
		assert!(range.end == size_of::<Udon>());

		return match UdonPrecursor::build(buf, cigar, packed_query, mdstr) {
			(_, None) => None,
			(buf, Some(precursor)) => Some(Self::compose_box(buf, precursor))
		}
	}


	/*
	 * decoder implementations
	 */
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
	pub unsafe fn decode_core_block(dst: &mut [u8], op: u32, ins: u64) -> (usize, u64) {

		/* load constants; expelled out of the innermost loop when inlined */
		let del_mask      = _mm_load_si128(&Self::DEL_MASK.v as *const [u8; 16] as *const __m128i);
		let scatter_mask  = _mm_load_si128(&Self::SCATTER_MASK.v as *const [u8; 16] as *const __m128i);
		let is_del_thresh = _mm_load_si128(&Self::IS_DEL_THRESH.v as *const [u8; 16] as *const __m128i);

		/* load op onto ymm */
		let rop = op as u64;

		/* compute deletion mask */
		let xop = _mm_cvtsi64_si128(rop as i64);
		let xop_scattered = _mm_shuffle_epi8(xop, scatter_mask);
		let is_del = _mm_cmpgt_epi8(xop_scattered, is_del_thresh);		/* signed comparison */

		/* compute mismatch / insertion mask */
		let yop = if rop == 0 { ins } else { rop>>5 };
		let ins_mismatch = _mm_cvtsi64_si128(yop as i64);

		/* merge deletion / insertion-mismatch vector */
		let merged = _mm_blendv_epi8(ins_mismatch, del_mask, is_del);
		_mm256_storeu_si256(dst as *mut [u8] as *mut __m256i, _mm256_castsi128_si256(merged));

		/*
		 * compute forward length; 31 is "continuous marker"
		 * critical path length is 6; rop -> rop
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

	fn check_span(&self, range: &Range<usize>) -> Option<bool> {
		let span = range.end - range.start;

		if span > self.ref_span { return None; }
		return Some(true);
	}

	fn scan_block(&self, pos: usize) -> (&[u8], usize) {

		let block_index = pos / BLOCK_PITCH;
		let block_rem   = pos & (BLOCK_PITCH - 1);

		let block = &self.block[block_index];
		let op_offset = block.op_offset() as usize;
		let op_skip   = block.op_skip() as usize;

		let ops = &self.op[op_offset ..];
		return (ops, block_rem + op_skip);
	}

	fn scan_op_array(ops: &[u8], rem: usize) -> (&[u8], usize) {

		let mut ops = ops.iter();
		let ofs = (&mut ops).peek_fold(0, |a, &x| {
			let len = x as usize & 0x1f;
			let len = a + len - (if len == 31 { 1 } else { 0 });

			return (len, len < rem);
		});
		return (ops.as_slice(), rem - ofs);		/* is this sound? */
	}

	pub fn decode_into(&self, dst: &mut Vec<u8>, ref_range: Range<usize>) -> Option<usize> {

		self.check_span(&ref_range)?;

		let (ops, rem) = self.scan_block(ref_range.start);
		let (ops, rem) = Self::scan_op_array(ops, rem);

		let used = Self::decode_core(dst, ops, rem, ref_range.end - ref_range.start)?;
		Some(used)
	}

	pub fn decode(&self, ref_range: Range<usize>) -> Option<Vec<u8>> {

		let size = ref_range.end - ref_range.start;
		let mut buf = Vec::<u8>::with_capacity(size);

		let used = self.decode_into(&mut buf, ref_range)?;
		buf.resize(used, 0);

		return Some(buf);
	}

	pub fn decode_scaled_into(&self, dst: &mut Vec<u8>, ref_range: Range<usize>, scale: f64) -> Option<usize> {
		None
	}

	pub fn decode_scaled(&self, ref_range: Range<usize>, scale: f64) -> Option<Vec<u8>> {
		None
	}

	/* ins */
	fn scan_ins(&self, ops: &[u8], rem: usize) -> (&[u8], usize) {
		unimplemented!();
	}

	pub fn get_ins_into(&self, dst: &mut Vec<u8>, pos: usize) -> Option<usize> {
		None
	}

	pub fn get_ins(&self, pos: usize) -> Option<Vec<u8>> {

		self.check_span(&Range::<usize> { start: 0, end: pos })?;

		let (ops, rem) = self.scan_block(pos);
		let (ops, rem) = Self::scan_op_array(ops, rem);

		// let mut buf = Vec::<u8>::with_capacity(10);
		// self.get_ins_into()
		None
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

	/*
	fn remove_clips(&mut self) -> Option<bool> {
		let c = self.peek_cigar()?;

		if c.op() == CigarOp::SoftClip as u32 {
			self.qofs += c.len() as usize;
			self.cigar.next();
			return Some(true);
		} else if c.op() == CigarOp::HardClip as u32 {
			self.cigar.next();
			return Some(true);
		}
		return Some(false);
	}

	fn finalize(self) -> Option<Vec<u8>> {
		return Some(self.buf);
	}
	*/



#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
