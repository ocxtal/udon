

use super::Udon;
use super::utils::{ SlicePrecursor, Writer, atoi_unchecked, transcode_base_unchecked };
use super::op::{ CigarOp, Cigar, CompressMark, IntoOpsIterator, op_len };
use super::index::{ Index, Block, BLOCK_PITCH };
use super::insbin::InsTracker;
use std::io::Write;
use std::mem::{ size_of, transmute, forget };
use std::ops::Range;
use std::ptr::copy_nonoverlapping;
use std::slice::{ Iter, IterMut };


/* builder APIs */
impl<'a, 'b> Index<'a> {

	pub(super) fn build(cigar: &'b [u32], packed_query: &'b [u8], mdstr: &'b [u8]) -> Option<Box<Udon<'a>>> {
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

	fn build_core(cigar: &'b [Cigar], packed_query: &'b [u8], mdstr: &'b [u8]) -> Option<Box<Udon<'a>>> {
		let mut buf = Vec::<u8>::new();

		/* header always at the head */
		let range = buf.reserve_to(1, |header: &mut [Udon], _: &[u8]| {
			header[0] = Default::default();
			1
		});
		assert!(range.start == 0);
		assert!(range.end == size_of::<Udon>());

		return match Precursor::build_core(buf, cigar, packed_query, mdstr) {
			(_, None) => None,
			(buf, Some(precursor)) => Some(Self::compose_box(buf, precursor))
		}
	}

	fn compose_box(buf: Vec<u8>, precursor: Precursor) -> Box<Udon<'a>> {
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
			let src = &header as *const Index<'a> as *const Udon<'a>;
			let dst = &mut udon as &mut Udon<'a> as *mut Udon<'a>;
			copy_nonoverlapping(src, dst, 1);

			udon
		};

		/* heap block inside buf was moved to udon so we have to release buf */
		forget(buf);		/* note: this allowed outside unsafe */

		udon
	}

	#[allow(dead_code)]
	pub(super) unsafe fn from_precursor(buf: &Vec<u8>, precursor: Precursor) -> Udon<'a> {
		let base: *const u8 = buf.as_ptr();

		Udon(Self::from_precursor_raw(base, &precursor))
	}

	unsafe fn from_precursor_raw(base: *const u8, precursor: &Precursor) -> Index<'a> {
		Index::<'a> {
			/* just copy */
			size:     precursor.size,
			ref_span: precursor.ref_span,

			/* adjusting pointer */
			op:    precursor.op.finalize_raw(base),
			block: precursor.block.finalize_raw(base),
			ins:   precursor.ins.finalize_raw(base)
		}
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


/* Precursor

Precursor is type that has some invalid pointers (slices) but the others are all sane.
The invalid pointer has an offset from a certain memory block, which may be reallocated
along with the object construction. The offset is converted to valid pointer, by adding
base pointer to the offset, at the very last of the object construction.
*/
#[repr(C)]
#[derive(Copy, Clone, Debug, Default)]
pub(super) struct Precursor {
	size: usize,
	ref_span: usize,
	op: SlicePrecursor<u8>,
	block: SlicePrecursor<Block>,
	ins: SlicePrecursor<u8>
}


/* build Precursor using Builder internally */
impl<'a> Precursor {
	fn build_core(buf: Vec<u8>, cigar: &'a [Cigar], packed_query: &'a [u8], mdstr: &'a [u8]) -> (Vec<u8>, Option<Precursor>) {
		/* save initial offset for unwinding */
		let base_offset = buf.len();

		/* compose working variables */
		let (cigar, qclip) = match strip_clips(cigar) {
			None => { return (buf, None); },		/* just return buffer (nothing happened) */
			Some((cigar, qclip)) => (cigar, qclip)
		};
		let mut state = Builder::<'a> {
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

	#[allow(dead_code)]
	pub unsafe fn build(buf: Vec<u8>, cigar: &'a [u32], packed_query: &'a [u8], mdstr: &'a [u8]) -> (Vec<u8>, Option<Precursor>) {
		
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
		let (buf, precursor) = unsafe { Precursor::build(self.buf, cigar, packed_query, mdstr) };

		/* put back */
		self.buf = buf;
		let precursor = precursor?;
		let ref_span = precursor.ref_span;
		self.precursors.push(precursor);

		return Some(ref_span);
	}

	unsafe fn from_precursor_vec(buf: Vec<u8>, precursors: Vec<Precursor>) -> Vec<Udon<'a>> {

		let base: *const u8 = buf.as_ptr();
		for mut precursor in precursors.iter_mut() {		/* consume? */

			/* compose udon header on stack */
			let header = Udon::<'a>::from_precursor_raw(base, &precursor);

			/* copy back udon on stack to vec */
			let src = &header as *const Udon<'a>;
			let dst = &mut precursor as &mut Precursor as *mut Precursor;

			let dst = transmute::<*mut Precursor, *mut Udon<'a>>(dst);
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
		let count = precursors.len() / size_of::<Precursor>();

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


/* transcoder state object */
pub struct Builder<'a> {
	buf: Vec<u8>,
	ins: Vec<u8>,
	cigar: Iter<'a, Cigar>,
	mdstr: Iter<'a, u8>,
	query: &'a [u8],
	qofs: usize
}

impl<'a, 'b> Builder<'a> {

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

	/* convert internal object (Builder) to Precursor */
	fn finalize(&self, ref_span: usize, op: &Range<usize>, ins: &Range<usize>, block: &Range<usize>) -> Precursor {
		assert!(op.end  <= ins.start);
		assert!(ins.end <= block.start);

		assert!(((block.end - block.start) % size_of::<Block>()) == 0);

		/* entire range on buffer for this object */
		let range = Range::<usize> {
			start: op.start,
			end:   ins.end
		};

		/* convert range to slice */
		Precursor {
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

