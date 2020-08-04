

// extern crate packed_struct;
// #[macro_use]
// extern crate packed_struct_codegen;
// use packed_struct::prelude::*;

#[macro_use]
extern crate bitfield;

use std::mem::size_of;
use std::slice::{ Iter, from_raw_parts, from_raw_parts_mut };
use std::ops::Range;


/*
 * bam cigar definition:
 * struct { uint32_t op : 4; uint32_t len : 28; }; as described in SAM/BAM spec pdf.
 */
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

bitfield!{
	#[derive(Copy, Clone, Default)]
	pub struct Cigar(u32);
	pub op, _:  4,  0;
	pub len, _: 32, 4;
}


/*
 * chunk for transcoded cigar.
 */
const CHUNK_PITCH: usize = 256;

bitfield!{
	#[derive(Copy, Clone, Default)]
	pub struct Chunk(u64);
	pub ins_offset, set_ins_offset: 29, 0;
	pub op_offset,  set_op_offset:  59, 29;
	pub op_skip,    set_op_skip:    64, 59;
}

/*
 * op definitions
 */
#[repr(u32)]
enum CompressMark {
	Ins      = 0x00,
	Mismatch = 0x04,
}

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

/*
#[derive(Copy, Clone, Default)]
struct MappingSpanPair {
	reference: u32,
	query: u32
}
*/

#[derive(Copy, Clone, Default)]
struct QueryClip {
	head: usize,
	tail: usize
}


/*
 * index object:
 * everything is laid out in a single flat memory block, thus constructed via some unsafe operations.
 */
#[derive(Copy, Clone, Default)]
pub struct Udon<'a> {
	size: usize,				/* object size for copying; for Clone */
	ref_span: usize,			/* reference side span */
	op: &'a [u8],
	chunk: &'a [Chunk],
	ins: &'a [u8]
}

#[repr(C)]
#[derive(Copy, Clone, Default)]
pub struct UdonPrecursor<'a>(Udon<'a>);



const HEADER_SIZE: usize = size_of::<Udon>();

trait Transcode {
	fn transcode<'a>(cigar: &[Cigar], packed_query: &[u8], query_length: usize, md: &[u8]) -> Option<Box<Udon<'a>>>;
}


/** transcoder helper functions **/

/*
 * writer
 */
trait Writer<T: Copy + Default> {
	type T;

	/*
	 * func supposed to return object count
	 */
	fn reserve_to<U, F>(&mut self, count: usize, func: F) -> Range<usize>
		where U: Copy,
		      F: FnOnce(&mut [U], &[T]) -> usize;
}

impl<T: Copy + Default> Writer<T> for Vec<T> {

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
		let work = unsafe { from_raw_parts_mut::<U>(ptr as *mut U, count) };

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


/*
 * pointer arithmetic for creating precursor and converting precursor to product
 */
fn compose_precursor<'a, T>(range: Range<usize>) -> &'a [T] {
	assert!(((range.end - range.start) % size_of::<T>()) == 0);

	unsafe {
		from_raw_parts(
			range.start as *const T,
			(range.end - range.start) / size_of::<T>()
		)
	}
}

fn finalize_precursor<'a, T>(precursor: &'a [T], base: *const u8) -> &'a [T] {
	let offset = precursor.as_ptr() as usize;
	unsafe {
		from_raw_parts(
			base.wrapping_add(offset) as *const T,
			precursor.len()
		)
	}
}


/*
 * convert ascii-printed number to u64, overflow is ignored.
 */
fn atoi_unchecked(v: &[u8]) -> (u64, &[u8]) {
	let mut n: u64 = 0;

	for (i, x) in v.iter().enumerate() {
		let m = (*x as u64) - ('0' as u64);

		/* return if non-numeric number found */
		if m >= 10 {
			return(n, v.split_at(i).1);
		}

		/* digit continues */
		n = n * 10 + m;
	}
	return(n, v.split_at(v.len()).1);
}

fn isnum(c: u8) -> bool {
	return (c as u64 - '0' as u64) < 10 as u64;
}

/*
 * transcode { 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't' } -> { 0x0, 0x1, 0x2, 0x3, 0x0, 0x1, 0x2, 0x3 }
 * and then mark mismatch.
 */
fn encode_base_unchecked(c: u8) -> u32 {
	let c = c as u32;
	let b2 = ((c>>1) - (c>>3)) & 0x03;
	return CompressMark::Mismatch as u32 + b2;
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


/* transcoder state object */
pub struct UdonBuilder<'a> {
	buf: Vec<u8>,
	ins: Vec<u8>,
	cigar: Iter<'a, Cigar>,
	mdstr: &'a [u8],
	query: &'a [u8],
	qofs: usize,
	rspan: usize
}

impl<'a, 'b> UdonBuilder<'a> {

	fn reserve_head(&mut self) -> Option<Range<usize>> {

		/* header always at the head */
		let range = self.buf.reserve_to(1, |header: &mut [Udon], _: &[u8]| {
			header[0] = Default::default();
			1
		});

		assert!(range.end - range.start == HEADER_SIZE);
		return Some(range);
	}

	/*
	 * input parsers
	 */
	fn next_base(&mut self) -> u32 {
		/*
		 * imitate iterator on 4-bit packed nucleotide sequence
		 */
		let ofs = self.qofs;
		self.qofs += 1;

		let c = self.query[ofs / 2];
		if (ofs & 0x01) == 0x01 {
			return encode_base_unchecked(c & 0x0f);
		}
		return encode_base_unchecked(c>>4);
	}

	fn peek_cigar(&self) -> Option<Cigar> {
		if self.cigar.as_slice().len() == 0 {
			return None;
		}

		return Some(self.cigar.as_slice()[0]);
	}

	fn is_double_mismatch(&self) -> bool {
		if self.mdstr.len() < 2 {
			return false;
		}
		return !isnum(self.mdstr[1]);
	}

	fn eat_md_eq(&mut self) -> usize {
		let (n, s) = atoi_unchecked(self.mdstr);

		self.mdstr = s;
		return n as usize;
	}

	/*
	 * writers, on top of impl Writer for Vec<u8>
	 */
	fn push_op(&mut self, match_len: usize, marker: u32) {
		assert!(marker    < 0x08);
		assert!(match_len < 0x20);

		self.buf.reserve_to(1, |arr: &mut [u8], _: &[u8]| -> usize {
			arr[0] = (marker<<5) as u8 | match_len as u8;
			1
		});
	}

	fn push_match(&mut self, match_len: usize, last_op: u32) {
		let mut op = last_op;
		let mut rem = match_len;

		while rem > 30 {
			self.push_op(31, op);
			op   = 0;
			rem -= 30;
		}
		self.push_op(rem, op);
	}

	fn push_ins(&mut self, qofs: usize, qlen: usize) {
		/* copy subsequence; 4-bit packed */
		let copy_len = (qlen + 1) / 2;
		let packet_size = copy_len + 1;

		let query = self.query;
		self.ins.reserve_to(packet_size, |arr: &mut [u8], _: &[u8]| -> usize {
			/* marker: valid bases never generate zero byte thus can be used for boundary marker */
			if let Some((head, payload)) = arr.split_first_mut() {
				*head = 0;

				/* payload */
				if (qofs & 0x01) == 0x01 {
					/* shift needed */
					payload.copy_from_slice(&query[qofs..(qofs + copy_len)]);
				} else {
					/* just copy except for the last column */
					payload.copy_from_slice(&query[qofs..(qofs + copy_len)]);
				}

				packet_size
			} else {
				0
			}
		});
	}

	fn eat_del(&mut self) -> Option<u32> {
		let c = self.cigar.next()?;
		assert!(c.op() == CigarOp::Del as u32);

		let len = c.len() as usize;
		if self.mdstr.len() < len + 1 {
			return None;
		}

		let mut rem = len;
		while rem > 3 {
			self.push_op(3, 3);
			rem -= 3;
		}

		assert!(rem > 0);
		return Some(rem as u32);
	}

	fn eat_ins(&mut self, xrem: usize) -> Option<usize> {
		assert!(self.qofs <= std::i32::MAX as usize);

		let c = self.cigar.next()?;
		assert!(c.op() == CigarOp::Ins as u32);

		let ofs = self.qofs;
		let len = c.len() as usize;

		self.qofs += len;
		self.push_ins(ofs, len);
		return Some(xrem - len);
	}

	fn eat_match(&mut self, xrem: usize, last_op: u32) -> Option<usize> {
		let c = self.cigar.next()?;

		let is_valid = c.op() != CigarOp::Ins as u32;
		if is_valid { self.cigar.next()?; }

		let mut crem = if is_valid { c.len() as usize } else { 0 } + xrem;
		let mut xrem = xrem;
		let mut op = last_op;

		while xrem < crem {
			self.push_match(xrem, op);
			crem      -= xrem;
			self.qofs += xrem;

			while self.is_double_mismatch() {
				let c = self.next_base();
				self.push_op(1, c);
				crem -= 1;
			}
			op = self.next_base();

			/* next match length */
			xrem = self.eat_md_eq();
		}

		self.push_match(crem, op);
		xrem      -= crem;
		self.qofs += crem;
		return Some(xrem);
	}

	fn parse_cigar(&mut self) -> Option<Range<usize>> {

		let base_offset = self.buf.len();
		let mut xrem = 0;

		let c = self.peek_cigar()?;
		if c.op() == CigarOp::Match as u32 {
			self.push_ins(0, 0);

			xrem = self.eat_md_eq();
			xrem = self.eat_match(xrem, CompressMark::Ins as u32)?;
		}

		loop {
			let c = self.peek_cigar()?;
			while c.op() == CigarOp::Del as u32 {
				if self.cigar.as_slice().len() < 2 { break; }

				let op = self.eat_del()?;
				xrem = self.eat_md_eq() + op as usize;
				xrem = self.eat_match(xrem, op)?;
			}

			let c = self.peek_cigar()?;
			if c.op() != CigarOp::Ins as u32 {
				return None;
			}

			if self.cigar.as_slice().len() < 2 { break; }

			xrem = self.eat_ins(xrem)?;
			xrem = self.eat_match(xrem, CompressMark::Ins as u32)?;
		}

		if self.cigar.as_slice().len() > 0 {
			let c = self.peek_cigar()?;

			if c.op() == CigarOp::Del as u32 {
				let op = self.eat_del()?;
				self.push_op(op as usize, op);
			} else if c.op() == CigarOp::Ins as u32 {
				self.eat_ins(xrem)?;
			}
		}

		return Some(Range::<usize> {
			start: base_offset,
			end:   self.buf.len()
		});
	}

	fn build_chunk(&mut self, op_range: &Range<usize>, query_length: usize) -> Option<Range<usize>> {

		/* chunk comes right after [header, op_array] */
		assert!(CHUNK_PITCH.next_power_of_two() == CHUNK_PITCH);

		/* we need this for storing return value */
		let mut ref_span = 0;

		/* rip buffer for disjoint ownership */
		let buf = &mut self.buf;

		/* FIXME: we need to cover entire query_length */
		let chunk_count = (query_length + CHUNK_PITCH - 1) / CHUNK_PITCH;
		let range = buf.reserve_to(chunk_count, |chunk: &mut [Chunk], base: &[u8]| {
			let ops = &base[op_range.start..op_range.end];

			let mut rpos: usize = 0;
			let mut iofs: usize = 0;
			let mut cofs: usize = 0;
			for (i, &op) in ops.iter().enumerate() {
				let next_boundary = (rpos | (CHUNK_PITCH - 1)) + 1;

				/* forward insertion array and reference-side offset */
				let marker = op>>5;
				let len = op as usize & 0x1f;

				iofs += if marker == 0 { 1 } else { 0 };
				rpos += if len == 31 { len - 1 } else { len };

				if rpos < next_boundary { continue; }

				/* save chunk info */
				chunk[cofs].set_ins_offset(iofs as u64);
				chunk[cofs].set_op_offset(i as u64);
				chunk[cofs].set_op_skip(rpos as u64 & 0x1f);
				cofs += 1;
			}

			/* save ref_span */
			ref_span = rpos;

			assert!(cofs == chunk_count);
			cofs
		});

		self.rspan = ref_span;
		return Some(range);
	}

	fn pack_ins(&mut self) -> Option<Range<usize>> {

		/* rip vectors from self for ownership */
		let buf = &mut self.buf;
		let ins = &self.ins;

		let range = buf.reserve_to(ins.len(), |arr: &mut [u8], _: &[u8]| {
			arr.copy_from_slice(ins.as_slice());
			ins.len()
		});
		return Some(range);
	}

	fn finalize(&self, header: Range<usize>, op: Range<usize>, chunk: Range<usize>, ins: Range<usize>) -> UdonPrecursor<'b> {
		
		assert!(header.end == op.start);
		assert!(op.end <= chunk.start);
		assert!(chunk.end <= ins.start);

		assert!(((chunk.end - chunk.start) % size_of::<Chunk>()) == 0);

		/* entire range on buffer for this object */
		let range = Range::<usize> {
			start: header.start,
			end: ins.end
		};

		/* convert range to slice */
		UdonPrecursor::<'b> (
			Udon::<'b> {
				/* just save */
				size: range.end - range.start,
				ref_span: self.rspan,

				/* compose fake slices (dereference causes SEGV) */
				op: compose_precursor(op),
				chunk: compose_precursor(chunk),
				ins: compose_precursor(ins)
			}
		)
	}

	/*
	 * APIs
	 */
	pub fn build_precursor_raw(buf: Vec<u8>, cigar: &'a [Cigar], packed_query: &'a [u8], query_length: usize, mdstr: &'a [u8])
		-> (Vec<u8>, Option<UdonPrecursor<'b>>)
	{
		/* save initial offset for unwinding */
		let base_offset = buf.len();

		/* compose working variables */
		let (cigar, qclip) = match strip_clips(cigar) {
			None => { return (buf, None); },
			Some((cigar, qclip)) => (cigar, qclip)
		};
		let mut state = Self {
			buf:   buf,				/* move */
			ins:   Vec::new(),
			cigar: cigar.iter(),
			mdstr: mdstr,
			query: packed_query,
			qofs:  qclip.head,
			rspan: 0				/* calcd in build_chunk */
		};

		/* if error detected, unwind destination vector and return it (so that the vector won't lost) */
		macro_rules! unwrap_or_return { ( $expr:expr ) => {
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
		let header = unwrap_or_return!(state.reserve_head());
		let op = unwrap_or_return!(state.parse_cigar());
		let chunk = unwrap_or_return!(state.build_chunk(&op, query_length));
		let ins = unwrap_or_return!(state.pack_ins());

		/* everything done; compose precursor */
		let precursor = state.finalize(header, op, chunk, ins);

		/* return back ownership of the buffer */
		let buf = state.buf;
		return (buf, Some(precursor));
	}

	pub fn build_precursor(cigar: &'a [Cigar], packed_query: &'a [u8], query_length: usize, mdstr: &'a [u8])
		-> Option<(Vec<u8>, UdonPrecursor<'b>)>
	{
		let buf = Vec::<u8>::new();

		match Self::build_precursor_raw(buf, cigar, packed_query, query_length, mdstr) {
			(_, None) => {
				return None;
			},
			(buf, Some(precursor)) => {
				return Some((buf, precursor));
			}
		}
	}
}

impl<'a> Udon<'a> {
	pub fn from_precursor(buf: &'a [u8], precursor: UdonPrecursor<'a>) -> Udon<'a> {
		let base: *const u8 = buf.as_ptr();

		Udon::<'a> {
			/* just copy */
			size: precursor.0.size,
			ref_span: precursor.0.ref_span,

			/* adjusting pointer */
			op: finalize_precursor(precursor.0.op, base),
			chunk: finalize_precursor(precursor.0.chunk, base),
			ins: finalize_precursor(precursor.0.ins, base)
		}
	}
}

/*
impl Transcode for Udon<'_> {
	fn transcode<'a>(cigar: &[Cigar], packed_query: &[u8], query_length: usize, md: &[u8]) -> Option<Box<Udon<'a>>> {

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
