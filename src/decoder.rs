

use super::UdonOp;
use super::utils::SimdAlignedU8;
use super::op::{ CompressMark, op_len, op_marker, op_is_cont };
use super::index::Index;
use super::scaler::Scaler;
use std::io::Write;
use std::ops::Range;


/* architecture-dependent stuffs */
#[cfg(all(target_arch = "x86_64", target_feature = "ssse3"))]
use core::arch::x86_64::*;

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
use core::arch::aarch64::*;


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


We also provide unsafe APIs, `Precursor::build` and `Udon::from_precursor`. This two
methods were originally implemented for use in `UdonVec::append` and `UdonVec::freeze` each.
The first function creates `Precursor` object, which is almost equivalent to `Udon`
except that pointers in slices are not valid. The invalid pointer in slice actually
retains an offset from the memory region. Being offset, not an absolute pointer, it allows
reallocation (on-the-fly extension) of the memory region without breaking consistency,
which is important for creating multiple objects in a single region one by one. The
offsets are converted to valid slices at once in the `freeze` procedure.

The raw unsafe APIs `Precursor::build` and `Udon::from_precursor` are useful if we
want to pack udon objects in a more complex way, not a flat vector. It is especially
useful when we design a data structure with additional metadata for visualization, by
keeping `Udon` object as one member of larger struct.

Please make sure the `Precursor::build` and `Udon::from_precursor` are unsafe APIs.
The unsafeness comes from the requirement that a vector precursor was built and a vector
the precursor is freezed are the same one, which can't be forced by the type system.
If the two vector don't match, the created `Udon` objects are definitely broken and
causes exception.
*/


/* decoder implementation */
impl<'a, 'b> Index<'a> {

	pub(super) fn decode_raw_into(&self, dst: &mut Vec<u8>, ref_span: &Range<usize>) -> Option<usize> {
		self.check_span(&ref_span)?;
		self.decode_core(dst, &ref_span)
	}

	pub(super) fn decode_raw(&self, ref_span: &Range<usize>) -> Option<Vec<u8>> {
		self.check_span(&ref_span)?;

		let size = ref_span.end - ref_span.start;
		let mut buf = Vec::<u8>::with_capacity(size);

		let used = self.decode_core(&mut buf, &ref_span)?;
		buf.resize(used, 0);

		return Some(buf);
	}

	pub(super) fn decode_scaled_into(&self, dst: &mut Vec<[[u8; 4]; 2]>, ref_span: &Range<usize>, offset_in_pixels: f64, scaler: &Scaler) -> Option<usize> {
		self.check_span(&ref_span)?;
		self.decode_scaled_core(dst, &ref_span, offset_in_pixels, &scaler)
	}

	pub(super) fn decode_scaled(&self, ref_span: &Range<usize>, offset_in_pixels: f64, scaler: &Scaler) -> Option<Vec<[[u8; 4]; 2]>> {
		self.check_span(&ref_span)?;

		let span = ref_span.end - ref_span.start;
		let size = scaler.expected_size(span);
		let mut buf = Vec::<[[u8; 4]; 2]>::with_capacity(size);
		// debug!("ref_span({:?}), span({}), size({})", ref_span, span, size);

		let used = self.decode_scaled_core(&mut buf, &ref_span, offset_in_pixels, &scaler)?;
		buf.resize(used, [[0; 4]; 2]);
		// debug!("used({})", used);

		return Some(buf);
	}

	const DEL_MASK: SimdAlignedU8 = {
		let mut x = [0u8; 16];
		x[0] = UdonOp::Del as u8;
		x[1] = UdonOp::Del as u8;
		x[2] = UdonOp::Del as u8;
		SimdAlignedU8::new(&x)
	};
	const SCATTER_MASK: SimdAlignedU8 = {
		let mut x = [0x80u8; 16];
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
		SimdAlignedU8::new(&x)
	};
	const IS_DEL_THRESH: SimdAlignedU8 = {
		let mut x = [0xffu8; 16];
		x[0] = 0x1f;
		x[1] = 0x3f;
		x[2] = 0x5f;
		SimdAlignedU8::new(&x)
	};

	#[cfg(all(target_arch = "x86_64", target_feature = "ssse3"))]
	unsafe fn decode_core_block(dst: &mut [u8], op: u8, ins: u32) -> (usize, u32) {

		/* load constants; expelled out of the innermost loop when inlined */
		let del_mask      = _mm_load_si128(Self::DEL_MASK.as_ref() as *const [u8; 16] as *const __m128i);
		let scatter_mask  = _mm_load_si128(Self::SCATTER_MASK.as_ref() as *const [u8; 16] as *const __m128i);
		let is_del_thresh = _mm_load_si128(Self::IS_DEL_THRESH.as_ref() as *const [u8; 16] as *const __m128i);

		/* compute deletion mask */
		let xop = _mm_cvtsi64_si128(op as i64);					/* copy op to xmm */
		let xop = _mm_shuffle_epi8(xop, scatter_mask);			/* [op, op, op, 0, ...] */
		let is_del = _mm_cmpgt_epi8(xop, is_del_thresh);		/* signed comparison */

		/* compute mismatch / insertion mask */
		let marker = if op_marker(op) == CompressMark::Ins as u32 { ins } else { op_marker(op) };
		let ins_mismatch_mask = _mm_cvtsi64_si128(marker as i64);

		/* merge deletion / insertion-mismatch vector */
		let merged = _mm_blendv_epi8(ins_mismatch_mask, del_mask, is_del);

		_mm_storeu_si128(&mut dst[0] as *mut u8 as *mut __m128i, merged);
		_mm_storeu_si128(&mut dst[16] as *mut u8 as *mut __m128i, _mm_setzero_si128());

		/*
		compute forward length; 31 is "continuous marker"
		rop-to-rop critical path length is 6
		*/
		let next_ins     = if op_is_cont(op) { 0 } else { UdonOp::Ins as u32 };	/* next ins will be masked if 0x1f */
		let adjusted_len = op_len(op);

		// debug!("{:#x}, {:#x}, {:#x}, {:#x}", op>>5, ins, marker, op_len(op));
		// debug!("{}, {:?}", adjusted_len, transmute::<__m128i, [u8; 16]>(merged));

		(adjusted_len, next_ins)
	}

	#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
	unsafe fn decode_core_block(dst: &mut [u8], op: u8, ins: u32) -> (usize, u32) {
		/* the following code correspond one-to-one to the x86_64 implementation above */
		let del_mask      = vld1q_u8(Self::DEL_MASK.as_ref() as *const u8);
		let is_del_thresh = vld1q_s8(Self::IS_DEL_THRESH.as_ref() as *const u8 as *const i8);

		let xop = [op, op, op, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
		let xop = vld1q_s8(&xop as *const u8 as *const i8);
		let is_del = vcgtq_s8(xop, is_del_thresh);			/* signed comparison */

		let marker = if op_marker(op) == CompressMark::Ins as u32 { ins } else { op_marker(op) };
		let ins_mismatch_mask = vsetq_lane_u8(marker as u8, vmovq_n_u8(0), 0);

		let merged = vbslq_u8(is_del, del_mask, ins_mismatch_mask);
		let zero = vmovq_n_u8(0);
		vst1q_u8(&mut dst[0] as *mut u8, merged);
		vst1q_u8(&mut dst[16] as *mut u8, zero);

		let next_ins     = if op_is_cont(op) { 0 } else { UdonOp::Ins as u32 };	/* next ins will be masked if 0x1f */
		let adjusted_len = op_len(op);
		(adjusted_len, next_ins)
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

	fn decode_scaled_core(&self, dst: &mut Vec<[[u8; 4]; 2]>, ref_span: &Range<usize>, offset_in_pixels: f64, scaler: &Scaler) -> Option<usize> {

		/* states (working variables) */
		let (mut offset, margin) = scaler.init(offset_in_pixels);
		let mut dst = dst;

		// println!("offset({}), margin({})", offset, margin);

		/* buffer */
		let bulk_size: usize = 16 * 1024;		/* 16KB */
		let mut buf = Vec::<u8>::with_capacity(bulk_size + margin);
		for _ in 0 .. margin { buf.push(0); }

		for spos in ref_span.clone().step_by(bulk_size) {
			/* decode a block */
			let decoded = self.decode_core(&mut buf, &Range::<usize> {
				start: spos,
				end:   ref_span.end.min(spos + bulk_size)
			})?;
			if decoded < bulk_size {
				for _ in 0 .. margin + 1 { buf.push(0); }
			}

			/* rescale to dst array, forward offset */
			let (next_offset, consumed) = scaler.scale(&mut dst, &buf, offset)?;
			offset = next_offset;
			if buf.len() < consumed || consumed < buf.len() - consumed { continue; }

			let (base, tail) = buf.split_at_mut(consumed);
			let (base, _) = base.split_at_mut(tail.len());
			base.copy_from_slice(tail);

			let next_margin = base.len();
			buf.resize(next_margin, 0);

			// println!("margin({}), len({}), decoded({}), consumed({}), next_margin({})", margin, buf.len(), decoded, consumed, next_margin);
		}
		Some(dst.len())
	}
}

