

use super::utils::{ Writer, PeekFold };
use super::op::{ CompressMark, op_len, op_marker, op_is_cont };
use super::index::{ Index, BLOCK_PITCH };
use std::ops::Range;
use std::slice::Iter;


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
pub(super) struct InsTracker<'a> {
	ignore_next: bool,
	len: usize,
	ins: Iter<'a, u8>
}

impl<'a> InsTracker<'a> {
	pub(super) fn new(last_op: u8, ins: &'a [u8]) -> Self {
		InsTracker::<'a> {
			ignore_next: op_is_cont(last_op),
			len: ins.len(),
			ins: ins.iter()
		}
	}

	pub(super) fn forward(&mut self, op: u8) -> bool {
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

	pub(super) fn get_offset(&self) -> usize {
		self.len - self.ins.as_slice().len()	/* is there better way for tracking the offset? */
	}

	pub(super) fn as_slice(&self) -> &'a [u8] {
		self.ins.as_slice()
	}
}


impl<'a, 'b> Index<'a> {

	fn fetch_ins_block(&self, pos: usize) -> &[u8] {

		let block_index = pos / BLOCK_PITCH;
		let block = &self.block[block_index];

		let ins_offset = block.ins_offset() as usize;
		return &self.ins[ins_offset as usize ..];
	}

	/* ins */
	fn scan_ins_array(&self, pos: usize) -> Option<&[u8]> {

		let (last_op, ops, rem) = self.fetch_ops_block(pos);
		let ins = self.fetch_ins_block(pos);

		/* linear polling on op array */
		let mut ops = ops.iter();
		let mut ins = InsTracker::new(last_op, ins);
		let len = (&mut ops).peek_fold(0, |a, &x| {
			let len = a + op_len(x);
			if len > rem { return None; }
			ins.forward(x);
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

	pub(super) fn get_ins_into(&self, dst: &mut Vec<u8>, pos: usize) -> Option<usize> {
		self.check_span(&Range::<usize> { start: 0, end: pos })?;

		let ins = self.scan_ins_array(pos)?;
		let len = Self::get_ins_core(dst, ins);

		Some(len)
	}

	pub(super) fn get_ins(&self, pos: usize) -> Option<Vec<u8>> {
		self.check_span(&Range::<usize> { start: 0, end: pos })?;

		/* fetch ins vector */
		let ins = self.scan_ins_array(pos)?;

		/* expand packed nucleotide to Vec */
		let mut v = Vec::with_capacity(ins.len() * 2);
		Self::get_ins_core(&mut v, ins);

		Some(v)
	}
}





