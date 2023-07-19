use super::op::op_len;
use super::utils::PeekFold;
// use std::alloc::{ Layout, dealloc };
// use std::mem::align_of;
use std::ops::Range;

/* Transcoded CIGAR and index object container

everything is laid out in a single flat memory block, thus constructed via
some unsafe operations. see `Precursor` and `SlicePrecursor`.

Impls are found in `decoder.rs` and `insbin.rs`. Builder implmentation is in `builder.rs`/
*/

/* block pitch, see above */
pub(super) const BLOCK_PITCH: usize = 256;

/* block object, see above */
bitfield! {
    #[derive(Copy, Clone, Debug, Default)]
    pub struct Block(u64);
    pub ins_offset, set_ins_offset: 28, 0;
    pub op_offset,  set_op_offset:  58, 29;
    pub op_skip,    set_op_skip:    63, 59;
}

#[derive(Default)]
pub(super) struct Index<'a> {
    pub(super) size: usize,     /* object size for copying */
    pub(super) ref_span: usize, /* reference side span */
    pub(super) op: &'a [u8],
    pub(super) block: &'a [Block],
    pub(super) ins: &'a [u8],
}

/* common; builders and decoders are found in `builder.rs`, `decoder.rs`, and `insbin.rs`. */
impl<'a> Index<'a> {
    /* API redirected */
    pub(super) fn reference_span(&self) -> usize {
        self.ref_span
    }

    /* Check sanity of the span. If queried span (range) is out of the indexed one, return None */
    pub(super) fn check_span(&self, range: &Range<usize>) -> Option<()> {
        // debug!("span({}, {}), ref_span({})", range.start, range.end, self.ref_span);

        if range.end < range.start {
            return None;
        }
        if range.end > self.ref_span {
            return None;
        }
        Some(())
    }

    /* fetch block head for this pos. */
    pub(super) fn fetch_ops_block(&self, pos: usize) -> (u8, &[u8], usize) {
        let block_index = pos / BLOCK_PITCH;
        let block_rem = pos & (BLOCK_PITCH - 1);
        let block = &self.block[block_index];

        let op_offset = block.op_offset() as usize;
        let op_skip = block.op_skip() as usize;
        // let ins_offset = block.ins_offset() as usize;

        let ops = &self.op[op_offset..];
        let last_op = if op_offset == 0 {
            0
        } else {
            self.op[op_offset - 1]
        };

        // debug!("pos({}), rem({}), skip({})", pos, block_rem, op_skip);
        (last_op, ops, block_rem + op_skip)
    }

    pub(super) fn scan_op_array(&self, pos: usize) -> (&[u8], usize) {
        /* get block head for this pos */
        let (_, ops, rem) = self.fetch_ops_block(pos);
        // debug!("rem({}), ops({:?})", rem, ops);

        /* linear polling */
        let mut ops = ops.iter();
        let ofs = ops.peek_fold(0, |a, &x| {
            let len = a + op_len(x);

            /* continue if at least one column remaining */
            if len >= rem {
                return None;
            }
            Some(len)
        });

        // debug!("rem({}), ofs({})", rem, ofs);
        return (ops.as_slice(), rem - ofs); /* is this sound? (I'm not sure...) */
    }
}
