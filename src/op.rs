

use std::slice::Iter;


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
pub(super) enum CigarOp {
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


/* op trailing events */
#[repr(u32)]
pub(super) enum CompressMark {
	Ins      = 0x00,

	/* mismatch base for 'A': 0x04, .., 'T': 0x07 */
	Mismatch = 0x04,
}

/*
op -> len conversion
*/
pub(super) fn op_len(x: u8) -> usize {
	let len = x as usize & 0x1f;
	len - (len == 31) as usize
}

pub(super) fn op_marker(x: u8) -> u32 {
	(x>>5) as u32
}

pub(super) fn op_is_cont(x: u8) -> bool {
	(x & 0x1f) == 0x1f
}


/* OpsIterator

Iterates over op array (&[u8]), accumulating reference-side offset.
*/
pub(super) struct OpsIter<'a> {
	iter: Iter<'a, u8>,
	rofs: usize
}

pub(super) trait IntoOpsIterator<'a> {
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

	fn next(&mut self) -> Option<Self::Item> {
		let x    = self.iter.next()?;
		self.rofs += op_len(*x);
		Some((*x, self.rofs))
	}
}



