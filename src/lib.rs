
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

Note: This file provides only APIs. See `builder.rs` and `decoder.rs` for the internals.
*/

mod utils;
mod op;
mod index;
mod builder;
mod insbin;
mod decoder;
mod scaler;

use self::{
	index::Index,
	scaler::Scaler
};

use std::ops::Range;


#[macro_use]
extern crate bitfield;

/* logging */
#[macro_use]
extern crate log;

/* Event for output (expanded) array. */
#[repr(u8)]
pub enum UdonOp {
	MisA = 0x04 | 0x00,
	MisC = 0x04 | 0x01,
	MisG = 0x04 | 0x02,
	MisT = 0x04 | 0x03,
	Del  = 0x08,
	Ins  = 0x10				/* may be OR-ed with the others */
}


/* aliasing type for convenience */
pub type UdonColorPair = [[u8; 4]; 2];


/* udon object */
pub struct Udon<'o>(Index<'o>);


/* wrapping is unnecessary; to keep this file looks simple... */
impl<'i, 'o> Udon<'o> {

	/* builders */
	pub fn build(cigar: &'i [u32], packed_query: &'i [u8], mdstr: &'i [u8]) -> Option<Box<Udon<'o>>> {
		let index = Index::build(cigar, packed_query, mdstr)?;
		let udon = unsafe {
			Box::<Udon<'o>>::from_raw(Box::into_raw(index) as *mut Udon)
		};
		Some(udon)
	}

	/* UdonPrecursor is not exported for now.
	pub unsafe fn from_precursor(buf: &Vec<u8>, precursor: UdonPrecursor) -> Udon<'a> {
		Index::from_precursor(buf, precursor)
	}
	*/


	/* reference side span */
	pub fn reference_span(&self) -> usize {
		self.0.reference_span()
	}


	/* decoders, output is actually Vec<Op> */
	pub fn decode_raw_into(&self, dst: &mut Vec<u8>, ref_span: &Range<usize>) -> Option<usize> {
		self.0.decode_raw_into(dst, ref_span)
	}

	pub fn decode_raw(&self, ref_span: &Range<usize>) -> Option<Vec<u8>> {
		self.0.decode_raw(ref_span)
	}

	pub fn decode_scaled_into(&self, dst: &mut Vec<UdonColorPair>, ref_span: &Range<usize>, offset_in_pixels: f64, scaler: &UdonScaler) -> Option<usize> {
		self.0.decode_scaled_into(dst, ref_span, offset_in_pixels, &scaler.0)
	}

	pub fn decode_scaled(&self, ref_span: &Range<usize>, offset_in_pixels: f64, scaler: &UdonScaler) -> Option<Vec<UdonColorPair>> {
		self.0.decode_scaled(ref_span, offset_in_pixels, &scaler.0)
	}


	/* fetch insertion sequence at a specific position */
	pub fn get_ins_into(&self, dst: &mut Vec<u8>, pos: usize) -> Option<usize> {
		self.0.get_ins_into(dst, pos)
	}

	pub fn get_ins(&self, pos: usize) -> Option<Vec<u8>> {
		self.0.get_ins(pos)
	}
}


/* UdonPalette

Color table for scaled decoder. Required to build `UdonScaler`.
`UdonColorPair` is defined as a pair of RGBa8. In the default color scheme,
the first channel is used for normal ribbon, and the second channel is for
insertion markers. The two channels are treated without any difference,
they can be used for arbitrary purposes.

See also: `UdonColorPair`, `UdonScaler`, and `Udon::decode_scaled`.
*/
#[derive(Copy, Clone, Debug)]
pub struct UdonPalette {
	/* all in [(r, g, b, alpha); 2] form; the larger alpha value, the more transparent */
	background:  UdonColorPair,
	del: UdonColorPair,
	ins: UdonColorPair,
	mismatch: [UdonColorPair; 4]
}

impl Default for UdonPalette {
	fn default() -> UdonPalette {
		UdonPalette {
			background: [[255, 255, 255, 255], [255, 255, 255, 255]],
			del: [[255, 255, 255,   0], [255, 255, 255, 255]],		/* white (transparent) */
			ins: [[255, 255, 255, 255], [153,   0, 153,   0]],		/* black (transparent) */
			mismatch: [
				[[  3, 175,  64, 0], [255, 255, 255, 255]],		/* green */
				[[  0,  90, 255, 0], [255, 255, 255, 255]],		/* blue */
				[[ 80,  32,  64, 0], [255, 255, 255, 255]],		/* dark brown */
				[[255,   0,   0, 0], [255, 255, 255, 255]]		/* yellow */
			]
		}
	}
}


/* UdonScaler

Holds constants for the scaled decoder. Built from `UdonPalette`.
*/
pub struct UdonScaler(Scaler);

impl UdonScaler {
	pub fn new(color: &UdonPalette, columns_per_pixel: f64) -> UdonScaler {
		UdonScaler(Scaler::new(color, columns_per_pixel))
	}
}


/* UdonUtils

Provides utilities on ribbon (&mut [u32]): blending and gamma correction
*/
pub trait UdonUtils {
	fn append_on_basecolor(&mut self, basecolor: &UdonColorPair) -> &mut Self;
	fn correct_gamma(&mut self) -> &mut Self;
}

impl UdonUtils for [UdonColorPair] {

	fn append_on_basecolor(&mut self, basecolor: &UdonColorPair) -> &mut Self {
		for x in self.iter_mut() {
			let alpha0 = x[0][3] as u16;
			let alpha1 = x[1][3] as u16;

			for i in 0 .. 4 {
				let base0 = basecolor[0][i] as u16;
				let base1 = basecolor[1][i] as u16;
				let color0 = (255 - x[0][i]) as u16;
				let color1 = (255 - x[1][i]) as u16;

				let color0 = base0 * (256 - alpha0) + color0 * alpha0;
				let color1 = base1 * (256 - alpha1) + color1 * alpha1;

				x[0][i] = (color0>>8) as u8;
				x[1][i] = (color1>>8) as u8;
			}
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
			for i in 0 .. 4 {
				x[0][i] = GAMMA[x[0][i] as usize];
				x[1][i] = GAMMA[x[1][i] as usize];
			}
		}
		self
	}
}


#[cfg(test)]
mod test {
	use std::ops::Range;
	use std::str::from_utf8;

	#[allow(unused_imports)]
	use super::{ Udon, UdonColorPair, UdonPalette, UdonScaler };
	use super::utils::{ encode_base_unchecked, decode_base_unchecked };
	use super::op::CigarOp;
	use super::index::BLOCK_PITCH;

	macro_rules! cigar {
		[ $( ( $op: ident, $len: expr ) ),* ] => ({
			vec![ $( CigarOp::$op as u32 | (( $len as u32 )<<4) ),* ]
		});
	}

	macro_rules! nucl {
		( $st: expr ) => ({
			let s = $st;
			// println!("{:?}", s);

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
				r.end = u.reference_span();
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

	const BG: [[u8; 4]; 2]   = [[0xff, 0xff, 0xff, 0xff], [0xff, 0xff, 0xff, 0xff]];
	const DEL: [[u8; 4]; 2]  = [[0x00, 0x00, 0xff, 0x00], [0x00, 0x00, 0xff, 0x00]];
	const INS: [[u8; 4]; 2]  = [[0xff, 0x00, 0xff, 0x00], [0xff, 0x00, 0xff, 0x00]];
	const MISA: [[u8; 4]; 2] = [[0x7f, 0x1f, 0x1f, 0x00], [0x7f, 0x1f, 0x1f, 0x00]];
	const MISC: [[u8; 4]; 2] = [[0x00, 0x00, 0xff, 0x00], [0x00, 0x00, 0xff, 0x00]];
	const MISG: [[u8; 4]; 2] = [[0x00, 0xff, 0x00, 0x00], [0x00, 0xff, 0x00, 0x00]];
	const MIST: [[u8; 4]; 2] = [[0xff, 0x00, 0x00, 0x00], [0xff, 0x00, 0x00, 0x00]];

	macro_rules! compare_color {
		( $cigar: expr, $nucl: expr, $mdstr: expr, $range: expr, $offset: expr, $scale: expr, $ribbon: expr, $color_factor: expr, $alpha_factor: expr ) => ({
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
				r.end = u.reference_span();
			}

			let s = UdonScaler::new(
				&UdonPalette {
					background: BG,
					del: DEL,
					ins: INS,
					mismatch: [ MISA, MISC, MISG, MIST ]
				},
				$scale
			);
			let b = match u.decode_scaled(&r, $offset, &s) {
				None    => Vec::<UdonColorPair>::new(),
				Some(v) => v
			};
			let n: Vec::<UdonColorPair> = $ribbon.iter().map(|x| {
				let mut x = *x;
				for i in 0 .. 3 {
					x[0][i] = ((0xff - x[0][i]) as f64 * $color_factor) as u8;
					x[1][i] = ((0xff - x[1][i]) as f64 * $color_factor) as u8;
				}
				for i in 3 .. 4 {
					x[0][i] = ((0xff - x[0][i]) as f64 * $alpha_factor) as u8;
					x[1][i] = ((0xff - x[1][i]) as f64 * $alpha_factor) as u8;
				}

				// u32::from_le_bytes(x)
				x
			}).collect();
			// println!("{:?}, {:?}", b, n);
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
	fn test_udon_build_extended() {
		compare!(
			cigar![(Eq, 4)],
			nucl!("ACGT"),
			"4",
			Range { start: 0, end: 0 },
			"MMMM",
			"----"
		);
		compare!(
			cigar![(Eq, 30)],
			nucl!("ACGTACGTACGTACGTACGTACGTACGTAC"),
			"30",
			Range { start: 0, end: 0 },
			"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
			"------------------------------"
		);
		compare!(
			cigar![(Eq, 4), (Mismatch, 1), (Eq, 4)],
			nucl!("ACGTACGTA"),
			"4T4",
			Range { start: 0, end: 0 },
			"MMMMAMMMM",
			"---------"
		);
	}

	#[test]
	fn test_udon_build_squash() {
		compare!(
			cigar![(Match, 2), (Match, 2)],
			nucl!("ACGT"),
			"4",
			Range { start: 0, end: 0 },
			"MMMM",
			"----"
		);
		compare!(
			cigar![(Eq, 5), (Mismatch, 1), (Match, 1), (Eq, 2)],
			nucl!("ACGTACGTA"),
			"5T3",
			Range { start: 0, end: 0 },
			"MMMMMCMMM",
			"---------"
		);
		compare!(
			cigar![(Ins, 2), (Ins, 2), (Match, 4)],
			nucl!("ACGTACGT"),
			"4",
			Range { start: 0, end: 0 },
			"MMMM",
			"I---"
		);
		compare!(
			cigar![(Del, 2), (Del, 2), (Match, 4)],
			nucl!("ACGT"),
			"^ACGT4",
			Range { start: 0, end: 0 },
			"DDDDMMMM",
			"--------"
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
			vec![BG, BG, BG, BG, DEL, BG, BG, BG, BG],
			1.0 / (1.0f64.log(3.5).max(1.0) + 1.0f64 / 10.0),
			1.0 / (1.0f64.log(2.5).max(1.0) + 1.0f64 / 5.0)
		);

		compare_color!(
			cigar![(Match, 4), (Del, 1), (Match, 4)],
			nucl!("ACGTACGT"),
			"4^A4",
			Range { start: 0, end: 0 },
			0.0, 3.0,
			vec![BG, DEL, BG],
			1.0 / (3.0f64.log(3.5).max(1.0) + 3.0f64 / 10.0),
			1.0 / (3.0f64.log(2.5).max(1.0) + 3.0f64 / 5.0)
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





