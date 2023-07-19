use super::{UdonOp, UdonPalette};
use std::f64::consts::PI;
use std::ops::{Add, AddAssign, Mul, MulAssign, Range, Sub};

/* Scaler

Scales decoded op array into a color bar. Taking scaling factor (`columns_per_pixel`) and
output vector of `Udon::decode`, it outputs Vec<u32> whose length is `input.len() / columns_per_pixel`.
Each element of the output vector is negated 8-bit RGB color where R is placed at the least
significant byte. The negated color is then overlaid on base color (typically pink for forward
and light blue for reverse) by `UdonUtils::append_on_basecolor` then applied gamma correction by
`UdonUtils::correct_gamma`.

Scaler implements approximated Lanczos filter (sinc interpolation) for arbitrary scaling
factor. It precomputes interpolation coefficient table on `new()` for sampling input columns for
an output bin. Sixteen tables for different column offset fractions for arbitrary combination of
scaling factor and drawing offset.
*/
#[derive(Copy, Clone, Debug, Default)]
struct Color {
    v: [i32; 8],
}

impl From<&[[u8; 4]; 2]> for Color {
    fn from(val: &[[u8; 4]; 2]) -> Color {
        let mut x = Color::default();
        for i in 0..4 {
            x.v[i] = val[0][i] as i32;
        }
        for i in 0..4 {
            x.v[i + 4] = val[1][i] as i32;
        }
        x
    }
}

impl From<&Color> for [[u8; 4]; 2] {
    fn from(val: &Color) -> [[u8; 4]; 2] {
        let mut x: [[u8; 4]; 2] = Default::default();
        for i in 0..4 {
            x[0][i] = val.v[i].min(255).max(0) as u8;
            x[1][i] = val.v[i + 4].min(255).max(0) as u8;
        }
        x
    }
}
/*
impl From<&Color> for (u32, u32) {
    fn from(val: &Color) -> (u32, u32) {
        let mut x: [[u8; 4]; 2] = Default::default();
        for i in 0 .. 4 {
            x[0][i] = val.v[i    ].min(255).max(0) as u8;
            x[1][i] = val.v[i + 4].min(255).max(0) as u8;
        }
        (u32::from_le_bytes(x[0]), u32::from_le_bytes(x[1]))
    }
}
*/

impl Add<Color> for Color {
    type Output = Color;
    fn add(self, other: Color) -> Color {
        let mut x = self;
        for i in 0..8 {
            x.v[i] += other.v[i];
        }
        x
    }
}
impl Sub<Color> for Color {
    type Output = Color;
    fn sub(self, other: Color) -> Color {
        let mut x = self;
        for i in 0..8 {
            x.v[i] -= other.v[i];
        }
        x
    }
}
impl Mul<Color> for Color {
    type Output = Color;
    fn mul(self, other: Color) -> Color {
        let mut x = self;
        for i in 0..8 {
            /* upcast, multiply, then downcast. (expect pmuludq) */
            let n = x.v[i] as i64 * other.v[i] as i64;
            x.v[i] = (n >> 24) as i32; /* I expect this won't overflow but no guarantee */
        }
        x
    }
}

impl Add<i32> for Color {
    type Output = Color;
    fn add(self, other: i32) -> Color {
        let mut x = self;
        for i in 0..8 {
            x.v[i] += other;
        }
        x
    }
}
impl Mul<i32> for Color {
    type Output = Color;
    fn mul(self, other: i32) -> Color {
        let mut x = self;
        for i in 0..8 {
            let n = x.v[i] as i64 * other as i64;
            x.v[i] = (n >> 24) as i32;
        }
        x
    }
}

/* arithmetic assign */
impl AddAssign<Color> for Color {
    fn add_assign(&mut self, other: Color) {
        *self = self.add(other);
    }
}
impl MulAssign<Color> for Color {
    fn mul_assign(&mut self, other: Color) {
        *self = self.mul(other);
    }
}

/* Scaler second impl */
#[derive(Default)]
pub(super) struct Scaler {
    columns_per_pixel: f64,
    window: f64,
    offset: f64,
    normalizer: Color,
    color: [Color; 12],
    table: [Vec<i32>; 17],
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
        if x < window { 0.0 } else { x - window }
    } else if x > -window { 0.0 } else { x + window }
}

impl Scaler {
    /* column -> color index */
    const INDEX: [u8; 32] = {
        let mut index = [0; 32];
        index[UdonOp::MisA as usize] = 1;
        index[UdonOp::MisC as usize] = 2;
        index[UdonOp::MisG as usize] = 3;
        index[UdonOp::MisT as usize] = 4;
        index[UdonOp::Del as usize] = 5;
        index[UdonOp::Ins as usize | UdonOp::MisA as usize] = 6;
        index[UdonOp::Ins as usize | UdonOp::MisC as usize] = 7;
        index[UdonOp::Ins as usize | UdonOp::MisG as usize] = 8;
        index[UdonOp::Ins as usize | UdonOp::MisT as usize] = 9;
        index[UdonOp::Ins as usize | UdonOp::Del as usize] = 10;
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
        let ff = Color::from(&[[0xff, 0xff, 0xff, 0xff], [0xff, 0xff, 0xff, 0xff]]);
        let mismatch: [Color; 4] = [
            ff - Color::from(&color.mismatch[0]),
            ff - Color::from(&color.mismatch[1]),
            ff - Color::from(&color.mismatch[2]),
            ff - Color::from(&color.mismatch[3]),
        ];
        let del = ff - Color::from(&color.del);
        let ins = ff - Color::from(&color.ins);
        let bg = ff - Color::from(&color.background);

        let mut x = [bg; 12];
        x[Self::index(UdonOp::MisA as u8)] = mismatch[0];
        x[Self::index(UdonOp::MisC as u8)] = mismatch[1];
        x[Self::index(UdonOp::MisG as u8)] = mismatch[2];
        x[Self::index(UdonOp::MisT as u8)] = mismatch[3];
        x[Self::index(UdonOp::Del as u8)] = del;
        x[Self::index(UdonOp::Ins as u8 + UdonOp::MisA as u8)] =
            ins * 0x800000 + mismatch[0] * 0x800000;
        x[Self::index(UdonOp::Ins as u8 + UdonOp::MisC as u8)] =
            ins * 0x800000 + mismatch[1] * 0x800000;
        x[Self::index(UdonOp::Ins as u8 + UdonOp::MisG as u8)] =
            ins * 0x800000 + mismatch[2] * 0x800000;
        x[Self::index(UdonOp::Ins as u8 + UdonOp::MisT as u8)] =
            ins * 0x800000 + mismatch[3] * 0x800000;
        x[Self::index(UdonOp::Ins as u8 + UdonOp::Del as u8)] = ins * 0x800000 + del * 0x800000;

        x
    }

    const WINDOW: f64 = 1.0;

    fn build_coef_table(v: &mut Vec<i32>, i: usize, scale: f64, pitch: f64, width: f64) {
        let span = 2 * ((0.5 * scale).ceil() as usize);
        let offset = (i as f64 - 8.0) / 16.0;

        /* FIXME */
        let center = pitch * pitch * offset + span as f64 / 2.0;

        for j in 0..=span {
            let dist = center - j as f64;
            let coef = sincx(clip(dist, width) / pitch, Self::WINDOW);

            /* FIXME */
            let coef = coef.max(0.0);
            v.push((0x01000000 as f64 * coef) as i32);
        }
    }

    pub(super) fn new(color: &UdonPalette, columns_per_pixel: f64) -> Scaler {
        let scale = columns_per_pixel.max(1.0);
        let pitch = columns_per_pixel / scale;
        let width = 0.5 * (scale - 1.0);
        let color_coef = 1.0 / (columns_per_pixel.log(3.5).max(1.0) + columns_per_pixel / 10.0);
        let alpha_coef = 1.0 / (columns_per_pixel.log(2.5).max(1.0) + columns_per_pixel / 5.0);

        let mut x = Scaler {
            columns_per_pixel,
            window: scale.ceil() + 1.0,
            offset: (scale.ceil() + 1.0) / 2.0,
            normalizer: Color {
                v: [
                    (0x01000000 as f64 * color_coef) as i32,
                    (0x01000000 as f64 * color_coef) as i32,
                    (0x01000000 as f64 * color_coef) as i32,
                    (0x01000000 as f64 * alpha_coef) as i32,
                    (0x01000000 as f64 * color_coef) as i32,
                    (0x01000000 as f64 * color_coef) as i32,
                    (0x01000000 as f64 * color_coef) as i32,
                    (0x01000000 as f64 * alpha_coef) as i32,
                ],
            },
            color: Self::build_color_table(color),
            table: Default::default(),
        };

        for i in 0..17 {
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

    pub fn new(color: &UdonPalette, columns_per_pixel: f64) -> Scaler {
        let scale = columns_per_pixel.max(1.0);
        let magnifier = scale / columns_per_pixel;

        let mut x = Scaler {
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

    pub(super) fn expected_size(&self, span: usize) -> usize {
        (span as f64 / self.columns_per_pixel) as usize + 2
    }

    pub(super) fn init(&self, offset_in_pixels: f64) -> (f64, usize) {
        let offset = (offset_in_pixels + 0.5) * self.columns_per_pixel;
        let margin = (/* offset + */self.offset) as usize;
        // println!("init, offset({}, {}), margin({})", offset, self.offset, margin);

        (offset, margin)
    }

    pub(super) fn scale(
        &self,
        dst: &mut Vec<[[u8; 4]; 2]>,
        src: &[u8],
        offset: f64,
    ) -> Option<(f64, usize)> {
        // println!("scale, offset({})", offset);
        for i in 0.. {
            /*
            offset := offset_in_columns
            */
            let base = offset + (i as f64 * self.columns_per_pixel);
            let range = Range::<usize> {
                start: base as usize,
                end: (base + self.window + 1.0).ceil() as usize,
            };
            // println!("base({}), range({:?})", base, range);

            if range.end > src.len() {
                return Some((base.fract(), range.start));
            }

            let table = &self.table[(base.fract() * 16.0) as usize];
            // println!("frac({}), {:?}", base.fract(), table);

            let mut a = Color::default();
            for (&coef, &column) in table.iter().zip(src[range].iter()) {
                // println!("col({}), color({:?}), coef({})", column, self.pick_color(column), coef);
                a += self.pick_color(column) * coef;
            }

            a *= self.normalizer;
            // println!("acc({:?})", <[[u8; 4]; 2]>::from(&a));
            dst.push(<[[u8; 4]; 2]>::from(&a));
        }

        None
    }
}

/*
/* Scaler and its impl */
struct Scaler {
    accum: Color,
    prev:  Color,
    color: [Color; 12],
    normalizer: u32
}

impl Scaler {

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

        Scaler {
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
