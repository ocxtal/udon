
/**
@file ribbon.rs
@brief example for udon, creates simple pileup from command line.

@author Hajime Suzuki
@license MIT
*/

use std::fs::File;
use std::io::Write;
use std::ops::Range;
use std::path::PathBuf;
use bam::{ BamReader, Record, RecordReader };
use bam::record::tags::{ TagValue };
use image::ColorType;
use image::png::PngEncoder;
use udon::{ Udon, UdonScaler, UdonPalette, UdonUtils };


/* argument parsing */
extern crate structopt;
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
#[structopt(name = "ribbon", about = "udon example")]
struct RibbonOpt {

	/* plot range; (chrom name, spos, epos), half-inclusive */
	#[structopt(short, long, default_value = "")]
	reference: String,
	#[structopt(short, long, default_value = "0")]
	start: usize,
	#[structopt(short, long, default_value = "1000000000000")]
	end: usize,

	/* image dimensions */
	#[structopt(short, long, default_value = "640")]
	width: usize,
	#[structopt(short, long, default_value = "480")]
	height: usize,
	#[structopt(short, long, default_value = "10")]
	margin: usize,

	/* output PNG filename */
	#[structopt(short, long, default_value = "")]
	output: PathBuf,

	/* input BAM filename (single positional argument) */
	input: Option<PathBuf>
}


/* utilities on ranges */
trait RangeUtils where Self: Sized {
	fn has_overlap(&self, query: &Self) -> bool;
	fn clip(&self, query: &Self) -> Option<Self>;
	fn scale(&self, divisor: f64) -> (Self, f64);
}

impl RangeUtils for Range<usize> {
	fn has_overlap(&self, query: &Range<usize>) -> bool {
		if query.start > self.end { return false; }
		if self.start > query.end { return false; }
		return true;
	}

	fn clip(&self, query: &Range<usize>) -> Option<Range<usize>> {
		/* clip query range by self (window). returns local range in the window */
		if query.start > self.end { return None; }
		if self.start > query.end { return None; }

		Some(Range::<usize> {
			start: query.start.saturating_sub(self.start),
			end:   query.end.min(self.end).saturating_sub(self.start)
		})
	}

	fn scale(&self, divisor: f64) -> (Range<usize>, f64) {
		let start  = self.start as f64 / divisor;
		let offset = start.fract();

		let range = Range::<usize> {
			start: start as usize,
			end: (self.end as f64 / divisor).ceil() as usize
		};

		(range, offset)
	}
}


/* piling up alignment ribbons */
#[derive(Copy, Clone, Debug)]
struct Dimension {
	x: usize,
	y: usize
}

#[derive(Copy, Clone, Debug)]
struct Border {
	thickness: usize,
	color: [u8; 4]
}

#[derive(Copy, Clone, Debug)]
struct RibbonAttributes {
	height: usize,
	border: Border
}

#[derive(Copy, Clone, Debug)]
struct PileupParams {
	window: Dimension,
	margin: Dimension,
	border: Border,
	background: [u8; 4],
	fontsize: f64,
	ribbon: RibbonAttributes
}

impl Default for PileupParams {
	fn default() -> PileupParams {
		PileupParams {
			window: Dimension { x: 640, y: 480 },
			margin: Dimension { x: 5, y: 5 },
			border: Border {
				thickness: 1,
				color: [128, 128, 128, 0]
			},
			background: [255, 255, 255, 0],
			fontsize: 9.0,
			ribbon: RibbonAttributes {
				height: 5,
				border: Border {
					thickness: 1,
					color: [32, 32, 32, 0]
				}
			}
		}
	}
}

struct Pileup {
	buf: Vec<u8>,
	height: usize,
	params: PileupParams
}

impl Pileup {
	fn new(params: &PileupParams) -> Self {
		let mut this = Pileup {
			buf: Vec::<u8>::new(),
			height: 0,
			params: *params
		};

		for _ in 0 .. this.params.margin.y {
			this.append_margin_row();
		}
		this.append_border_row();
		this
	}

	fn push(&mut self, ribbon: &[[[u8; 4]; 2]], horizontal_offset: usize) -> Option<()> {

		let background = self.params.background;
		let border = self.params.border.color;

		let left_blank  = horizontal_offset;
		let right_blank = self.params.window.x.saturating_sub(ribbon.len() + horizontal_offset);
		let ribbon_len  = self.params.window.x - (left_blank + right_blank);

		// println!("{:?}, {:?}, {:?}", left_blank, right_blank, ribbon_len);

		for i in 0 .. self.params.ribbon.height {
			let idx = if i == self.params.ribbon.height / 2 { 1 } else { 0 };

			/* left margin */
			self.fill(&background, self.params.margin.x);
			self.fill(&border, self.params.border.thickness);
			self.fill(&background, left_blank);

			/* body */
			for &x in &ribbon[.. ribbon_len] {
				self.buf.write(&x[idx][.. 3]).ok()?;
			}

			/* right margin */
			self.fill(&background, right_blank);
			self.fill(&border, self.params.border.thickness);
			self.fill(&background, self.params.margin.x);
			self.height += 1;
		}
		Some(())
	}

	fn finalize(&mut self) -> Option<()> {

		/* fill blanks */
		for _ in self.height .. self.params.window.y {
			self.append_blank_row();
		}

		/* bottom border and margin */
		self.append_border_row()?;
		for _ in 0 .. self.params.margin.y {
			self.append_margin_row()?;
		}
		Some(())
	}

	fn render(&self, filename: &PathBuf) -> Option<()> {
		let output = File::create(filename).ok()?;
		let encoder = PngEncoder::new(output);

		// debug!("{:?}, {:?}, {:?}, {:?}", self.total_width(), self.total_height(), self.height, self.buf.len());

		encoder.encode(&self.buf[.. 3 * self.total_width() * self.total_height()],
			self.total_width() as u32,
			self.total_height() as u32,
			ColorType::Rgb8
		).ok()?;
		Some(())
	}


	/* internal utilities */
	fn total_width(&self) -> usize {
		self.params.window.x + 2 * (self.params.margin.x + self.params.border.thickness)
	}

	fn total_height(&self) -> usize {
		self.params.window.y + 2 * (self.params.margin.y + self.params.border.thickness)
	}

	fn fill(&mut self, color: &[u8; 4], size: usize) -> Option<usize> {
		for _ in 0 .. size {
			self.buf.write(&color[.. 3]).ok()?;
		}
		Some(size)
	}

	fn append_margin_row(&mut self) -> Option<usize> {
		let len = self.total_width();
		let background = self.params.background;
		self.fill(&background, len);
		Some(len)
	}

	fn append_blank_row(&mut self) -> Option<usize> {
		let background = self.params.background;
		let border = self.params.border.color;

		self.fill(&background, self.params.margin.x);
		self.fill(&border, self.params.border.thickness);
		self.fill(&background, self.params.window.x);
		self.fill(&border, self.params.border.thickness);
		self.fill(&background, self.params.margin.x);
		Some(self.total_width())
	}

	fn append_border_row(&mut self) -> Option<usize> {
		let background = self.params.background;
		let border = self.params.border.color;

		self.fill(&background, self.params.margin.x);
		self.fill(&border, self.params.window.x + 2 * self.params.border.thickness);
		self.fill(&background, self.params.margin.x);
		Some(self.total_width())
	}

}


fn main() {

	env_logger::init();

	/* parse args */
	let opt = RibbonOpt::from_args();
	let filename = opt.input.unwrap();

	/* then open file */
	let mut reader = BamReader::from_path(&filename, 2).expect(
		format!("Failed to open file `{:?}'. Please check the file exists.", &filename).as_str()
	);


	/* get reference sequence id, then extract valid range */
	let id = reader.header().reference_id(&opt.reference).unwrap_or(0);
	let window = Range::<usize> {
		start: opt.start,
		end:   opt.end.min(reader.header().reference_len(id).unwrap() as usize)
	};


	/* prepare ribbon scaler and color */
	let columns_per_pixel = window.len() as f64 / opt.width as f64;
	let scaler = UdonScaler::new(&UdonPalette::default(), columns_per_pixel);
	let base_color: [[[u8; 4]; 2]; 2] = [
		[[255, 202, 191, 255], [255, 202, 191, 255]],
		[[191, 228, 255, 255], [191, 228, 255, 255]]
	];


	/* everything successful; create PNG buffer */
	let mut pileup = Pileup::new(
		&PileupParams {
			window: Dimension { x: opt.width, y: opt.height },
			.. PileupParams::default()
		}
	);


	/* for each alignment... */
	let mut record = Record::new();
	while let Ok(true) = reader.read_into(&mut record) {
		if !record.flag().is_mapped() { continue; }

		/* construct indexed ribbon (udon) */
		let udon = Udon::build(
			record.cigar().raw(),
			record.sequence().raw(),
			if let Some(TagValue::String(s, _)) = record.tags().get(b"MD") { s } else {
				panic!("Each BAM record must have MD string. Inspect `samtools calmd` for restoring missing MD strings.")
			}
		).expect(&format!("Failed to create udon index. Would be a bug. ({:?})", &record.name()));

		/* compose span, skip if out of the window */
		let range = Range::<usize> {
			start: record.start() as usize,
			end:   record.start() as usize + udon.reference_span()
		};
		if !window.has_overlap(&range) { continue; }

		/* compute local ranges */
		let udon_range   = range.clip(&window).unwrap();
		let window_range = window.clip(&range).unwrap();
		if 3 * window_range.len() < window.len() { continue; }

		let (window_range, offset_in_pixel) = window_range.scale(columns_per_pixel);

		/* slice ribbon scaled */
		let mut ribbon = udon.decode_scaled(
			&udon_range,
			offset_in_pixel,
			&scaler
		).expect(&format!("Failed to decode udon ribbon. Would be a bug. ({:?})", &record.name()));

		/* put forward / reverse color then do gamma correction */
		ribbon.append_on_basecolor(&base_color[record.flag().is_reverse_strand() as usize]).correct_gamma();

		/* then pileup; break when buffer is full */
		pileup.push(&ribbon, window_range.start);
		// println!("{:?}, {:?}, {}", udon_range, window_range, offset_in_pixel);
	}

	/* done!!! */
	pileup.finalize();
	pileup.render(&opt.output).expect(
		format!("failed to dump image to `{:?}'", &opt.output).as_str()
	);
}

