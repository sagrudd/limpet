//! Streaming record sampler (`sample`).
//!
//! Randomly pick *n* **raw records** from an input **FASTA** or **FASTQ** (optionally `.gz`) using
//! **reservoir sampling**. Output records are written **unmodified** in the **same logical format**.
//! If the output filename ends with `.gz`, the output is gzipped.
//!
//! ### Why reservoir sampling?
//! Reservoir sampling uses *O(n)* memory (for your requested sample size) and *O(1)* extra work per record,
//! enabling fair sampling without a prior pass to count records.
//!
//! ### Example
//! ```text
//! limpet sample --input reads.fastq.gz --n 10000 --output subset.fastq.gz --seed 123
//! ```

use anyhow::{anyhow, Context, Result};
use clap::Args;
use flate2::write::GzEncoder;
use flate2::Compression;
use rand::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Args, Debug, Clone)]
pub struct SampleArgs {
    /// Input file (FASTA/FASTQ; optionally .gz)
    #[arg(short = 'i', long = "input", value_name = "INPUT")]
    pub input: PathBuf,

    /// Number of sequences to sample
    #[arg(short = 'n', long = "n", value_name = "INT")]
    pub n: usize,

    /// Output file; format will match the input (FASTA vs FASTQ). If the name ends with .gz, output will be gzipped.
    #[arg(short = 'o', long = "output", value_name = "OUTPUT")]
    pub output: PathBuf,

    /// Optional RNG seed for reproducibility
    #[arg(long = "seed", value_name = "INT")]
    pub seed: Option<u64>,
}

enum Format {
    Fasta,
    Fastq,
}

fn is_gz(path: &Path) -> bool {
    path.extension().map(|e| e.eq_ignore_ascii_case("gz")).unwrap_or(false)
}

fn open_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    use flate2::read::MultiGzDecoder;
    let f = File::open(path).with_context(|| format!("Failed to open {}", path.display()))?;
    if is_gz(path) {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(f))))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

fn detect_format(path: &Path) -> Result<Format> {
    let mut rdr = open_reader(path)?;
    let mut line = String::new();
    loop {
        line.clear();
        let bytes = rdr.read_line(&mut line)?;
        if bytes == 0 { break; }
        let s = line.trim_end();
        if s.is_empty() { continue; }
        let first = s.as_bytes()[0];
        return match first {
            b'>' => Ok(Format::Fasta),
            b'@' => Ok(Format::Fastq),
            _ => Err(anyhow!("Unrecognized input format: first non-empty line did not start with '>' (FASTA) or '@' (FASTQ)")),
        };
    }
    Err(anyhow!("Input appears empty: {}", path.display()))
}

/// Read next FASTA record as a raw string (including trailing newline). Returns None on EOF.
fn read_fasta_record<R: BufRead>(rdr: &mut R, first_header: Option<String>) -> Result<Option<String>> {
    let header = match first_header {
        Some(h) => h,
        None => {
            let mut line = String::new();
            loop {
                line.clear();
                if rdr.read_line(&mut line)? == 0 { return Ok(None); }
                if line.starts_with('>') { break; }
            }
            line
        }
    };

    let mut raw = String::new();
    raw.push_str(&header);

    let mut line = String::new();
    loop {
        line.clear();
        let bytes = rdr.read_line(&mut line)?;
        if bytes == 0 { break; }
        if line.starts_with('>') {
            // push header back by returning it as first_header in the next call
            return Ok(Some(raw));
        } else {
            raw.push_str(&line);
        }
    }
    Ok(Some(raw))
}

/// Read next FASTQ record as a raw string (including trailing newline). Returns None on EOF.
fn read_fastq_record<R: BufRead>(rdr: &mut R) -> Result<Option<String>> {
    let mut header = String::new();
    loop {
        header.clear();
        if rdr.read_line(&mut header)? == 0 { return Ok(None); }
        if header.trim_end().is_empty() { continue; }
        if !header.starts_with('@') {
            return Err(anyhow!("FASTQ header line must start with '@'"));
        }
        break;
    }

    let mut raw = String::new();
    raw.push_str(&header);

    // sequence lines until '+'
    let mut seq_len: usize = 0;
    let mut line = String::new();
    loop {
        line.clear();
        if rdr.read_line(&mut line)? == 0 {
            return Err(anyhow!("Unexpected EOF while reading FASTQ sequence"));
        }
        let s = line.trim_end();
        if s.starts_with('+') {
            raw.push_str(&line); // include '+' line
            break;
        } else {
            seq_len += s.as_bytes().len();
            raw.push_str(&line);
        }
    }
    // quality lines until length >= seq_len
    let mut qlen: usize = 0;
    loop {
        line.clear();
        if rdr.read_line(&mut line)? == 0 {
            return Err(anyhow!("Unexpected EOF while reading FASTQ qualities"));
        }
        let s = line.trim_end();
        qlen += s.as_bytes().len();
        raw.push_str(&line);
        if qlen >= seq_len { break; }
    }
    Ok(Some(raw))
}

/// Execute the `sample` subcommand.
/// Streams input, performs reservoir sampling, and writes output in matching format.
pub fn run(args: SampleArgs) -> Result<()> {
    if args.n == 0 {
        return Err(anyhow!("--n must be greater than 0"));
    }

    let input_fmt = detect_format(&args.input)?;
    let mut rdr = open_reader(&args.input)?;

    // Reservoir sample of raw records
    let mut rng: StdRng = match args.seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_rng(thread_rng()).context("Failed to initialize RNG")?,
    };
    let mut reservoir: Vec<String> = Vec::with_capacity(args.n);
    let mut seen: usize = 0;

    match input_fmt {
        Format::Fastq => {
            while let Some(rec) = read_fastq_record(&mut rdr)? {
                seen += 1;
                if reservoir.len() < args.n {
                    reservoir.push(rec);
                } else {
                    let j = rng.gen_range(0..seen);
                    if j < args.n {
                        reservoir[j] = rec;
                    }
                }
            }
        }
        Format::Fasta => {
            // FASTA: we need to manage lookahead of the next '>' header
            let mut buf = String::new();
            // Read first header
            loop {
                buf.clear();
                if rdr.read_line(&mut buf)? == 0 { break; }
                if buf.starts_with('>') { break; }
            }
            if !buf.is_empty() {
                let mut next_header = Some(buf.clone());
                loop {
                    if let Some(rec) = read_fasta_record(&mut rdr, next_header.take())? {
                        seen += 1;
                        if reservoir.len() < args.n {
                            reservoir.push(rec);
                        } else {
                            let j = rng.gen_range(0..seen);
                            if j < args.n {
                                reservoir[j] = rec;
                            }
                        }
                        // Now we need to peek if next char is '>' — read_fasta_record stops before reading next header
                        // We'll attempt to read the next header line here
                        let mut h = String::new();
                        loop {
                            h.clear();
                            let bytes = rdr.read_line(&mut h)?;
                            if bytes == 0 { break; }
                            if h.starts_with('>') { next_header = Some(h); break; }
                        }
                        if next_header.is_none() { break; }
                    } else {
                        break;
                    }
                }
            }
        }
    }

    if reservoir.is_empty() {
        return Err(anyhow!("No records found in {}", args.input.display()));
    }

    // Shuffle selected to randomize order
    reservoir.shuffle(&mut rng);

    // Open output writer, gz if .gz
    let f = File::create(&args.output).with_context(|| format!("Failed to create {}", args.output.display()))?;
    let gz = args.output.extension().map(|e| e.eq_ignore_ascii_case("gz")).unwrap_or(false);
    if gz {
        let enc = GzEncoder::new(f, Compression::default());
        let mut w = BufWriter::new(enc);
        for rec in &reservoir {
            w.write_all(rec.as_bytes())?;
        }
        w.flush()?;
    } else {
        let mut w = BufWriter::new(f);
        for rec in &reservoir {
            w.write_all(rec.as_bytes())?;
        }
        w.flush()?;
    }

    eprintln!("Sampled {} records (from ~{} seen) into {}", reservoir.len(), seen, args.output.display());
    Ok(())
}
