//! Sequence scrambler (`scramble`).
//!
//! Reads multiple **FASTA/FASTQ** files (plain or `.gz`), **loads all sequences into memory**, shuffles the global order,
//! and writes a single FASTA. Each output header begins with a new sequential accession (`scramble_00001`), followed by
//! `src=<original_accession>` and `file=<source_file>`, and finally the original header text.
//!
//! ### Memory considerations
//! This command is explicitly **in‑memory**. Handling ≳1 Gbp of sequence is reasonable on a modern laptop, but for very
//! large datasets consider running in batches.
//!
//! ### Example
//! ```text
//! limpet scramble input1.fa input2.fq.gz -o scrambled.fa --seed 42
//! ```

use crate::seqio::{read_sequences, write_fasta, FastaRecord};
use anyhow::{anyhow, Context, Result};
use clap::Args;
use rand::prelude::*;
use std::path::PathBuf;

/// Scramble: read multiple inputs (FASTA/FASTQ and .gz variants), shuffle all records, and write a single FASTA.
#[derive(Args, Debug, Clone)]
pub struct ScrambleArgs {
    /// One or more input files (FASTA/FASTQ/FASTA.GZ/FASTQ.GZ)
    #[arg(value_name = "INPUT", required = true)]
    pub inputs: Vec<PathBuf>,

    /// Output FASTA path
    #[arg(short = 'o', long = "output", value_name = "FASTA")]
    pub output: PathBuf,

    /// Optional RNG seed for reproducibility
    #[arg(long = "seed", value_name = "INT")]
    pub seed: Option<u64>,
}

/// Execute the `scramble` subcommand.
/// Loads all inputs, shuffles records, rewrites headers, and writes FASTA output.
pub fn run(args: ScrambleArgs) -> Result<()> {
    if args.inputs.is_empty() {
        return Err(anyhow!("Provide at least one input file."));
    }
    // Load all sequences (+ provenance) into memory
    let mut all: Vec<(String, String, Vec<u8>, String)> = Vec::new(); // (orig_name, header_full, seq, file_base)
    for path in &args.inputs {
        let recs = read_sequences(path)
            .with_context(|| format!("Failed to read input {}", path.display()))?;
        let file_base = path.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_else(|| path.display().to_string());
        for c in recs {
            all.push((c.name, c.header, c.seq, file_base.clone()));
        }
    }
    if all.is_empty() {
        return Err(anyhow!("No sequences found in provided inputs."));
    }

    // Shuffle globally
    let mut rng: StdRng = match args.seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_rng(thread_rng()).context("Failed to initialize RNG")?,
    };
    all.shuffle(&mut rng);

    // Build new headers: scramble_00001..N + source provenance + original header
    let mut out: Vec<FastaRecord<'_>> = Vec::with_capacity(all.len());
    let mut owned_headers: Vec<String> = Vec::with_capacity(all.len());
    let mut owned_seqs: Vec<Vec<u8>> = Vec::with_capacity(all.len());

    for (i, (orig_name, header_full, seq, file_base)) in all.into_iter().enumerate() {
        let new_name = format!("scramble_{:05}", i + 1);
        let hdr = format!("{} src={} file={} | {}", new_name, orig_name, file_base, header_full);
        owned_headers.push(hdr);
        owned_seqs.push(seq);
    }

    for i in 0..owned_headers.len() {
        let hdr = owned_headers[i].clone();
        let seq_ref = owned_seqs[i].as_slice();
        out.push(FastaRecord { header: hdr, seq: seq_ref });
    }

    write_fasta(&out, &args.output, 80)?;
    eprintln!("Wrote {} sequences to {}", out.len(), args.output.display());
    Ok(())
}
