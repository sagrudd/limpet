//! Header reducer (`strip`).
//!
//! Accepts a **single FASTA** (plain or `.gz`) and writes a FASTA where each header is reduced to just
//! the **accession** (first whitespace‑separated token). Sequence content is unchanged.
//!
//! ### Example
//! ```text
//! limpet strip --input reference.fa.gz --output reference_accessions.fa
//! ```

use crate::seqio::{read_sequences, write_fasta, FastaRecord};
use anyhow::{anyhow, Context, Result};
use clap::Args;
use std::path::PathBuf;

/// Strip FASTA headers down to just the accession (first token), preserving sequences.
#[derive(Args, Debug, Clone)]
pub struct StripArgs {
    /// Input FASTA (optionally gzipped). FASTQ is not intended for this command.
    #[arg(short = 'i', long = "input", value_name = "FASTA")]
    pub input: PathBuf,

    /// Output FASTA path
    #[arg(short = 'o', long = "output", value_name = "FASTA")]
    pub output: PathBuf,
}

/// Execute the `strip` subcommand.
/// Reduces headers to accession tokens and writes FASTA.
pub fn run(args: StripArgs) -> Result<()> {
    // Load sequences (FASTA or FASTA.GZ). read_sequences will also parse FASTQ, but this
    // command is intended for FASTA; we simply use the accession token `name` for headers.
    let records = read_sequences(&args.input)
        .with_context(|| format!("Failed to read {}", args.input.display()))?;

    if records.is_empty() {
        return Err(anyhow!("No sequences found in {}", args.input.display()));
    }

    let out: Vec<FastaRecord<'_>> = records
        .iter()
        .map(|c| FastaRecord { header: c.name.clone(), seq: c.seq.as_slice() })
        .collect();

    write_fasta(&out, &args.output, 80)?;
    eprintln!("Wrote {} sequences to {}", out.len(), args.output.display());
    Ok(())
}
