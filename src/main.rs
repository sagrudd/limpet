//! # limpet — bioinformatics CLI
//!
//! `limpet` is a teaching-friendly and lab-friendly command‑line toolkit for DNA/RNA sequence
//! manipulation written in Rust. It emphasizes **clarity**, **safety**, and **reproducibility**.
//!
//! ## Subcommands (overview)
//! - **`seq_sample`** — sample *n* random genomic intervals from a reference and write FASTA.
//! - **`scramble`** — ingest many FASTA/FASTQ (plain or `.gz`), randomize global order, write one FASTA;
//!   headers are rewritten to `scramble_00001` with provenance retained.
//! - **`strip`** — reduce FASTA headers to the accession token only (first whitespace‑separated token).
//! - **`sample`** — randomly pick *n* raw records from the input (FASTA or FASTQ) and write them **unmodified**,
//!   preserving the file format; gzip if output ends with `.gz`.
//!
//! ## Installation
//! ```bash
//! cargo build --release
//! ```
//!
//! ## Reproducibility
//! Where appropriate, commands expose a `--seed` option for deterministic behavior.
//!
//! ## Safety & Scope
//! `limpet` is intended for **in silico** education and analysis. It does not interact with lab equipment
//! and does not attempt to evaluate biological risk. Prefer non‑pathogenic and openly available reference data,
//! and follow your institution's biosafety and data governance policies.

//!
//! Subcommands:
//! - `seq_sample`: sample random sequences from a reference FASTA

mod seqio;
mod seq_sample;
mod scramble;
mod sample;
mod strip;

use anyhow::Result;
use clap::{Parser, Subcommand};

/// limpet CLI
#[derive(Parser, Debug)]
#[command(name = "limpet")]
#[command(author, version, about = "Bioinformatics utilities in Rust", long_about = None)]
struct Cli {
    /// Subcommands
    #[command(subcommand)]
    command: Commands,
}

/// Top-level subcommands
#[derive(Subcommand, Debug)]
enum Commands {
    /// Sample random sequences from a reference FASTA
    SeqSample(seq_sample::SeqSampleArgs),
    /// Randomly sample N records from an input, keeping original format
    Sample(sample::SampleArgs),
    /// Strip FASTA headers to accession-only
    Strip(strip::StripArgs),
    /// Scramble sequences from multiple inputs into one FASTA
    Scramble(scramble::ScrambleArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::SeqSample(args) => seq_sample::run(args)?,
        Commands::Scramble(args) => scramble::run(args)?,
        Commands::Strip(args) => strip::run(args)?,
        Commands::Sample(args) => sample::run(args)?,
    }
    Ok(())
}
