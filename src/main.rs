//! limpet: bioinformatics CLI
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
