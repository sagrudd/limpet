use crate::seqio::{read_sequences as read_fasta, write_fasta, Contig, FastaRecord};
use anyhow::{anyhow, Context, Result};
use clap::Args;
use rand::prelude::*;
use std::path::PathBuf;

/// Arguments for `limpet seq_sample`
#[derive(Args, Debug, Clone)]
pub struct SeqSampleArgs {
    /// Reference (FASTA/FASTA.GZ/FASTQ/FASTQ.GZ)
    #[arg(short = 'r', long = "reference", value_name = "INPUT")]
    pub reference: PathBuf,

    /// Number of sequences to sample
    #[arg(short = 'n', long = "n", value_name = "INT")]
    pub n: usize,

    /// Minimum length (inclusive)
    #[arg(long = "min", value_name = "INT")]
    pub min: usize,

    /// Maximum length (inclusive)
    #[arg(long = "max", value_name = "INT")]
    pub max: usize,

    /// Output FASTA path
    #[arg(short = 'o', long = "output", value_name = "FASTA")]
    pub output: PathBuf,

    /// Optional RNG seed for reproducibility
    #[arg(long = "seed", value_name = "INT")]
    pub seed: Option<u64>,
}

fn has_long_n_run(seq: &[u8], max_run: usize) -> bool {
    let mut run = 0usize;
    for &b in seq {
        if b == b'N' {
            run += 1;
            if run > max_run {
                return true;
            }
        } else {
            run = 0;
        }
    }
    false
}

pub fn run(args: SeqSampleArgs) -> Result<()> {
    if args.n == 0 {
        return Err(anyhow!("--n must be greater than 0"));
    }
    if args.min == 0 {
        return Err(anyhow!("--min must be greater than 0"));
    }
    if args.min > args.max {
        return Err(anyhow!("--min must be <= --max"));
    }

    let contigs = read_fasta(&args.reference)?;
    if !contigs.iter().any(|c| c.seq.len() >= args.min) {
        return Err(anyhow!(
            "No sequences in {} are at least {} bp long.",
            args.reference.display(),
            args.min
        ));
    }

    let mut rng: StdRng = match args.seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_rng(thread_rng()).context("Failed to initialize RNG")?,
    };

    let mut out_records: Vec<(String, Vec<u8>)> = Vec::with_capacity(args.n);

    while out_records.len() < args.n {
        let len = rng.gen_range(args.min..=args.max);

        // Compute weights = available start positions per contig
        let mut weights: Vec<u64> = Vec::with_capacity(contigs.len());
        let mut total: u128 = 0;
        for c in &contigs {
            if c.seq.len() >= len {
                let w = (c.seq.len() - len + 1) as u64;
                weights.push(w);
                total += w as u128;
            } else {
                weights.push(0);
            }
        }
        if total == 0 {
            // No contig can fit this length; try another length
            continue;
        }

        // Weighted choose contig
        let mut pick = rng.gen_range(0..total) as u128;
        let mut chosen_idx = 0usize;
        for (i, &w) in weights.iter().enumerate() {
            if w == 0 { continue; }
            if pick < w as u128 {
                chosen_idx = i;
                break;
            }
            pick -= w as u128;
        }

        let c: &Contig = &contigs[chosen_idx];
        let max_start = c.seq.len() - len;
        let start = rng.gen_range(0..=max_start);
        let end = start + len;
        let slice = &c.seq[start..end];

        // Reject sequences with long runs of 'N' (>2)
        if has_long_n_run(slice, 2) {
            continue;
        }

        // Build header: use 1-based inclusive coordinates for human-friendliness
        let header = format!(
            "seq{:06} src={} range={}..{} len={}",
            out_records.len() + 1,
            c.name,
            start + 1,
            end,
            len
        );
        out_records.push((header, slice.to_vec()));
    }

    // convert to records for writing
    let records: Vec<_> = out_records
        .iter()
        .map(|(h, s)| FastaRecord { header: h.clone(), seq: s.as_slice() })
        .collect();

    write_fasta(&records, &args.output, 80)?;
    eprintln!("Wrote {} sequences to {}", records.len(), args.output.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::{self, File};
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn rejects_long_n_runs() {
        let dir = tempdir().unwrap();
        let ref_path = dir.path().join("ref.fa");
        let out_path = dir.path().join("out.fa");

        let mut f = File::create(&ref_path).unwrap();
        writeln!(f, ">chrA").unwrap();
        writeln!(f, "ACGTNNNACGTACGT").unwrap(); // contains NNN
        writeln!(f, ">chrB").unwrap();
        writeln!(f, "ACGTACGTACGTACGT").unwrap();

        let args = SeqSampleArgs {
            reference: ref_path.clone(),
            n: 3,
            min: 4,
            max: 6,
            output: out_path.clone(),
            seed: Some(123),
        };
        run(args).unwrap();

        let out = fs::read_to_string(out_path).unwrap();
        // Should have samples, and none should contain NNN
        for line in out.lines() {
            if !line.starts_with('>') {
                assert!(!line.contains("NNN"));
            }
        }
    }
}
