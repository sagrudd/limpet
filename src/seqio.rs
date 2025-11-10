use anyhow::{anyhow, Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A single contig / record
#[derive(Debug, Clone)]
pub struct Contig {
    /// First token of the header (accession)
    pub name: String,
    /// Full header text without the leading '>' or '@'
    pub header: String,
    /// Uppercase sequence
    pub seq: Vec<u8>,
}

enum Format {
    Fasta,
    Fastq,
}

fn is_gz(path: &Path) -> bool {
    path.extension().map(|e| e.eq_ignore_ascii_case("gz")).unwrap_or(false)
}

fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>> {
    let f = File::open(path)
        .with_context(|| format!("Failed to open input: {}", path.display()))?;
    if is_gz(path) {
        let gz = MultiGzDecoder::new(f);
        Ok(Box::new(BufReader::new(gz)))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

fn detect_format(path: &Path) -> Result<Format> {
    let mut rdr = open_maybe_gz(path)?;
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

fn parse_fasta<R: BufRead>(reader: R) -> Result<Vec<Contig>> {
    let mut contigs: Vec<Contig> = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_header: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line_res in reader.lines() {
        let line = line_res?;
        if line.is_empty() { continue; }
        if line.starts_with('>') {
            // flush previous
            if let Some(name) = current_name.take() {
                let header = current_header.take().unwrap_or_else(|| name.clone());
                contigs.push(Contig { name, header, seq: current_seq.clone() });
                current_seq.clear();
            }
            // capture full header and name token
            let header_full = line[1..].trim().to_string();
            let name = header_full.split_whitespace().next().unwrap_or(header_full.as_str()).to_string();
            current_name = Some(name);
            current_header = Some(header_full);
        } else {
            for b in line.bytes() {
                let c = b as char;
                if c.is_ascii_alphabetic() {
                    current_seq.push(c.to_ascii_uppercase() as u8);
                }
            }
        }
    }
    if let Some(name) = current_name.take() {
        let header = current_header.take().unwrap_or_else(|| name.clone());
        contigs.push(Contig { name, header, seq: current_seq });
    }
    if contigs.is_empty() {
        return Err(anyhow!("No sequences found in FASTA."));
    }
    Ok(contigs)
}

fn parse_fastq<R: BufRead>(mut reader: R) -> Result<Vec<Contig>> {
    // Robust FASTQ parser supporting wrapped sequence/quality.
    let mut contigs: Vec<Contig> = Vec::new();
    let mut line = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break; // EOF
        }
        let header = line.trim_end().to_string();
        if header.is_empty() { continue; }
        if !header.starts_with('@') {
            return Err(anyhow!("FASTQ record does not start with '@' header"));
        }
        let header_full = header[1..].trim().to_string();
        let name = header_full.split_whitespace().next().unwrap_or("").to_string();

        // Read sequence lines until '+' line
        let mut seq_buf: Vec<u8> = Vec::new();
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                return Err(anyhow!("Unexpected EOF while reading FASTQ sequence"));
            }
            let s = line.trim_end();
            if s.starts_with('+') { break; } // next stage
            for b in s.bytes() {
                let c = b as char;
                if c.is_ascii_alphabetic() {
                    seq_buf.push(c.to_ascii_uppercase() as u8);
                }
            }
        }

        // Read quality lines until we have as many quality chars as sequence length
        let mut qlen: usize = 0;
        while qlen < seq_buf.len() {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                return Err(anyhow!("Unexpected EOF while reading FASTQ quality"));
            }
            let s = line.trim_end();
            qlen += s.as_bytes().len();
        }

        contigs.push(Contig { name, header: header_full, seq: seq_buf });
    }

    if contigs.is_empty() {
        return Err(anyhow!("No sequences found in FASTQ."));
    }
    Ok(contigs)
}

/// Read a reference/input file that may be FASTA/FASTQ and optionally gzipped.
pub fn read_sequences<P: AsRef<Path>>(path: P) -> Result<Vec<Contig>> {
    let path_ref: &Path = path.as_ref();
    let fmt = detect_format(path_ref)?;
    let rdr = open_maybe_gz(path_ref)?;
    let contigs = match fmt {
        Format::Fasta => parse_fasta(rdr)?,
        Format::Fastq => parse_fastq(rdr)?,
    };
    Ok(contigs)
}

/// A small record for writing to FASTA
pub struct FastaRecord<'a> {
    pub header: String,
    pub seq: &'a [u8],
}

/// Write records to a FASTA file (wrapped to `line_width` chars).
pub fn write_fasta<P: AsRef<Path>>(records: &[FastaRecord<'_>], path: P, line_width: usize) -> Result<()> {
    use std::io::Write;
    let mut fh = std::fs::File::create(&path)
        .with_context(|| format!("Failed to create output FASTA: {}", path.as_ref().display()))?;
    let lw = if line_width == 0 { usize::MAX } else { line_width };

    for rec in records {
        writeln!(fh, ">{}", rec.header)?;
        let mut start = 0usize;
        while start < rec.seq.len() {
            let end = (start + lw).min(rec.seq.len());
            fh.write_all(&rec.seq[start..end])?;
            writeln!(fh)?;
            start = end;
        }
    }
    Ok(())
}
