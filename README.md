# limpet

Bioinformatics utilities in Rust. First tool: `seq_sample` samples random sequences from a reference in FASTA/FASTA.GZ/FASTQ/FASTQ.GZ. Output is always FASTA.

## Install & Build

```bash
cd limpet
cargo build --release
```

## Usage

```bash
limpet seq_sample   --reference reference.fa   --n 100   --min 50   --max 120   --output samples.fa   --seed 12345   # optional for reproducibility
```

**Notes**

- The tool treats multi-FASTA inputs correctly.
- Length is chosen uniformly in [min, max]. For a chosen length `L`, contigs are selected *weighted by the number of valid start positions* (`len(contig) - L + 1`), so positions across the genome are close to uniformly likely.
- Output headers include the source contig and 1-based inclusive coordinates, e.g.:
  `>seq000123 src=chr1 range=10001..10120 len=120`
- Sequences are wrapped at 80 columns.
- If no contig is at least `min` bp long, the command exits with an error.

## Example

```bash
limpet seq_sample -r ref.fa -n 10 --min 100 --max 250 -o sampled.fa
```

## License

MIT OR Apache-2.0


### Accepted input formats
- FASTA (`.fa`, `.fasta`) and `.fa.gz` / `.fasta.gz`
- FASTQ (`.fq`, `.fastq`) and `.fq.gz` / `.fastq.gz`

> Format is auto-detected by the first non-empty line: `>` for FASTA, `@` for FASTQ (extension not required).

###  Filtering of N-runs
Candidates containing `NNN` (or longer) are rejected before writing to the output FASTA.

> The `--tui` dashboard opens immediately and updates ~60fps while sampling; it closes when the run completes.


**Tip:** Add `--tui-wait` to keep the dashboard open at the end until you press a key.


## `scramble` — merge & randomize multiple inputs
Reads one or more inputs (FASTA/FASTQ and `.gz` variants), loads all sequences into memory, **shuffles the global order**, ensures **unique accessions** (appends `-dupN` if duplicates), and writes a single FASTA.

```bash
limpet scramble \  inputs1.fa inputs2.fq.gz inputs3.fa.gz \  -o scrambled.fa \  --seed 1234   # optional
```


**scramble header format:** first token is a new sequential accession `scramble_00001`, then provenance fields `src=<original_accession>` and `file=<source_file>`, followed by the original header text. Example:

```
>scramble_00001 src=seqA file=reads_1.fq.gz | seqA some description here
```


## `strip` — keep only accessions
Accepts a single FASTA (optionally `.gz`) and writes a FASTA with headers reduced to the accession (first token) only.

```bash
limpet strip \  --input reference.fa.gz \  --output reference_accessions.fa
```


## `sample` — reservoir sample FASTA/FASTQ
Randomly picks `n` sequences from the input file and writes them **unmodified** in the **same format** (FASTA vs FASTQ). If the output filename ends with `.gz`, the output is gzipped.

```bash
# Sample 10k reads from a gzipped FASTQ, keep FASTQ formatting + gzip
limpet sample \  --input reads_R1.fastq.gz \  --n 10000 \  --output subset_R1.fastq.gz \  --seed 42
```
