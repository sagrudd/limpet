# limpet

**A small, sharp, teaching‑first bioinformatics CLI in Rust.**

`limpet` provides a handful of carefully designed commands for working with DNA/RNA sequences,
optimized for clarity, correctness, and reproducibility. It is ideal for workshops, high‑school
labs, and university teaching, while remaining useful for everyday research tasks.

---

## Quick install

```bash
# Build from source
cargo build --release
# Run the binary
./target/release/limpet -h
```

---

## Subcommands at a glance

- `seq_sample` — sample *n* random genomic intervals from a reference FASTA. Output FASTA.
- `scramble` — load many FASTA/FASTQ (plain or `.gz`), shuffle *all* sequences into a single FASTA with provenance‑rich headers.
- `strip` — reduce FASTA headers to accession tokens only.
- `sample` — randomly pick *n* raw records from a FASTA/FASTQ (optionally `.gz`) and write them unmodified; output format matches input.

Each command supports `-h/--help` for usage details.

---

## `seq_sample` — sample random genomic intervals

**Goal:** create synthetic fragments from a reference genome.

```bash
limpet seq_sample   --reference genome.fa.gz   --n 1000   --min 100   --max 300   --output fragments.fa   --seed 1
```

**How it works:** length `L` is uniform in `[min, max]`; for that `L`, contigs are weighted by `(len(contig) - L + 1)`.
This approximates a *uniform* distribution over reference coordinates. Candidates with runs of `N` longer than 2 are rejected.

---

## `scramble` — merge & randomize multiple inputs

**Goal:** create a shuffled corpus of sequences from many inputs (FASTA/FASTQ, gzipped or not).

```bash
limpet scramble   genomeA.fa genomeB.fq.gz genomeC.fa.gz   -o scrambled.fa   --seed 42
```

**Provenance‑rich headers:** each output record begins with `scramble_00001` (sequential),
then `src=<original_accession> file=<source_file> | <full_original_header>`.

**Memory note:** `scramble` is intentionally in‑memory. 1+ Gbp is fine on modern laptops; scale accordingly.

---

## `strip` — keep only accessions

**Goal:** normalize FASTA headers.

```bash
limpet strip   --input genome.fa.gz   --output genome_accessions.fa
```

Resulting headers are the first token (accession) from each original header.

---

## `sample` — reservoir sample FASTA/FASTQ

**Goal:** take a small, representative subset of a very large file without loading it all.

```bash
limpet sample   --input reads_R1.fastq.gz   --n 10000   --output subset_R1.fastq.gz   --seed 123
```

Uses **reservoir sampling**: keeps only `n` records in memory, scans the input once, and preserves records *exactly*,
including FASTQ qualities and line wrapping.

---

## Building “exotic” metagenomic‑type datasets for the classroom

Use `limpet` to craft controlled mixtures of sequences from discrete, safe genomes to simulate real‑world metagenomes:

### 1) Start from complete genomes
Pick several non‑pathogenic reference genomes (e.g., model bacteria or phage) from public databases.
Download their FASTA references.

### 2) Fragment each genome
Use `seq_sample` to fragment each genome at realistic lengths:

```bash
limpet seq_sample -r ecoli.fa -n 20000 --min 150 --max 300 -o ecoli_frag.fa --seed 11
limpet seq_sample -r bacillus.fa -n 12000 --min 150 --max 300 -o bacillus_frag.fa --seed 12
```

### 3) Combine and randomize
Scramble everything into a single “metagenomic” pool:

```bash
limpet scramble ecoli_frag.fa bacillus_frag.fa -o classroom_metagenome.fa --seed 99
```

### 4) Optional: down‑sample reads
If you already have read data (FASTQ), use `sample` to draw subsets per file to simulate varying species abundances:

```bash
limpet sample --input speciesA.fastq.gz --n 5000 --output speciesA_sub.fastq.gz --seed 7
```

### 5) Analyze with your favorite tools
The resulting FASTA/FASTQ can be used to practice mapping, assembly, k‑mer profiling, or taxonomic classification.
Because you defined the mixture, you know the ground truth for class discussions.

> **Teaching tip:** assign different mixtures to student groups and ask them to infer community composition,
> compare pipelines, and discuss sources of bias (GC content, repeats, N regions, read length).

---

## Formats & conventions

- **FASTA** headers written by limpet never include spaces before the accession token; additional metadata follows as `key=value` pairs or free text.
- **Coordinates** reported by `seq_sample` are **1‑based inclusive** (`range=start..end`).
- **Gzip**: input is auto‑detected; output is gzipped only when the output filename ends with `.gz`.

---

## Performance & resource planning

- `seq_sample` loads the reference in memory (O(genome_size)). Use gzip to reduce disk I/O.
- `scramble` loads all sequences in memory. Expect roughly ~1 byte per base plus overhead for headers and vectors.
- `sample` streams the input; memory is O(n) where `n` is your requested sample size.

---

## Reproducibility

Most commands expose a `--seed`. Store the seed, command line, and limpet version alongside your results.

---

## Biosafety & ethics

Use only **non‑pathogenic** genomes from reputable sources for classroom activities and follow your institution’s policies.
`limpet` operates on public sequence files and does not generate biological materials.

---

## License

Mozilla Public License 2.0 (MPL-2.0)

