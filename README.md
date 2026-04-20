# maxbin-rs

A Rust reimplementation of [MaxBin2](https://sourceforge.net/projects/maxbin2/)
by [Yu-Wei Wu, Blake A. Simmons, and Steven W. Singer](https://doi.org/10.1093/bioinformatics/btv638) —
a metagenome binning tool that clusters assembled contigs into individual
genomes using tetranucleotide frequency and abundance information.

> Wu Y-W, Simmons BA, Singer SW. **MaxBin 2.0: an automated binning algorithm to
> recover genomes from multiple metagenomic datasets.** *Bioinformatics.*
> 2016;32(4):605–607. doi:[10.1093/bioinformatics/btv638](https://doi.org/10.1093/bioinformatics/btv638)

This project reimplements their algorithm while addressing long-standing
packaging and reliability issues. Not yet tested inside nf-core/mag or
other pipeline managers — if you try it, please report issues.

If referring specifically to this reimplementation, please also cite:

> Kroneman W. **maxbin-rs: a Rust reimplementation of MaxBin2.** 2026.
> https://github.com/werner291/maxbin-rs

Always cite the original MaxBin2 first — the algorithm is theirs.

## Quick start

```bash
# Nix (builds everything from source, pinned deps)
nix run github:werner291/maxbin-rs -- \
  --contig contigs.fa.gz --reads reads.fastq.gz --out my_bins/

# Docker (for pipeline integration — pin to a version tag)
docker run --rm -v /path/to/data:/data \
  ghcr.io/werner291/maxbin-rs:v0.3.0 \
  --contig /data/contigs.fa.gz --reads /data/reads.fastq.gz --out /data/bins
```

Both bundle HMMER, Bowtie2, and FragGeneScan. A bioconda package is
planned. You can also `cargo build --release` if you have these tools
on PATH.

For a drop-in replacement that matches the original MaxBin2 CLI
exactly, use a v0.1.x release (e.g. `ghcr.io/werner291/maxbin-rs:v0.1.3`).

### Inputs

- **`--contig`** — assembled contigs (FASTA, optionally gzipped).
- **`--abund`** — per-contig read depth (tab-separated: contig name, depth).
  Pass once per sample:
  `--abund sample1.txt --abund sample2.txt --abund sample3.txt`
- **`--reads`** — raw reads (FASTQ, optionally gzipped). Abundance is
  computed via Bowtie2 mapping. Pass once per sample:
  `--reads sample1.fq.gz --reads sample2.fq.gz`

Use `--abund` or `--reads`, not both. Default minimum contig length is
1000 bp. Run `maxbin-rs --help` for all options.

### Output

```
my_bins/
├── 001.fasta    # bin 1
├── 002.fasta    # bin 2
├── ...
├── summary      # per-bin completeness, genome size, GC
├── noclass      # contigs not assigned to any bin
├── marker       # marker gene counts per bin
└── log          # tool version and run metadata
```

### EM subcommand

The EM algorithm is the core of MaxBin. If you already have abundance
data and seed contigs (e.g. from your own marker gene pipeline), you
can run it directly without any external tools:

```bash
maxbin-rs em --contig filtered.fa --abund depth.txt --seed seeds.txt --out result
```

The default mode (`maxbin-rs pipeline`) wraps this with gene calling
(FragGeneScanRs), marker detection (HMMER), and read mapping (Bowtie2) to
produce seeds and abundance automatically. You can also supply your own
intermediates via `--faa` (protein FASTA) or `--hmmout` (HMMER output)
to substitute individual stages. See `maxbin-rs --help` for all options.

## Migrating from MaxBin2

```bash
# before
run_MaxBin.pl -contig contigs.fa -reads reads.fq -out results/my_sample
# produces: results/my_sample.001.fasta, results/my_sample.summary, ...

# after
maxbin-rs --contig contigs.fa --reads reads.fq --out results/my_sample/
# produces: results/my_sample/001.fasta, results/my_sample/summary, ...
```

What to update in your scripts:

- Binary: `maxbin-rs` instead of `run_MaxBin.pl`
- Flags use double dashes: `--contig`, `--reads`, `--out`
- `--out` is a directory, not a filename prefix
- Output files: `001.fasta` instead of `prefix.001.fasta`
- Version string: `maxbin-rs X.Y.Z` instead of `MaxBin 2.2.7`
- Progress goes to stderr, not stdout
- Exit code `1` on error (original uses `-1` / 255)
- Not yet produced: per-bin marker tarball, multi-sample `.abundance`
  file (both `optional: true` in nf-core/mag)

If you can't change your scripts at all, **v0.1.x** uses the original
single-dash flags and prefix-based `--out` — it's a drop-in replacement
for `run_MaxBin.pl`. Bug fixes and performance improvements are
backported as long as they don't change output.

See [VERSIONING.md](VERSIONING.md) for the full version policy.

## Verification

The v0.1.x series was verified to produce identical output to MaxBin2
2.2.7. Later versions build on that as a baseline, adding bug fixes,
performance improvements, and interface changes. The evidence is
described in detail in [PAPER.md](PAPER.md). In short:

- **Unit/property tests** compare each Rust function against the
  original C++ via FFI on randomized inputs (`cargo nextest run`).
- **End-to-end tests** run the full pipeline (including recursive
  binning) on CAMI I High: all bins byte-identical.
- The EM core is approximately **8x faster** than the original C++
  (see PAPER.md for methodology).

These results have not been independently verified. If you do run your
own comparison, we'd welcome the feedback.

