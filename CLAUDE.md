# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

maxbin-rs is a Rust reimplementation of MaxBin2, a metagenome binning tool. It clusters assembled contigs into individual genomes using an EM algorithm over tetranucleotide frequency and abundance data. Status: early planning — no Rust code exists yet.

The project follows the [rewrites.bio](https://rewrites.bio) principles: faithful algorithm preservation while fixing packaging, performance, and reliability problems. See TODO.md for the full complaint catalog, roadmap, and algorithm summary.

## Development Environment

Uses a Nix flake devshell (enter via `direnv allow` or `nix develop`). The shell provides:
- Rust stable toolchain with rust-analyzer
- cargo-nextest (preferred test runner)
- hmmer and bowtie2 (reference tool dependencies)

## Build & Test Commands

```bash
cargo build              # build
cargo nextest run        # run all tests (preferred over cargo test)
cargo nextest run <name> # run a single test by name
cargo clippy             # lint
cargo fmt --check        # check formatting
```

## Algorithm Overview

1. Detect single-copy marker genes (107 bacterial / 40 universal) via HMMER + gene caller
2. Seed initial bins from marker gene clusters
3. Iteratively assign contigs to bins based on tetranucleotide frequency + abundance probability (EM)
4. Refine until convergence
5. Contigs below confidence threshold remain unclassified

Inputs: assembled contigs (FASTA) + per-sample read depth (abundance files or BAM).
Outputs: one FASTA per bin + summary statistics.

## Equivalence Testing

Correctness is verified through **equivalence tests**: run a component from the original MaxBin2 and the corresponding Rust component on the same input, then assert identical output — including reproducing known bugs. The original MaxBin2 is available in the devshell via `run_MaxBin.pl`.

- Test at the component level (EM core, FASTA parsing, abundance loading, etc.), not just end-to-end.
- When an equivalence test reveals a bug in the original, **document it** (in TODO.md's complaint catalog and/or inline) but still reproduce the buggy behavior in the default code path.
- Bug fixes come later as intentional, separately-tested changes — never mixed with the rewrite.
- The Perl orchestration script (`run_MaxBin.pl`) does significant preprocessing
  (contig filtering, header munging, abundance generation, output assembly) that
  must also be equivalence-tested — not just the C++ core.

## Key Design Goals

- Single static binary — no Perl, no `setting` file, no runtime path resolution
- Replace FragGeneScan with Prodigal, adopt GTDB markers (per maxbin2_custom fork approach)
- Drop-in replacement for nf-core/mag pipeline
- Fix known correctness bugs: `prob_threshold` default mismatch, duplicate contigs across bins, compressed input detection
