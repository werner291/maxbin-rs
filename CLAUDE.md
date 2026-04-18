# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

maxbin-rs is a Rust reimplementation of MaxBin2, a metagenome binning tool. It clusters assembled contigs into individual genomes using an EM algorithm over tetranucleotide frequency and abundance data. Currently at v0.2.0: correctness fixes on top of the equivalence-tested rewrite.

The project follows the [rewrites.bio](https://rewrites.bio) principles: faithful algorithm preservation while fixing packaging, performance, and reliability problems. See TODO.md for the full complaint catalog, roadmap, and algorithm summary.

## Development Environment

Uses a Nix flake devshell (enter via `direnv allow` or `nix develop`). The shell provides:
- Rust stable toolchain (edition 2024) with rust-analyzer
- cargo-nextest (preferred test runner)
- hmmer and bowtie2 (reference tool dependencies)
- Original MaxBin2 available as `run_MaxBin.pl` for equivalence testing

## Build & Test Commands

```bash
cargo build              # build
cargo nextest run        # run all tests (preferred over cargo test)
cargo nextest run <name> # run a single test by name
cargo clippy             # lint
cargo fmt --check        # check formatting
```

Shell-based equivalence and integration tests live in `tests/`:
- `tests/equivalence-e2e.sh` — end-to-end comparison against original MaxBin2
- `tests/pipeline-stages.sh` — individual subcommand tests
- `tests/cli-equivalence.sh` / `tests/cli-integration.sh` — CLI behavior

## Architecture

### CLI (clap derive, `src/cli.rs`)

Uses standard double-dash flags (`--contig`, `--abund`, etc.). Single-dash compatibility was dropped in v0.2.0. Subcommands expose individual pipeline stages:

- `maxbin-rs pipeline` — full pipeline (default)
- `maxbin-rs filter` — filter contigs by minimum length
- `maxbin-rs seeds` — generate seed file from HMMER marker gene hits
- `maxbin-rs em` — run the Rust EM algorithm
- `maxbin-rs cpp-em` — run the original C++ EM via FFI (equivalence testing only)
- `maxbin-rs sam-to-abund` — compute abundance from a SAM file

### Key source files

- `pipeline.rs` — orchestration for each subcommand
- `emanager.rs` — EM algorithm core (E-step, M-step, convergence)
- `distance.rs` — tetranucleotide distance and abundance probability
- `kmer_map.rs` — tetranucleotide frequency computation
- `abundance.rs` — abundance file parsing
- `fasta.rs` — FASTA reading/writing
- `original_ffi.rs` — C++ FFI bridge for equivalence testing against the original EM
- `external.rs` — shelling out to HMMER, Bowtie2, gene caller
- `profiler.rs` — timing instrumentation

### Parallelism

Both the inner abundance-file loop and the outer contig loop in the EM E-step are parallelized with Rayon, controlled by `--thread`.

## Algorithm Overview

1. Detect single-copy marker genes (107 bacterial / 40 universal) via HMMER + gene caller
2. Seed initial bins from marker gene clusters
3. Iteratively assign contigs to bins based on tetranucleotide frequency + abundance probability (EM)
4. Refine until convergence
5. Contigs below confidence threshold remain unclassified

Inputs: assembled contigs (FASTA) + per-sample read depth (abundance files or BAM).
Outputs: one FASTA per bin + summary statistics + `.log` file.

## Testing Strategy

- **Equivalence tests** compare Rust output against original MaxBin2 on the same input, including reproducing known bugs. The C++ EM is available via FFI (`cpp-em` subcommand) for component-level comparison.
- **Property tests** (proptest) cover parsing, distance computation, EM numerics, and sorting.
- **Benchmarks** via `cargo nextest run --cargo-profile bench bench_components -- --ignored`.
- When an equivalence test reveals a bug in the original, **document it** (in TODO.md and/or inline) but still reproduce the buggy behavior in the v0.1.x code path.
- Bug fixes are intentional, separately-tested changes — never mixed with the rewrite.

## Versioning

See [VERSIONING.md](VERSIONING.md). Key rule: **"breaking change" = changes binning
output on the same input.** v0.1.x is bug-for-bug compatible with original MaxBin2.
v0.2.x introduces correctness fixes that change output (recoverable via explicit flags).
Bug fixes (crashes, packaging, diagnostics) never need a minor version bump.

## Remaining Design Goals

- **Decouple from gene calling** — the tool no longer embeds a gene caller. The `--faa` flag accepts pre-computed protein FASTA from any gene caller (Prodigal, MetaGeneMark, etc.). The `maxbin-rs-em` flake output provides the EM and utility subcommands without external tool dependencies.
- Drop-in replacement for nf-core/mag pipeline
- Per-bin marker tarball (`marker_of_each_bin.tar.gz`) for full nf-core/mag compatibility
