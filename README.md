# maxbin-rs

A Rust reimplementation of [MaxBin2](https://sourceforge.net/projects/maxbin2/),
a metagenome binning tool that clusters assembled contigs into individual genomes
using tetranucleotide frequency and abundance information.

MaxBin2 was created by
[Yu-Wei Wu, Blake A. Simmons, and Steven W. Singer](https://doi.org/10.1093/bioinformatics/btv638).
This project faithfully preserves their algorithm while addressing long-standing
packaging and reliability issues. It is not a fork — it is a clean
reimplementation guided by the [rewrites.bio](https://rewrites.bio) principles,
and would not exist without the original MaxBin2.

For a detailed account of how this rewrite was done and verified, see
[PAPER.md](PAPER.md).

**Status:** Pipeline stage equivalence verified on four datasets (B. fragilis,
minigut, CAPES_S7, CAMI I High). Component-level equivalence verified via
proptest + FFI against the original C++ (109 tests). See [Verifying the claims](#verifying-the-claims)
below.

## Quick start

```bash
# Run directly via Nix (bundles HMMER, Bowtie2, FragGeneScan automatically)
nix run github:werner291/maxbin-rs -- \
  -contig contigs.fa.gz \
  -reads reads.fastq.gz \
  -out my_bins

# Or build and install
nix build github:werner291/maxbin-rs
./result/bin/maxbin-rs --help
```

**Note:** maxbin-rs shells out to HMMER (`hmmsearch`), Bowtie2, and a gene
caller (FragGeneScan or Prodigal) at runtime. The Nix package bundles these
automatically. If you build with `cargo build` instead, you must have these
tools on your PATH.

### Adding to a NixOS / home-manager config

```nix
# flake.nix inputs
inputs.maxbin-rs.url = "github:werner291/maxbin-rs";

# Then in your config:
environment.systemPackages = [ inputs.maxbin-rs.packages.${system}.default ];
# or in home-manager:
home.packages = [ inputs.maxbin-rs.packages.${system}.default ];
```

## Verifying the claims

This project claims equivalence with the original MaxBin2 — byte-identical
for the EM core, within floating-point tolerance for abundance, identical
as sets for seed selection.
Here is how to check, in order of increasing time commitment:

### 1. Component-level tests (~1 minute)

Proptest-based equivalence: each Rust function is tested against the original
C++ via FFI on randomized inputs. 109 tests covering all major components.

```bash
nix develop
cargo nextest run
```

See: `tests/proptest_*.rs`, `tests/emanager_equivalence.rs`

### 2. Pipeline stage tests (~1–10 minutes per dataset)

The primary verification. Tests each pipeline stage independently using
pre-computed intermediates from the original MaxBin2. No Bowtie2/HMMER runs
during the test itself — the slow parts are cached by Nix.

```bash
# B. fragilis — small, fast (~1 min)
nix run .#test-pipeline-stages

# minigut — multi-organism (~2 min)
nix run .#test-pipeline-stages-minigut

# CAPES_S7 — 25K contigs, ~2.5 GB download (~10 min)
nix run .#test-pipeline-stages-capes
```

See: `tests/pipeline-stages.sh` for detailed documentation of what each
stage tests and example output.

### 3. Component performance comparison

Per-function Rust vs C++ throughput comparison:

```bash
nix run .#bench-components
```

See: `tests/bench_components.rs`

## What's different from the original

**Same**: algorithm, output format, known bugs (all reproduced intentionally).
This is a faithful reimplementation, not an improvement — yet.

### Known behavioral difference: seed ordering

The original `_getmarker.pl` iterates marker gene clusters with Perl's
`keys %hash`, which produces a **non-deterministic order** (randomized since
Perl 5.18). This means the original MaxBin2 can produce different bin counts
and assignments on the same input across runs — the EM algorithm converges to
different local optima depending on which contigs are used as initial seeds and
in what order.

By default, maxbin-rs **shuffles seeds randomly** at startup, matching the
original's effectively-random behavior. Setting `MAXBIN_RS_DETERMINISTIC=1`
sorts seeds alphabetically instead, producing reproducible output. The
pipeline stage tests use this mode, together with a patch to the original
Perl, so both tools produce comparable output.

**I lack the domain expertise to determine whether the original's random
seed ordering serves as a form of regularization.** If you are a metagenomics
researcher and have insight on this, I would welcome your input.

### Other differences

**Currently addressed**:
- Installability: Nix flake bundles the binary with all runtime dependencies
  (HMMER, Bowtie2, FragGeneScan). No `setting` file, no CPAN, no manual path
  configuration. Note: FragGeneScan's `run_FragGeneScan.pl` still requires
  a Perl interpreter (provided by Nix).
- Reproducibility: Nix flake builds everything from source with pinned deps
- Transparent provenance: original C++ and Perl are patched via `.patch` files
  (see `nix/`), not vendored — every change from upstream is auditable
- No temp files next to input: the original shells out to `gunzip -c` and
  writes a decompressed copy next to the input contig file, which fails on
  read-only filesystems (e.g. Nix store). maxbin-rs decompresses via
  streaming gzip instead, avoiding the temp file at the cost of some
  additional memory usage during decompression.

**Observed**: On CAMI I High (36K contigs × 577 bins), the Rust EM is
approximately 8x faster than the original C++ EM (see
[PAPER.md](PAPER.md) for methodology and caveats).

**Not yet addressed** (see `TODO.md`):
- Correctness bugs (prob_threshold default, duplicate contigs, etc.)
- FragGeneScan → Prodigal switch
- GTDB marker gene sets
- nf-core/mag drop-in module

## Citation

If you use this tool, please cite the original work:

> Wu Y-W, Simmons BA, Singer SW. **MaxBin 2.0: an automated binning algorithm to
> recover genomes from multiple metagenomic datasets.** *Bioinformatics.*
> 2016;32(4):605–607. doi:[10.1093/bioinformatics/btv638](https://doi.org/10.1093/bioinformatics/btv638)

This reimplementation follows the [rewrites.bio](https://rewrites.bio) manifesto
for responsible rewrites of bioinformatics tools, and would not exist without
the original MaxBin2.
