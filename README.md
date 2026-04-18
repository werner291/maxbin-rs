# maxbin-rs

**Status: Equivalence-tested, not yet pipeline-tested.** At the same float
width (`double`/`f64`), all output is bit-for-bit identical to MaxBin2 2.2.7
through 5 levels of recursive binning. The only divergence is the original's
use of 80-bit `long double` on x86-64 Linux (~4 contigs on CAMI I High).
~8x faster EM core. Not yet tested inside nf-core/mag or other pipelines —
if you try it, please report issues.

---

A Rust reimplementation of [MaxBin2](https://sourceforge.net/projects/maxbin2/),
a metagenome binning tool that clusters assembled contigs into individual genomes
using tetranucleotide frequency and abundance information.

MaxBin2 was created by
[Yu-Wei Wu, Blake A. Simmons, and Steven W. Singer](https://doi.org/10.1093/bioinformatics/btv638).
This project reimplements their algorithm while addressing long-standing
packaging and reliability issues, guided by the
[rewrites.bio](https://rewrites.bio) principles. It would not exist
without the original MaxBin2.

For a detailed account of how this rewrite was done and verified, see
[PAPER.md](PAPER.md).


## Quick start

```bash
# Full pipeline — bundles HMMER, Bowtie2, FragGeneScan via Nix
nix run github:werner291/maxbin-rs -- \
  --contig contigs.fa.gz \
  --reads reads.fastq.gz \
  --out my_bins
```

### Skip the pipeline, just run the EM

If you already have seeds (e.g. from your own marker gene detection),
you can skip the gene caller, HMMER, and Bowtie2 entirely:

```bash
# EM only — no external tools needed
maxbin-rs em --contig filtered.fa --abund depth.txt --seed seeds.txt --out result
```

This is the core algorithm. Everything else in the pipeline is
orchestration to produce the seeds and abundance data.

### Bring your own intermediates

The pipeline also accepts pre-computed intermediates, letting you skip
expensive stages you've already run:

```bash
# Skip gene calling + HMMER (bring your own HMMER output)
maxbin-rs --contig contigs.fa --abund depth.txt --hmmout hits.txt --out result

# Skip read alignment (bring your own abundance file)
maxbin-rs --contig contigs.fa --abund depth.txt --out result
```

**Note:** the full pipeline shells out to HMMER, Bowtie2, and
FragGeneScan at runtime. The Nix package bundles these automatically.
If you build with `cargo build` instead, you must have these tools on
your PATH.

### Docker

```bash
docker pull ghcr.io/werner291/maxbin-rs:v0.1.1
docker run --rm \
  -v /path/to/data:/data \
  ghcr.io/werner291/maxbin-rs:v0.1.1 \
  -contig /data/contigs.fa.gz \
  -reads /data/reads.fastq.gz \
  -out /data/output/bins
```

The image is built by Nix and includes all runtime dependencies. Published
automatically on version tags.

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

See: `tests/proptest_emanager.rs`, `tests/proptest_fasta.rs`, `tests/emanager_equivalence.rs`, and other `tests/proptest_*.rs` files

### 2. Pipeline stage tests (~1–10 minutes per dataset)

The primary verification. Tests each pipeline stage independently using
pre-computed intermediates from the original MaxBin2. No Bowtie2/HMMER runs
during the test itself — the slow parts are cached by Nix.

```bash
# B. fragilis — small, fast (~1 min, runs in CI)
nix build .#test-pipeline-stages

# minigut — multi-organism (~2 min)
nix build .#test-pipeline-stages-minigut

# CAPES_S7 — 25K contigs, ~2.5 GB download (~10 min)
nix build .#test-pipeline-stages-capes

# CAMI I High — 36K contigs, 240 bins (~50 min)
nix build .#test-pipeline-stages-cami
```

These run in the Nix build sandbox — isolated filesystem, no network, no
host environment. The only inputs are the explicit Nix store paths.

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

maxbin-rs always sorts seeds deterministically. The equivalence tests
patch the original Perl to do the same (via sorted hash iteration), so
both tools produce identical seed ordering and comparable output.

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
- Per-bin marker gene tarball (`marker_of_each_bin.tar.gz`)
- Multi-sample `.abundance` output file

## Compatibility

### What works

All CLI flags accepted by the original: `-contig`, `-abund` (including
`-abund2`, `-abund3` etc), `-reads`, `-thread`, `-out`,
`-min_contig_length`, `-max_iteration`, `-prob_threshold`, `-markerset`,
`-verbose`, `-plotmarker`. Output files: `*.fasta` bins, `.noclass`,
`.tooshort`, `.summary`, `.marker`, `.log`.

### Known differences from the original

These are intentional or unavoidable and should not affect correctness:

- **Binary name**: `maxbin-rs`, not `run_MaxBin.pl`. Pipeline modules
  that hardcode the binary name need updating.
- **Version output**: `maxbin-rs 0.1.3` (via `-v` or `--version`),
  not `MaxBin 2.2.7`. Pipelines that parse version strings (e.g.
  nf-core/mag's `sed 's/MaxBin //'`) need a one-line fix.
- **Progress output goes to stderr**, not stdout. The original writes
  to both stdout and the log file. Pipelines capturing stdout will see
  less output but this shouldn't affect behavior.
- **Exit code**: `1` on error (the original uses `-1` / 255 for some
  errors). Pipelines checking for non-zero should be unaffected.
- **No temp files next to input**: the original writes decompressed
  contigs and other temp files next to the input file, which fails on
  read-only filesystems. We write temp files relative to the output
  prefix instead.

### Missing output files

- Per-bin marker tarball (`marker_of_each_bin.tar.gz`)
- Multi-sample `.abundance` file
- `--verbose` and `--plotmarker` are accepted but silently ignored
- `--gene-caller prodigal` is accepted but **untested** — the code
  path exists but prodigal is not bundled in the Nix flake and has
  never been run end-to-end. Do not use in production.

### nf-core/mag

The goal is a drop-in replacement for the MaxBin2 module in
[nf-core/mag](https://nf-co.re/mag). The missing output files above
are `optional: true` in the nf-core module and won't crash the pipeline.
**Not tested in a real nf-core/mag run yet.** If you try this, please
report any issues.

## Versioning

See [VERSIONING.md](VERSIONING.md) for the full policy. In short:

- **v0.1.x** — bug-for-bug compatible with the original MaxBin2.
  Output is bit-identical at the same float width. Intended as a
  drop-in replacement, but not yet tested inside real pipelines.
- **v0.2.x** — correctness fixes that change output (e.g.,
  `prob_threshold` default, duplicate contigs).
- **v0.3+** — new features (Prodigal, GTDB markers, parallelized
  outer loop, nf-core/mag module).

"Breaking change" means "changes binning output on the same input." Crash fixes,
better errors, and performance improvements are never breaking.

## Citation

If you use this tool, please cite the original work:

> Wu Y-W, Simmons BA, Singer SW. **MaxBin 2.0: an automated binning algorithm to
> recover genomes from multiple metagenomic datasets.** *Bioinformatics.*
> 2016;32(4):605–607. doi:[10.1093/bioinformatics/btv638](https://doi.org/10.1093/bioinformatics/btv638)

This reimplementation follows the [rewrites.bio](https://rewrites.bio) manifesto
for responsible rewrites of bioinformatics tools, and would not exist without
the original MaxBin2.
