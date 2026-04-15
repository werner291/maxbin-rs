# Acceptance Criteria — Drop-in Replaceability of MaxBin2

## Claim

**maxbin-rs produces identical output to MaxBin2 2.2.7 when given identical input,
accepts the same command-line interface, and can replace it in the nf-core/mag
pipeline without any downstream change.**

## How this is verified

We verify at three levels, each building on the one below:

1. **End-to-end**: run both tools on the same contigs + reads with the same flags.
   Diff every output file. If outputs match, the tools are interchangeable.
2. **Component**: link the original C++ code at compile time. For each function
   (parser, sorter, profiler, EM core, ...) feed random inputs to both
   implementations and assert identical output. This catches bugs that
   end-to-end tests might miss due to limited input diversity.
3. **Interface**: systematically exercise every CLI flag, input mode, and error
   path on both tools and compare behavior (exit codes, files produced, errors).

The evidence for each level:

| Level | Test | Dataset | Run time | Invocation |
|-------|------|---------|----------|------------|
| Pipeline stages | B. fragilis (tiny) | ~1 MB contigs + reads | ~1 min | `nix run .#test-pipeline-stages` |
| Pipeline stages | test_minigut (small) | ~10 MB contigs + reads | ~2 min | `nix run .#test-pipeline-stages-minigut` |
| Pipeline stages | CAPES_S7 (real) | 25K contigs + 2.5 GB reads | ~10 min | `nix run .#test-pipeline-stages-capes` |
| Component | proptest + FFI | random inputs, 200-500 cases each | seconds | `cargo nextest run` |
| Interface | CLI flag matrix | synthetic args | seconds | `cargo nextest run` (unit tests in cli.rs) |

All E2E tests run both `run_MaxBin.pl` and `maxbin-rs` with **identical flags**
and compare:

| Output file | Comparison method | Why not byte-for-byte |
|-------------|-------------------|-----------------------|
| `*.NNN.fasta` (bins) | sha256 of each bin, compared as unordered set | Perl hash ordering randomizes bin numbering |
| `*.noclass` | byte-for-byte | Deterministic |
| `*.tooshort` | byte-for-byte | Deterministic |
| `*.summary` | sorted rows, numeric tolerance 1e-6 | Float formatting may differ (`{:.4}` vs `%0.4Lf`) |
| `*.marker` | sorted rows | Bin order may differ |
| `*.abund*` | numeric tolerance 1e-6 | Float formatting |
| `*.log` | existence + non-empty | Free-text format differs |

---

## 1. End-to-end output identity

On each dataset (B. fragilis, test_minigut, CAPES_S7):

- [ ] **1.1** Bin count matches.
- [ ] **1.2** Bin contents match (same contigs in same bins, unordered set comparison).
- [ ] **1.3** Unclassified contigs (`*.noclass`) byte-for-byte identical.
- [ ] **1.4** Too-short contigs (`*.tooshort`) byte-for-byte identical.
- [ ] **1.5** Summary data numerically identical within 1e-6.
- [ ] **1.6** Marker file identical after sorting.
- [ ] **1.7** Abundance files (when generated from reads) numerically identical within 1e-6.

## 2. CLI interface identity

### 2.A Every valid invocation of the original works identically

- [ ] **2.1** Reads-only mode: `-contig C -reads R -out O` → both succeed, abundance files generated.
- [ ] **2.2** Abundance-only mode: `-contig C -abund A -out O` → both succeed, no Bowtie2.
- [ ] **2.3** Multiple reads: `-reads R1 -reads2 R2` → both succeed, per-sample abundance.
- [ ] **2.4** Multiple abundance: `-abund A1 -abund2 A2` → both succeed.
- [ ] **2.5** File-of-filenames: `-reads_list L` and `-abund_list L` → both succeed.
- [ ] **2.6** Single-dash and double-dash syntax both accepted: `-contig` = `--contig`.
- [ ] **2.7** All optional parameters (`-min_contig_length`, `-max_iteration`,
      `-prob_threshold`, `-thread`, `-markerset`, `-plotmarker`, `-verbose`,
      `-preserve_intermediate`) accepted with same defaults.

### 2.B Every invalid invocation fails the same way

- [ ] **2.8** Missing `-contig`: both exit non-zero.
- [ ] **2.9** Missing `-out`: both exit non-zero.
- [ ] **2.10** No reads or abundance provided: both exit non-zero.
- [ ] **2.11** Nonexistent input file: both exit non-zero (not segfault).
- [ ] **2.12** Unknown flag: both exit non-zero.

### 2.C Defaults match

- [ ] **2.13** `-min_contig_length` = 1000, `-max_iteration` = 50,
      `-prob_threshold` = 0.5 (not 0.9 — known bug, reproduced),
      `-thread` = 1, `-markerset` = 107.

## 3. Component-level equivalence (proptest + FFI)

Each Rust function is tested against the original C++ compiled into a static
library (`libmaxbin2_ffi.a`). Random inputs are generated via proptest; both
implementations receive identical input and must produce identical output.

- [ ] **3.1** FASTA parser: 200+ random files, headers/sequences/record count match.
- [ ] **3.2** Abundance parser: 200+ random files with all separator types (tab, space, comma, semicolon).
- [ ] **3.3** Quicksort: 500+ random arrays including duplicates, sorted array + index match.
- [ ] **3.4** K-mer map: kmerlen 2-5, symmetric/non-symmetric, full table + per-kmer lookups match.
- [ ] **3.5** Profiler: 100+ random DNA sequences including N characters, frequency profiles match within 1e-12.
- [ ] **3.6** Normal distribution: 200+ random (mean, std, x) triples, PDF values match within 1e-12.
- [ ] **3.7** Distance metrics: 100+ cases each for Euclidean and Spearman (sequence + profile variants).
- [ ] **3.8** EM probability functions: formula cross-checked against C++ NormalDistribution FFI.
- [ ] **3.9** EM pipeline: full cycle on synthetic + real data produces byte-for-byte identical bins.

## 4. Known bugs reproduced

These bugs exist in MaxBin2 2.2.7. For equivalence, we reproduce them exactly.

- [ ] **4.1** `prob_threshold` defaults to 0.5 (help text says 0.9).
- [ ] **4.2** Abundance separator `&&` bug: consecutive separators not skipped.

## 5. No side effects

- [ ] **5.1** No files created in CWD (the original creates `<random>.log`).
- [ ] **5.2** No files created next to input contigs (the original writes there, breaking read-only mounts).
- [ ] **5.3** All output under the `-out` prefix only.

## 6. Performance (not a gate, but tracked)

- [ ] **6.1** Wall-clock time on CAPES_S7 ≤ 1.5× the original.
- [ ] **6.2** Peak RSS on CAPES_S7 ≤ 2× the original.
      Check: `nix run .#bench-capes-s7`.

## 7. Reproducibility

- [ ] **7.1** `cargo nextest run` passes in the devshell.
- [ ] **7.2** `nix run .#test-pipeline-stages` passes (B. fragilis).
- [ ] **7.3** `nix run .#test-pipeline-stages-minigut` passes.
- [ ] **7.4** `nix run .#test-pipeline-stages-capes` passes (CAPES_S7).
