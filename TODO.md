# maxbin-rs

Rust replacement for MaxBin2, a metagenome binning tool developed by
[Yu-Wei Wu, Blake A. Simmons, and Steven W. Singer](https://doi.org/10.1093/bioinformatics/btv638).

MaxBin2 is widely used (default binner in nf-core/mag) and has not seen
updates since approximately 2021. The issues cataloged below are about
packaging, maintenance, and edge-case correctness — the underlying
algorithm is sound and this project faithfully preserves it. Many of
these issues are common in academic software of this age and reflect the
lack of ongoing maintenance rather than poor original engineering.

## Known issues (with sources)

### Installation and dependency management

- [ ] Perl `@INC` path breakage in bioconda — recurring, never permanently fixed
  [bioconda-recipes #24707](https://github.com/bioconda/bioconda-recipes/issues/24707),
  [atlas #328](https://github.com/metagenome-atlas/atlas/issues/328),
  [metaWRAP #207](https://github.com/bxlab/metaWRAP/issues/207)
- [ ] `setting` file requires manual absolute paths to Bowtie2, FragGeneScan, HMMER3, IDBA-UD
  [Biostars #9473674](https://www.biostars.org/p/9473674/)
- [ ] FragGeneScan segfaults on macOS
  [SourceForge ticket #2](https://sourceforge.net/p/maxbin2/tickets/2/)
- [ ] C++11 compiler requirement breaks on old HPC clusters
  [SourceForge discussion](https://sourceforge.net/p/maxbin2/discussion/general/thread/7e55b8f0/)
- [ ] External deps: Bowtie2, FragGeneScan, HMMER3, IDBA-UD, Perl 5 + CPAN, R + gplots
- [ ] `use LWP::Simple` imported in `run_MaxBin.pl` but never used — pulls in the
  entire libwww-perl dependency tree unnecessarily
  (observed directly: had to patch it out to avoid pulling in HTTP::Status etc.)
- [ ] FragGeneScan 1.30 (bundled version) doesn't compile with modern GCC — missing
  forward declarations (`print_usage`, `get_rc_dna_indel`) are errors since GCC 14
  (observed directly: needed `-Wno-error=implicit-function-declaration` to build)

### Performance

- [ ] Single-threaded bottleneck despite `-thread` flag — one user ran 77 days
  [SqueezeMeta #69](https://github.com/jtamames/SqueezeMeta/issues/69),
  [SourceForge ticket #13](https://sourceforge.net/p/maxbin2/tickets/13/),
  [metaWRAP #516](https://github.com/bxlab/metaWRAP/issues/516)
- [ ] FragGeneScan initialization is slow — processes all marker genes including
  unreliable short sequences
  [metaWRAP #516](https://github.com/bxlab/metaWRAP/issues/516)

### File system side effects

- [ ] Creates a temporary log file (`<random_number>.log`) in the current working directory
  at startup — writes to CWD without the user's request, which fails on read-only
  directories and clutters shared workspaces
  (observed directly: `run_MaxBin.pl` line ~118–126)
- [ ] Writes temp files next to the input contig file (`$contig.fa`, etc.) —
  fails on read-only filesystems (e.g. Nix store, network mounts).
  nf-core works around this by copying contigs into a local `input/` directory.
  (observed directly: fails with "Read-only file system" when contig is in `/nix/store`)

### Crashes

- [ ] Segfault on missing contig in abundance file — no error message
  [SourceForge ticket #1](https://sourceforge.net/p/maxbin2/tickets/1/)
- [ ] Misparses decimal abundance values on large datasets (2.3M lines) — treats
  `0.0825536598789213` as a line number
  [SourceForge ticket #12](https://sourceforge.net/p/maxbin2/tickets/12/)
- [ ] Silent C++ core crash with empty log file
  [SourceForge ticket #5](https://sourceforge.net/p/maxbin2/tickets/5/),
  [SourceForge ticket #8](https://sourceforge.net/p/maxbin2/tickets/8/)
- [ ] Resource exhaustion → `readline() on closed filehandle` on large datasets —
  16 comments, never resolved
  [metaWRAP #26](https://github.com/bxlab/metaWRAP/issues/26)
- [ ] Non-deterministic failures finding its own dependencies — workaround: `--restart-times 4`
  [atlas #314](https://github.com/metagenome-atlas/atlas/issues/314)
- [ ] Compressed input (.fa.gz) not detected, produces misleading errors
  [metaWRAP #469](https://github.com/bxlab/metaWRAP/issues/469)

### Correctness

- [ ] `prob_threshold` default mismatch: help says 0.9, code uses 0.5
- [ ] AbundanceLoader separator-skipping loop uses `&&` instead of `||` (line 110) —
  a char can't equal `\t` AND ` ` AND `,` AND `;` simultaneously, so multiple
  consecutive separators are never skipped. Harmless in practice since abundance
  files typically use a single tab, but would misparse files with e.g. `\t\t`.
  (observed directly in `AbundanceLoader.cpp`)
  [SourceForge ticket #7](https://sourceforge.net/p/maxbin2/tickets/7/)
- [ ] Same contigs in multiple bins
  [SourceForge ticket #10](https://sourceforge.net/p/maxbin2/tickets/10/)
- [ ] Marker gene detection fails on valid data (CheckM finds 568 genes, MaxBin2 finds <=1)
  [SourceForge ticket #6](https://sourceforge.net/p/maxbin2/tickets/6/)
- [ ] Drops all viral contigs (by design — uses bacterial/archaeal markers only)
  [Biostars #371607](https://www.biostars.org/p/371607/)
- [ ] nf-core/mag was feeding averaged depth instead of per-sample depth in co-assembly
  [nf-core/mag #690](https://github.com/nf-core/mag/issues/690)

### Upstream context

- [ ] All 13 SourceForge tickets open, most unanswered
- [ ] Decluttered fork [mruehlemann/maxbin2_custom](https://github.com/mruehlemann/maxbin2_custom)
  ripped out FragGeneScan, adopted GTDB markers, takes Prodigal output instead
- [ ] 2025 benchmark (Nature Comms, PMC11933696): MetaBAT2 competitive with deep learning
  binners while finishing in minutes; MaxBin2 not separately benchmarked but known slower

### Pipeline equivalence gaps (found via CLI surface comparison)

Remaining:

- [ ] **`long double` vs `f64` precision gap** — the original uses 80-bit
  `long double` on x86-64 Linux (but 64-bit on macOS/Windows/ARM). This
  is an inherent ~2.5 decimal digit precision difference that cannot be
  closed without 80-bit floats in Rust. Affects a small number of contigs
  at decision boundaries where probability ≈ 0.5 ± 10⁻¹⁶.
- [ ] **Summary "Bins without sequences" counter bug** — C++ reuses loop variable `j`
  in write_result (EManager.cpp:1424-1438), causing wrong numbering when ab_num > 1.
  Our Rust uses a separate counter and numbers correctly. Not reproduced — unlikely
  anyone parses this section.
- [ ] **Summary sprintf rounding** — C++ uses `%0.4Lf` (long double) vs Rust `{:.4}`
  (f64). Rounding at the 4th decimal boundary could differ due to different input
  precision. Cosmetic — affects summary file, not bins.
- [ ] **Seed ordering always deterministic** — seeds are always sorted alphabetically.
  The `MAXBIN_RS_DETERMINISTIC` env var documented in README is never read; the random
  shuffle path (matching the original's Perl hash iteration randomness) is not implemented.
  README claims the env var controls this. Either implement the shuffle or update docs.
- [ ] **Abundance file header stripping** — original strips leading `>` from abundance
  headers (run_MaxBin.pl:1505-1532 `checkAbundFile`). We don't, so abundance files
  with `>contig_name\tabundance` format will fail to match FASTA headers.

Moderate — missing features:

- [ ] `-verbose` flag accepted but silently ignored (original passes to C++ core)
- [ ] `-plotmarker` flag accepted but silently ignored (original generates heatmap PDF via R)
- [ ] No `.log` file output (original writes structured log to `$out_f.log`)
- [ ] No `.abundance` file for multi-sample runs (run_MaxBin.pl:841-851)
- [ ] No per-bin marker tarball (`marker_of_each_bin.tar.gz`)
- [ ] `-reassembly` flag not accepted (original does IDBA-UD reassembly per bin)

Minor:

- [ ] Reads format detection: only checks extension, no file-content fallback, no bzip2
  (original reads first line and supports `.q.gz`, `.a.gz`, `.q.bz2`, `.a.bz2`)
- [ ] Incomplete intermediate cleanup (SAM files, `.bt2` index files, reads-derived
  abundance files not deleted)

### Code quality

- [ ] Replace hand-rolled CLI parsing (~1000 lines in `cli.rs`) with clap derive.
  Clap doesn't support MaxBin2's single-dash long flags (`-contig`), so preprocess
  argv to rewrite `-contig` → `--contig` before clap sees it. Keeps backwards
  compatibility, eliminates most of the parsing code.

### Testing


## Plan

1. **Nix flake** — package the original MaxBin2 + all deps reproducibly. Useful immediately,
   forces understanding of the dependency graph, provides a reference implementation to test against.
2. **Characterization tests** — run original MaxBin2 on known datasets, capture outputs.
   Can't rewrite what you can't verify.
3. **Rust core** — replace the C++ EM algorithm. Eliminate segfaults, fix parsing bugs,
   proper error handling.
   - **Done**: The E-step abundance probability computation is now parallelized over
     abundance files using Rayon, matching the C++ `ThreadPool` + `threadfunc_E` behavior.
     Controlled by the `-thread` flag (via `EmParams::thread_num`). Output is bit-for-bit
     identical to the sequential version (the per-file computations are independent).
   - **Future**: Parallelize the *outer* contig loop (which the C++ does NOT do).
     This is where the real speedup is — 25K contigs × 10 bins = 250K independent
     distance computations per EM iteration, embarrassingly parallel with Rayon.
4. **Replace Perl orchestration** — rewrite the orchestration scripts in Rust. Eliminates
   the entire class of `@INC` / `readline()` / path resolution bugs.
5. **Adopt maxbin2_custom approach** — Prodigal instead of FragGeneScan, GTDB markers.
   - FragGeneScan is optimized for short reads but MaxBin2 feeds it assembled contigs —
     a mismatch that Prodigal handles better for this use case. Prodigal is faster,
     better at gene boundary prediction on contigs, and actively maintained.
   - Implement as a CLI option (e.g. `--gene-caller prodigal|fraggenescan`) so users
     can switch. Default to Prodigal, keep FragGeneScan for backwards compatibility.
   - When shelling out to external tools (HMMER, Bowtie2, gene caller), add comments
     citing the known issues: FragGeneScan segfaults (ticket #2), gene calling on
     assembled contigs is a misuse of FGS, and the `setting` file path resolution
     mess. These are the problems the rewrite eliminates.
   - Bug fixes (prob_threshold default, separator parsing, etc.) are out of scope for
     the initial rewrite — equivalence first, fixes as separate tracked changes.
6. **nf-core/mag module** — drop-in replacement.

## Algorithm summary

EM algorithm for metagenome binning:
1. Detect single-copy marker genes (107 bacterial / 40 universal) via HMMER + gene caller
2. Seed initial bins from marker gene clusters
3. Iteratively assign contigs to bins based on tetranucleotide frequency + abundance probability
4. Refine until convergence
5. Contigs below confidence threshold remain unclassified

Inputs: assembled contigs (FASTA) + per-sample read depth (abundance files or BAM)
Outputs: binned contigs (one FASTA per bin) + summary statistics

## References

- [MaxBin2 paper](https://pubmed.ncbi.nlm.nih.gov/26515820/)
- [MaxBin2 SourceForge](https://sourceforge.net/projects/maxbin2/)
- [maxbin2_custom fork](https://github.com/mruehlemann/maxbin2_custom)
- [nf-core/mag pipeline](https://nf-co.re/mag/latest/docs/usage/)
- [rewrites.bio manifesto](https://rewrites.bio)
- [2025 binning benchmark](https://www.nature.com/articles/s41467-025-57957-6)
