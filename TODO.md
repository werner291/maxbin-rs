# maxbin-rs

Rust replacement for MaxBin2, a metagenome binning tool.

MaxBin2 is popular but unmaintained (last release ~2021, SourceForge, single
maintainer, all 13 tickets unanswered). It is a default binner in nf-core/mag.
The nf-core/mag team flagged it as a pain point (2026-04-13 Slack conversation
with James Fellows Yates and Jim Downie).

## Known complaints (with sources)

### Install / dependency hell

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

### Performance

- [ ] Single-threaded bottleneck despite `-thread` flag — one user ran 77 days
  [SqueezeMeta #69](https://github.com/jtamames/SqueezeMeta/issues/69),
  [SourceForge ticket #13](https://sourceforge.net/p/maxbin2/tickets/13/),
  [metaWRAP #516](https://github.com/bxlab/metaWRAP/issues/516)
- [ ] FragGeneScan initialization is slow — processes all marker genes including
  unreliable short sequences
  [metaWRAP #516](https://github.com/bxlab/metaWRAP/issues/516)

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

## Plan

1. **Nix flake** — package the original MaxBin2 + all deps reproducibly. Useful immediately,
   forces understanding of the dependency graph, provides a reference implementation to test against.
2. **Characterization tests** — run original MaxBin2 on known datasets, capture outputs.
   Can't rewrite what you can't verify.
3. **Rust core** — replace the C++ EM algorithm. Eliminate segfaults, fix parsing bugs,
   proper error handling.
4. **Kill Perl** — replace orchestration scripts with Rust CLI. Eliminate the entire class
   of `@INC` / `readline()` / path resolution bugs.
5. **Adopt maxbin2_custom approach** — Prodigal instead of FragGeneScan, GTDB markers.
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
