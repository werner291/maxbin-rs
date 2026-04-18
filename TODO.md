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

- FragGeneScan segfaults on macOS
  [SourceForge ticket #2](https://sourceforge.net/p/maxbin2/tickets/2/)
- Still shells out to Bowtie2, FragGeneScan, HMMER3 (bundled via Nix)

### Performance

- Outer contig loop not yet parallelized (embarrassingly parallel, biggest remaining win)
- FragGeneScan initialization is slow — future: replace with Prodigal

### Correctness

- [ ] Marker gene detection fails on valid data (CheckM finds 568 genes, MaxBin2 finds <=1)
  [SourceForge ticket #6](https://sourceforge.net/p/maxbin2/tickets/6/)
- [ ] Drops all viral contigs (by design — uses bacterial/archaeal markers only)
  [Biostars #371607](https://www.biostars.org/p/371607/)
- [ ] nf-core/mag was feeding averaged depth instead of per-sample depth in co-assembly
  [nf-core/mag #690](https://github.com/nf-core/mag/issues/690)

### Upstream context

- Decluttered fork [mruehlemann/maxbin2_custom](https://github.com/mruehlemann/maxbin2_custom)
  ripped out FragGeneScan, adopted GTDB markers, takes Prodigal output instead
- 2025 benchmark (Nature Comms, PMC11933696): MetaBAT2 competitive with deep learning
  binners while finishing in minutes; MaxBin2 not separately benchmarked but known slower

### Known gaps

- **`long double` vs `f64` precision gap** — inherent ~2.5 decimal digit
  difference. Affects a small number of contigs at decision boundaries.
  Cannot fix without 80-bit floats in Rust.
- **Summary formatting** — minor differences in "Bins without sequences"
  counter and `%0.4Lf` vs `{:.4}` rounding. Cosmetic, doesn't affect bins.
- `-verbose` and `-plotmarker` accepted but ignored
- No `.abundance` file for multi-sample runs
- No per-bin marker tarball (`marker_of_each_bin.tar.gz`) — requires gene
  caller output. Optional in nf-core/mag.
- `-reassembly` not implemented (IDBA-UD reassembly per bin)
- Reads format detection: extension-only, no bzip2
- Incomplete intermediate cleanup (SAM, `.bt2` index files)

### Code quality

- Replace hand-rolled CLI parsing (~1000 lines in `cli.rs`) with clap derive.
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
