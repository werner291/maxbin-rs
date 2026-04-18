# TODO

## Known issues (with sources)

### Installation and dependency management

- Full pipeline mode still shells out to Bowtie2, FragGeneScan, HMMER3 (bundled via Nix).
  The `--faa` flag and `maxbin-rs-em` package avoid the gene caller dependency entirely.
- FragGeneScan segfaults on macOS — moot when using `--faa`
  [SourceForge ticket #2](https://sourceforge.net/p/maxbin2/tickets/2/)

### Correctness

- [ ] Marker gene detection fails on valid data (CheckM finds 568 genes, MaxBin2 finds <=1)
  [SourceForge ticket #6](https://sourceforge.net/p/maxbin2/tickets/6/)
- [ ] Drops all viral contigs (by design — uses bacterial/archaeal markers only)
  [Biostars #371607](https://www.biostars.org/p/371607/)
- [ ] nf-core/mag was feeding averaged depth instead of per-sample depth in co-assembly
  [nf-core/mag #690](https://github.com/nf-core/mag/issues/690)

### Known gaps

- **Noclass file assembly order** — on full CAMI I High (42k contigs),
  the noclass file differs from the original despite containing the same
  contigs. All 253 bins are byte-identical. The difference is in how
  noclass files from recursive EM rounds are merged. No pipeline parses
  noclass internals.
- **Summary formatting** — minor differences in "Bins without sequences"
  counter and `%0.4Lf` vs `{:.4}` rounding. Cosmetic, doesn't affect bins.
- `-verbose` and `-plotmarker` accepted but ignored
- No `.abundance` file for multi-sample runs
- No per-bin marker tarball (`marker_of_each_bin.tar.gz`) — requires gene
  caller output. Optional in nf-core/mag.
- `-reassembly` not implemented (IDBA-UD reassembly per bin)
- Reads format detection: extension-only, no bzip2
- **Cleanup deletes user-provided `--hmmout` file** — pipeline.rs cleanup
  doesn't distinguish generated files from user inputs. Live bug, fixed by v0.3
  work-dir approach.

## Plan

1. **v0.3: File I/O rework** — clean up how the pipeline handles intermediate and output files.
   - `--out <dir>` becomes an output directory (default `./output/`). Bins are
     `001.fasta`, `002.fasta`, etc. inside it. Errors if directory exists and is
     non-empty, unless `--force-overwrite`.
   - Intermediates (filtered contigs, SAM, .bt2, .faa, .hmmout, seeds, abundance)
     go in a tempdir (`$TMPDIR/maxbin-rs-<id>/`), auto-deleted on success.
   - `--keep-intermediates` keeps the work dir; defaults to
     `./intermediates/maxbin-rs-<id>/` when no `--work-dir` is set.
   - `--work-dir <path>` overrides intermediate location (always creates
     `<path>/maxbin-rs-<id>/` subdirectory).
   - On error: always preserve intermediates for debugging.
   - Gitignore `intermediates/`.
   - Fixes the bug where cleanup deletes user-provided input files — only the
     work dir (which only contains generated files) is ever deleted.
   - Pipe bowtie2 SAM output to stdout where possible to avoid landing SAM files.
   - Keep filtered contigs and seeds in memory instead of writing to disk.
2. **nf-core/mag module** — drop-in replacement.
