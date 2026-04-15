# MaxBin2 Original Codebase Analysis

Research notes for the maxbin-rs reimplementation.

## Codebase Size: ~7,300 lines

| Language | Lines | Role |
|----------|-------|------|
| C++ (16 .cpp + 14 .h) | ~4,500 | EM binning core |
| Perl (4 scripts) | ~2,260 | Pipeline orchestration + helpers |
| Perl (buildapp) | ~261 | Dependency bootstrapper |
| Shell (autobuild_auxiliary) | ~48 | Auxiliary build |
| R (heatmap.r) | ~19 | Marker gene heatmap |

### C++ Breakdown

| File | Lines | Purpose |
|------|-------|---------|
| `EManager.cpp` | ~1,850 | **EM algorithm core** — the largest file by far |
| `fastaReader.cpp` | ~646 | FASTA I/O |
| `kmerMap.cpp` | ~286 | Tetranucleotide frequency profiling |
| `SpeciesLoader.cpp` | ~262 | Marker gene / species loading |
| `AbundanceLoader.cpp` | ~187 | Abundance file parsing |
| `quickSort.cpp` | ~124 | Sorting utility |
| `main.cpp` | ~107 | Entry point for the C++ binary |
| `Profiler.cpp` | ~104 | K-mer profiling |
| `SpearmanDist.cpp` | ~92 | Spearman distance metric |
| `AbstractDist.cpp` | ~75 | Base distance class |
| `ThreadPool.cpp` | ~74 | Thread pool for parallelism |
| `ManhattanDist.cpp` | ~73 | Manhattan distance metric |
| `logger.cpp` | ~57 | Logging |
| `EucDist.cpp` | ~47 | Euclidean distance metric |
| `NormalDistribution.cpp` | ~18 | Normal distribution |
| Distance variants (Kendall, MinManhattan, RatioManhattan) | ~200 | Additional distance metrics |
| 14 header files | ~500-700 | Declarations (EManager.h alone is 102) |

### Perl Breakdown

| File | Lines | Purpose |
|------|-------|---------|
| `run_MaxBin.pl` | ~1,300 | **Main orchestration** — the entry point users call |
| `_getmarker.pl` | ~551 | HMM marker gene extraction |
| `_sepReads.pl` | ~291 | Read separation into bins |
| `_getabund.pl` | ~119 | Abundance calculation from SAM files |

### Data Files

- `marker.hmm` — 107 bacterial single-copy marker gene HMM profiles
- `bacar_marker.hmm` — 40 bacterial+archaeal universal marker HMM profiles

## CLI Interface

### Usage

```
run_MaxBin.pl
  -contig (contig file)                          # REQUIRED
  -out (output file prefix)                      # REQUIRED

  # Abundance input (at least one required)
  [-abund (abundance file)]
  [-abund2 (abundfile) -abund3 (abundfile) ...]
  [-abund_list (list of abundance files)]

  # Reads input (alternative to abundance — MaxBin maps via Bowtie2)
  [-reads (reads file)]
  [-reads2 (readsfile) -reads3 (readsfile) ...]
  [-reads_list (list of reads files)]

  # Parameters
  [-min_contig_length (minimum contig length. Default 1000)]
  [-max_iteration (maximum EM iteration number. Default 50)]
  [-thread (thread num; default 1)]
  [-prob_threshold (probability threshold. Default 0.9)]
  [-markerset (107 or 40. Default 107)]
  [-plotmarker]

  # Debug
  [-version] [-v]
  [-verbose]
  [-preserve_intermediate]
```

### Required Arguments

| Flag | Description |
|------|-------------|
| `-contig` | Assembled contigs/scaffolds (FASTA). Gzip poorly detected. |
| `-out` | Output file prefix. All outputs named from this. |
| `-reads`/`-abund` | At least one abundance source required. |

### Optional Arguments

| Flag | Type | Default | Notes |
|------|------|---------|-------|
| `-min_contig_length` | int | 1000 | Contigs shorter than this are filtered to `.tooshort` |
| `-max_iteration` | int | 50 | EM iteration cap |
| `-thread` | int | 1 | Actual parallelism is limited (known bottleneck) |
| `-prob_threshold` | float | 0.9 | **Known bug:** help says 0.9, code uses 0.5 |
| `-markerset` | int | 107 | `107` = bacterial, `40` = universal |
| `-plotmarker` | flag | off | Generates R heatmap |
| `-verbose` | flag | off | |
| `-preserve_intermediate` | flag | off | Keep temp files |
| `-version` / `-v` | flag | | Print version |

### Input File Formats

**Contig file (-contig):** Standard FASTA.

**Abundance file (-abund):** Tab-delimited, two columns:
```
contig_name\tabundance_value
```
If first line starts with `>`, the `>` is stripped. One file per sample.

**Reads file (-reads):** FASTQ or FASTA. MaxBin2 maps these via Bowtie2 internally.

### The `setting` File

Located in MaxBin2 install directory. Format:
```
[HMMER3] /path/to/hmmer3/bin
[BOWTIE2] /path/to/bowtie2
[FRAGGENESCAN] /path/to/FragGeneScan
[IDBA_UD] /path/to/idba_ud
```
Only `[HMMER3]` is actually used in practice — the rest are also searched via `$PATH`.

### Output Files

| File | Description |
|------|-------------|
| `{out}.summary` | Tab-delimited: bin name, completeness, genome size, GC content |
| `{out}.001.fasta`, `.002.fasta`, ... | One FASTA per bin (zero-padded 3 digits, 4 if >1000 bins) |
| `{out}.noclass` | Contigs below probability threshold (unclassified) |
| `{out}.tooshort` | Contigs filtered for being below `-min_contig_length` |
| `{out}.marker` | Marker gene counts per bin |
| `{out}.marker_of_each_bin.tar.gz` | Per-bin marker gene details |
| `{out}.abundance` | Multi-sample abundance matrix (only with >1 sample) |
| `{out}.log` | Execution log |

### Hardcoded Constants

- Version: `2.2.7`
- Min bin size: 100,000 bp
- K-mer length: 55
- Max recursion depth: 5
- Default marker set: 107

### Example Invocations

```bash
# With pre-computed abundance (most common)
run_MaxBin.pl -contig assembly.fasta -abund sample1.abund -abund2 sample2.abund -out mybin

# With reads (MaxBin2 maps internally via Bowtie2)
run_MaxBin.pl -contig assembly.fasta -reads sample1.fastq -reads2 sample2.fastq -out mybin -thread 8

# With abundance list file
run_MaxBin.pl -contig assembly.fasta -abund_list abundance_files.txt -out mybin
```

### External Dependencies Called at Runtime

| Program | Purpose |
|---------|---------|
| `hmmsearch` (HMMER3) | Search contigs against marker gene HMM profiles |
| `bowtie2` + `bowtie2-build` | Map reads to contigs for abundance (reads workflow only) |
| `FragGeneScan` | Gene calling on contigs (original; Prodigal in custom fork) |
| `idba_ud` | Reassembly (referenced but appears disabled) |
| `Rscript` | Marker heatmap (only with `-plotmarker`) |
| `MaxBin` (C++ binary) | The actual EM algorithm |
| `tar` | Archive per-bin marker files |

## Pipeline Flow

What `run_MaxBin.pl` actually does, in order:

1. Parse arguments, validate inputs
2. Filter contigs below `-min_contig_length` → write `.tooshort`
3. **If reads provided:** build Bowtie2 index, map each reads file, compute abundance
4. Run FragGeneScan (or Prodigal) on contigs for gene prediction
5. Run `hmmsearch` against marker HMM profiles on predicted genes
6. Parse HMM results → identify marker genes per contig (`_getmarker.pl`)
7. Seed initial bins from marker gene clusters
8. Invoke the C++ `MaxBin` binary with:
   - Contig file
   - Abundance file(s)
   - K-mer/tetranucleotide profiles
   - Marker gene assignments
9. C++ binary runs EM iterations → outputs bin assignments
10. Split contigs into per-bin FASTA files
11. Write `.summary`, `.marker`, `.noclass`, `.log`
12. If `-plotmarker`: generate R heatmap
13. Clean up intermediate files (unless `-preserve_intermediate`)

## maxbin2_custom Fork Differences

The [mruehlemann/maxbin2_custom](https://github.com/mruehlemann/maxbin2_custom) fork:
- Replaces FragGeneScan with **Prodigal** (adds `-prodigal` flag, requires pre-computed amino acid FASTA)
- Removes `-reads` support entirely (abundance-only)
- Adds GTDB marker sets: `120` (GTDB bacteria), `122` (GTDB archaea)
- Version string: `2.2.7MR`

This fork's approach aligns with our plan (step 5 in TODO.md).

## Sources

- GitHub mirror: https://github.com/galsang-git/maxbin2 (full source)
- Custom fork: https://github.com/mruehlemann/maxbin2_custom
- SourceForge: https://sourceforge.net/projects/maxbin2/
- Bioconda recipe metadata
