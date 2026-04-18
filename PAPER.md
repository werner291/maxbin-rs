# Rewriting MaxBin2 in Rust: A Case Study in Equivalence-First Bioinformatics Rewrites

**Werner Kroneman, April 2026**\
**Applies to: v0.2.0**

## Disclosure

This paper has two things to disclose upfront.

First: the Rust code and much of this text were produced in collaboration with
an LLM (Claude, by Anthropic). I directed the design, made architectural
decisions, and verified results; the LLM did most of the typing. I iterated
extensively on its output — rejecting drafts, correcting tone, reshaping
structure — but the literal words on this page were largely produced by the
model. It would be dishonest to present this as solely human work.

Second: I am not a biologist. My background is in computer science
(MSc, TU Eindhoven) and motion planning for autonomous drones. I have no training in metagenomics, microbial ecology,
or marker gene biology. My understanding of biology does not extend much
beyond high school level and watching science videos online. To be quite
honest, I have a superficial understanding of the EM algorithm at the core,
but not of the upstream stages — the gene calling, marker detection, and
abundance estimation that feed into it. I acknowledge the
apparent absurdity of attempting to rewrite a biology tool under these
conditions. Nevertheless, I claim to have done so successfully.

This is less absurd than it sounds. At the level of implementation,
MaxBin2 is not a biology problem — it is matrix arithmetic, probability
formulas, and file format wrangling. I do not need to understand the
biology to verify that two implementations produce identical output.

I state these upfront because they matter for evaluating this work.
Readers who consider either point disqualifying are welcome to stop here.

## Abstract

I present a Rust reimplementation of MaxBin2, a widely-used metagenome
binning tool. This case study touches on three topics: equivalence-first
rewriting as a methodology, Nix-based reproducibility, and LLM-assisted
development in practice. At the same float width (`double`/`f64`), all pipeline stages and the
full recursive EM produce bit-for-bit identical output. The only
divergence is between the original's `long double` (80-bit, x86-64
Linux only) and Rust's `f64` (64-bit): on the full CAMI I High
benchmark (42,038 contigs), ~4 contigs are classified differently at a
decision boundary where the classifier has no confidence (probability
≈ 0.5 ± 10⁻¹⁶). Known bugs are reproduced intentionally.

I extend the rewrites.bio approach by introducing Nix as a packaging and
reproducibility layer. The original MaxBin2 can be installed via Conda,
which handles much of the dependency pain — but Nix goes further: the Nix flake in this repository builds the tool
and all its dependencies from source with pinned inputs, and also fetches
and builds the original MaxBin2 and test data for equivalence testing. Computational reproducibility — historically a
substantial problem in bioinformatics — becomes straightforward, at least on the computational side. (Wet lab work is regrettably not yet declarative.) The main caveat is a non-deterministic
component in the original that I preserve as-is.

The accompanying repository contains the full implementation, equivalence
test suite, and a Nix-packaged copy of the original MaxBin2 for
independent verification.

## Introduction

This paper documents a Rust reimplementation of MaxBin2, a metagenome
binning tool. The reimplementation preserves the original algorithm and
reproduces its output — including known bugs — while replacing the
packaging and build infrastructure. It is not a fork; no original source
code was copied. The following sections cover the background and
motivation, the equivalence-first methodology used, and the results
obtained on four test datasets.

### Original work

MaxBin2 was created by Yu-Wei Wu, Blake A. Simmons, and Steven W.
Singer. This reimplementation would not exist without their work, and
the algorithm, design, and scientific contribution are entirely theirs.
If you are referring to MaxBin2 as a tool or algorithm — rather than to
this case study specifically — please cite the original:

> Wu Y-W, Simmons BA, Singer SW. **MaxBin 2.0: an automated binning
> algorithm to recover genomes from multiple metagenomic datasets.**
> *Bioinformatics.* 2016;32(4):605–607.
> doi:[10.1093/bioinformatics/btv638](https://doi.org/10.1093/bioinformatics/btv638)

## Background

Metagenome binning is the problem of taking a mixed bag of DNA sequences
from an environmental sample — soil, gut, ocean water — and sorting them
into groups that (hopefully) correspond to individual organisms. The input
is a set of assembled contigs; the output is bins, each ideally containing
the genome of one species. I will not pretend to explain the biology
further than that.

[MaxBin2](https://doi.org/10.1093/bioinformatics/btv638) (Wu, Simmons,
and Singer, 2016) is one of several tools that does this. It uses
tetranucleotide frequencies and read abundance across samples to cluster
contigs, seeded by single-copy marker genes found via
[HMMER](http://hmmer.org/). It is one of the binners included in the
[nf-core/mag](https://nf-co.re/mag) metagenomics pipeline and is widely
cited.

The algorithm is sound. The packaging has not kept pace.

MaxBin2 is implemented as a Perl orchestration script (`run_MaxBin.pl`)
driving a C++ core, with runtime dependencies on
[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/),
[HMMER3](http://hmmer.org/),
[FragGeneScan](https://sourceforge.net/projects/fraggenescan/), and
[IDBA-UD](https://github.com/loneknightpy/idba) (referenced in the
configuration but apparently unused in practice). These are located via a
hand-edited `setting` file containing absolute paths. Its
[bioconda](https://bioconda.github.io/) package has recurring Perl `@INC`
path breakage
([bioconda-recipes #24707](https://github.com/bioconda/bioconda-recipes/issues/24707),
[atlas #328](https://github.com/metagenome-atlas/atlas/issues/328),
[metaWRAP #207](https://github.com/bxlab/metaWRAP/issues/207)).
FragGeneScan [segfaults on macOS](https://sourceforge.net/p/maxbin2/tickets/2/)
and its bundled version (1.30) does not compile with GCC 14+ without
suppressing implicit function declaration errors (observed directly during
this project's Nix packaging effort). The Perl script imports
`LWP::Simple` but never uses it, pulling in the entire `libwww-perl`
dependency tree unnecessarily.
All [13 SourceForge tickets](https://sourceforge.net/p/maxbin2/tickets/)
are open; 3 received maintainer replies (2017–2018), the rest are
unanswered. The last commit was around 2021.

These are not engineering failures — they are the normal consequences of
academic software outliving its maintenance window. The `setting` file,
the Perl orchestration, the dependency management were all reasonable
choices in 2014–2016.

This project eliminates the Perl layer, the setting file, and the
unnecessary dependencies. The tool still shells out to FragGeneScan,
Bowtie2, and HMMER — Nix bundles these, which tames the packaging, but
does not fix upstream issues like the FragGeneScan macOS segfault.
Replacing FragGeneScan with
[Prodigal](https://github.com/hyattpd/Prodigal) (as in
[maxbin2_custom](https://github.com/mruehlemann/maxbin2_custom)) is
planned but not yet done.

## Related Work

MaxBin2 is one of many metagenome binners — others include
[MetaBAT2](https://bitbucket.org/berkeleylab/metabat),
[CONCOCT](https://github.com/BinPro/CONCOCT), and more recent deep
learning-based approaches. Benchmarks exist (e.g.
[Döring et al., 2025](https://doi.org/10.1038/s41467-025-57957-6)) but
are beyond my expertise to evaluate. MaxBin2 remains part of the
nf-core/mag pipeline and users depend on it. The question this
project asks is not "which binner is best" but "can a working tool be
repackaged without breaking it."

Two prior efforts are relevant. Rühlemann's
[maxbin2_custom](https://github.com/mruehlemann/maxbin2_custom) fork
strips out FragGeneScan and IDBA-UD, adopts
[GTDB](https://gtdb.ecogenomic.org/) marker gene sets, and takes
[Prodigal](https://github.com/hyattpd/Prodigal) output instead. This
addresses several of the same packaging complaints but as a fork — it
inherits the Perl orchestration and C++ core. I plan to adopt the same
Prodigal/GTDB approach in a later phase.

The [rewrites.bio](https://rewrites.bio) manifesto advocates for
rewrites that preserve the original algorithm exactly, verify
equivalence through testing, and improve only the engineering. This
project follows those principles and extends them with Nix-based
reproducibility — which rewrites.bio recommends in spirit (document
exact versions) but does not prescribe a mechanism for.

## Method

### Tools Used

There is an inherent tension in this project's toolchain. On one end of
the spectrum: [Nix](https://nixos.org/), a build system that aspires to
make every build perfectly deterministic and reproducible. On the other
end: an LLM, which is by nature stochastic and produces different output
on every run. These work together: Nix pins every dependency and build input, so when
the LLM generates something broken, you can at least reproduce the
breakage exactly.

[Rust](https://www.rust-lang.org/). Chosen for memory safety, strong
typing, and good C FFI — the last point critical for calling original C++
functions directly from test code. (Rust's ecosystem also makes
improvements easy — e.g. replacing the original's C++ thread pool with
[Rayon](https://github.com/rayon-rs/rayon) for the EM E-step was a
one-line iterator change, safe by construction. The outer contig loop —
the bigger optimization opportunity — remains single-threaded in both
implementations, waiting until equivalence is fully proven.)

[Nix](https://nixos.org/)
([Dolstra et al., 2004](https://edolstra.github.io/pubs/nspfssd-lisa2004-final.pdf)).
Packaging, reproducibility, and testing. The entire dependency graph —
Rust toolchain, HMMER, Bowtie2, FragGeneScan, the original MaxBin2, and
test datasets — is built from source with pinned inputs. Anyone with Nix
can build exactly what I built and run exactly the tests I ran. The
equivalence tests themselves are also pinned: if a dependency version
drifts, hashes change and the tests catch it. Nix also packages the original
MaxBin2, which was itself nontrivial (FragGeneScan doesn't compile with
modern GCC without patches, the Perl script imports unused modules, etc.).

[Claude™](https://claude.ai/) (Anthropic). Opus 4.6 for most code and
prose, Sonnet for web searches. The use of LLMs is controversial — the
environmental costs are real and the economic disruption is not
hypothetical. This rewrite would not have happened without one.

To mitigate LLM confirmation bias (see [Discussion](#discussion)), I used subagents
with fresh context for self-review — both for code and for this paper
(see [docs/code-review-prompts.md](docs/code-review-prompts.md) and
[docs/review-prompts.md](docs/review-prompts.md)).

[Proptest](https://github.com/proptest-rs/proptest). Property-based
testing framework for Rust. Each Rust function is tested against the
original C++ via FFI on randomized inputs, checking that both produce
identical output across thousands of cases.

### Equivalence Testing

Correctness in this project means "produces the same output as the
original." Not correct by theory, not correct by inspection — identical
output on the same input. This includes reproducing known bugs. If the original C++ computes a
wrong value for a degenerate input, the Rust reimplementation computes
the same wrong value, and a test verifies this. Bug fixes come later as
separate, documented changes.

This follows from [rewrites.bio](https://rewrites.bio)'s "emulate
exactly" principle (2.2), and it is what made the project tractable. With an oracle (the original C++ and
Perl, callable via FFI and Nix), every function has a testable answer.

I test at two levels.

#### Component-level

The build script ([build.rs](build.rs)) compiles the original C++ into a static library,
and [src/original_ffi.rs](src/original_ffi.rs) exposes it to Rust via
FFI. The [proptest](https://github.com/proptest-rs/proptest) framework
then generates random inputs, feeds them to both the Rust function and
the C++ function, and asserts identical output. 109 tests cover all major
components: FASTA parsing, abundance loading, kmer frequencies, distance
metrics, EM probability functions, profiling, normal distribution, and
quicksort. This layer caught real divergences: C's `atof` skips leading
whitespace while Rust's `parse::<f64>()` does not, and the original C++
computes a nonsensical `percent_N` (NaN or -0.0) for sequences shorter
than the kmer length due to signed integer underflow.
Both are now reproduced exactly. (See [tests/proptest_*.rs](tests/).)

#### Pipeline-level

The MaxBin2 pipeline has four stages: contig filtering, read mapping +
abundance, marker gene detection + seed selection, and EM binning.
Intermediate data (SAM alignments, HMMER output, abundance tables, seed
files) is produced by running the original MaxBin2 once per dataset,
cached as a Nix derivation ([nix/intermediates.nix](nix/intermediates.nix)) — the slow Bowtie2 and HMMER runs happen
once and are never repeated (both implementations call the same external
tools, so their output is identical). Each stage is then tested by feeding the
original's intermediate output to the Rust reimplementation and checking
that it produces matching results. The pipeline stage test
([tests/pipeline-stages.sh](tests/pipeline-stages.sh)) compares
at each interface: byte-exact for text files and abundance values,
and as unordered sets for bin contents.
(See the script header for detailed documentation and example output.)

Adding a new dataset requires one function call in the Nix configuration
([nix/datasets.nix](nix/datasets.nix)). I currently test on four
datasets, with a fifth
([MetaHIT](https://doi.org/10.7717/peerj.1165), 60K contigs) wired up
but not yet run:

- **B. fragilis** — nf-core/modules MaxBin2 test data, small and fast.
- **minigut** — nf-core/mag test profile, multi-organism.
- **CAPES_S7** — 25,244 contigs, ~2.5 GB reads from ENA.
- **CAMI I High** — 36,116 contigs, 577 bins. From the
  [CAMI challenge](https://doi.org/10.1038/nmeth.4458) ([Sczyrba et al., 2017](https://doi.org/10.1038/nmeth.4458)).

A caveat: the original MaxBin2 is nondeterministic. Perl 5.18+
randomizes hash key iteration order, and `_getmarker.pl` iterates marker
gene clusters with `keys %hash`. This means the seed ordering — and
therefore the EM convergence point — can differ between runs of the
original itself. I patch the original Perl to sort seeds alphabetically
for testing (`sort keys %hash` — see
[nix/maxbin2-deterministic.patch](nix/maxbin2-deterministic.patch)),
and maxbin-rs does the same when `MAXBIN_RS_DETERMINISTIC=1` is set.
There are further sources of nondeterminism and a float precision
difference (`long double` vs `f64`) that I expand on in the [Discussion](#nondeterminism-and-float-precision).

### Source Traceability

The Rust code is annotated with
comments mapping to the original source (e.g. "Matches
kmerMap.cpp:219-244" throughout [src/kmer_map.rs](src/kmer_map.rs), and
"Matches run_MaxBin.pl:366" in [src/pipeline.rs](src/pipeline.rs)),
allowing a reviewer to verify structural correspondence by reading both
implementations side by side.

The upstream source is fetched by hash from SourceForge
([flake.nix](flake.nix), [nix/maxbin2.nix](nix/maxbin2.nix)) and
modified only via `.patch` files
([nix/maxbin2-cpp-ffi.patch](nix/maxbin2-cpp-ffi.patch),
[nix/maxbin2-deterministic.patch](nix/maxbin2-deterministic.patch)),
so every change from the original is auditable.

## Results

### Equivalence

All pipeline stage tests pass on all four tested datasets (B. fragilis,
minigut, CAPES_S7, CAMI I High). At the component level, 109 tests
pass — proptest equivalence, inline FFI tests, and CLI parsing tests.

To reproduce:

```bash
# Component-level tests (~1 min)
nix develop -c cargo nextest run

# Pipeline stage tests per dataset
nix run .#test-pipeline-stages         # B. fragilis (~1 min)
nix run .#test-pipeline-stages-minigut # minigut (~2 min)
nix run .#test-pipeline-stages-capes   # CAPES_S7 (~10 min)
nix run .#test-pipeline-stages-cami    # CAMI I High (~50 min)
```

The pipeline stage tests verify:

| Stage | What is compared | Result |
|---|---|---|
| Contig filtering | tooshort reject file | Byte-identical |
| SAM → abundance | Per-contig depth values | Byte-identical |
| HMM → seeds | Seed contig names | Identical as sets |
| EM binning | Bin FASTA files, noclass, tooshort | Byte-identical |

Note: the EM stage comparison runs the original C++ EM core via FFI from
the Rust binary, not by invoking the full original `run_MaxBin.pl`. This
verifies that the algorithm implementation matches, but does not test
the Perl orchestration layer end-to-end. The earlier stages (filtering,
abundance, seeds) cover the orchestration logic independently.

CAMI I High is the largest dataset tested (36,863 filtered contigs,
577 seeds, 240 bins). Real-world metagenomics datasets can be
substantially larger — broader testing is needed.

The equivalence testing process uncovered several bugs in the original,
most harmless in practice. The most user-visible is a `prob_threshold`
default mismatch (help says 0.9, code uses 0.5). All are reproduced in
the reimplementation and documented in [TODO.md](TODO.md) for later
fixing. The reimplementation is approximately 4,400 lines of Rust
across 14 source files, with 1,500 lines of test code (the original is
roughly 2,500 lines of C++ and 1,500 lines of Perl).

CLI argument parsing is tested for `-reads_list`, `-abund_list`,
`-min_contig_length`, `-markerset`, and double-dash normalization
(see [tests/cli_list_files.rs](tests/cli_list_files.rs)). However,
the pipeline stage tests only exercise single-reads and pre-computed
abundance input modes on real data. End-to-end testing with
`-reads_list` on a multi-sample dataset has not been done.

### Performance

Performance was not a goal of this rewrite, but the results are worth
reporting.

**Component-level benchmarks**
([tests/bench_components.rs](tests/bench_components.rs)) compare
individual Rust functions against their C++ counterparts via FFI. Most
components are 1.2–4.5x faster in Rust. The one exception (`kmermap
lookup`, 0.65x) has not been investigated; proptest equivalence passes,
suggesting a performance issue rather than a correctness one.

**EM algorithm on CAMI I High** (36,863 contigs, 577 bins, 50
iterations). All measurements on a Framework 16 laptop (AMD Ryzen 9
7940HS, Zen 4, 60 GB RAM) running NixOS 26.05. Observed across multiple runs during development — C++ EM
iterations consistently take 30–40s each; Rust completes the full 50
iterations in under 5 minutes.
Reproducible via `nix run .#test-pipeline-stages-cami`.

Early measurements were confounded by IO overhead: the Rust code had
unbuffered file writes (since fixed), and the C++ had `sprintf` with
`%Lf` formatting inside the innermost loop — per contig, per seed,
per iteration — plus per-call file write and flush in the logger
(see [nix/maxbin2-cpp-ffi.patch](nix/maxbin2-cpp-ffi.patch)). After
stripping all logging and IO from both sides (both compiled at `-O3`,
timing via `clock_gettime` on the C++ side and `std::time::Instant`
on the Rust side):

| Implementation | EM loop | classify | write | total |
|---|---|---|---|---|
| C++ with `long double` (original) | ~2200s | <1s | ~4s | ~37 min |
| C++ with `double` (experiment) | 1886s | <1s | ~4s | ~31 min |
| Rust with `f64` | 231s | <1s | <1s | ~3.9 min |

The C++ uses `long double` (80-bit x87 on this platform) while Rust
uses `f64` (64-bit IEEE 754). At the same float width, output is
bit-identical; the `long double`-vs-`f64` difference affects ~4 contigs
on CAMI I High (see [Discussion](#nondeterminism-and-float-precision)).

I did not expect an ~8x gap and spent considerable effort trying to
close it from the C++ side. I stripped all `sprintf` and logging calls
from the inner loop
([nix/maxbin2-cpp-ffi.patch](nix/maxbin2-cpp-ffi.patch)), switched
from `long double` to `double`
([nix/maxbin2-cpp-ffi-f64.patch](nix/maxbin2-cpp-ffi-f64.patch)) —
~15% improvement. I suspected the remaining gap was due to
cross-compilation-unit calls preventing inlining — the C++ EM inner
loop calls into `EucDist.cpp`, `Profiler.cpp`, and `fastaReader.cpp`,
which GCC cannot inline across separate `.cpp` files without link-time
optimization. To test this, I compiled all C++ into a single translation
unit (unity build). It made no difference.

Disassembly of both sides (`objdump -d`, reproducible via
`nix build .#disasm-em`) explains why: neither is vectorized, but the
C++ `run_EM` makes 38 function calls — `malloc`/`free`, non-inlined
helpers, even calls to logging functions whose bodies are empty — while
Rust's `run_em` makes 39 calls over 990 instructions because LLVM
inlined the helpers into the function body. GCC does not eliminate the
empty-body calls even in the unity build. Closing the gap would require
rewriting the C++ inner loop, not changing compiler flags.

The EM is one stage of the pipeline — Bowtie2, HMMER, and FragGeneScan
also contribute runtime. The ~8x EM speedup does not translate to ~8x
end-to-end, but ~37 minutes of C++ EM time is not negligible.

## Discussion

The rewrites.bio equivalence-first principle made the project tractable.
"Correct means same output as the original" gives you an oracle for
every function — no ambiguity, no design decisions about what "better"
means, just match the output and move on.

The initial reimplementation took approximately two days of part-time
work — the EM core, the largest component, reached equivalence on the
first day. Testing on larger datasets, investigating the performance
results, and writing this paper took considerably longer.

### Methodology

The combination of an equivalence target, proptest, and C++ FFI meant
that every function had a testable answer from the start. Nix eliminated
"works on my machine" problems entirely — the original MaxBin2, its
dependencies, the test data, and the equivalence tests are all
reproducible from a single configuration.

In hindsight, many of the improvements did not require a rewrite —
buffered IO, Nix packaging, removing unused imports — these are
straightforward patches to the existing code.
The case for a Rust rewrite is weaker than "rewrite it in Rust" memes
suggest. Whether a rewrite is justified over direct maintenance is a
judgment call. This project exists because an LLM made the rewrite
cheap enough to attempt as a side project. If the manual effort had
been the full cost, patching the original would likely have been the
better investment.

### Technical surprises

None of the difficult parts were where I expected them. The algorithm
was straightforward to translate. The biology never came up.

The nondeterminism in seed ordering (see [Method](#equivalence-testing))
cost hours of debugging before I realized both implementations were
correct, just differently ordered. Nondeterminism is pervasive in the
stack — both Perl (since 5.18) and Rust randomize the iteration order
of their standard dictionary types (`%hash`, `HashMap`) to prevent
hash-flooding DoS attacks
([Klink and Wälde, 28C3 2011](https://events.ccc.de/congress/2011/Fahrplan/events/4680.en.html)),
and threaded tools like Bowtie2 and HMMER may have their own ordering
dependencies. This is why I
test at pipeline stage interfaces rather than end-to-end.

Matching C's `%.15g` float formatting in Rust was attempted but
ultimately abandoned — no built-in equivalent exists, and the effort
was disproportionate. The equivalence tests accept floating-point
tolerance instead.

Floating-point equivalence across language boundaries is harder than it
looks. The original C++ uses `long double` (80-bit extended precision on
x86-64 Linux, but 64-bit on macOS, Windows, and ARM — the original is
not even consistent with itself across platforms). Rust has no `long
double`; it uses `f64` (64-bit IEEE 754). Through recursive binning,
this precision difference produced different bin assignments on a small
number of contigs.

To isolate the cause, I patched the C++ to use `double` instead of
`long double` — matching the Rust float width exactly. This changes the original's behavior, but since `long double` is
already 64-bit on macOS, Windows, and ARM, the patched version is
equivalent to how the original behaves on those platforms. Even
with matching float widths, the output still diverged. Per-iteration hex
dumps of all EM state showed the first difference at iteration 7: every
individual math function (`lgamma`, `exp`, `sqrt`, normal distribution
PDF) was bit-identical between C++ and Rust, but the M-step
accumulation was not. The cause: IEEE 754 multiplication is not
associative, and the C++ evaluates `abundance * length * probability`
left-to-right while the Rust code had pre-computed
`weight = length * probability` and then multiplied
`abundance * weight`. The 1 ULP difference compounded over iterations,
eventually flipping a threshold decision at probability ≈ 0.5.

After matching the evaluation order, all output is bit-for-bit
identical at the same float width — verified across all EM iterations
and all 10 component-level FFI equivalence tests.

The remaining divergence is `long double` (80-bit) vs `f64` (64-bit):
a ~2.5 decimal digit precision gap. On CAMI I High (42,038 contigs)
this affects ~4 contigs — cases where the classifier converges to
0.5 ± 10⁻¹⁶ and has no confidence in either assignment.

### Working with LLMs

The LLM was productive at mechanical translation (C++ or Perl to
equivalent Rust), source annotation, and test scaffolding. It also
wrote tautological tests — comparing output to a reimplementation of
the same formula, which proves nothing. The equivalence target largely
prevents this: it is hard to write a tautological test when the oracle
is a separate C++ binary called via FFI. LLMs will also agree that
broken code works if presented confidently; fresh-context subagents
(separate invocations that see only the code, not the conversation)
helped catch this. Throughout, the LLM had to be repeatedly redirected
away from "improving" the original rather than reproducing it — the
core [rewrites.bio](https://rewrites.bio) principle that a rewrite
which changes output is a different tool.

A subtler pattern: when investigating the float precision divergence,
the LLM repeatedly suggested accepting the difference and moving on —
framing it as an inherent `long double` vs `f64` gap, proposing
tolerance-based tests, and characterizing the affected contigs as
"meaningless edge cases." This was reasonable advice on the surface,
but it would have left a multiplication-ordering discrepancy unresolved.
The Rust code computed `abundance * (length * probability)` where the
C++ evaluated `abundance * length * probability` left-to-right — both
are mathematically equivalent and neither is wrong, but IEEE 754
multiplication is not associative and the difference compounds over EM
iterations. This is not a bug in either implementation; it is an
incidental property of the original's source code that must be matched
for bit-identical output. The root cause was only found by insisting on
a bit-level investigation that the LLM initially framed as unnecessary.

Whether this tendency to satisfice is a net positive (preventing rabbit
holes) or a net negative (missing reproducibility requirements) likely
depends on context; in equivalence work, the instinct to accept "close
enough" is dangerous.

The most significant failure: the original MaxBin2's Perl orchestration
runs the EM algorithm *recursively* — each output bin is re-seeded from
the HMMER results and re-run through the EM, up to 5 levels deep. This
recursive splitting is the core of MaxBin2's bin refinement. The LLM
translated the C++ EM core faithfully, but never implemented the Perl
loop that calls it repeatedly. I did not catch this either. The
component-level equivalence tests all passed — they test each stage in
isolation, not the orchestration between stages. The gap was only
discovered the next day when end-to-end CLI equivalence tests were added
that run both tools on the same input and compare output. By that point
I had already announced the tool on the nf-core/mag Slack channel.

In hindsight, the pre-processing stages (contig filtering, abundance
parsing, seed selection) had enough subtle nondeterminism and formatting
differences to make component-level equivalence feel like the hard
problem. It was not. The hard problem was the orchestration, and
component-level tests create a false sense of completeness. End-to-end
tests should come first, not last.

LLM reviewer agents are sensitive to prompt framing — adversarial
prompts find overclaims, neutral prompts find the same gaps but frame
them as fixable. This paper was reviewed using multiple framings
([docs/review-prompts.md](docs/review-prompts.md)). More fundamentally,
I have not read every line of the generated code with the depth of
someone who wrote it by hand. The original MaxBin2 authors *wrote*
their code — they understand it in a way I cannot claim to. The FFI
equivalence tests verify that output matches, but not that I understand
*why* it matches.

### Limitations

This work has been tested on five datasets, up to the full CAMI I High
benchmark (42,038 contigs). I do not know how well it generalizes to
larger or more complex metagenomes. It has not been tested in a real
pipeline (nf-core/mag). It still shells out to HMMER, Bowtie2, and
FragGeneScan (which itself still requires Perl). No domain expert has
reviewed the reimplementation. Bin quality has not been assessed with
tools like [CheckM](https://github.com/Ecogenomics/CheckM) or
[BUSCO](https://busco.ezlab.org/) — the equivalence tests verify
identical output, not whether that output is biologically meaningful.

## What's Next

Now that equivalence is established, bug fixes can be applied as
separate, documented changes — each with its own test showing the old
(buggy) output and the new (correct) output. The `prob_threshold`
default mismatch is the most user-visible fix. Integration with the
[nf-core/mag](https://nf-co.re/mag) pipeline — the real goal — is the
next step. Longer-term, replacing FragGeneScan with
[Prodigal](https://github.com/hyattpd/Prodigal) (following
[maxbin2_custom](https://github.com/mruehlemann/maxbin2_custom)) and
adopting GTDB marker gene sets would address the upstream dependency
issues that Nix currently papers over. The EM inner loop is
embarrassingly parallel and neither implementation exploits this yet —
it is the biggest remaining performance opportunity.

More broadly, the original MaxBin2 has not seen a release since 2019.
This rewrite is, in practice, a maintenance succession: bugs can be
fixed, dependencies updated, and the tool kept working as the ecosystem
around it evolves.

## Conclusion

At the same float width, the reimplementation produces bit-for-bit
identical output to the original MaxBin2.
The only divergence is between the original's 80-bit `long double` and
Rust's 64-bit `f64`: on the full CAMI I High benchmark (42,038 contigs),
~4 contigs differ at a meaningless decision boundary. The EM core is
~8x faster and the tool is installable with a single `nix run` command. The equivalence-first methodology —
match the original's output, bugs included, then fix things separately
— made the project tractable without domain expertise. Whether that
generalizes to other bioinformatics tools in similar situations remains
to be seen. The code, tests, and this paper are available for
independent verification.

