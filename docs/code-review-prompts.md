# Code Review Agent Prompts

These prompts were used during development to review code using LLM
subagents with fresh context. Each agent saw only the prompt and the
repository — no conversation history from the writing process.

## 1. Source-Tracing Comments

Add comments mapping every piece of Rust logic back to the original
C++/Perl source it was derived from.

> Your task: add comments to every Rust source file that trace each piece of logic back to the original source code it was derived from.
>
> For every non-trivial function, branch, or computation in the Rust code, add a comment like:
> - `// Matches fastaReader.cpp:516-587 (parse())`
> - `// Matches run_MaxBin.pl:366 (checkContig call)`
> - `// Matches _getmarker.pl:94-96 (coverage filter)`
>
> Rules:
> - Read BOTH the Rust file and its corresponding C++/Perl file
> - Only add comments, do not change any logic
> - Be precise with line numbers — verify them by reading the source
> - Don't over-annotate trivial things, focus on logic
> - For KNOWN BUGS being reproduced, keep existing annotations and add the source line reference if missing

## 2. Equivalence Test Coverage Audit

Audit every public function for FFI equivalence test coverage and add
proptest tests for gaps.

> Your task: audit every public function in the Rust codebase, check whether it has an FFI equivalence test (calling the original C++ via FFI and comparing outputs), and add proptest equivalence tests for any that are missing.
>
> For each public function, determine:
> - Does it have an FFI equivalence test? (either inline or in proptest)
> - If not, can we add one? (is there a corresponding FFI wrapper?)
>
> Add proptest equivalence tests for any gaps. Focus on:
> - Profiler: add_profile() and calc_profile()
> - Distance: profile-based variants with proptest
> - emanager: get_prob_dist() and get_prob_abund()
> - Any kmer_map methods not yet proptested
>
> Important:
> - Do NOT modify any existing implementation files — only add tests
> - Use proptest, not just hardcoded test cases

## 3. Proptest Edge Case Audit (First Pass)

Exhaustive audit of proptest generators for missed edge cases and
untested code paths.

> You are a meticulous test quality auditor. Your job is to find every gap, every missed edge case, every weak generator, every untested code path in the proptest equivalence tests. Do not give anything the benefit of the doubt.
>
> For each proptest file, cross-reference against the source code to find:
> 1. Code paths never exercised by any generated input
> 2. Input domain gaps (boundary values, special characters, empty inputs, NaN, infinity, etc.)
> 3. Assertion gaps (output fields not checked, tolerances too loose)
> 4. Generator biases (clustering in a boring subspace)
> 5. Missing equivalence checks for public functions
> 6. Interaction effects (bugs that only appear when components chain)
>
> Be brutal. Be thorough. Read every line of every file. Do not summarize — enumerate.
> Rank findings: CRITICAL / MODERATE / LOW.

## 4. Proptest Second-Pass Audit

After fixes from the first audit, a fresh agent verified the fixes and
checked for regressions.

> You are doing a SECOND pass. The first audit found 6 CRITICAL and 12 MODERATE issues. Most have been fixed. Verify the fixes and find anything missed or introduced.
>
> Previous CRITICAL findings were:
> 1. FASTA: headers with descriptions never generated (header truncation untested)
> 2. Abundance: headers never compared
> 3. Abundance: consecutive separators never equivalence-tested
> 4. Abundance: atof vs f64::parse divergence untested
> 5. Profiler: N characters never in proptest
> 6. emanager: prob_abund_formula_consistency was tautological
>
> Check for NEW issues: generator changes that exclude valid inputs, assertions too lenient or strict, tests that look like they test something but don't, local reimplementations tested instead of the source.

## 5. README Accuracy Audit

Cross-reference every claim in the README against the actual codebase.

> Flag anything that is overstated, no longer accurate, missing caveats, or promising things that don't exist yet.
>
> Check specifically:
> 1. "single binary" — does it actually work standalone or need HMMER/Bowtie2/FGS on PATH?
> 2. "no Perl" — true?
> 3. "no setting file" — true?
> 4. CLI flags — actually flag-compatible with the original?
> 5. The nix run command — does it work?
> 6. The NixOS config snippet — package path correct?

## 6. Nondeterminism Investigation

Investigate why the original and reimplementation produce different bin
counts on the same dataset. This agent discovered the Perl hash
iteration ordering problem.

> I need to understand exactly how seed selection works in both the original MaxBin2 and the Rust reimplementation. The issue is that the original produces 13 bins and Rust produces 10 on the same CAPES_S7 dataset, likely due to different seed selection.
>
> Key questions:
> - Where exactly in the pipeline does seed selection happen?
> - What data structure holds the marker→contig mapping?
> - Is there a sort or ordering step that could differ?
> - Could we make the Perl deterministic by sorting its hash keys?
> - Does the Rust side use deterministic ordering?

## 7. Code Quality and Organization Audit

Review for dead code, AI-generated smell, and organizational
inconsistencies. Instructed to preserve source-mapping comments.

> Phase 1: Research. Look for organizational inconsistencies, dead code,
> AI-generated smell (verbose comments restating code, unnecessary
> abstractions, copy-pasted blocks), inconsistent style, anything
> confusing.
>
> Phase 2: Conservative fixes.
> - Source-mapping comments ("Matches kmerMap.cpp:219-244") are CRITICALLY IMPORTANT — preserve them.
> - Prefer keeping code over cleaning.
> - Do not change logic or behavior.
> - Run cargo build and cargo nextest run to verify.
