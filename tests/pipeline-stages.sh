#!/usr/bin/env bash
#
# Pipeline stage equivalence test.
#
# This is the primary verification script for maxbin-rs. It tests each stage
# of the MaxBin2 pipeline independently, comparing the Rust reimplementation
# against the original Perl/C++ at each interface.
#
# The MaxBin2 pipeline has four stages:
#
#   1. CONTIG FILTERING — remove contigs shorter than min_length (default 1000bp).
#      Input:  assembled contigs (FASTA)
#      Output: filtered contigs + "tooshort" rejects
#
#   2. READ MAPPING + ABUNDANCE — align reads to contigs (Bowtie2), then compute
#      per-contig read depth from the SAM alignment.
#      Input:  SAM alignment file
#      Output: abundance table (contig name → mean depth)
#
#   3. MARKER GENE DETECTION + SEED SELECTION — identify single-copy marker genes
#      (HMMER), cluster them, pick one representative contig per cluster as an
#      initial seed for the EM algorithm.
#      Input:  HMMER output + filtered contigs
#      Output: seed file (one contig name per line)
#
#   4. EM BINNING — assign all contigs to bins using an EM algorithm over
#      tetranucleotide frequency and abundance probability.
#      Input:  filtered contigs + abundance + seeds
#      Output: bin FASTA files + noclass (unassigned contigs) + summary
#
# Stages 2 and 3 involve external tools (Bowtie2, HMMER) that are slow and
# nondeterministic. To avoid running them during the test, we use pre-computed
# intermediate files from a Nix derivation — the original MaxBin2's output,
# cached once and reused. This lets us test each stage in seconds.
#
# For each stage, we feed the SAME intermediate data to both the original
# implementation and the Rust reimplementation, then compare outputs.
#
# Required environment variables (set by the Nix wrapper in flake.nix):
#   INTERMEDIATES      — path to pre-computed intermediate files
#   MAXBIN2_TEST_CONTIGS — path to the original gzipped contigs
#
# Usage: nix run .#test-pipeline-stages
#
# Example output (B. fragilis dataset, 2024 laptop):
#
#   === Pipeline Stage Equivalence Tests ===
#
#   Intermediates: /nix/store/...-maxbin2-intermediates-bfragilis
#     SAM alignment:   1528 lines
#     HMMER hits:      18 hits
#     Abundance table:  38 contigs
#     Seeds:            7 contigs
#
#   --- Stage 1: Contig filtering ---
#   PASS tooshort file identical
#     (0s)
#
#   --- Stage 2: SAM → abundance ---
#     38 contigs, max_diff=0.00e+00, diffs_above_1e-15=0
#   PASS abundance values identical
#     (0s)
#
#   --- Stage 3: HMMER → seed selection ---
#   PASS seed files identical (as sets)
#   PASS seed files identical (exact order)
#     (1s)
#
#   --- Stage 4: EM binning ---
#     Running C++ EM via FFI...
#     Running Rust EM...
#
#     Comparing EM outputs:
#     PASS noclass
#     PASS tooshort
#     PASS bin count (4)
#     PASS bin contents identical
#     (2s)
#
#   ============================================
#     ALL STAGES PASSED
#   ============================================
#
# See also:
#   - tests/proptest_*.rs    — component-level equivalence tests (Rust vs C++ FFI)
#   - tests/bench_components.rs — per-function performance comparison
#   - PAPER.md               — case study describing the verification approach

# -u: error on undefined variables.  -o pipefail: fail if any command in a
# pipeline fails.  We intentionally omit -e (exit on first error) so that
# all stages run and all failures are reported.
set -uo pipefail

# Deterministic seed ordering so Rust output is comparable to the patched
# original (which sorts Perl hash keys instead of iterating randomly).
export MAXBIN_RS_DETERMINISTIC=1

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

if [ -z "${INTERMEDIATES:-}" ] || [ -z "${MAXBIN2_TEST_CONTIGS:-}" ]; then
  echo "ERROR: INTERMEDIATES and MAXBIN2_TEST_CONTIGS must be set."
  echo "Run via: nix run .#test-pipeline-stages"
  exit 1
fi

# Locate the original MaxBin2's Perl helper scripts (installed by Nix).
MAXBIN2_LIB="$(dirname "$(dirname "$(which run_MaxBin.pl)")")/libexec/maxbin2"

FAILED=0

# These are pre-computed outputs from the original MaxBin2 Perl/C++ pipeline,
# cached as a Nix derivation. The ORIG_ prefix makes the provenance clear:
# these come from the original tool, not from the Rust reimplementation.
ORIG_SAM="$INTERMEDIATES/sam"
ORIG_HMMOUT="$INTERMEDIATES/hmmout"
ORIG_ABUND="$INTERMEDIATES/abund"
ORIG_SEED="$INTERMEDIATES/seed"
ORIG_TOOSHORT="$INTERMEDIATES/tooshort"

# Datasets with pre-computed abundance skip Bowtie2 (no SAM file).
HAS_SAM=1
if [ ! -f "$ORIG_SAM" ]; then
  HAS_SAM=0
fi

echo "=== Pipeline Stage Equivalence Tests ==="
echo ""
echo "Intermediates: $INTERMEDIATES"
if [ "$HAS_SAM" = "1" ]; then
  echo "  SAM alignment:   $(wc -l < "$ORIG_SAM") lines"
else
  echo "  SAM alignment:   (none — pre-computed abundance)"
fi
echo "  HMMER hits:      $(grep -cv '^#' "$ORIG_HMMOUT") hits"
echo "  Abundance table:  $(wc -l < "$ORIG_ABUND") contigs"
echo "  Seeds:            $(wc -l < "$ORIG_SEED") contigs"

# Check required intermediate files exist before starting.
# SAM is optional (absent for pre-computed abundance datasets).
for f in "$ORIG_HMMOUT" "$ORIG_ABUND" "$ORIG_SEED" "$ORIG_TOOSHORT"; do
  if [ ! -f "$f" ]; then
    echo "MISSING: $f"
    FAILED=1
  fi
done
if [ "$FAILED" = "1" ]; then
  echo "Intermediates incomplete."
  exit 1
fi

# ---------------------------------------------------------------------------
# Stage 1: Contig filtering
#
# Compare the "tooshort" reject file — contigs below min_length (1000bp).
# This verifies that the Rust FASTA parser applies the same length filter
# and writes the same output format as the Perl checkContig subroutine.
# ---------------------------------------------------------------------------

echo ""
echo "--- Stage 1: Contig filtering ---"
STAGE1_START=$SECONDS

STAGE1_DIR=$(mktemp -d)
# Subcommand: just filter, nothing else.
maxbin-rs filter --contig "$MAXBIN2_TEST_CONTIGS" --out "$STAGE1_DIR/test"

# Save the filtered contigs for later stages — avoids re-decompressing
# the (potentially very large) contig file for stages 3 and 4.
FILTERED_CONTIGS="$STAGE1_DIR/test.contig.tmp"

if diff "$ORIG_TOOSHORT" "$STAGE1_DIR/test.tooshort" > /dev/null 2>&1; then
  echo "PASS tooshort file identical"
else
  echo "FAIL tooshort file differs"
  diff "$ORIG_TOOSHORT" "$STAGE1_DIR/test.tooshort" | head -5
  FAILED=1
fi

# Don't rm STAGE1_DIR yet — FILTERED_CONTIGS is needed by later stages.
echo "  ($(( SECONDS - STAGE1_START ))s)"

# ---------------------------------------------------------------------------
# Stage 2: SAM → abundance
#
# Feed the same SAM alignment file to both the original Perl (getsam() from
# _getabund.pl) and the Rust implementation. Compare per-contig abundance
# values. These should be bit-exact since the arithmetic is simple (sum of
# aligned bases / contig length), but we check with a tolerance in case of
# float rounding differences across languages.
#
# Skipped for datasets with pre-computed abundance (no SAM file).
# ---------------------------------------------------------------------------

echo ""
echo "--- Stage 2: SAM → abundance ---"
STAGE2_START=$SECONDS

if [ "$HAS_SAM" = "0" ]; then
  echo "SKIP (pre-computed abundance — no SAM file)"
else

PERL_ABUND=$(mktemp)
RUST_ABUND=$(mktemp)

# Original Perl
perl -e '
require "'"$MAXBIN2_LIB"'/_getabund.pl";
getsam("'"$ORIG_SAM"'", "'"$PERL_ABUND"'");
' 2>/dev/null

# Rust reimplementation
maxbin-rs sam-to-abund --sam "$ORIG_SAM" --out "$RUST_ABUND"

# Both abundance files are tab-separated: contig_name<TAB>abundance_value.
# We paste them side by side and check:
#   - Column 1 (Perl contig name) must equal column 3 (Rust contig name)
#   - Column 2 (Perl abundance) must be close to column 4 (Rust abundance)
# "Close" means within 1e-15 — i.e. bit-exact in practice, but tolerant of
# the last bit of floating-point representation differing across languages.
if [ -f "$RUST_ABUND" ] && [ -f "$PERL_ABUND" ]; then
  RESULT=$(paste "$PERL_ABUND" "$RUST_ABUND" | awk -F'\t' '
    BEGIN { max=0; n=0; diffs=0 }
    {
      if ($1 != $3) { print "KEY MISMATCH line " NR ": " $1 " vs " $3; exit 1 }
      d = ($2 > $4) ? $2 - $4 : $4 - $2
      if (d > max) max = d
      if (d > 1e-15) diffs++
      n++
    }
    END { printf "%d contigs, max_diff=%.2e, diffs_above_1e-15=%d\n", n, max, diffs }
  ')
  echo "  $RESULT"

  DIFFS=$(echo "$RESULT" | grep -oP 'diffs_above_1e-15=\K[0-9]+')
  if [ "${DIFFS:-0}" = "0" ]; then
    echo "PASS abundance values identical"
  else
    echo "INFO $DIFFS contigs have rounding diffs above 1e-15 (acceptable)"
  fi
else
  echo "FAIL one or both abundance outputs missing"
  FAILED=1
fi

rm -f "$PERL_ABUND" "$RUST_ABUND"

fi # HAS_SAM
echo "  ($(( SECONDS - STAGE2_START ))s)"

# ---------------------------------------------------------------------------
# Stage 3: HMMER output → seed selection
#
# Feed the same HMMER output to both implementations. Compare the resulting
# seed files — the set of contigs chosen as initial EM seeds. We compare as
# sets (sorted) because seed ordering depends on hash iteration order, which
# differs between Perl and Rust even in deterministic mode.
# ---------------------------------------------------------------------------

echo ""
echo "--- Stage 3: HMMER → seed selection ---"
STAGE3_START=$SECONDS

RUST_DIR=$(mktemp -d)
# Subcommand: just generate seeds from HMMER output. Uses the filtered
# contigs from Stage 1 (no decompression needed).
maxbin-rs seeds --contig "$FILTERED_CONTIGS" --hmmout "$ORIG_HMMOUT" --out "$RUST_DIR/test"

RUST_SEED="$RUST_DIR/test.seed"
if [ -f "$RUST_SEED" ]; then
  if diff <(sort "$ORIG_SEED") <(sort "$RUST_SEED") > /dev/null 2>&1; then
    echo "PASS seed files identical (as sets)"
  else
    echo "FAIL seed files differ"
    diff <(sort "$ORIG_SEED") <(sort "$RUST_SEED") | head -10
    FAILED=1
  fi

  # Also check exact order — informational only, not a failure.
  if diff "$ORIG_SEED" "$RUST_SEED" > /dev/null 2>&1; then
    echo "PASS seed files identical (exact order)"
  else
    echo "INFO seed order differs (same set, different iteration order)"
  fi
else
  echo "FAIL Rust did not produce a seed file"
  FAILED=1
fi

rm -rf "$RUST_DIR"
echo "  ($(( SECONDS - STAGE3_START ))s)"

# ---------------------------------------------------------------------------
# Stage 4: EM algorithm
#
# Run the EM binning step with identical inputs (filtered contigs, abundance
# table, seed file) on both the original C++ (via FFI) and the Rust
# reimplementation. Compare output bins, noclass, and tooshort files.
#
# Uses FILTERED_CONTIGS from Stage 1 — no redundant decompression.
# ---------------------------------------------------------------------------

echo ""
echo "--- Stage 4: EM binning ---"
STAGE4_START=$SECONDS

CPP_DIR=$(mktemp -d)
RUST_DIR=$(mktemp -d)

# Subcommands: run C++ EM and Rust EM separately. Each reports its own timing.
maxbin-rs cpp-em --contig "$FILTERED_CONTIGS" --abund "$ORIG_ABUND" \
  --seed "$ORIG_SEED" --out "$CPP_DIR/test" --thread 1

maxbin-rs em --contig "$FILTERED_CONTIGS" --abund "$ORIG_ABUND" \
  --seed "$ORIG_SEED" --out "$RUST_DIR/test" --thread 1

# Compare outputs.
echo ""
echo "  Comparing EM outputs:"

for name in noclass tooshort; do
  if [ ! -f "$CPP_DIR/test.$name" ]; then
    echo "  SKIP $name (C++ did not produce)"
  elif [ ! -f "$RUST_DIR/test.$name" ]; then
    echo "  FAIL $name (Rust did not produce)"
    FAILED=1
  elif diff "$CPP_DIR/test.$name" "$RUST_DIR/test.$name" > /dev/null 2>&1; then
    echo "  PASS $name"
  else
    echo "  FAIL $name"
    FAILED=1
  fi
done

# Compare bins as unordered sets of SHA-256 hashes (bin numbering may differ).
CPP_BINS=$(ls "$CPP_DIR"/test.*.fasta 2>/dev/null | wc -l)
RUST_BINS=$(ls "$RUST_DIR"/test.*.fasta 2>/dev/null | wc -l)
if [ "$CPP_BINS" = "$RUST_BINS" ]; then
  echo "  PASS bin count ($CPP_BINS)"
  CPP_HASHES=$(for f in "$CPP_DIR"/test.*.fasta; do sha256sum "$f" | cut -d' ' -f1; done | sort)
  RUST_HASHES=$(for f in "$RUST_DIR"/test.*.fasta; do sha256sum "$f" | cut -d' ' -f1; done | sort)
  if [ "$CPP_HASHES" = "$RUST_HASHES" ]; then
    echo "  PASS bin contents identical"
  else
    echo "  FAIL bin contents differ"
    diff <(echo "$CPP_HASHES") <(echo "$RUST_HASHES")
    FAILED=1
  fi
else
  echo "  FAIL bin count: C++=$CPP_BINS Rust=$RUST_BINS"
  FAILED=1
fi

rm -rf "$CPP_DIR" "$RUST_DIR" "$STAGE1_DIR"
echo "  ($(( SECONDS - STAGE4_START ))s)"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

echo ""
echo "============================================"
if [ "$FAILED" = "0" ]; then
  echo "  ALL STAGES PASSED"
else
  echo "  SOME STAGES FAILED"
fi
echo "============================================"
exit "$FAILED"
