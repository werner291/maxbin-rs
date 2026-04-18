#!/usr/bin/env bash
#
# CLI-level equivalence tests.
#
# Runs the original MaxBin2 (run_MaxBin.pl) and maxbin-rs with the same
# flags, then compares output bins. Tests that the Rust tool is a drop-in
# replacement at the CLI level, not just at the component level.
#
# Required environment variables (set by the Nix wrapper):
#   CONTIGS        — path to raw contigs (FASTA, possibly gzipped)
#   READS1         — path to reads file (optional, for reads-based tests)
#   INTERMEDIATES  — path to pre-computed intermediates (abund)
#
# Usage: nix build .#test-cli-equivalence

set -uo pipefail

FAILED=0
PASSED=0
SKIPPED=0
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

export MAXBIN_RS_DETERMINISTIC=1

ABUND="$INTERMEDIATES/abund"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

pass() { echo "  PASS $1"; PASSED=$((PASSED + 1)); }
fail() { echo "  FAIL $1"; FAILED=$((FAILED + 1)); }
skip() { echo "  SKIP $1"; SKIPPED=$((SKIPPED + 1)); }

# Compare bin FASTA files as unordered sets of content hashes.
# Bin numbering may differ, so we compare content, not filenames.
compare_bins() {
  local dir1="$1" dir2="$2" label="$3"

  local bins1 bins2
  bins1=$(ls "$dir1"/*.fasta 2>/dev/null | wc -l)
  bins2=$(ls "$dir2"/*.fasta 2>/dev/null | wc -l)

  if [ "$bins1" != "$bins2" ]; then
    fail "$label bin count (original=$bins1, rust=$bins2)"
    return
  fi

  if [ "$bins1" -eq 0 ]; then
    pass "$label (both produced 0 bins)"
    return
  fi

  local hashes1 hashes2
  hashes1=$(for f in "$dir1"/*.fasta; do sha256sum "$f" | cut -d' ' -f1; done | sort)
  hashes2=$(for f in "$dir2"/*.fasta; do sha256sum "$f" | cut -d' ' -f1; done | sort)

  if [ "$hashes1" = "$hashes2" ]; then
    pass "$label ($bins1 bins, identical)"
  else
    fail "$label (bin contents differ)"
    echo "    original hashes: $(echo "$hashes1" | tr '\n' ' ')"
    echo "    rust hashes:     $(echo "$hashes2" | tr '\n' ' ')"
  fi
}

# Compare two text files, ignoring trailing whitespace differences.
compare_files() {
  local f1="$1" f2="$2" label="$3"
  if diff <(sed 's/[[:space:]]*$//' "$f1") <(sed 's/[[:space:]]*$//' "$f2") > /dev/null 2>&1; then
    pass "$label"
  else
    fail "$label"
    diff <(sed 's/[[:space:]]*$//' "$f1") <(sed 's/[[:space:]]*$//' "$f2") | head -10
  fi
}

echo "=== CLI-Level Equivalence Tests ==="
echo ""

# =========================================================================
# 1. Full pipeline with pre-computed abundance
# =========================================================================
# Both tools run with -prob_threshold 0.9 (the documented default).
# The original's code defaults to 0.5 but we use 0.9; passing it
# explicitly ensures both tools use the same threshold.
echo "--- 1. Full pipeline, -prob_threshold 0.9, -abund ---"

ORIG="$WORK/1_orig"
RUST="$WORK/1_rust"
mkdir -p "$ORIG" "$RUST"

# The original MaxBin2 needs contigs copied to a writable location
# (it writes temp files next to the input).
cp "$CONTIGS" "$WORK/contigs_input.fa.gz"

echo "  Running original MaxBin2..."
run_MaxBin.pl -contig "$WORK/contigs_input.fa.gz" \
  -abund "$ABUND" -out "$ORIG/test" -thread 1 \
  -prob_threshold 0.9 \
  > "$ORIG/stdout.log" 2>&1 || true

echo "  Running maxbin-rs..."
maxbin-rs --contig "$CONTIGS" \
  --abund "$ABUND" -out "$RUST/test" -thread 1 \
  2> "$RUST/stderr.log" || true

compare_bins "$ORIG" "$RUST" "1.1 bins"

# Compare noclass
if [ -f "$ORIG/test.noclass" ] && [ -f "$RUST/test.noclass" ]; then
  compare_files "$ORIG/test.noclass" "$RUST/test.noclass" "1.2 noclass"
elif [ ! -f "$ORIG/test.noclass" ] && [ ! -f "$RUST/test.noclass" ]; then
  pass "1.2 noclass (both absent)"
else
  fail "1.2 noclass (one present, one absent)"
fi

# Compare tooshort
if [ -f "$ORIG/test.tooshort" ] && [ -f "$RUST/test.tooshort" ]; then
  compare_files "$ORIG/test.tooshort" "$RUST/test.tooshort" "1.3 tooshort"
elif [ ! -f "$ORIG/test.tooshort" ] && [ ! -f "$RUST/test.tooshort" ]; then
  pass "1.3 tooshort (both absent)"
else
  fail "1.3 tooshort (one present, one absent)"
fi

echo ""

# =========================================================================
# 2. Different min_contig_length
# =========================================================================
echo "--- 2. min_contig_length=500 ---"

ORIG="$WORK/2_orig"
RUST="$WORK/2_rust"
mkdir -p "$ORIG" "$RUST"

echo "  Running original MaxBin2..."
run_MaxBin.pl -contig "$WORK/contigs_input.fa.gz" \
  -abund "$ABUND" -out "$ORIG/test" -thread 1 \
  -min_contig_length 500 -prob_threshold 0.9 \
  > "$ORIG/stdout.log" 2>&1 || true

echo "  Running maxbin-rs..."
maxbin-rs --contig "$CONTIGS" \
  --abund "$ABUND" -out "$RUST/test" -thread 1 \
  --min_contig_length 500 \
  2> "$RUST/stderr.log" || true

compare_bins "$ORIG" "$RUST" "2.1 bins (min_contig_length=500)"

echo ""

# =========================================================================
# 3. prob_threshold=0.5 (original's code default, backward compat)
# =========================================================================
echo "--- 3. prob_threshold=0.5 (original's buggy default) ---"

ORIG="$WORK/3_orig"
RUST="$WORK/3_rust"
mkdir -p "$ORIG" "$RUST"

echo "  Running original MaxBin2..."
run_MaxBin.pl -contig "$WORK/contigs_input.fa.gz" \
  -abund "$ABUND" -out "$ORIG/test" -thread 1 \
  -prob_threshold 0.5 \
  > "$ORIG/stdout.log" 2>&1 || true

echo "  Running maxbin-rs..."
maxbin-rs --contig "$CONTIGS" \
  --abund "$ABUND" -out "$RUST/test" -thread 1 \
  --prob_threshold 0.5 \
  2> "$RUST/stderr.log" || true

compare_bins "$ORIG" "$RUST" "3.1 bins (prob_threshold=0.5)"

if [ -f "$ORIG/test.noclass" ] && [ -f "$RUST/test.noclass" ]; then
  compare_files "$ORIG/test.noclass" "$RUST/test.noclass" "3.2 noclass (prob_threshold=0.5)"
fi

echo ""

# =========================================================================
# 4. markerset=40
# =========================================================================
echo "--- 4. markerset=40 ---"

ORIG="$WORK/4_orig"
RUST="$WORK/4_rust"
mkdir -p "$ORIG" "$RUST"

echo "  Running original MaxBin2..."
run_MaxBin.pl -contig "$WORK/contigs_input.fa.gz" \
  -abund "$ABUND" -out "$ORIG/test" -thread 1 \
  -markerset 40 -prob_threshold 0.9 \
  > "$ORIG/stdout.log" 2>&1 || true

echo "  Running maxbin-rs..."
maxbin-rs --contig "$CONTIGS" \
  --abund "$ABUND" -out "$RUST/test" -thread 1 \
  --markerset 40 \
  2>&1 || true
echo "  maxbin-rs output files: $(ls "$RUST"/test.* 2>/dev/null || echo '(none)')"

compare_bins "$ORIG" "$RUST" "4.1 bins (markerset=40)"

echo ""

# =========================================================================
# 5. With reads (Bowtie2 mapping)
# =========================================================================
echo "--- 5. With reads (full pipeline including Bowtie2) ---"

if [ -n "${READS1:-}" ]; then
  ORIG="$WORK/5_orig"
  RUST="$WORK/5_rust"
  mkdir -p "$ORIG" "$RUST"

  # Copy reads to writable location for the original
  cp "$READS1" "$WORK/reads1_input.fastq.gz"

  echo "  Running original MaxBin2..."
  run_MaxBin.pl -contig "$WORK/contigs_input.fa.gz" \
    --reads "$WORK/reads1_input.fastq.gz" \
    --out "$ORIG/test" -thread 1 -prob_threshold 0.9 \
    > "$ORIG/stdout.log" 2>&1 || true

  echo "  Running maxbin-rs..."
  maxbin-rs --contig "$CONTIGS" \
    --reads "$READS1" \
    --out "$RUST/test" -thread 1 \
    2> "$RUST/stderr.log" || true

  compare_bins "$ORIG" "$RUST" "5.1 bins (reads mode)"
else
  skip "5.x reads mode (no READS1)"
fi

echo ""

# =========================================================================
# 6. max_iteration=1
# =========================================================================
echo "--- 6. max_iteration=1 ---"

ORIG="$WORK/6_orig"
RUST="$WORK/6_rust"
mkdir -p "$ORIG" "$RUST"

echo "  Running original MaxBin2..."
run_MaxBin.pl -contig "$WORK/contigs_input.fa.gz" \
  -abund "$ABUND" -out "$ORIG/test" -thread 1 \
  -max_iteration 1 -prob_threshold 0.9 \
  > "$ORIG/stdout.log" 2>&1 || true

echo "  Running maxbin-rs..."
maxbin-rs --contig "$CONTIGS" \
  --abund "$ABUND" -out "$RUST/test" -thread 1 \
  --max_iteration 1 \
  2> "$RUST/stderr.log" || true

compare_bins "$ORIG" "$RUST" "6.1 bins (max_iteration=1)"

echo ""

# =========================================================================
# Summary
# =========================================================================
echo "============================================"
echo "  PASSED: $PASSED"
echo "  FAILED: $FAILED"
echo "  SKIPPED: $SKIPPED"
echo "============================================"

if [ "$FAILED" -gt 0 ]; then
  exit 1
fi
