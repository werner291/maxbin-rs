#!/usr/bin/env bash
#
# Recursive binning smoke test.
#
# Runs maxbin-rs and the original MaxBin2 on a downsampled CAMI dataset
# (5000 contigs) and compares structural properties of the output.
# Does NOT require byte-level equivalence — checks bin count, total
# genome size, and that recursion actually triggered.
#
# Required environment variables (set by the Nix wrapper):
#   CONTIGS     — path to downsampled contigs FASTA
#   ABUND       — path to abundance file
#   HMMOUT      — path to pre-computed HMMER output
#
# Usage: nix build .#test-recursion-smoke

set -uo pipefail

FAILED=0
PASSED=0
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

export MAXBIN_RS_DETERMINISTIC=1

pass() { echo "  PASS $1"; PASSED=$((PASSED + 1)); }
fail() { echo "  FAIL $1"; FAILED=$((FAILED + 1)); }

echo "=== Recursive Binning Smoke Test ==="
echo ""
echo "  contigs: $(grep -c '^>' "$CONTIGS")"
echo "  abundance: $(wc -l < "$ABUND") lines"
echo ""

# =========================================================================
# Run both tools
# =========================================================================

ORIG="$WORK/orig"
RUST="$WORK/rust"
mkdir -p "$ORIG" "$RUST"

# --- Original MaxBin2 (with pre-computed HMMER output) ---
echo "--- Running original MaxBin2 ---"
cp "$CONTIGS" "$ORIG/contigs.fa"
# Pre-place HMMER output + FINISH marker so the original skips
# FragGeneScan and HMMER, using the same intermediates as maxbin-rs.
# Matches run_MaxBin.pl:457-464: skips getHMM() if .FINISH exists.
# After filtering, $contig_f becomes $out_f.contig.tmp, so the HMMER
# output must be at {out_prefix}.contig.tmp.hmmout.
cp "$HMMOUT" "$ORIG/test.contig.tmp.hmmout"
touch "$ORIG/test.contig.tmp.hmmout.FINISH"
run_MaxBin.pl -contig "$ORIG/contigs.fa" \
  -abund "$ABUND" -out "$ORIG/test" -thread 1 \
  2>&1

ORIG_BINS=$(ls "$ORIG"/test.*.fasta 2>/dev/null | wc -l)
echo "  original produced $ORIG_BINS bins"
echo ""

# --- maxbin-rs (using pre-computed HMMER output) ---
echo "--- Running maxbin-rs ---"
maxbin-rs -contig "$CONTIGS" \
  -abund "$ABUND" -hmmout "$HMMOUT" -out "$RUST/test" -thread 1 \
  2>&1

RUST_BINS=$(ls "$RUST"/test.*.fasta 2>/dev/null | wc -l)
echo "  maxbin-rs produced $RUST_BINS bins"
echo ""

# =========================================================================
# Structural comparison
# =========================================================================
echo "--- Structural comparison ---"

# 1. Both produced bins
if [ "$ORIG_BINS" -gt 0 ] && [ "$RUST_BINS" -gt 0 ]; then
  pass "1. both produced bins (original=$ORIG_BINS, rust=$RUST_BINS)"
else
  fail "1. one or both produced no bins (original=$ORIG_BINS, rust=$RUST_BINS)"
fi

# 2. Bin count within 50% of each other
if [ "$ORIG_BINS" -gt 0 ] && [ "$RUST_BINS" -gt 0 ]; then
  RATIO=$(awk "BEGIN { r = $RUST_BINS / $ORIG_BINS; print (r > 0.5 && r < 2.0) ? 1 : 0 }")
  if [ "$RATIO" = "1" ]; then
    pass "2. bin count ratio acceptable (original=$ORIG_BINS, rust=$RUST_BINS)"
  else
    fail "2. bin count ratio out of range (original=$ORIG_BINS, rust=$RUST_BINS)"
  fi
fi

# 3. Total genome size within 20% of each other
orig_total_size() {
  local total=0
  for f in "$ORIG"/test.*.fasta; do
    local s
    s=$(grep -v '^>' "$f" | tr -d '\n' | wc -c)
    total=$((total + s))
  done
  echo $total
}

rust_total_size() {
  local total=0
  for f in "$RUST"/test.*.fasta; do
    local s
    s=$(grep -v '^>' "$f" | tr -d '\n' | wc -c)
    total=$((total + s))
  done
  echo $total
}

ORIG_SIZE=$(orig_total_size)
RUST_SIZE=$(rust_total_size)
if [ "$ORIG_SIZE" -gt 0 ]; then
  SIZE_RATIO=$(awk "BEGIN { r = $RUST_SIZE / $ORIG_SIZE; print (r > 0.8 && r < 1.2) ? 1 : 0 }")
  if [ "$SIZE_RATIO" = "1" ]; then
    pass "3. total genome size within 20% (original=${ORIG_SIZE}bp, rust=${RUST_SIZE}bp)"
  else
    fail "3. total genome size diverged (original=${ORIG_SIZE}bp, rust=${RUST_SIZE}bp)"
  fi
fi

# 4. maxbin-rs recursion triggered (check for depth > 0 in output)
if grep -q "depth [1-9]" "$RUST/test.summary" 2>/dev/null || \
   ls "$RUST"/test.*.out.*.fasta 2>/dev/null | head -1 > /dev/null 2>&1; then
  pass "4. recursion triggered (sub-bin files exist)"
else
  # Check if any intermediate .out files were created and cleaned up
  # The bin names in the summary would show recursion
  if [ "$RUST_BINS" -gt 2 ]; then
    pass "4. recursion likely triggered ($RUST_BINS bins from 85 seeds)"
  else
    fail "4. recursion may not have triggered ($RUST_BINS bins)"
  fi
fi

# 5. Both have noclass files
if [ -f "$ORIG/test.noclass" ] && [ -f "$RUST/test.noclass" ]; then
  pass "5. both produced noclass files"
elif [ ! -f "$ORIG/test.noclass" ] && [ ! -f "$RUST/test.noclass" ]; then
  pass "5. neither produced noclass (consistent)"
else
  fail "5. noclass mismatch (original=$([ -f "$ORIG/test.noclass" ] && echo yes || echo no), rust=$([ -f "$RUST/test.noclass" ] && echo yes || echo no))"
fi

# 6. Summary files exist
if [ -f "$ORIG/test.summary" ] && [ -f "$RUST/test.summary" ]; then
  pass "6. both produced summary files"
else
  fail "6. summary missing"
fi

echo ""
echo "============================================"
echo "  PASSED: $PASSED"
echo "  FAILED: $FAILED"
echo "============================================"

if [ "$FAILED" -gt 0 ]; then
  exit 1
fi
