#!/usr/bin/env bash
#
# Full pipeline trace with recursion.
#
# Runs maxbin-rs and the original MaxBin2 on the same input through the
# full recursive pipeline. Prints a structural trace of bin sizes and
# wall-clock timing for manual inspection. No pass/fail verdict.
#
# Usage: pipeline-trace.sh <contigs> <abund> <hmmout>

set -uo pipefail

if [ $# -lt 3 ]; then
  echo "Usage: $0 <contigs.fa> <abund> <hmmout>" >&2
  exit 1
fi

CONTIGS="$1"
ABUND="$2"
HMMOUT="$3"

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

export MAXBIN_RS_DETERMINISTIC=1

# Print a table of final bins: name, contigs, total bp
bin_table() {
  local dir="$1" prefix="$2"
  printf "  %-40s %8s %12s\n" "bin" "contigs" "bp"
  printf "  %-40s %8s %12s\n" "---" "-------" "--"
  local total_contigs=0 total_bp=0 n_bins=0
  for f in "$dir"/${prefix}.*.fasta; do
    [ -f "$f" ] || continue
    local name
    name=$(basename "$f")
    local contigs bp
    contigs=$(grep -c '^>' "$f")
    bp=$(grep -v '^>' "$f" | tr -d '\n' | wc -c)
    printf "  %-40s %8d %12d\n" "$name" "$contigs" "$bp"
    total_contigs=$((total_contigs + contigs))
    total_bp=$((total_bp + bp))
    n_bins=$((n_bins + 1))
  done
  echo ""
  printf "  total: %d bins, %d contigs, %d bp\n" "$n_bins" "$total_contigs" "$total_bp"

  if [ -f "$dir/${prefix}.noclass" ]; then
    local nc_contigs nc_bp
    nc_contigs=$(grep -c '^>' "$dir/${prefix}.noclass")
    nc_bp=$(grep -v '^>' "$dir/${prefix}.noclass" | tr -d '\n' | wc -c)
    printf "  noclass: %d contigs, %d bp\n" "$nc_contigs" "$nc_bp"
  fi

  if [ -f "$dir/${prefix}.tooshort" ]; then
    local ts_contigs ts_bp
    ts_contigs=$(grep -c '^>' "$dir/${prefix}.tooshort")
    ts_bp=$(grep -v '^>' "$dir/${prefix}.tooshort" | tr -d '\n' | wc -c)
    printf "  tooshort: %d contigs, %d bp\n" "$ts_contigs" "$ts_bp"
  fi
}

echo "=== Pipeline Trace ==="
echo ""
N_CONTIGS=$(zgrep -c '^>' "$CONTIGS" 2>/dev/null || grep -c '^>' "$CONTIGS")
echo "  contigs: $CONTIGS ($N_CONTIGS sequences)"
echo "  abundance: $ABUND ($(wc -l < "$ABUND") lines)"
echo "  hmmout: $HMMOUT"
echo ""

# =========================================================================
# Run original MaxBin2
# =========================================================================
ORIG="$WORK/orig"
mkdir -p "$ORIG"

echo "--- Running original MaxBin2 ---"
# Decompress contigs if gzipped — MaxBin2 detects by extension, not magic bytes
if gzip -t "$CONTIGS" 2>/dev/null; then
  zcat "$CONTIGS" > "$ORIG/contigs.fa"
else
  cp "$CONTIGS" "$ORIG/contigs.fa"
fi
cp "$HMMOUT" "$ORIG/test.contig.tmp.hmmout"
touch "$ORIG/test.contig.tmp.hmmout.FINISH"
ORIG_START=$SECONDS
run_MaxBin.pl -contig "$ORIG/contigs.fa" \
  -abund "$ABUND" -out "$ORIG/test" -thread 1 -prob_threshold 0.9 \
  2>&1 | grep --line-buffered -E '^\[trace\]|MaxBin 2\.|^Iteration' | grep --line-buffered -v 'skipped (no seeds)' | sed -u 's/^/  orig│ /'
ORIG_ELAPSED=$((SECONDS - ORIG_START))
echo ""
printf "  wall time: %dm%02ds\n\n" $((ORIG_ELAPSED / 60)) $((ORIG_ELAPSED % 60))

# =========================================================================
# Run maxbin-rs
# =========================================================================
RUST="$WORK/rust"
mkdir -p "$RUST"

echo "--- Running maxbin-rs ---"
RUST_START=$SECONDS
maxbin-rs --contig "$CONTIGS" \
  --abund "$ABUND" -hmmout "$HMMOUT" -out "$RUST/test" -thread 1 \
  2>&1 | grep --line-buffered -v -E 'candidate: marker=|^  seed: |EM iteration [0-9]|no seeds .* skipping|median ≤ 1' | sed -u 's/^/  rust│ /'
RUST_ELAPSED=$((SECONDS - RUST_START))
echo ""
printf "  wall time: %dm%02ds\n\n" $((RUST_ELAPSED / 60)) $((RUST_ELAPSED % 60))

# =========================================================================
# Output comparison
# =========================================================================
echo "==========================================="
echo "  Original MaxBin2 — final bins"
echo "==========================================="
bin_table "$ORIG" "test"

echo ""
echo "==========================================="
echo "  maxbin-rs — final bins"
echo "==========================================="
bin_table "$RUST" "test"

echo ""

for tool_dir in "$ORIG" "$RUST"; do
  label="original"
  [ "$tool_dir" = "$RUST" ] && label="maxbin-rs"
  if [ -f "$tool_dir/test.summary" ]; then
    echo "--- $label summary ---"
    sed 's/^/  /' "$tool_dir/test.summary"
    echo ""
  fi
done

echo "--- timing ---"
printf "  original: %dm%02ds\n" $((ORIG_ELAPSED / 60)) $((ORIG_ELAPSED % 60))
printf "  maxbin-rs: %dm%02ds\n" $((RUST_ELAPSED / 60)) $((RUST_ELAPSED % 60))

# =========================================================================
# Verdict: byte-identical bins?
# =========================================================================
echo ""
VERDICT="PASS"

# Check each bin by sorted hash (bin numbering may differ)
ORIG_HASHES=$(for f in "$ORIG"/test.*.fasta; do [ -f "$f" ] && sha256sum "$f" | cut -d' ' -f1; done | sort)
RUST_HASHES=$(for f in "$RUST"/test.*.fasta; do [ -f "$f" ] && sha256sum "$f" | cut -d' ' -f1; done | sort)

if [ "$ORIG_HASHES" != "$RUST_HASHES" ]; then
  echo "FAIL: bin hashes differ"
  # Find which bins don't have a matching hash on the other side
  ORIG_ONLY=$(comm -23 <(echo "$ORIG_HASHES") <(echo "$RUST_HASHES"))
  RUST_ONLY=$(comm -13 <(echo "$ORIG_HASHES") <(echo "$RUST_HASHES"))
  N_ORIG=$(echo "$ORIG_ONLY" | grep -c .)
  N_RUST=$(echo "$RUST_ONLY" | grep -c .)
  echo "  $N_ORIG bins unique to original, $N_RUST bins unique to rust"

  # For each unmatched original bin, find which rust bin has the closest
  # contig count and diff the headers
  for f in "$ORIG"/test.*.fasta; do
    h=$(sha256sum "$f" | cut -d' ' -f1)
    if echo "$RUST_ONLY" | grep -q . && ! echo "$RUST_HASHES" | grep -q "$h"; then
      bname=$(basename "$f")
      rf="$RUST/$bname"
      if [ -f "$rf" ]; then
        echo "  --- $bname headers diff ---"
        diff <(grep '^>' "$f" | sort) <(grep '^>' "$rf" | sort) | head -10
      fi
    fi
  done
  VERDICT="FAIL"
fi

# Check noclass
if ! diff -q "$ORIG/test.noclass" "$RUST/test.noclass" > /dev/null 2>&1; then
  echo "FAIL: noclass differs"
  echo "  --- noclass headers diff ---"
  diff <(grep '^>' "$ORIG/test.noclass" | sort) <(grep '^>' "$RUST/test.noclass" | sort) | head -20
  VERDICT="FAIL"
fi

if [ "$VERDICT" = "PASS" ]; then
  echo "PASS: all bins and noclass byte-identical"
else
  exit 1
fi
