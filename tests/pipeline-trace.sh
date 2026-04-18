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
echo "  contigs: $CONTIGS ($(grep -c '^>' "$CONTIGS") sequences)"
echo "  abundance: $ABUND ($(wc -l < "$ABUND") lines)"
echo "  hmmout: $HMMOUT"
echo ""

# =========================================================================
# Run original MaxBin2
# =========================================================================
ORIG="$WORK/orig"
mkdir -p "$ORIG"

echo "--- Running original MaxBin2 ---"
cp "$CONTIGS" "$ORIG/contigs.fa"
cp "$HMMOUT" "$ORIG/test.contig.tmp.hmmout"
touch "$ORIG/test.contig.tmp.hmmout.FINISH"
ORIG_START=$SECONDS
run_MaxBin.pl -contig "$ORIG/contigs.fa" \
  -abund "$ABUND" -out "$ORIG/test" -thread 1 \
  2>&1 | sed 's/^/  orig│ /'
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
maxbin-rs -contig "$CONTIGS" \
  -abund "$ABUND" -hmmout "$HMMOUT" -out "$RUST/test" -thread 1 \
  2>&1 | sed 's/^/  rust│ /'
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
