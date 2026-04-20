#!/usr/bin/env bash
# Benchmark: FragGeneScan (C) vs FragGeneScanRs (Rust) on assembled contigs.
#
# Compares wall-clock time, output size, and predicted gene counts.
# Both tools are run with -complete=0 -train=complete (the MaxBin2 defaults).
#
# Environment variables (set by the Nix wrapper):
#   CONTIGS           — input contig FASTA file
#   TRAIN_DIR         — path to FragGeneScanRs training data directory
#   THREADS           — thread count (default: 1)
#
# Both run_FragGeneScan.pl and FragGeneScanRs must be on $PATH.

set -euo pipefail

THREADS="${THREADS:-1}"

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# Decompress if gzipped — gene callers need raw FASTA.
if [[ "$CONTIGS" == *.gz ]]; then
    gzip -dc "$CONTIGS" > "$WORK/contigs.fa"
    CONTIGS="$WORK/contigs.fa"
fi

echo "=== Gene Caller Benchmark ==="
echo ""
echo "  contigs:  $CONTIGS"
echo "  threads:  $THREADS"
echo "  contigs size: $(du -h "$CONTIGS" | cut -f1)"
num_contigs=$(grep -c '^>' "$CONTIGS")
echo "  contigs count: $num_contigs"
echo ""

# --- FragGeneScan (original C + Perl wrapper) ---
echo "--- FragGeneScan (C, v1.30) ---"
FGS_OUT="$WORK/fgs"

# FragGeneScan needs to run from its own directory (resolves paths relative to $0).
FGS_START=$(date +%s%N)
run_FragGeneScan.pl \
    "-genome=$CONTIGS" \
    "-out=$FGS_OUT" \
    "-complete=0" \
    "-train=complete" \
    "-thread=$THREADS"
FGS_END=$(date +%s%N)

FGS_MS=$(( (FGS_END - FGS_START) / 1000000 ))
FGS_FAA="$FGS_OUT.faa"
fgs_genes=$(grep -c '^>' "$FGS_FAA")
fgs_size=$(du -h "$FGS_FAA" | cut -f1)

echo "  time:  ${FGS_MS} ms"
echo "  genes: $fgs_genes"
echo "  .faa:  $fgs_size"
echo ""

# --- FragGeneScanRs (Rust) ---
echo "--- FragGeneScanRs (Rust, v1.1.0) ---"
FGSRS_OUT="$WORK/fgsrs"

FGSRS_START=$(date +%s%N)
FragGeneScanRs \
    -s "$CONTIGS" \
    -o "$FGSRS_OUT" \
    -w 0 \
    -t complete \
    -r "$TRAIN_DIR" \
    -p "$THREADS"
FGSRS_END=$(date +%s%N)

FGSRS_MS=$(( (FGSRS_END - FGSRS_START) / 1000000 ))
FGSRS_FAA="$FGSRS_OUT.faa"
fgsrs_genes=$(grep -c '^>' "$FGSRS_FAA")
fgsrs_size=$(du -h "$FGSRS_FAA" | cut -f1)

echo "  time:  ${FGSRS_MS} ms"
echo "  genes: $fgsrs_genes"
echo "  .faa:  $fgsrs_size"
echo ""

# --- Comparison ---
echo "=== Summary ==="
if [ "$FGS_MS" -gt 0 ]; then
    speedup=$(awk "BEGIN { printf \"%.2f\", $FGS_MS / $FGSRS_MS }")
    echo "  speedup: ${speedup}x"
fi
gene_diff=$(( fgsrs_genes - fgs_genes ))
echo "  gene count diff: $gene_diff (FGS=$fgs_genes, FGSrs=$fgsrs_genes)"
echo ""

# Check if .faa files have the same gene IDs (order may differ).
# Extract contig names from gene IDs and compare sets.
extract_contigs() {
    grep '^>' "$1" | sed 's/^>//' | sort
}
extract_contigs "$FGS_FAA" > "$WORK/fgs_ids.txt"
extract_contigs "$FGSRS_FAA" > "$WORK/fgsrs_ids.txt"

if diff -q "$WORK/fgs_ids.txt" "$WORK/fgsrs_ids.txt" > /dev/null 2>&1; then
    echo "  gene IDs: IDENTICAL"
else
    only_fgs=$(comm -23 "$WORK/fgs_ids.txt" "$WORK/fgsrs_ids.txt" | wc -l)
    only_fgsrs=$(comm -13 "$WORK/fgs_ids.txt" "$WORK/fgsrs_ids.txt" | wc -l)
    echo "  gene IDs: DIFFER (only in FGS: $only_fgs, only in FGSrs: $only_fgsrs)"
fi

echo ""
echo "=== Done ==="
