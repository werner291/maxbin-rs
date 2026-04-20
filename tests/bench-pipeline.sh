#!/usr/bin/env bash
# Benchmark: full maxbin-rs pipeline with per-stage timing.
#
# Runs the pipeline on a dataset with pre-computed abundance. The pipeline
# itself prints per-stage timing to stderr — this script captures that
# and adds wall-clock timing around the whole run.
#
# Environment variables (set by the Nix wrapper):
#   CONTIGS    — input contig FASTA file (may be gzipped)
#   ABUND      — pre-computed abundance file
#   THREADS    — thread count (default: 1)

set -euo pipefail

THREADS="${THREADS:-1}"

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# Decompress if gzipped.
if [[ "$CONTIGS" == *.gz ]]; then
    gzip -dc "$CONTIGS" > "$WORK/contigs.fa"
    CONTIGS="$WORK/contigs.fa"
fi

num_contigs=$(grep -c '^>' "$CONTIGS")

echo "=== maxbin-rs Pipeline Benchmark ==="
echo "  contigs: $num_contigs"
echo "  threads: $THREADS"
echo ""

export MAXBIN_RS_DETERMINISTIC=1

START=$(date +%s%N)
maxbin-rs pipeline \
    --contig "$CONTIGS" \
    --abund "$ABUND" \
    --out "$WORK/out" \
    --thread "$THREADS" 2>&1
END=$(date +%s%N)

TOTAL_MS=$(( (END - START) / 1000000 ))

echo ""
BIN_COUNT=$(find "$WORK/out" -name '*.fasta' | wc -l)
echo "  total wall-clock: ${TOTAL_MS} ms"
echo "  bins produced:    $BIN_COUNT"
echo ""
echo "=== Done ==="
