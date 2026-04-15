#!/usr/bin/env bash
# End-to-end equivalence test: run original MaxBin2 and maxbin-rs on the same
# input with the same prefix, compare outputs byte-for-byte.
# Intended to be run via `nix run .#test-equivalence`.
set -uo pipefail

# Enable deterministic seed ordering so maxbin-rs output matches the
# patched original (which also sorts hash keys).
export MAXBIN_RS_DETERMINISTIC=1

if [ -z "${MAXBIN2_TEST_CONTIGS:-}" ] || [ -z "${MAXBIN2_TEST_READS1:-}" ]; then
  echo "ERROR: MAXBIN2_TEST_CONTIGS and MAXBIN2_TEST_READS1 must be set."
  exit 1
fi

TOTAL_FAILED=0

# Run a single equivalence scenario.
# Usage: run_scenario "name" "orig_args" "rust_args"
run_scenario() {
  local SCENARIO="$1"
  local ORIG_ARGS="$2"
  local RUST_ARGS="$3"

  echo ""
  echo "============================================"
  echo "  Scenario: $SCENARIO"
  echo "============================================"

  local ORIG_DIR RUST_DIR PREFIX
  ORIG_DIR=$(mktemp -d)
  RUST_DIR=$(mktemp -d)
  PREFIX="test1"

  # Copy inputs to both dirs
  cp "$MAXBIN2_TEST_CONTIGS" "$ORIG_DIR/contigs.fa.gz"
  cp "$MAXBIN2_TEST_CONTIGS" "$RUST_DIR/contigs.fa.gz"
  chmod u+w "$ORIG_DIR"/* "$RUST_DIR"/*

  # Copy reads if they exist
  for var in MAXBIN2_TEST_READS1 MAXBIN2_TEST_READS2; do
    if [ -n "${!var:-}" ]; then
      local basename
      basename=$(echo "$var" | sed 's/MAXBIN2_TEST_//' | tr '[:upper:]' '[:lower:]')
      cp "${!var}" "$ORIG_DIR/${basename}.fastq.gz"
      cp "${!var}" "$RUST_DIR/${basename}.fastq.gz"
      chmod u+w "$ORIG_DIR/${basename}.fastq.gz" "$RUST_DIR/${basename}.fastq.gz"
    fi
  done

  # HMM marker files for Rust
  if [ -n "${MAXBIN2_MARKER_HMM:-}" ]; then
    ln -sf "$MAXBIN2_MARKER_HMM" "$RUST_DIR/marker.hmm"
  fi
  if [ -n "${MAXBIN2_BACAR_MARKER_HMM:-}" ]; then
    ln -sf "$MAXBIN2_BACAR_MARKER_HMM" "$RUST_DIR/bacar_marker.hmm"
  fi

  echo "--- Running original MaxBin2 ---"
  (cd "$ORIG_DIR" && eval run_MaxBin.pl "$ORIG_ARGS" 2>&1) | tail -3

  echo "--- Running maxbin-rs ---"
  local RUST_EXIT=0
  (cd "$RUST_DIR" && eval maxbin-rs "$RUST_ARGS" 2>&1) | tail -3 || RUST_EXIT=$?
  if [ "$RUST_EXIT" != "0" ]; then
    echo "WARNING: maxbin-rs exited with code $RUST_EXIT"
  fi

  echo "--- Comparing outputs ---"
  local FAILED=0

  # Non-binned outputs: must match exactly
  for name in noclass tooshort; do
    local orig="$ORIG_DIR/$PREFIX.$name"
    local rust="$RUST_DIR/$PREFIX.$name"
    if [ ! -f "$orig" ]; then
      echo "SKIP $name"
    elif [ ! -f "$rust" ]; then
      echo "FAIL $name (not produced)"
      FAILED=1
    elif [ "$(sha256sum "$orig" | cut -d' ' -f1)" = "$(sha256sum "$rust" | cut -d' ' -f1)" ]; then
      echo "PASS $name"
    else
      echo "FAIL $name"
      diff "$orig" "$rust" | head -10
      FAILED=1
    fi
  done

  # Bin contents: compare as a set (Perl hash ordering is non-deterministic, see
  # _getmarker.pl keys %tmphash — Perl 5.18+ randomizes hash iteration order,
  # so bin numbering varies between runs of the original itself).
  local ORIG_BINS RUST_BINS
  ORIG_BINS=$(ls "$ORIG_DIR/$PREFIX".*.fasta 2>/dev/null | wc -l)
  RUST_BINS=$(ls "$RUST_DIR/$PREFIX".*.fasta 2>/dev/null | wc -l)
  if [ "$ORIG_BINS" != "$RUST_BINS" ]; then
    echo "FAIL bin count: original=$ORIG_BINS maxbin-rs=$RUST_BINS"
    FAILED=1
  else
    echo "PASS bin count ($ORIG_BINS bins)"
    local ORIG_HASHES RUST_HASHES
    ORIG_HASHES=$(for f in "$ORIG_DIR/$PREFIX".*.fasta; do sha256sum "$f" | cut -d' ' -f1; done | sort)
    RUST_HASHES=$(for f in "$RUST_DIR/$PREFIX".*.fasta; do sha256sum "$f" | cut -d' ' -f1; done | sort)
    if [ "$ORIG_HASHES" = "$RUST_HASHES" ]; then
      echo "PASS bin contents (all bin hashes match as a set)"
    else
      echo "FAIL bin contents"
      diff <(echo "$ORIG_HASHES") <(echo "$RUST_HASHES")
      FAILED=1
    fi
  fi

  # Summary and marker: compare with sorted data lines (bin order may differ)
  for name in summary marker; do
    local orig="$ORIG_DIR/$PREFIX.$name"
    local rust="$RUST_DIR/$PREFIX.$name"
    if [ ! -f "$orig" ]; then
      echo "SKIP $name"
    elif [ ! -f "$rust" ]; then
      echo "FAIL $name (not produced)"
      FAILED=1
    else
      local osorted rsorted
      osorted=$(head -1 "$orig"; tail -n+2 "$orig" | sort)
      rsorted=$(head -1 "$rust"; tail -n+2 "$rust" | sort)
      if [ "$osorted" = "$rsorted" ]; then
        echo "PASS $name"
      else
        echo "FAIL $name"
        diff <(echo "$osorted") <(echo "$rsorted") | head -10
        FAILED=1
      fi
    fi
  done

  rm -rf "$ORIG_DIR" "$RUST_DIR"

  if [ "$FAILED" = "0" ]; then
    echo "=== SCENARIO PASSED ==="
  else
    echo "=== SCENARIO FAILED ==="
    TOTAL_FAILED=$((TOTAL_FAILED + 1))
  fi
}

# ---- Scenarios ----

echo "=== End-to-end equivalence tests ==="

# Both binaries are called with IDENTICAL flags — proving backwards compatibility.
# maxbin-rs accepts the original single-dash flag syntax.

# Scenario 1: single reads file (the baseline test)
COMMON_ARGS="-contig contigs.fa.gz -reads reads1.fastq.gz -thread 2 -out test1"
run_scenario "single reads file" "$COMMON_ARGS" "$COMMON_ARGS"

# Scenario 2: two reads files (multi-sample abundance)
if [ -n "${MAXBIN2_TEST_READS2:-}" ]; then
  COMMON_ARGS="-contig contigs.fa.gz -reads reads1.fastq.gz -reads2 reads2.fastq.gz -thread 2 -out test1"
  run_scenario "two reads files (multi-sample)" "$COMMON_ARGS" "$COMMON_ARGS"
fi

# ---- Summary ----
echo ""
echo "============================================"
if [ "$TOTAL_FAILED" = "0" ]; then
  echo "  ALL SCENARIOS PASSED"
  echo "============================================"
  exit 0
else
  echo "  $TOTAL_FAILED SCENARIO(S) FAILED"
  echo "============================================"
  exit 1
fi
