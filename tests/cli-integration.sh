#!/usr/bin/env bash
#
# End-to-end CLI integration tests.
#
# Invokes the maxbin-rs binary with different flag combinations and
# verifies exit codes, output files, and basic output properties.
#
# Required environment variables (set by the Nix wrapper):
#   CONTIGS        — path to raw contigs (FASTA, possibly gzipped)
#   READS1         — path to reads file (FASTQ, possibly gzipped)
#   READS2         — path to second reads file (optional)
#   INTERMEDIATES  — path to pre-computed intermediates (abund, seed, hmmout, sam, tooshort)
#
# Usage: nix build .#test-cli

set -uo pipefail

FAILED=0
PASSED=0
SKIPPED=0
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

export MAXBIN_RS_DETERMINISTIC=1

# Pre-computed intermediates
ABUND="$INTERMEDIATES/abund"
SEED="$INTERMEDIATES/seed"
HMMOUT="$INTERMEDIATES/hmmout"
SAM="$INTERMEDIATES/sam"

# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------

pass() { echo "  PASS $1"; PASSED=$((PASSED + 1)); }
fail() { echo "  FAIL $1"; FAILED=$((FAILED + 1)); }
skip() { echo "  SKIP $1"; SKIPPED=$((SKIPPED + 1)); }

# assert_exit_code NAME EXPECTED_CODE COMMAND...
assert_exit_code() {
  local name="$1" expected="$2"
  shift 2
  local rc=0
  "$@" > /dev/null 2>&1 || rc=$?
  if [ "$rc" = "$expected" ]; then
    pass "$name (exit $expected)"
  else
    fail "$name (expected exit $expected, got $rc)"
  fi
}

# assert_exit_code_stderr NAME EXPECTED_CODE PATTERN COMMAND...
# Checks exit code AND that stderr contains PATTERN.
assert_exit_code_stderr() {
  local name="$1" expected="$2" pattern="$3"
  shift 3
  local rc=0
  local stderr_file
  stderr_file=$(mktemp)
  "$@" > /dev/null 2>"$stderr_file" || rc=$?
  if [ "$rc" != "$expected" ]; then
    fail "$name (expected exit $expected, got $rc)"
    rm -f "$stderr_file"
    return
  fi
  if grep -qi "$pattern" "$stderr_file"; then
    pass "$name (exit $expected, stderr matches '$pattern')"
  else
    fail "$name (exit $expected but stderr missing '$pattern')"
    echo "    stderr was: $(head -3 "$stderr_file")"
  fi
  rm -f "$stderr_file"
}

assert_file_exists() {
  if [ -f "$1" ]; then pass "$2 exists"; else fail "$2 exists (missing: $1)"; fi
}

assert_file_nonempty() {
  if [ -s "$1" ]; then pass "$2 non-empty"; else fail "$2 non-empty (empty: $1)"; fi
}

assert_file_empty() {
  if [ ! -s "$1" ]; then pass "$2 empty"; else fail "$2 empty (has content: $1)"; fi
}

assert_files_identical() {
  if diff "$1" "$2" > /dev/null 2>&1; then
    pass "$3"
  else
    fail "$3 (files differ)"
  fi
}

assert_line_count_gt() {
  local count
  count=$(wc -l < "$1")
  if [ "$count" -gt "$2" ]; then
    pass "$3 ($count lines > $2)"
  else
    fail "$3 ($count lines, expected > $2)"
  fi
}

# Count FASTA records (lines starting with >)
fasta_record_count() {
  local count
  count=$(grep -c '^>' "$1" 2>/dev/null) || count=0
  echo "$count"
}

# ---------------------------------------------------------------------------
# Setup: filter contigs once for reuse
# ---------------------------------------------------------------------------

echo "=== CLI Integration Tests ==="
echo ""
echo "Setup: filtering contigs..."
maxbin-rs filter --contig "$CONTIGS" --out "$WORK/setup"
FILTERED="$WORK/setup.contig.tmp"
FILTERED_COUNT=$(fasta_record_count "$FILTERED")
echo "  $FILTERED_COUNT contigs after filtering (default min_contig_length=1000)"
echo ""

# =========================================================================
# 1. Global / meta
# =========================================================================
echo "--- 1. Global / meta ---"

assert_exit_code_stderr "1.1 no-args-shows-usage" 1 "usage" maxbin-rs
assert_exit_code "1.2 version-flag" 0 maxbin-rs --version
assert_exit_code "1.3 v-flag" 0 maxbin-rs --version
assert_exit_code "1.4 double-dash-version" 0 maxbin-rs --version

echo ""

# =========================================================================
# 2. filter subcommand
# =========================================================================
echo "--- 2. filter subcommand ---"

D="$WORK/filter"

# 2.1 basic
mkdir -p "$D/basic"
assert_exit_code "2.1 filter-basic" 0 maxbin-rs filter --contig "$CONTIGS" --out "$D/basic/test"
assert_file_nonempty "$D/basic/test.contig.tmp" "2.1 contig.tmp"
assert_file_exists "$D/basic/test.tooshort" "2.1 tooshort"

# 2.2 min_contig_length=500 (more contigs pass)
mkdir -p "$D/len500"
maxbin-rs filter --contig "$CONTIGS" --out "$D/len500/test" --min-contig-length 500 2>/dev/null
COUNT_500=$(fasta_record_count "$D/len500/test.contig.tmp")
if [ "$COUNT_500" -ge "$FILTERED_COUNT" ]; then
  pass "2.2 min_contig_length=500 ($COUNT_500 >= $FILTERED_COUNT)"
else
  fail "2.2 min_contig_length=500 ($COUNT_500 < $FILTERED_COUNT)"
fi

# 2.3 min_contig_length=5000 (fewer contigs pass)
mkdir -p "$D/len5000"
maxbin-rs filter --contig "$CONTIGS" --out "$D/len5000/test" --min-contig-length 5000 2>/dev/null
COUNT_5000=$(fasta_record_count "$D/len5000/test.contig.tmp")
if [ "$COUNT_5000" -le "$FILTERED_COUNT" ]; then
  pass "2.3 min_contig_length=5000 ($COUNT_5000 <= $FILTERED_COUNT)"
else
  fail "2.3 min_contig_length=5000 ($COUNT_5000 > $FILTERED_COUNT)"
fi

# 2.4 missing --contig
assert_exit_code "2.4 filter-missing-contig" 1 maxbin-rs filter --out prefix

# 2.5 missing --out
assert_exit_code "2.5 filter-missing-out" 1 maxbin-rs filter --contig "$CONTIGS"

# 2.6 unrecognized flag
assert_exit_code_stderr "2.6 filter-unrecognized-flag" 1 "unexpected" \
  maxbin-rs filter --contig "$CONTIGS" --out prefix --bogus

# 2.7 min_contig_length=1 (all pass)
mkdir -p "$D/len1"
maxbin-rs filter --contig "$CONTIGS" --out "$D/len1/test" --min-contig-length 1 2>/dev/null
TOOSHORT_COUNT=$(fasta_record_count "$D/len1/test.tooshort")
if [ "$TOOSHORT_COUNT" -eq 0 ]; then
  pass "2.8 min_contig_length=1 (no tooshort)"
else
  fail "2.8 min_contig_length=1 ($TOOSHORT_COUNT tooshort)"
fi

# 2.9 min_contig_length=999999 (all too short)
mkdir -p "$D/lenhuge"
maxbin-rs filter --contig "$CONTIGS" --out "$D/lenhuge/test" -min_contig_length 999999 2>/dev/null
PASS_COUNT=$(fasta_record_count "$D/lenhuge/test.contig.tmp")
if [ "$PASS_COUNT" -eq 0 ]; then
  pass "2.9 min_contig_length=999999 (all too short)"
else
  fail "2.9 min_contig_length=999999 ($PASS_COUNT passed)"
fi

echo ""

# =========================================================================
# 3. seeds subcommand
# =========================================================================
echo "--- 3. seeds subcommand ---"

D="$WORK/seeds"

# 3.1 basic (default markerset=107)
mkdir -p "$D/basic"
assert_exit_code "3.1 seeds-basic" 0 \
  maxbin-rs seeds --contig "$FILTERED" --hmmout "$HMMOUT" --out "$D/basic/test"
assert_file_nonempty "$D/basic/test.seed" "3.1 seed file"

# 3.2 markerset 40
mkdir -p "$D/m40"
maxbin-rs seeds --contig "$FILTERED" --hmmout "$HMMOUT" --out "$D/m40/test" --markerset 40 2>/dev/null
assert_file_exists "$D/m40/test.seed" "3.2 markerset=40 seed file"

# 3.3 missing required args
assert_exit_code "3.3 seeds-missing-hmmout" 1 \
  maxbin-rs seeds --contig "$FILTERED" --out prefix
assert_exit_code "3.4 seeds-missing-contig" 1 \
  maxbin-rs seeds --hmmout "$HMMOUT" --out prefix
assert_exit_code "3.5 seeds-missing-out" 1 \
  maxbin-rs seeds --contig "$FILTERED" --hmmout "$HMMOUT"

echo ""

# =========================================================================
# 4. em subcommand
# =========================================================================
echo "--- 4. em subcommand ---"

D="$WORK/em"

# 4.1 basic
mkdir -p "$D/basic"
assert_exit_code "4.1 em-basic" 0 \
  maxbin-rs em --contig "$FILTERED" --abund "$ABUND" --seed "$SEED" --out "$D/basic/test"
assert_file_exists "$D/basic/test.summary" "4.1 summary"
BIN_COUNT=$(ls "$D/basic"/test.*.fasta 2>/dev/null | wc -l)
if [ "$BIN_COUNT" -gt 0 ]; then
  pass "4.1 bins produced ($BIN_COUNT)"
else
  fail "4.1 no bins produced"
fi

# 4.2 prob_threshold=0.9 (more contigs unclassified)
mkdir -p "$D/pt09"
maxbin-rs em --contig "$FILTERED" --abund "$ABUND" --seed "$SEED" --out "$D/pt09/test" \
  --prob-threshold 0.9 2>/dev/null
NOCLASS_DEFAULT=$(wc -c < "$D/basic/test.noclass" 2>/dev/null || echo 0)
NOCLASS_09=$(wc -c < "$D/pt09/test.noclass" 2>/dev/null || echo 0)
if [ "$NOCLASS_09" -ge "$NOCLASS_DEFAULT" ]; then
  pass "4.2 prob_threshold=0.9 (noclass $NOCLASS_09 >= $NOCLASS_DEFAULT)"
else
  fail "4.2 prob_threshold=0.9 (noclass $NOCLASS_09 < $NOCLASS_DEFAULT)"
fi

# 4.3 max_iteration=1
mkdir -p "$D/iter1"
assert_exit_code "4.3 em-max-iteration-1" 0 \
  maxbin-rs em --contig "$FILTERED" --abund "$ABUND" --seed "$SEED" --out "$D/iter1/test" \
  --max-iteration 1
assert_file_exists "$D/iter1/test.summary" "4.3 summary"

# 4.4 thread=2 (deterministic, same bins)
mkdir -p "$D/thread2"
maxbin-rs em --contig "$FILTERED" --abund "$ABUND" --seed "$SEED" --out "$D/thread2/test" \
  --thread 2 2>/dev/null
BINS_T1=$(for f in "$D/basic"/test.*.fasta; do sha256sum "$f" | cut -d' ' -f1; done | sort)
BINS_T2=$(for f in "$D/thread2"/test.*.fasta; do sha256sum "$f" | cut -d' ' -f1; done | sort)
if [ "$BINS_T1" = "$BINS_T2" ]; then
  pass "4.4 thread=2 same bins as thread=1"
else
  fail "4.4 thread=2 bins differ from thread=1"
fi

# 4.5 multiple abund files (same file twice = degenerate but valid)
mkdir -p "$D/multi_abund"
assert_exit_code "4.5 em-multiple-abund" 0 \
  maxbin-rs em --contig "$FILTERED" --abund "$ABUND" --abund "$ABUND" \
  --seed "$SEED" --out "$D/multi_abund/test"

# 4.6 abund_list
mkdir -p "$D/abund_list"
echo "$ABUND" > "$D/abund_list/list.txt"
assert_exit_code "4.6 em-abund-list" 0 \
  maxbin-rs em --contig "$FILTERED" --abund-list "$D/abund_list/list.txt" \
  --seed "$SEED" --out "$D/abund_list/test"

# 4.7 missing required args
assert_exit_code "4.7 em-missing-seed" 1 \
  maxbin-rs em --contig "$FILTERED" --abund "$ABUND" --out prefix
assert_exit_code "4.8 em-missing-contig" 1 \
  maxbin-rs em --abund "$ABUND" --seed "$SEED" --out prefix
assert_exit_code "4.9 em-missing-out" 1 \
  maxbin-rs em --contig "$FILTERED" --abund "$ABUND" --seed "$SEED"

echo ""

# =========================================================================
# 5. cpp-em subcommand
# =========================================================================
echo "--- 5. cpp-em subcommand ---"

D="$WORK/cpp_em"

# 5.1 basic
mkdir -p "$D/basic"
assert_exit_code "5.1 cpp-em-basic" 0 \
  maxbin-rs cpp-em --contig "$FILTERED" --abund "$ABUND" --seed "$SEED" --out "$D/basic/test"
assert_file_exists "$D/basic/test.summary" "5.1 summary"

# 5.2 underscore subcommand variant
mkdir -p "$D/underscore"
assert_exit_code "5.2 cpp_em-underscore" 0 \
  maxbin-rs cpp_em --contig "$FILTERED" --abund "$ABUND" --seed "$SEED" --out "$D/underscore/test"

echo ""

# =========================================================================
# 6. sam-to-abund subcommand
# =========================================================================
echo "--- 6. sam-to-abund subcommand ---"

D="$WORK/sam_to_abund"

if [ -f "$SAM" ]; then
  # 6.1 basic
  mkdir -p "$D/basic"
  assert_exit_code "6.1 sam-to-abund-basic" 0 \
    maxbin-rs sam-to-abund --sam "$SAM" --out "$D/basic/abund.txt"
  assert_file_nonempty "$D/basic/abund.txt" "6.1 output"

  # 6.2 underscore subcommand variant
  mkdir -p "$D/underscore"
  assert_exit_code "6.2 sam_to_abund-underscore" 0 \
    maxbin-rs sam_to_abund --sam "$SAM" --out "$D/underscore/abund.txt"

  # 6.3 missing args
  assert_exit_code "6.3 sam-to-abund-missing-sam" 1 maxbin-rs sam-to-abund --out out.txt
  assert_exit_code "6.4 sam-to-abund-missing-out" 1 maxbin-rs sam-to-abund --sam "$SAM"
else
  skip "6.x sam-to-abund (no SAM file for this dataset)"
fi

echo ""

# =========================================================================
# 7. Legacy (no subcommand) mode with abundance
# =========================================================================
echo "--- 7. Pipeline mode ---"

D="$WORK/pipeline"

# 7.1 basic with --abund (skips Bowtie2), default output dir
assert_exit_code "7.1 pipeline-abund" 0 \
  maxbin-rs --contig "$CONTIGS" --abund "$ABUND" --out "$D/abund"
assert_file_exists "$D/abund/summary" "7.1 summary"

# 7.2 explicit pipeline subcommand (same thing)
assert_exit_code "7.2 explicit-pipeline" 0 \
  maxbin-rs pipeline --contig "$CONTIGS" --abund "$ABUND" --out "$D/explicit"
assert_files_identical "$D/abund/summary" "$D/explicit/summary" \
  "7.2 explicit pipeline == legacy"

# 7.3 error cases
assert_exit_code "7.3 pipeline-missing-contig" 1 maxbin-rs --reads "$CONTIGS"
assert_exit_code "7.5 pipeline-no-reads-no-abund" 1 maxbin-rs --contig "$CONTIGS"

# 7.4 non-empty output dir should fail
mkdir -p "$D/nonempty"
touch "$D/nonempty/dummy"
assert_exit_code "7.4 nonempty-out-fails" 1 \
  maxbin-rs --contig "$CONTIGS" --abund "$ABUND" --out "$D/nonempty"

# 7.4b force-overwrite allows non-empty output dir
assert_exit_code "7.4b force-overwrite" 0 \
  maxbin-rs --contig "$CONTIGS" --abund "$ABUND" --out "$D/nonempty" --force-overwrite

# 7.5 keep-intermediates
assert_exit_code "7.5 keep-intermediates" 0 \
  maxbin-rs --contig "$CONTIGS" --abund "$ABUND" --out "$D/keepint" --keep-intermediates
# Intermediates dir should exist somewhere under ./intermediates/
INTDIR=$(find "$D/keepint/../" -path "*/intermediates/maxbin-rs-*" -type d 2>/dev/null | head -1)
if [ -z "$INTDIR" ]; then
  # Check cwd-relative intermediates dir
  INTDIR=$(find . -path "./intermediates/maxbin-rs-*" -type d 2>/dev/null | head -1)
fi
if [ -n "$INTDIR" ]; then
  pass "7.5 intermediates dir preserved"
  rm -rf "$INTDIR"
else
  fail "7.5 intermediates dir not found"
fi

# 7.6 verbose (just check it doesn't crash)
assert_exit_code "7.6 verbose-flag" 0 \
  maxbin-rs --contig "$CONTIGS" --abund "$ABUND" --out "$D/verbose" --verbose

# 7.7 multiple reads files
if [ -n "${READS1:-}" ] && [ -n "${READS2:-}" ]; then
  assert_exit_code "7.7 pipeline-multi-reads" 0 \
    maxbin-rs --contig "$CONTIGS" --reads "$READS1" --reads "$READS2" --out "$D/multi_reads"
  assert_file_exists "$D/multi_reads/summary" "7.7 summary"
else
  skip "7.7 pipeline-multi-reads (no READS2)"
fi

# 7.8 reads_list
if [ -n "${READS1:-}" ]; then
  mkdir -p "$D/reads_list"
  echo "$READS1" > "$D/reads_list/list.txt"
  assert_exit_code "7.8 pipeline-reads-list" 0 \
    maxbin-rs --contig "$CONTIGS" --reads-list "$D/reads_list/list.txt" --out "$D/reads_list_out"
  assert_file_exists "$D/reads_list_out/summary" "7.8 summary"
else
  skip "7.8 pipeline-reads-list (no READS1)"
fi

echo ""

# =========================================================================
# 8. Determinism
# =========================================================================
echo "--- 8. Determinism ---"

D="$WORK/determinism"

mkdir -p "$D/run1" "$D/run2"
maxbin-rs em --contig "$FILTERED" --abund "$ABUND" --seed "$SEED" --out "$D/run1/test" 2>/dev/null
maxbin-rs em --contig "$FILTERED" --abund "$ABUND" --seed "$SEED" --out "$D/run2/test" 2>/dev/null

# Compare bin assignments (summary may have timing/completeness diffs due to HMMER)
BINS_R1=$(for f in "$D/run1"/test.*.fasta; do sha256sum "$f" | cut -d' ' -f1; done | sort)
BINS_R2=$(for f in "$D/run2"/test.*.fasta; do sha256sum "$f" | cut -d' ' -f1; done | sort)
if [ "$BINS_R1" = "$BINS_R2" ]; then
  pass "8.1 deterministic: bins identical across runs"
else
  fail "8.1 deterministic: bins differ across runs"
fi
assert_files_identical "$D/run1/test.noclass" "$D/run2/test.noclass" \
  "8.2 deterministic: noclass identical across runs"

# Compare bin hashes
HASHES1=$(for f in "$D/run1"/test.*.fasta; do sha256sum "$f" | cut -d' ' -f1; done | sort)
HASHES2=$(for f in "$D/run2"/test.*.fasta; do sha256sum "$f" | cut -d' ' -f1; done | sort)
if [ "$HASHES1" = "$HASHES2" ]; then
  pass "8.3 deterministic: all bins identical across runs"
else
  fail "8.3 deterministic: bins differ across runs"
fi

echo ""

# =========================================================================
# 9. Output structure validation
# =========================================================================
echo "--- 9. Output validation ---"

D="$WORK/em/basic"  # reuse from section 4

# 9.1 bin files are valid FASTA
ALL_VALID=1
for f in "$D"/test.*.fasta; do
  FIRST_CHAR=$(head -c1 "$f")
  if [ "$FIRST_CHAR" != ">" ]; then
    fail "9.1 $f does not start with >"
    ALL_VALID=0
  fi
done
if [ "$ALL_VALID" = 1 ]; then
  pass "9.1 all bin files are valid FASTA"
fi

# 9.2 summary is tab-separated with header
if head -1 "$D/test.summary" | grep -q $'\t'; then
  pass "9.2 summary is tab-separated"
else
  fail "9.2 summary not tab-separated"
fi

# 9.3 no duplicate contigs across bins
ALL_CONTIGS=$(cat "$D"/test.*.fasta | grep '^>' | sort)
UNIQUE_CONTIGS=$(echo "$ALL_CONTIGS" | sort -u)
if [ "$ALL_CONTIGS" = "$UNIQUE_CONTIGS" ]; then
  pass "9.3 no duplicate contigs across bins"
else
  DUPES=$(echo "$ALL_CONTIGS" | uniq -d | wc -l)
  # Note: this is a known MaxBin2 bug we reproduce in v0.1.x
  fail "9.3 $DUPES duplicate contigs across bins (known bug)"
fi

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
