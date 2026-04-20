# Test and benchmark wrappers.
#
# Each entry here creates a small shell script that:
#   1. Sets environment variables pointing to Nix-cached test data
#   2. Runs the actual test script from tests/
#
# This indirection exists because the test scripts are plain bash (readable
# without knowing Nix), but they need paths to data in the Nix store —
# which only Nix knows. These wrappers bridge that gap.
#
# Component benchmark:
#   nix run .#bench-components              — per-function Rust vs C++ timing

{
  writeShellApplication,
  runCommand,
  binutils-unwrapped,
  maxbin2,
  maxbin2-f64,
  maxbin2-trace,
  maxbin2-f64-trace,
  maxbin-rs,
  maxbin-rs-lto,
  fraggenescan,
  fraggenescan-rs,
  rust,
  cargo-nextest,
  datasets,
  intermediates,
}:

let
  # Helper: create a sandboxed CLI integration test for a given dataset.
  # Runs tests/cli-integration.sh with the binary on PATH and test data
  # available. Tests flag parsing, subcommands, error cases, and output
  # validation end-to-end.
  mkCliTest =
    {
      name,
      dataset,
      intermediates',
    }:
    runCommand "test-cli-${name}"
      {
        nativeBuildInputs = [
          maxbin2
          maxbin-rs
        ];
      }
      ''
        export CONTIGS="${dataset.contigs}"
        export READS1="${dataset.reads1}"
        ${if dataset ? reads2 then ''export READS2="${dataset.reads2}"'' else ""}
        export INTERMEDIATES="${intermediates'}"
        export MAXBIN_RS_DETERMINISTIC=1
        export HOME=$(mktemp -d)

        bash ${../tests/cli-integration.sh}

        touch $out
      '';

  # Helper: create a sandboxed CLI equivalence test for a given dataset.
  # Runs both maxbin-rs and the original MaxBin2 with the same flags,
  # then compares output bins. Tests drop-in replacement at the CLI level.
  mkCliEquivTest =
    {
      name,
      dataset,
      intermediates',
    }:
    runCommand "test-cli-equivalence-${name}"
      {
        nativeBuildInputs = [
          maxbin2
          maxbin-rs
        ];
      }
      ''
        export CONTIGS="${dataset.contigs}"
        ${if dataset ? reads1 then ''export READS1="${dataset.reads1}"'' else ""}
        export INTERMEDIATES="${intermediates'}"
        export MAXBIN_RS_DETERMINISTIC=1
        export HOME=$(mktemp -d)

        bash ${../tests/cli-equivalence.sh}

        touch $out
      '';

  # Pipeline trace — runs both maxbin-rs and original MaxBin2 on the same
  # input through the full recursive pipeline, and prints a structural
  # trace of bin sizes across recursion depths. No pass/fail verdict.
  #
  # Factory: mkPipelineTrace { name, contigs, intermediates', maxbin2' }
  # Defaults to the f64-patched (double) MaxBin2 with trace logging — this is
  # the canonical comparison since Rust f64 matches C++ double exactly.
  # Use maxbin2-trace for the long-double comparison (inherent ~4 contig gap).
  #
  #   nix run .#trace-cami-small           (downsampled CAMI, C++ double)
  #   nix run .#trace-cami-small-ldouble   (same, C++ long double)
  #   nix run .#trace-cami                 (full CAMI I High, C++ double)
  mkPipelineTrace =
    {
      name,
      contigs,
      intermediates',
      maxbin2' ? maxbin2-f64-trace,
    }:
    writeShellApplication {
      name = "trace-${name}";
      runtimeInputs = [
        maxbin2'
        maxbin-rs
      ];
      text = ''
        bash ${../tests/pipeline-trace.sh} \
          "${contigs}" \
          "${intermediates'}/abund" \
          "${intermediates'}/hmmout"
      '';
    };
in
{
  # CLI integration tests — end-to-end flag parsing, subcommands, error cases.
  test-cli-bfragilis = mkCliTest {
    name = "bfragilis";
    dataset = datasets.bfragilis;
    intermediates' = intermediates.bfragilis;
  };

  # CLI equivalence tests — run both tools, compare output.
  test-cli-equivalence-bfragilis = mkCliEquivTest {
    name = "bfragilis";
    dataset = datasets.bfragilis;
    intermediates' = intermediates.bfragilis;
  };

  test-cli-equivalence-capes = mkCliEquivTest {
    name = "capes";
    dataset = datasets.capes-s7;
    intermediates' = intermediates.capes;
  };

  trace-cami-small = mkPipelineTrace {
    name = "cami-small";
    contigs = "${intermediates.cami-small}/contigs.fa";
    intermediates' = intermediates.cami-small;
  };

  trace-cami-small-ldouble = mkPipelineTrace {
    name = "cami-small-ldouble";
    contigs = "${intermediates.cami-small}/contigs.fa";
    intermediates' = intermediates.cami-small;
    maxbin2' = maxbin2-trace;
  };

  trace-cami = mkPipelineTrace {
    name = "cami";
    contigs = "${datasets.cami-i-high.contigs}";
    intermediates' = intermediates.cami;
  };

  # Precision divergence test — runs the Rust EM and C++ EM (via FFI)
  # on a 14-contig sub-bin from CAMI I High where the long double vs f64
  # precision difference causes different classification. Documents and
  # asserts the known divergence.
  #   nix build .#test-precision-divergence
  test-precision-divergence =
    runCommand "test-precision-divergence"
      {
        nativeBuildInputs = [
          rust
          cargo-nextest
        ];
      }
      ''
        export MAXBIN2_SRC_TARBALL="${maxbin-rs.MAXBIN2_SRC_TARBALL}"
        cp -r ${../.}/* .
        chmod -R u+w .
        cargo nextest run emanager_precision_divergence --no-fail-fast
        touch $out
      '';

  # LTO A/B benchmark: runs the CAMI I High EM three ways (C++ baseline,
  # C++ with -flto, Rust) to test whether link-time optimization explains
  # the ~8x performance gap.
  bench-cpp-lto = writeShellApplication {
    name = "bench-cpp-lto";
    runtimeInputs = [
      maxbin2
      maxbin-rs
      maxbin-rs-lto
    ];
    text = ''
      echo "=== C++ LTO Benchmark: CAMI I High EM ==="
      echo ""

      INTERMEDIATES="${intermediates.cami}"
      CONTIGS="${datasets.cami-i-high.contigs}"

      WORK=$(mktemp -d)
      trap 'rm -rf "$WORK"' EXIT

      echo "--- Filtering contigs ---"
      "${maxbin-rs}/bin/maxbin-rs" filter -contig "$CONTIGS" -out "$WORK/test"
      FILTERED="$WORK/test.contig.tmp"

      ABUND="$INTERMEDIATES/abund"
      SEED="$INTERMEDIATES/seed"

      echo ""
      echo "--- C++ EM (baseline, no LTO) ---"
      CPP_DIR=$(mktemp -d)
      MAXBIN_RS_DETERMINISTIC=1 "${maxbin-rs}/bin/maxbin-rs" cpp-em \
        -contig "$FILTERED" -abund "$ABUND" -seed "$SEED" \
        -out "$CPP_DIR/test" -thread 1
      rm -rf "$CPP_DIR"

      echo ""
      echo "--- C++ EM (with -flto) ---"
      LTO_DIR=$(mktemp -d)
      MAXBIN_RS_DETERMINISTIC=1 "${maxbin-rs-lto}/bin/maxbin-rs" cpp-em \
        -contig "$FILTERED" -abund "$ABUND" -seed "$SEED" \
        -out "$LTO_DIR/test" -thread 1
      rm -rf "$LTO_DIR"

      echo ""
      echo "--- Rust EM ---"
      RUST_DIR=$(mktemp -d)
      MAXBIN_RS_DETERMINISTIC=1 "${maxbin-rs}/bin/maxbin-rs" em \
        -contig "$FILTERED" -abund "$ABUND" -seed "$SEED" \
        -out "$RUST_DIR/test" -thread 1
      rm -rf "$RUST_DIR"

      echo ""
      echo "=== Done ==="
    '';
  };

  # Per-function performance comparison: runs each component (profiler,
  # distance metrics, sort, etc.) in a timed loop, measuring Rust throughput
  # against the original C++ via FFI. Produces a comparison table.
  bench-components = writeShellApplication {
    name = "bench-components";
    runtimeInputs = [
      rust
      cargo-nextest
    ];
    text = ''
      export MAXBIN2_TEST_CONTIGS="${datasets.bfragilis.contigs}"
      # Copy source to a writable temp dir (Nix store is read-only).
      WORK=$(mktemp -d)
      cp -r ${../.}/* "$WORK/"
      chmod -R u+w "$WORK"
      cd "$WORK"
      cargo test --profile bench bench_components -- --ignored --nocapture
      rm -rf "$WORK"
    '';
  };

  # =====================================================================
  # Disassembly of EM hot functions — reproducible assembly comparison
  # =====================================================================
  #
  # Extracts disassembly of the EM hot-path functions from both the
  # standard and unity (LTO) builds, plus instruction count summaries.
  #
  #   nix build .#disasm-em
  #
  # Output ($out/):
  #   standard/cpp-run_EM.s        — C++ EManager::run_EM(int)
  #   standard/cpp-get_prob_abund.s — C++ EManager::get_prob_abund(double, double)
  #   standard/rust-run_em.s       — Rust maxbin_rs::emanager::run_em
  #   standard/rust-compute_abund.s — Rust …::compute_abund_prob_for_contig
  #   standard/summary.txt         — call counts, vector vs scalar instruction counts
  #   unity/…                      — same layout for unity (MAXBIN2_CPP_LTO=1) build
  #   comparison.txt               — side-by-side summary of both builds

  disasm-em =
    let
      nm = "${binutils-unwrapped}/bin/nm";
      objdump = "${binutils-unwrapped}/bin/objdump";

      # Shell function library shared by both variants.
      # Uses nm -n to find symbol start addresses and compute size from the
      # next symbol, then objdump -d to disassemble exactly that range.
      disasmLib = ''
        # disasm_function BINARY SYMBOL_PATTERN OUTPUT_FILE LABEL
        # Finds the function matching SYMBOL_PATTERN in BINARY using nm,
        # computes its address range, and disassembles it.
        disasm_function() {
          local binary="$1" pattern="$2" outfile="$3" label="$4"

          # nm -n -C: numeric sort, demangled. Find the line matching pattern,
          # then the next line gives us the end address.
          local match
          match=$(${nm} -n -C "$binary" | grep -m1 -E "$pattern")
          if [ -z "$match" ]; then
            echo "WARNING: symbol matching '$pattern' not found in $binary" > "$outfile"
            echo "  (label: $label)" >> "$outfile"
            return
          fi

          local addr
          addr=$(echo "$match" | awk '{print $1}')
          local sym_name
          sym_name=$(echo "$match" | awk '{$1=""; $2=""; print}' | sed 's/^ *//')

          # Find the next symbol after this address to determine the end.
          local next_addr
          next_addr=$(${nm} -n "$binary" | awk -v a="$addr" 'found {print $1; exit} $1 == a {found=1}')

          if [ -z "$next_addr" ]; then
            echo "WARNING: could not determine end address for $sym_name" > "$outfile"
            return
          fi

          echo "; $label" > "$outfile"
          echo "; Symbol: $sym_name" >> "$outfile"
          echo "; Address range: 0x$addr — 0x$next_addr" >> "$outfile"
          echo "" >> "$outfile"

          ${objdump} -d --no-show-raw-insn \
            --start-address="0x$addr" --stop-address="0x$next_addr" \
            "$binary" >> "$outfile"
        }

        # summarize_asm FILE
        # Count call instructions, vector ops (vmulpd/vaddpd/vfmadd*), and
        # scalar ops (mulsd/addsd) in a disassembly file.
        summarize_asm() {
          local file="$1" label="$2"

          calls=$(grep -c 'call' "$file" || true)
          vector=$(grep -ciE '(vmulpd|vaddpd|vmulsd|vaddsd|vfmadd|vfmsub|vfnmadd|vfnmsub)' "$file" || true)
          scalar=$(grep -ciE '(mulsd|addsd|divsd|subsd)' "$file" || true)
          total=$(wc -l < "$file")

          echo "  $label: $total lines, $calls calls, $vector vector, $scalar scalar"
        }

        # disasm_variant BINARY OUTDIR VARIANT_LABEL
        # Disassemble all four EM hot functions from BINARY into OUTDIR.
        disasm_variant() {
          local binary="$1" outdir="$2" label="$3"
          mkdir -p "$outdir"

          disasm_function "$binary" \
            'T EManager::run_EM\(int\)$' \
            "$outdir/cpp-run_EM.s" \
            "$label — C++ EManager::run_EM(int)"

          disasm_function "$binary" \
            'EManager::get_prob_abund\(double,? ?double\)' \
            "$outdir/cpp-get_prob_abund.s" \
            "$label — C++ EManager::get_prob_abund(double, double)"

          disasm_function "$binary" \
            'T maxbin_rs.*emanager.*run_em$' \
            "$outdir/rust-run_em.s" \
            "$label — Rust run_em"

          disasm_function "$binary" \
            't maxbin_rs.*emanager.*compute_abund_prob_for_contig$' \
            "$outdir/rust-compute_abund.s" \
            "$label — Rust compute_abund_prob_for_contig"

          # Per-function summary
          {
            echo "=== $label ==="
            echo ""
            summarize_asm "$outdir/cpp-run_EM.s"        "C++ EManager::run_EM(int)"
            summarize_asm "$outdir/cpp-get_prob_abund.s" "C++ EManager::get_prob_abund(double, double)"
            summarize_asm "$outdir/rust-run_em.s"        "Rust run_em"
            summarize_asm "$outdir/rust-compute_abund.s" "Rust compute_abund_prob_for_contig"
            echo ""
          } > "$outdir/summary.txt"
        }
      '';

      # Resolve the unwrapped binary (skip the shell wrapper).
      binaryOf = pkg: "${pkg}/bin/.maxbin-rs-wrapped";
      # Fallback if the wrapper layout differs.
      binaryFallback = pkg: "${pkg}/bin/maxbin-rs";

    in
    runCommand "disasm-em"
      {
        nativeBuildInputs = [ binutils-unwrapped ];
      }
      ''
        ${disasmLib}

        # Find the actual ELF binary (wrapProgram creates a shell wrapper).
        find_binary() {
          local pkg="$1"
          if [ -f "$pkg/bin/.maxbin-rs-wrapped" ]; then
            echo "$pkg/bin/.maxbin-rs-wrapped"
          else
            echo "$pkg/bin/maxbin-rs"
          fi
        }

        STD_BIN=$(find_binary "${maxbin-rs}")
        LTO_BIN=$(find_binary "${maxbin-rs-lto}")

        mkdir -p $out/standard $out/unity

        echo "Disassembling standard build: $STD_BIN"
        disasm_variant "$STD_BIN" "$out/standard" "standard"

        echo "Disassembling unity/LTO build: $LTO_BIN"
        disasm_variant "$LTO_BIN" "$out/unity" "unity (MAXBIN2_CPP_LTO=1)"

        # Side-by-side comparison
        {
          echo "=== EM Hot Function Disassembly Comparison ==="
          echo ""
          echo "Binary (standard): $STD_BIN"
          echo "Binary (unity):    $LTO_BIN"
          echo ""
          cat "$out/standard/summary.txt"
          echo ""
          cat "$out/unity/summary.txt"
        } > $out/comparison.txt

        echo ""
        echo "=== Results ==="
        cat $out/comparison.txt
      '';

  # Gene caller benchmark: FragGeneScan (C) vs FragGeneScanRs (Rust).
  # Compares wall-clock time and output on assembled contigs.
  #   nix run .#bench-genecaller-bfragilis    — B. fragilis (small, seconds)
  #   nix run .#bench-genecaller-capes        — CAPES_S7 (25K contigs, minutes)
  bench-genecaller-bfragilis = writeShellApplication {
    name = "bench-genecaller-bfragilis";
    runtimeInputs = [
      fraggenescan-rs
    ];
    text = ''
      export CONTIGS="${datasets.bfragilis.contigs}"
      export TRAIN_DIR="${fraggenescan-rs}/share/FragGeneScanRs/train"
      export PATH="${fraggenescan}/libexec/FragGeneScan:$PATH"
      bash ${../tests/bench-genecaller.sh}
    '';
  };

  bench-genecaller-capes = writeShellApplication {
    name = "bench-genecaller-capes";
    runtimeInputs = [
      fraggenescan-rs
    ];
    text = ''
      export CONTIGS="${datasets.capes-s7.contigs}"
      export TRAIN_DIR="${fraggenescan-rs}/share/FragGeneScanRs/train"
      export PATH="${fraggenescan}/libexec/FragGeneScan:$PATH"
      bash ${../tests/bench-genecaller.sh}
    '';
  };
}
