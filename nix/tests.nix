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
# Pipeline stage tests (primary verification):
#   nix run .#test-pipeline-stages          — B. fragilis (~1 min)
#   nix run .#test-pipeline-stages-minigut  — minigut (~2 min)
#   nix run .#test-pipeline-stages-capes    — CAPES_S7 (~10 min, 2.5 GB download)
#   nix run .#test-pipeline-stages-cami     — CAMI I High (~1-3 hours, 779 MB contigs)
#   nix run .#test-pipeline-stages-metahit  — MetaHIT (~3 hours, 275 MB contigs)
#
# Component benchmark:
#   nix run .#bench-components              — per-function Rust vs C++ timing

{
  writeShellApplication,
  runCommand,
  binutils-unwrapped,
  maxbin2,
  maxbin-rs,
  maxbin-rs-lto,
  rust,
  cargo-nextest,
  datasets,
  intermediates,
}:

let
  # Helper: create a pipeline stage test for a given dataset.
  # All three tests run the same script (tests/pipeline-stages.sh) with
  # different data — the script doesn't know or care which dataset it gets.
  mkStageTest =
    {
      name,
      dataset,
      intermediates',
    }:
    writeShellApplication {
      name = "test-pipeline-stages-${name}";
      # Both the original MaxBin2 and maxbin-rs are on $PATH during the test,
      # so the script can invoke either.
      runtimeInputs = [
        maxbin2
        maxbin-rs
      ];
      text = ''
        # INTERMEDIATES: path to pre-computed outputs from the original MaxBin2
        export INTERMEDIATES="${intermediates'}"
        # MAXBIN2_TEST_CONTIGS: the raw input contigs for this dataset
        export MAXBIN2_TEST_CONTIGS="${dataset.contigs}"
        exec ${../tests/pipeline-stages.sh}
      '';
    };
in
{
  test-pipeline-stages = mkStageTest {
    name = "bfragilis";
    dataset = datasets.bfragilis;
    intermediates' = intermediates.bfragilis;
  };

  test-pipeline-stages-minigut = mkStageTest {
    name = "minigut";
    dataset = datasets.minigut;
    intermediates' = intermediates.minigut;
  };

  test-pipeline-stages-capes = mkStageTest {
    name = "capes";
    dataset = datasets.capes-s7;
    intermediates' = intermediates.capes;
  };

  test-pipeline-stages-cami = mkStageTest {
    name = "cami";
    dataset = datasets.cami-i-high;
    intermediates' = intermediates.cami;
  };

  test-pipeline-stages-metahit = mkStageTest {
    name = "metahit";
    dataset = datasets.metahit;
    intermediates' = intermediates.metahit;
  };

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
}
