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

{ writeShellApplication, maxbin2, maxbin-rs, rust, cargo-nextest,
  datasets, intermediates }:

let
  # Helper: create a pipeline stage test for a given dataset.
  # All three tests run the same script (tests/pipeline-stages.sh) with
  # different data — the script doesn't know or care which dataset it gets.
  mkStageTest = { name, dataset, intermediates' }:
    writeShellApplication {
      name = "test-pipeline-stages-${name}";
      # Both the original MaxBin2 and maxbin-rs are on $PATH during the test,
      # so the script can invoke either.
      runtimeInputs = [ maxbin2 maxbin-rs ];
      text = ''
        # INTERMEDIATES: path to pre-computed outputs from the original MaxBin2
        export INTERMEDIATES="${intermediates'}"
        # MAXBIN2_TEST_CONTIGS: the raw input contigs for this dataset
        export MAXBIN2_TEST_CONTIGS="${dataset.contigs}"
        exec ${../tests/pipeline-stages.sh}
      '';
    };
in {
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

  # Per-function performance comparison: runs each component (profiler,
  # distance metrics, sort, etc.) in a timed loop, measuring Rust throughput
  # against the original C++ via FFI. Produces a comparison table.
  bench-components = writeShellApplication {
    name = "bench-components";
    runtimeInputs = [ rust cargo-nextest ];
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
}
