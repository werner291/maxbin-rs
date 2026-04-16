# Nix flake — the top-level build and test configuration.
#
# If you are not familiar with Nix: a flake is a declarative build file
# (similar in role to a Makefile or CMakeLists.txt) that pins every
# dependency to an exact version. Running `nix build` or `nix run` on
# any machine with Nix installed will produce the same binary from the
# same source, with the same compiler, libraries, and tools.
#
# This file declares:
#   - How to build maxbin-rs and the original MaxBin2
#   - Where to find test datasets (see nix/datasets.nix)
#   - How to generate pre-computed intermediates (see nix/intermediates.nix)
#   - How to run equivalence tests (see nix/tests.nix)
#
# Quick reference:
#   nix build                               — build maxbin-rs
#   nix run .#test-pipeline-stages          — run B. fragilis stage tests
#   nix run .#test-pipeline-stages-minigut  — run minigut stage tests
#   nix run .#test-pipeline-stages-capes    — run CAPES_S7 stage tests
#   nix run .#test-pipeline-stages-cami     — run CAMI I High stage tests
#   nix run .#test-pipeline-stages-metahit  — run MetaHIT stage tests
#   nix build .#disasm-em                   — disassemble EM hot functions
#   nix develop                             — enter development shell

{
  description = "maxbin-rs — Rust reimplementation of MaxBin2";

  # Pinned dependencies. The flake.lock file records exact commit hashes
  # for each input, ensuring reproducible builds.
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  # `outputs` defines what this flake produces. The `eachDefaultSystem` call
  # generates packages for Linux and macOS automatically.
  outputs = { self, nixpkgs, flake-utils, rust-overlay }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        overlays = [ (import rust-overlay) ];
        pkgs = import nixpkgs {
          inherit system overlays;
          # FragGeneScan has no explicit license file — allow it here.
          config.allowUnfreePredicate = pkg:
            builtins.elem (nixpkgs.lib.getName pkg) [ "fraggenescan" ];
        };
        rust = pkgs.rust-bin.stable.latest.default.override {
          extensions = [ "rust-src" "rust-analyzer" ];
        };

        # =====================================================================
        # Packages
        # =====================================================================

        # Gene caller used by MaxBin2 (see nix/fraggenescan.nix for build details).
        fraggenescan = pkgs.callPackage ./nix/fraggenescan.nix { };

        # Upstream MaxBin2 source tarball — used in two places:
        #   1. nix/maxbin2.nix builds the original MaxBin2 from it
        #   2. build.rs extracts C++ source from it for FFI equivalence testing
        maxbin2-src-tarball = pkgs.fetchurl {
          url = "https://downloads.sourceforge.net/project/maxbin2/MaxBin-2.2.7.tar.gz";
          hash = "sha256-y2Qp6FcoDCt1gjyM1VBY7Rack7xweka94MQ4Pyv/4J4=";
        };

        # The original MaxBin2, built from source (see nix/maxbin2.nix).
        maxbin2 = pkgs.callPackage ./nix/maxbin2.nix {
          inherit fraggenescan;
          src = maxbin2-src-tarball;
        };

        # The Rust reimplementation — this is the main output of this project.
        maxbin-rs = pkgs.rustPlatform.buildRustPackage {
          pname = "maxbin-rs";
          version = "0.1.0";
          src = ./.;
          cargoHash = "sha256-5KrBRpVM8yM8QViX0/G/IjmFpankZL8AeCTD+pTf/zw=";
          # build.rs extracts C++ source from this tarball for FFI testing.
          MAXBIN2_SRC_TARBALL = "${maxbin2-src-tarball}";
          # Wrap the binary so HMMER, Bowtie2, and FragGeneScan are on $PATH
          # at runtime — no manual 'setting' file needed.
          postInstall = ''
            wrapProgram $out/bin/maxbin-rs \
              --prefix PATH : "${pkgs.lib.makeBinPath [ pkgs.hmmer pkgs.bowtie2 ]}:${fraggenescan}/libexec/FragGeneScan"
          '';
          nativeBuildInputs = [ pkgs.makeWrapper ];
        };

        # Same as maxbin-rs but with C++ LTO enabled for the FFI library.
        # Used for A/B benchmarking to test whether -flto closes the ~8x gap.
        maxbin-rs-lto = pkgs.rustPlatform.buildRustPackage {
          pname = "maxbin-rs-lto";
          version = "0.1.0";
          src = ./.;
          cargoHash = "sha256-5KrBRpVM8yM8QViX0/G/IjmFpankZL8AeCTD+pTf/zw=";
          MAXBIN2_SRC_TARBALL = "${maxbin2-src-tarball}";
          MAXBIN2_CPP_LTO = "1";
          postInstall = ''
            wrapProgram $out/bin/maxbin-rs \
              --prefix PATH : "${pkgs.lib.makeBinPath [ pkgs.hmmer pkgs.bowtie2 ]}:${fraggenescan}/libexec/FragGeneScan"
          '';
          nativeBuildInputs = [ pkgs.makeWrapper ];
        };

        # =====================================================================
        # Test data and pre-computed intermediates
        # =====================================================================

        # Pinned test datasets — see nix/datasets.nix for URLs and hashes.
        datasets = import ./nix/datasets.nix { inherit (pkgs) fetchurl; };

        # Factory for generating cached pipeline intermediates — see
        # nix/intermediates.nix for what this produces and why.
        intermediatesLib = import ./nix/intermediates.nix {
          inherit (pkgs) runCommand perl gawk;
          inherit maxbin2;
        };

        # One set of cached intermediates per dataset. Each is computed once
        # (by running the original MaxBin2) and cached in the Nix store.
        intermediates = {
          bfragilis = intermediatesLib.mkPipelineIntermediates {
            name = "bfragilis";
            contigs = datasets.bfragilis.contigs;
            reads = datasets.bfragilis.reads1;
          };
          minigut = intermediatesLib.mkPipelineIntermediates {
            name = "minigut";
            contigs = datasets.minigut.contigs;
            reads = datasets.minigut.reads1;
          };
          capes = intermediatesLib.mkPipelineIntermediates {
            name = "capes";
            contigs = datasets.capes-s7.contigs;
            reads = datasets.capes-s7.reads1;
          };
          # Large datasets with pre-computed abundance (no raw reads).
          # These skip Bowtie2 and use the provided depth files directly.
          cami = intermediatesLib.mkPipelineIntermediatesFromAbund {
            name = "cami-i-high";
            contigs = datasets.cami-i-high.contigs;
            depth = datasets.cami-i-high.depth;
            depthFormat = "metabat";  # NERSC MetaBAT depth format
          };
          metahit = intermediatesLib.mkPipelineIntermediatesFromAbund {
            name = "metahit";
            contigs = datasets.metahit.contigs;
            depth = datasets.metahit.depth;
            depthFormat = "maxbin";   # Already in MaxBin format
          };
        };

        # =====================================================================
        # Tests and benchmarks — see nix/tests.nix
        # =====================================================================

        tests = import ./nix/tests.nix {
          inherit (pkgs) writeShellApplication runCommand binutils-unwrapped;
          inherit maxbin2 maxbin-rs maxbin-rs-lto rust;
          inherit (pkgs) cargo-nextest;
          inherit datasets intermediates;
        };

      in {
        # Everything listed here can be built with `nix build .#<name>`
        # or run with `nix run .#<name>`.
        packages = {
          default = maxbin-rs;
          inherit maxbin-rs maxbin2 fraggenescan;
          inherit (intermediates) bfragilis minigut capes cami metahit;
          inherit (tests)
            test-pipeline-stages test-pipeline-stages-minigut
            test-pipeline-stages-capes test-pipeline-stages-cami
            test-pipeline-stages-metahit bench-components bench-cpp-lto
            disasm-em;
        };

        # `nix develop` drops you into this shell with all tools available.
        devShells.default = pkgs.mkShell {
          buildInputs = [ rust pkgs.cargo-nextest maxbin2 ];
          # Environment variables available inside the devshell:
          MAXBIN2_TEST_CONTIGS = "${datasets.bfragilis.contigs}";
          MAXBIN2_TEST_READS1 = "${datasets.bfragilis.reads1}";
          MAXBIN2_SRC_TARBALL = "${maxbin2-src-tarball}";
          shellHook = ''
            export PATH="${fraggenescan}/libexec/FragGeneScan:$PATH"
            echo "maxbin-rs devshell"
            echo "  maxbin2 (original): $(run_MaxBin.pl -v 2>&1 | head -1 || echo 'available')"
            echo "  test contigs: $MAXBIN2_TEST_CONTIGS"
            echo "  test reads:   $MAXBIN2_TEST_READS1"
          '';
        };
      });
}
