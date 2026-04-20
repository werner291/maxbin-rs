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
#   nix run .#test-pipeline-stages-bfragilis — run B. fragilis stage tests
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
    crane.url = "github:ipetkov/crane";
  };

  # `outputs` defines what this flake produces. Currently x86_64-linux only;
  # add more systems to the list below when they are tested.
  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      rust-overlay,
      crane,
    }:
    flake-utils.lib.eachSystem [ "x86_64-linux" ] (
      system:
      let
        overlays = [ (import rust-overlay) ];
        pkgs = import nixpkgs {
          inherit system overlays;
          # FragGeneScan has no explicit license file — allow it here.
          config.allowUnfreePredicate = pkg: builtins.elem (nixpkgs.lib.getName pkg) [ "fraggenescan" ];
        };
        rust = pkgs.rust-bin.stable.latest.default.override {
          extensions = [
            "rust-src"
            "rust-analyzer"
            "clippy"
          ];
        };
        rustPlatform = pkgs.makeRustPlatform {
          cargo = rust;
          rustc = rust;
        };

        # Crane — builds Rust with cached dependency artifacts.
        # Dependencies are built once (craneLib.buildDepsOnly), then source
        # changes only rebuild our code, not all of crates.io.
        craneLib = (crane.mkLib pkgs).overrideToolchain rust;
        # Include Rust source + files needed by build.rs (FFI wrappers,
        # patches) and tests (fixtures, test scripts).
        src = pkgs.lib.cleanSourceWith {
          src = ./.;
          filter =
            path: type:
            (craneLib.filterCargoSources path type) || builtins.match ".*/(nix|vendor|tests)/.*" path != null;
        };
        commonArgs = {
          inherit src;
          pname = "maxbin-rs";
          version = "0.2.0";
          MAXBIN2_SRC_TARBALL = "${maxbin2-src-tarball}";
          nativeBuildInputs = [ pkgs.makeWrapper ];
        };
        cargoArtifacts = craneLib.buildDepsOnly commonArgs;

        # =====================================================================
        # Packages
        # =====================================================================

        # Gene caller used by MaxBin2 (see nix/fraggenescan.nix for build details).
        fraggenescan = pkgs.callPackage ./nix/fraggenescan.nix { };

        # Rust reimplementation of FragGeneScan (see nix/fraggenescan-rs.nix).
        fraggenescan-rs = pkgs.callPackage ./nix/fraggenescan-rs.nix {
          inherit rustPlatform;
        };

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

        # MaxBin2 with long double → double patch for precision testing.
        maxbin2-f64 = pkgs.callPackage ./nix/maxbin2.nix {
          inherit fraggenescan;
          src = maxbin2-src-tarball;
          extraPatches = [ ./nix/maxbin2-cpp-ffi-f64.patch ];
        };

        # MaxBin2 variants with trace logging for pipeline-trace.sh.
        # Adds [trace] lines to stderr showing recursion depth, bin names,
        # and sub-bin counts at each step — observation only, no behavior change.
        tracePatch = ./nix/maxbin2-trace-logging.patch;
        maxbin2-trace = maxbin2.overrideAttrs (old: {
          patches = (old.patches or [ ]) ++ [ tracePatch ];
        });
        maxbin2-f64-trace = maxbin2-f64.overrideAttrs (old: {
          patches = (old.patches or [ ]) ++ [ tracePatch ];
        });

        # EM-only binary — no external tool dependencies.
        # Supports `em`, `filter`, and `seeds` subcommands only.
        # Pipeline mode errors out telling the user to use the full package.
        maxbin-rs-em = craneLib.buildPackage (
          commonArgs
          // {
            inherit cargoArtifacts;
            pname = "maxbin-rs-em";
          }
        );

        # The Rust reimplementation — this is the main output of this project.
        maxbin-rs = craneLib.buildPackage (
          commonArgs
          // {
            inherit cargoArtifacts;
            # Wrap the binary so HMMER, Bowtie2, and FragGeneScan are on $PATH
            # at runtime — no manual 'setting' file needed.
            postInstall = ''
              cp ${maxbin2}/libexec/maxbin2/marker.hmm $out/bin/
              cp ${maxbin2}/libexec/maxbin2/bacar_marker.hmm $out/bin/

              # FragGeneScanRs training data — same files as the original FragGeneScan.
              # Resolved by the binary via ../share/FragGeneScanRs/train/ relative
              # to the bin directory.
              mkdir -p $out/share/FragGeneScanRs
              ln -s ${fraggenescan}/libexec/FragGeneScan/train $out/share/FragGeneScanRs/train

              wrapProgram $out/bin/maxbin-rs \
                --prefix PATH : "${
                  pkgs.lib.makeBinPath [
                    pkgs.hmmer
                    pkgs.bowtie2
                  ]
                }:${fraggenescan}/libexec/FragGeneScan"
            '';
          }
        );

        # Same as maxbin-rs but with C++ LTO enabled for the FFI library.
        # Used for A/B benchmarking to test whether -flto closes the ~8x gap.
        maxbin-rs-lto = craneLib.buildPackage (
          commonArgs
          // {
            inherit cargoArtifacts;
            pname = "maxbin-rs-lto";
            MAXBIN2_CPP_LTO = "1";
            postInstall = ''
              wrapProgram $out/bin/maxbin-rs \
                --prefix PATH : "${
                  pkgs.lib.makeBinPath [
                    pkgs.hmmer
                    pkgs.bowtie2
                  ]
                }:${fraggenescan}/libexec/FragGeneScan"
            '';
          }
        );

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
            depthFormat = "metabat"; # NERSC MetaBAT depth format
          };
          # Downsampled CAMI I High: first 5000 filtered contigs with matching
          # abundance and HMMER hits. Small enough to run the full recursive
          # pipeline in ~7 minutes, large enough to trigger depth-5 recursion.
          cami-small =
            pkgs.runCommand "cami-small-intermediates"
              {
                nativeBuildInputs = [
                  maxbin-rs
                  pkgs.gawk
                ];
              }
              ''
                # Filter full CAMI contigs
                maxbin-rs filter --contig ${datasets.cami-i-high.contigs} --out $TMPDIR/full 2>/dev/null

                # Take first 5000 filtered contigs
                awk '/^>/{n++} n>5000{exit} {print}' $TMPDIR/full.contig.tmp > $TMPDIR/small.fa

                # Extract contig names
                grep '^>' $TMPDIR/small.fa | sed 's/^>//;s/ .*//' > $TMPDIR/names.txt

                # Subset abundance (exact match on first column)
                awk -F'\t' 'NR==FNR{names[$1]; next} $1 in names' \
                  $TMPDIR/names.txt ${intermediates.cami}/abund > $TMPDIR/abund.txt

                # Subset HMMER hits (match contig name from gene ID)
                grep '^#' ${intermediates.cami}/hmmout > $TMPDIR/hmmout
                awk 'NR==FNR{names[$1]; next} /^#/{next} {
                  n=split($1,a,"_"); contig="";
                  for(i=1;i<=n-3;i++){if(i>1)contig=contig"_"; contig=contig a[i]}
                  if(contig in names) print
                }' $TMPDIR/names.txt ${intermediates.cami}/hmmout >> $TMPDIR/hmmout

                mkdir -p $out
                cp $TMPDIR/small.fa $out/contigs.fa
                cp $TMPDIR/abund.txt $out/abund
                cp $TMPDIR/hmmout $out/hmmout
                cp $TMPDIR/full.tooshort $out/tooshort

                # Generate seeds from subsetted HMMER output
                MAXBIN_RS_DETERMINISTIC=1 maxbin-rs seeds \
                  --contig $out/contigs.fa --hmmout $out/hmmout --out $TMPDIR/seeds 2>/dev/null
                cp $TMPDIR/seeds.seed $out/seed

                echo "=== cami-small summary ==="
                echo "  contigs: $(grep -c '^>' $out/contigs.fa)"
                echo "  abundance lines: $(wc -l < $out/abund)"
                echo "  HMMER hits: $(grep -cv '^#' $out/hmmout)"
                echo "  seeds: $(wc -l < $out/seed)"
              '';

          metahit = intermediatesLib.mkPipelineIntermediatesFromAbund {
            name = "metahit";
            contigs = datasets.metahit.contigs;
            depth = datasets.metahit.depth;
            depthFormat = "maxbin"; # Already in MaxBin format
          };
        };

        # =====================================================================
        # Tests and benchmarks — see nix/tests.nix
        # =====================================================================

        tests = import ./nix/tests.nix {
          inherit (pkgs) writeShellApplication runCommand binutils-unwrapped;
          inherit
            maxbin2
            maxbin2-f64
            maxbin2-trace
            maxbin2-f64-trace
            maxbin-rs
            maxbin-rs-lto
            fraggenescan
            fraggenescan-rs
            rust
            ;
          inherit (pkgs) cargo-nextest;
          inherit datasets intermediates;
        };

        # Docker image for ghcr.io — minimal, just the wrapped binary.
        dockerImage = pkgs.dockerTools.buildLayeredImage {
          name = "ghcr.io/werner291/maxbin-rs";
          tag = "latest";
          contents = [
            maxbin-rs
            pkgs.coreutils
            pkgs.bashInteractive
            pkgs.gawk
            # Bowtie2's Perl wrapper shells out via /bin/sh
            (pkgs.runCommand "sh-symlink" { } ''
              mkdir -p $out/bin
              ln -s ${pkgs.bashInteractive}/bin/bash $out/bin/sh
            '')
            # /tmp for intermediate files (Nix Docker images have no FHS skeleton)
            (pkgs.runCommand "tmp-dir" { } ''
              mkdir -p $out/tmp
            '')
          ];
          config = {
            Entrypoint = [ "${maxbin-rs}/bin/maxbin-rs" ];
            Labels = {
              "org.opencontainers.image.source" = "https://github.com/werner291/maxbin-rs";
            };
          };
        };

        # NixOS VM test — run the Docker image on real data.
        dockerTest = import ./nix/docker-test.nix {
          inherit pkgs dockerImage datasets;
        };

      in
      {
        # Everything listed here can be built with `nix build .#<name>`
        # or run with `nix run .#<name>`.
        packages = {
          default = maxbin-rs;
          inherit
            maxbin-rs
            maxbin2
            fraggenescan
            fraggenescan-rs
            dockerImage
            dockerTest
            ;
          inherit (intermediates)
            bfragilis
            capes
            cami
            metahit
            ;
          "cami-small" = intermediates.cami-small;
          inherit (tests)
            test-pipeline-stages-bfragilis
            test-pipeline-stages-capes
            test-pipeline-stages-cami
            test-pipeline-stages-metahit
            test-cli-bfragilis
            test-cli-equivalence-bfragilis
            test-cli-equivalence-capes
            trace-cami-small
            trace-cami-small-ldouble
            trace-cami
            test-precision-divergence
            bench-components
            bench-cpp-lto
            disasm-em
            bench-genecaller-bfragilis
            bench-genecaller-capes
            ;
        };

        # `nix develop` drops you into this shell with all tools available.
        devShells.default = pkgs.mkShell {
          buildInputs = [
            rust
            pkgs.cargo-nextest
            pkgs.gh
            pkgs.cachix
            pkgs.nixfmt-rfc-style
            maxbin2
          ];
          # Environment variables available inside the devshell:
          MAXBIN2_TEST_CONTIGS = "${datasets.bfragilis.contigs}";
          MAXBIN2_TEST_READS1 = "${datasets.bfragilis.reads1}";
          MAXBIN2_SRC_TARBALL = "${maxbin2-src-tarball}";
          shellHook = ''
            export PATH="${fraggenescan}/libexec/FragGeneScan:$PATH"

            # Auto-format pre-commit hook
            mkdir -p .git/hooks
            cat > .git/hooks/pre-commit <<'HOOK'
            #!/bin/sh
            nix develop -c cargo fmt
            git add -u -- '*.rs'
            nix run nixpkgs#nixfmt -- flake.nix nix/*.nix
            git add -u -- '*.nix'
            HOOK
            chmod +x .git/hooks/pre-commit

            # Pre-push hook: runs the same checks as CI
            cat > .git/hooks/pre-push <<'HOOK'
            #!/bin/sh
            echo "Running pre-push checks (nix flake check)..."
            nix flake check || { echo "Pre-push checks failed."; exit 1; }
            echo "Pre-push checks passed."
            HOOK
            chmod +x .git/hooks/pre-push

            echo "maxbin-rs devshell"
            echo "  maxbin2 (original): $(run_MaxBin.pl -v 2>&1 | head -1 || echo 'available')"
            echo "  test contigs: $MAXBIN2_TEST_CONTIGS"
            echo "  test reads:   $MAXBIN2_TEST_READS1"
          '';
        };

        # `nix flake check` runs these automatically.
        # The pre-push hook and CI both call `nix flake check`.
        checks = {
          fmt =
            pkgs.runCommand "check-fmt"
              {
                nativeBuildInputs = [
                  rust
                  pkgs.nixfmt-rfc-style
                ];
              }
              ''
                cd ${./.}
                cargo fmt --check
                nixfmt --check flake.nix nix/*.nix
                touch $out
              '';
          clippy = craneLib.cargoClippy (
            commonArgs
            // {
              inherit cargoArtifacts;
              cargoClippyExtraArgs = "--tests -- -D warnings";
            }
          );
          tests = craneLib.cargoNextest (
            commonArgs
            // {
              inherit cargoArtifacts;
              nativeBuildInputs = (commonArgs.nativeBuildInputs or [ ]) ++ [ pkgs.cargo-nextest ];
            }
          );
          inherit (tests) test-pipeline-stages-bfragilis test-cli-bfragilis test-cli-equivalence-bfragilis;
          inherit dockerTest;
          # End-to-end recursive equivalence: f64-patched C++ vs Rust on
          # downsampled CAMI (5000 contigs, depth-5 recursion). Asserts
          # byte-identical bin output via the trace script's verdict.
          recursive-equivalence =
            pkgs.runCommand "check-recursive-equivalence"
              {
                nativeBuildInputs = [
                  maxbin2-f64
                  maxbin-rs
                ];
              }
              ''
                export MAXBIN_RS_DETERMINISTIC=1
                export HOME=$(mktemp -d)
                bash ${./tests/pipeline-trace.sh} \
                  "${intermediates.cami-small}/contigs.fa" \
                  "${intermediates.cami-small}/abund" \
                  "${intermediates.cami-small}/hmmout"
                touch $out
              '';
        };
      }
    );
}
