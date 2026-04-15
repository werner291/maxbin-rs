# Build the original MaxBin2 2.2.7 from source.
#
# This packages the original tool (Perl scripts + C++ core) with all its
# runtime dependencies (HMMER, Bowtie2, FragGeneScan) wired in via Nix.
# It replaces the manual `setting` file and CPAN dependency management
# with deterministic Nix store paths.
#
# The deterministic patch (maxbin2-deterministic.patch) is applied to make
# output reproducible for equivalence testing — see that file for details.
#
# This package exists primarily as a reference implementation: the pipeline
# stage tests (tests/pipeline-stages.sh) run this original MaxBin2 to
# generate cached intermediate files, then compare the Rust reimplementation
# against those outputs.

{ lib, stdenv, perl, hmmer, bowtie2, fraggenescan, makeWrapper, src }:

stdenv.mkDerivation rec {
  pname = "maxbin2";
  version = "2.2.7";

  inherit src;
  sourceRoot = "MaxBin-${version}";

  nativeBuildInputs = [ makeWrapper ];
  buildInputs = [ perl ];

  # See nix/maxbin2-deterministic.patch for details:
  # 1. Remove unused `use LWP::Simple` (avoids pulling in libwww-perl)
  # 2. Sort all `foreach (keys %hash)` iterations for deterministic output
  #    (Perl 5.18+ randomizes hash key order, causing different seed selection
  #    and bin assignments on every run)
  patches = [ ./maxbin2-deterministic.patch ];

  buildPhase = ''
    # Only build the C++ core — we provide dependencies via Nix, not buildapp
    cd src
    make
    cd ..
  '';

  installPhase = ''
    mkdir -p $out/bin $out/libexec/maxbin2/src

    # Install the C++ binary — must be at src/MaxBin relative to the Perl
    # scripts because run_MaxBin.pl hardcodes: $MAXBIN = "$Bin/src/MaxBin"
    cp src/MaxBin $out/libexec/maxbin2/src/

    # Install Perl scripts
    cp run_MaxBin.pl $out/libexec/maxbin2/
    cp _getmarker.pl $out/libexec/maxbin2/
    cp _getabund.pl $out/libexec/maxbin2/
    cp _sepReads.pl $out/libexec/maxbin2/
    chmod +x $out/libexec/maxbin2/run_MaxBin.pl

    # Install HMM marker gene databases
    cp marker.hmm $out/libexec/maxbin2/
    cp bacar_marker.hmm $out/libexec/maxbin2/

    # Install R heatmap script
    cp heatmap.r $out/libexec/maxbin2/

    # Generate the setting file pointing to Nix store paths.
    # The regex in checkProgram() does NOT tolerate leading whitespace.
    cat > $out/libexec/maxbin2/setting <<EOF
[FragGeneScan] ${fraggenescan}/libexec/FragGeneScan
[Bowtie2] ${bowtie2}/bin
[HMMER3] ${hmmer}/bin
EOF

    # Create wrapper that sets up PATH so the Perl script can find all tools.
    # run_MaxBin.pl uses FindBin to locate sibling scripts and data files
    # relative to $Bin, so the wrapper must exec the real script (makeWrapper
    # does this correctly — $0 in Perl resolves through the wrapper).
    makeWrapper $out/libexec/maxbin2/run_MaxBin.pl $out/bin/run_MaxBin.pl \
      --prefix PATH : "${hmmer}/bin" \
      --prefix PATH : "${bowtie2}/bin" \
      --prefix PATH : "${fraggenescan}/libexec/FragGeneScan" \
      --prefix PATH : "$out/libexec/maxbin2" \
      --prefix PATH : "${perl}/bin"
  '';

  meta = with lib; {
    description = "Automated binning of assembled metagenomic sequences";
    homepage = "https://sourceforge.net/projects/maxbin2/";
    # BSD-3-Clause (DOE variant) — see LICENSE file in tarball
    license = licenses.bsd3;
    platforms = platforms.unix;
  };
}
