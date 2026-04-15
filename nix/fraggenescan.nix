# Build FragGeneScan 1.30, a gene prediction tool used by MaxBin2.
#
# FragGeneScan identifies protein-coding regions in short DNA sequences.
# MaxBin2 uses it to find open reading frames before searching for marker
# genes with HMMER. This is 2010-era C code that requires a compiler flag
# workaround to build with modern GCC (see NIX_CFLAGS_COMPILE below).

{ lib, stdenv, fetchurl, perl }:

stdenv.mkDerivation rec {
  pname = "fraggenescan";
  version = "1.30";

  src = fetchurl {
    url = "https://downloads.sourceforge.net/project/maxbin-auxiliary/FragGeneScan${version}.tar.gz";
    hash = "sha256-8tfw36Sl9LvqKV7YZdy/7fFslU6hU0wqh569z7hlDZU=";
  };

  sourceRoot = "FragGeneScan${version}";

  buildInputs = [ perl ];

  # Old C code (2010-era) with missing forward declarations that modern GCC
  # rejects. Suppress the specific error to build unmodified upstream source.
  env.NIX_CFLAGS_COMPILE = "-Wno-error=implicit-function-declaration";

  buildPhase = ''
    make
    make fgs
  '';

  # MaxBin2's setting file points at a directory containing run_FragGeneScan.pl,
  # the FragGeneScan binary, and the train/ data. The Perl wrapper resolves
  # paths relative to $0 with a hardcoded offset, so everything must be
  # co-located. We install into libexec/ and let MaxBin2 reference that path.
  installPhase = ''
    mkdir -p $out/libexec/FragGeneScan

    cp FragGeneScan $out/libexec/FragGeneScan/
    cp run_FragGeneScan.pl $out/libexec/FragGeneScan/
    chmod +x $out/libexec/FragGeneScan/run_FragGeneScan.pl
    cp -r train $out/libexec/FragGeneScan/
  '';

  meta = with lib; {
    description = "Application for finding (fragmented) genes in short reads";
    # NOTE: No license file found in the FragGeneScan 1.30 tarball.
    # It is freely distributed academic software from Indiana University.
    # The original paper is: Rho et al., Nucleic Acids Research, 2010.
    # Treat as unfree/unknown until clarified.
    license = licenses.unfreeRedistributable;
    platforms = platforms.unix;
  };
}
