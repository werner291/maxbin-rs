# Build FragGeneScanRs — a Rust reimplementation of FragGeneScan.
#
# Van der Jeugt et al., BMC Bioinformatics 2022.
# https://github.com/unipept/FragGeneScanRs
#
# Unlike the original C FragGeneScan (which requires a Perl wrapper and
# 2010-era compiler workarounds), this builds cleanly with a modern Rust
# toolchain. Training data ships with the source and is installed alongside
# the binary.
#
# The upstream repo does not commit a Cargo.lock, so we supply one
# (fraggenescan-rs-Cargo.lock) generated from the v1.1.0 tag.

{
  lib,
  rustPlatform,
  fetchFromGitHub,
}:

rustPlatform.buildRustPackage rec {
  pname = "fraggenescan-rs";
  version = "1.1.0";

  src = fetchFromGitHub {
    owner = "unipept";
    repo = "FragGeneScanRs";
    rev = "v${version}";
    hash = "sha256-YmOZ0wnjfAnhSd/aIMGCXGsR9TQUnPzTTiz9XfL5W9Y=";
  };

  cargoLock.lockFile = ./fraggenescan-rs-Cargo.lock;

  postPatch = ''
    ln -s ${./fraggenescan-rs-Cargo.lock} Cargo.lock
  '';

  postInstall = ''
    # Install training data alongside the binary so -r can find it.
    mkdir -p $out/share/FragGeneScanRs
    cp -r train $out/share/FragGeneScanRs/
  '';

  meta = with lib; {
    description = "Rust implementation of FragGeneScan gene prediction";
    homepage = "https://github.com/unipept/FragGeneScanRs";
    license = licenses.gpl3Plus;
    platforms = platforms.unix;
  };
}
