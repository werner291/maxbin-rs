//! FFI equivalence test for the FASTA parser.
//!
//! Compares `maxbin_rs::fasta` against the original MaxBin2 C++ fastaReader.
//! Requires nf-core test data via the `MAXBIN2_TEST_CONTIGS` env var (set by
//! the Nix devshell); skips otherwise. Previously an in-source test in
//! `fasta.rs`.
use maxbin_rs::fasta::parse_file;
use maxbin_rs_equivalence::original_ffi::OriginalFastaReader;

/// Decompress a .gz file to a temp path for the C++ parser (which can't read gzip).
fn decompress_gz(gz_path: &std::path::Path) -> std::path::PathBuf {
    use std::io::{Read, Write};
    let out_path = std::env::temp_dir().join("maxbin_rs_test_contigs.fa");
    let file = std::fs::File::open(gz_path).unwrap();
    let mut decoder = flate2::read::GzDecoder::new(file);
    let mut contents = Vec::new();
    decoder.read_to_end(&mut contents).unwrap();
    let mut out = std::fs::File::create(&out_path).unwrap();
    out.write_all(&contents).unwrap();
    out_path
}

#[test]
fn equivalence_with_original_fasta_parser() {
    // This test requires nf-core test data (only available in devshell).
    let gz_path = match std::env::var("MAXBIN2_TEST_CONTIGS") {
        Ok(p) => std::path::PathBuf::from(p),
        Err(_) => return,
    };

    // Parse with our Rust implementation (handles gzip natively)
    let rust_records = parse_file(&gz_path).expect("rust parser failed");

    // Decompress for the C++ parser, then parse with the original
    let plain_path = decompress_gz(&gz_path);
    let original = OriginalFastaReader::new(&plain_path);

    // Same number of records
    assert_eq!(
        rust_records.len(),
        original.num_records() as usize,
        "record count mismatch"
    );

    // Field-by-field comparison for every record
    for (i, rust_rec) in rust_records.iter().enumerate() {
        let idx = i as u32;

        assert_eq!(
            rust_rec.header,
            original.header(idx),
            "header mismatch at record {i}"
        );

        assert_eq!(
            rust_rec.seq,
            original.seq(idx),
            "sequence mismatch at record {i} (header: {})",
            rust_rec.header
        );

        assert_eq!(
            rust_rec.len(),
            original.seq_len(idx) as usize,
            "length mismatch at record {i} (header: {})",
            rust_rec.header
        );
    }

    std::fs::remove_file(&plain_path).ok();
}
