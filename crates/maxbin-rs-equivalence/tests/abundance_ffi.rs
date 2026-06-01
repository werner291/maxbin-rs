//! FFI equivalence test for the abundance-file parser.
//!
//! Compares `maxbin_rs::abundance` against the original MaxBin2 C++
//! AbundanceLoader. Previously an in-source test in `abundance.rs`.
use maxbin_rs::abundance::parse_file;
use maxbin_rs_equivalence::original_ffi::OriginalAbundanceLoader;

#[test]
fn equivalence_with_original() {
    // Writes a temp abundance file, parses with both Rust and C++ FFI,
    // compares field by field.
    use std::io::Write;

    let test_data = "contig_1\t10.5\ncontig_2\t0.00001\ncontig_3\t100.0\ncontig_4 0.5\n";

    let tmp_path = std::env::temp_dir().join("maxbin_rs_test_abund.txt");
    {
        let mut f = std::fs::File::create(&tmp_path).unwrap();
        f.write_all(test_data.as_bytes()).unwrap();
    }

    // Rust
    let rust_records = parse_file(&tmp_path).unwrap();

    // Original C++
    let original = OriginalAbundanceLoader::new(&tmp_path);

    assert_eq!(rust_records.len(), original.num_records() as usize);
    assert!(original.is_parse_success());

    for (i, rec) in rust_records.iter().enumerate() {
        let orig_abund = original.abundance_by_index(i as i32);
        assert!(
            (rec.abundance - orig_abund).abs() < 1e-10,
            "abundance mismatch at record {i} (header: {}): rust={} original={}",
            rec.header,
            rec.abundance,
            orig_abund,
        );
    }

    std::fs::remove_file(&tmp_path).ok();
}
