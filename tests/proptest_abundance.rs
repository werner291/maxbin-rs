/// Property-based equivalence tests for abundance.rs.
///
/// Coverage:
/// - parse(): abundance_equivalence (proptest, 200 cases)
///   Generators exercise: all four separator types (tab, space, comma, semicolon),
///   consecutive separators (the && bug), scientific notation, empty lines,
///   edge-case float strings. Both abundance values AND headers are compared.

use proptest::prelude::*;
use std::io::{Cursor, Write};

/// Generate an abundance line with varied separators and value formats.
fn arb_abundance_line() -> impl Strategy<Value = String> {
    let header = "[A-Za-z0-9_.-]{1,30}";
    // All four separators the parser accepts
    let separator = prop_oneof![
        5 => Just("\t".to_string()),
        2 => Just(" ".to_string()),
        1 => Just(",".to_string()),
        1 => Just(";".to_string()),
        // Consecutive separators — exercises the && bug (only first char consumed)
        1 => Just("\t\t".to_string()),
        1 => Just("\t ".to_string()),
        1 => Just("  ".to_string()),
    ];
    let value = prop_oneof![
        // Normal positive floats
        5 => (0.0001f64..1000.0).prop_map(|v| format!("{v}")),
        // Scientific notation (exercises atof vs f64::parse)
        2 => (0.0001f64..1000.0).prop_map(|v| format!("{v:e}")),
        // Very small values (exercises VERY_SMALL_NUM clamping)
        1 => Just("0.00001".to_string()),
        1 => Just("0.0".to_string()),
        // Negative values (should be clamped to VERY_SMALL_NUM by both)
        1 => Just("-1.0".to_string()),
    ];
    (header, separator, value).prop_map(|(h, s, v)| format!("{h}{s}{v}"))
}

fn arb_abundance_file() -> impl Strategy<Value = String> {
    prop::collection::vec(
        prop_oneof![
            9 => arb_abundance_line().prop_map(|l| l + "\n"),
            // Empty lines mid-file (the parser should skip them)
            1 => Just("\n".to_string()),
        ],
        1..50,
    )
    .prop_map(|lines| lines.concat())
    // Must have at least one data line
    .prop_filter("must have at least one data line", |s| {
        s.lines().any(|l| !l.trim().is_empty())
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    #[test]
    fn abundance_equivalence(input in arb_abundance_file()) {
        let tmp_path = std::env::temp_dir().join(format!(
            "maxbin_rs_proptest_abund_{}.txt",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&tmp_path).unwrap();
            f.write_all(input.as_bytes()).unwrap();
        }

        // Rust
        let rust_records = maxbin_rs::abundance::parse(Cursor::new(input.as_bytes())).unwrap();

        // C++
        let original = maxbin_rs::original_ffi::OriginalAbundanceLoader::new(&tmp_path);

        prop_assert!(original.is_parse_success(), "C++ parse failed");
        prop_assert_eq!(
            rust_records.len(),
            original.num_records() as usize,
            "record count mismatch"
        );

        // Collect headers to detect duplicates (abundance_by_header returns the
        // first match, so we can only cross-check unique headers).
        let mut seen_once = std::collections::HashSet::new();
        let mut seen_twice = std::collections::HashSet::new();
        for rec in &rust_records {
            if !seen_once.insert(rec.header.clone()) {
                seen_twice.insert(rec.header.clone());
            }
        }

        for (i, rec) in rust_records.iter().enumerate() {
            // Compare abundance by index
            let orig_abund = original.abundance_by_index(i as i32);
            let diff = (rec.abundance - orig_abund).abs();
            prop_assert!(
                diff < 1e-10,
                "abundance mismatch at record {} (header: {}): rust={} original={} diff={}",
                i, rec.header, rec.abundance, orig_abund, diff
            );

            // Cross-check headers via name lookup (only for unique headers,
            // since abundance_by_header returns the first match for duplicates).
            if !seen_twice.contains(&rec.header) {
                let orig_by_header = original.abundance_by_header(&rec.header);
                let header_diff = (rec.abundance - orig_by_header).abs();
                prop_assert!(
                    header_diff < 1e-10,
                    "header lookup mismatch at record {} (header '{}'): \
                     rust_abund={} cpp_by_header={} — headers likely parsed differently",
                    i, rec.header, rec.abundance, orig_by_header
                );
            }
        }

        std::fs::remove_file(&tmp_path).ok();
    }
}
