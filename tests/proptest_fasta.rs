/// Property-based equivalence tests for fasta.rs.
///
/// Coverage:
/// - parse(): fasta_equivalence (proptest, 200 cases)
///   Generators exercise: headers with descriptions (whitespace truncation),
///   empty sequences, N/IUPAC characters, lowercase, \r\n line endings,
///   junk lines before first header.
use proptest::prelude::*;
use std::io::{Cursor, Write};

/// Generate a FASTA record with realistic variation.
/// - Header may include a description after whitespace (tests truncation logic).
/// - Sequence may contain N and lowercase bases.
/// - Sequence may be empty (consecutive headers).
fn arb_fasta_record() -> impl Strategy<Value = (String, String)> {
    let name = "[A-Za-z0-9_.-]{1,30}";
    let description = prop::option::of("[A-Za-z0-9_ .;:=|/-]{1,60}");
    let header = (name, description).prop_map(|(n, d)| match d {
        Some(desc) => format!("{n} {desc}"),
        None => n,
    });
    let seq = prop_oneof![
        // Normal: 1-5 lines with ACGT, N, and lowercase
        8 => prop::collection::vec("[ACGTacgtNnRYSWKM]{1,80}", 1..5)
                .prop_map(|lines| lines.join("\n")),
        // Empty sequence (header immediately followed by next header or EOF)
        1 => Just(String::new()),
        // Single very short sequence (1-3 bases)
        1 => "[ACGTNn]{1,3}".prop_map(|s| s),
    ];
    (header, seq)
}

fn arb_fasta_file() -> impl Strategy<Value = String> {
    (
        // Optionally prepend junk lines before the first header
        prop::option::of("[A-Za-z0-9 ]{0,40}"),
        prop::collection::vec(arb_fasta_record(), 1..20),
        // Optionally use \r\n line endings
        proptest::bool::ANY,
    )
        .prop_map(|(junk, records, use_crlf)| {
            let newline = if use_crlf { "\r\n" } else { "\n" };
            let mut s = String::new();
            if let Some(junk_line) = junk {
                s.push_str(&junk_line);
                s.push_str(newline);
            }
            for (header, seq) in records {
                s.push('>');
                s.push_str(&header);
                s.push_str(newline);
                if !seq.is_empty() {
                    // Replace \n in multi-line sequences with the chosen line ending
                    s.push_str(&seq.replace('\n', newline));
                    s.push_str(newline);
                }
            }
            s
        })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    #[test]
    fn fasta_equivalence(input in arb_fasta_file()) {
        // Write to temp file for C++ parser
        let tmp_path = std::env::temp_dir().join(format!(
            "maxbin_rs_proptest_fasta_{}.fa",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&tmp_path).unwrap();
            f.write_all(input.as_bytes()).unwrap();
        }

        // Rust parser
        let rust_records = maxbin_rs::fasta::parse(Cursor::new(input.as_bytes()));

        // C++ parser
        let original = maxbin_rs::original_ffi::OriginalFastaReader::new(&tmp_path);

        prop_assert_eq!(
            rust_records.len(),
            original.num_records() as usize,
            "record count mismatch"
        );

        for i in 0..rust_records.len() {
            let idx = i as u32;
            prop_assert_eq!(
                &rust_records[i].header,
                &original.header(idx),
                "header mismatch at record {}", i
            );
            prop_assert_eq!(
                &rust_records[i].seq,
                &original.seq(idx),
                "seq mismatch at record {} (header: {})", i, rust_records[i].header
            );
            // Also check sequence length
            prop_assert_eq!(
                rust_records[i].len() as u32,
                original.seq_len(idx),
                "seq_len mismatch at record {} (header: {})", i, rust_records[i].header
            );
        }

        std::fs::remove_file(&tmp_path).ok();
    }
}
