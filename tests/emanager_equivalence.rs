/// EM algorithm equivalence test: run both the Rust and C++ EManager on
/// constructed input and compare output files byte-for-byte.
///
/// This test always runs (no external data required). It constructs a
/// minimal dataset with two clearly separable contig groups (AT-rich vs
/// GC-rich at different abundances), so the EM must produce two bins.
///
/// For equivalence testing on real metagenomics data, see the pipeline
/// stage tests: tests/pipeline-stages.sh (run via nix run .#test-pipeline-stages).
use std::io::Write;
use std::path::Path;

fn generate_seq(pattern: &[u8], len: usize) -> Vec<u8> {
    (0..len).map(|i| pattern[i % pattern.len()]).collect()
}

fn write_seed_file(seeds: &[String], path: &Path) {
    let mut f = std::fs::File::create(path).unwrap();
    for seed in seeds {
        writeln!(f, "{}", seed).unwrap();
    }
}

#[test]
fn emanager_cpp_vs_rust() {
    let test_dir = std::env::temp_dir().join("maxbin_rs_em_equiv");
    let _ = std::fs::remove_dir_all(&test_dir);
    std::fs::create_dir_all(&test_dir).unwrap();

    // Two groups of contigs with distinct tetranucleotide profiles and abundances.
    let contig_len = 2000;
    let groups: &[(&str, &[&[u8]], f64)] = &[
        ("A", &[b"AATTAATT", b"AAATATAT", b"ATATATAT"], 5.0),
        ("B", &[b"GGCCGGCC", b"GCGCGCGC", b"GGGCCCGC"], 10.0),
    ];

    let mut fasta = String::new();
    let mut abund_lines = String::new();
    for (group, patterns, base_abund) in groups {
        for (i, pat) in patterns.iter().enumerate() {
            let name = format!("contig{}_{}", group, i);
            let seq = String::from_utf8(generate_seq(pat, contig_len)).unwrap();
            fasta.push_str(&format!(">{}\n{}\n", name, seq));
            let abund = base_abund + i as f64 * 0.1;
            abund_lines.push_str(&format!("{}\t{:.6}\n", name, abund));
        }
    }
    // Short contig — should go to noclass.
    fasta.push_str(">short_contig\nACGTACGT\n");
    abund_lines.push_str("short_contig\t1.000000\n");

    let fasta_path = test_dir.join("test.fa");
    let abund_path = test_dir.join("test.abund");
    let seed_path = test_dir.join("seeds.txt");
    std::fs::write(&fasta_path, &fasta).unwrap();
    std::fs::write(&abund_path, &abund_lines).unwrap();

    let seeds = vec!["contigA_0".to_string(), "contigB_0".to_string()];
    write_seed_file(&seeds, &seed_path);

    // Run C++ EManager via FFI.
    let cpp_prefix = test_dir.join("cpp").to_str().unwrap().to_string();
    {
        let em =
            maxbin_rs::original_ffi::OriginalEManager::new(&fasta_path, &abund_path, &cpp_prefix);
        em.set_thread_num(1);
        assert_eq!(em.run(&seed_path), 0, "C++ EManager returned error");
    }

    // Run Rust EManager.
    let rust_prefix = test_dir.join("rust").to_str().unwrap().to_string();
    {
        let params = maxbin_rs::emanager::EmParams::default();
        let abund_ref: &Path = &abund_path;
        maxbin_rs::emanager::run_pipeline(&fasta_path, &[abund_ref], &seeds, &rust_prefix, &params);
    }

    // Compare output files byte-for-byte.
    for suffix in &["noclass", "0001.fasta", "0002.fasta"] {
        let rust_file = format!("{}.{}", rust_prefix, suffix);
        let cpp_file = format!("{}.{}", cpp_prefix, suffix);
        let rust_bytes = std::fs::read(&rust_file).unwrap_or_default();
        let cpp_bytes = std::fs::read(&cpp_file).unwrap_or_default();
        assert_eq!(
            rust_bytes,
            cpp_bytes,
            "{} differs (rust {} bytes, cpp {} bytes)",
            suffix,
            rust_bytes.len(),
            cpp_bytes.len()
        );
    }

    // Summary file embeds the output path, so compare after normalizing it.
    let rust_summary = std::fs::read_to_string(format!("{}.summary", rust_prefix)).unwrap();
    let cpp_summary = std::fs::read_to_string(format!("{}.summary", cpp_prefix)).unwrap();
    let rust_normalized = rust_summary.replace(&rust_prefix, "PREFIX");
    let cpp_normalized = cpp_summary.replace(&cpp_prefix, "PREFIX");
    assert_eq!(
        rust_normalized, cpp_normalized,
        "summary differs after normalizing paths"
    );

    let _ = std::fs::remove_dir_all(&test_dir);
}
