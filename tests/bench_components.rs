/// Component-level performance comparison: Rust vs. original C++ (FFI).
///
/// For each component, counts how many iterations complete in a fixed time window.
/// Reports ops/sec for both implementations and the ratio.
///
/// Run via:
///   cargo nextest run --cargo-profile bench bench_components -- --ignored
/// (Uses release optimizations for both Rust and C++. Without --cargo-profile bench,
/// both sides run at -O0, which is fair but not representative of real performance.)
///
/// Requires MAXBIN2_TEST_CONTIGS env var for FASTA/profiler/distance benchmarks.
use std::io::Cursor;
use std::time::{Duration, Instant};

const BENCH_DURATION: Duration = Duration::from_secs(3);

/// Run `f` repeatedly for `duration`, return the number of completed iterations.
fn count_ops(duration: Duration, mut f: impl FnMut()) -> u64 {
    let mut count = 0u64;
    let start = Instant::now();
    while start.elapsed() < duration {
        f();
        count += 1;
    }
    count
}

/// Format ops/sec with thousands separators.
fn fmt_ops(ops: f64) -> String {
    let s = format!("{:.0}", ops);
    let bytes: Vec<u8> = s.bytes().collect();
    let mut result = String::new();
    for (i, &b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i) % 3 == 0 {
            result.push(',');
        }
        result.push(b as char);
    }
    result
}

fn print_row(name: &str, cpp_ops: f64, rust_ops: f64) {
    let ratio = rust_ops / cpp_ops;
    eprintln!(
        "  {:<30} {:>14} {:>14}    {:.2}x",
        name,
        fmt_ops(cpp_ops),
        fmt_ops(rust_ops),
        ratio
    );
}

fn print_header() {
    eprintln!();
    eprintln!(
        "  {:<30} {:>14} {:>14}    {}",
        "Component", "C++ (ops/s)", "Rust (ops/s)", "Ratio"
    );
    eprintln!("  {:-<30} {:-<14} {:-<14}    {:-<6}", "", "", "", "");
}

// --- Benchmarks ---

fn bench_normal_distribution() -> (f64, f64) {
    let mean = 0.0676654;
    let std_dev = 0.03419337;
    let inputs: Vec<f64> = (0..1000).map(|i| i as f64 * 0.001).collect();

    let cpp_count = count_ops(BENCH_DURATION, || {
        let nd = maxbin_rs::original_ffi::OriginalNormalDistribution::new(mean, std_dev);
        for &x in &inputs {
            std::hint::black_box(nd.prob(x));
        }
    });

    let rust_count = count_ops(BENCH_DURATION, || {
        let nd = maxbin_rs::normal_distribution::NormalDistribution::new(mean, std_dev);
        for &x in &inputs {
            std::hint::black_box(nd.prob(x));
        }
    });

    let cpp_ops = cpp_count as f64 * inputs.len() as f64 / BENCH_DURATION.as_secs_f64();
    let rust_ops = rust_count as f64 * inputs.len() as f64 / BENCH_DURATION.as_secs_f64();
    (cpp_ops, rust_ops)
}

fn bench_quicksort() -> (f64, f64) {
    // Generate a realistic-sized array (136 entries = profile length for kmerlen=4 symmetric)
    let base: Vec<f64> = (0..136).map(|i| (i as f64 * 0.73).sin()).collect();

    let cpp_count = count_ops(BENCH_DURATION, || {
        let mut arr = base.clone();
        let mut idx: Vec<i32> = (0..arr.len() as i32).collect();
        maxbin_rs::original_ffi::original_quicksort(&mut arr, Some(&mut idx));
        std::hint::black_box(&arr);
    });

    let rust_count = count_ops(BENCH_DURATION, || {
        let mut arr = base.clone();
        let mut idx: Vec<i32> = (0..arr.len() as i32).collect();
        maxbin_rs::quicksort::sort_descending(&mut arr, Some(&mut idx));
        std::hint::black_box(&arr);
    });

    let cpp_ops = cpp_count as f64 / BENCH_DURATION.as_secs_f64();
    let rust_ops = rust_count as f64 / BENCH_DURATION.as_secs_f64();
    (cpp_ops, rust_ops)
}

fn bench_kmermap_lookup() -> (f64, f64) {
    let kmers: Vec<String> = {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut v = Vec::new();
        for &a in &bases {
            for &b in &bases {
                for &c in &bases {
                    for &d in &bases {
                        v.push(String::from_utf8(vec![a, b, c, d]).unwrap());
                    }
                }
            }
        }
        v
    };

    let cpp_count = count_ops(BENCH_DURATION, || {
        let km = maxbin_rs::original_ffi::OriginalKmerMap::new(4, true);
        for kmer in &kmers {
            std::hint::black_box(km.get_mapping(kmer));
        }
    });

    let rust_count = count_ops(BENCH_DURATION, || {
        let km = maxbin_rs::kmer_map::KmerMap::new(4, true);
        for kmer in &kmers {
            std::hint::black_box(km.get_mapping(kmer.as_bytes()));
        }
    });

    let cpp_ops = cpp_count as f64 * kmers.len() as f64 / BENCH_DURATION.as_secs_f64();
    let rust_ops = rust_count as f64 * kmers.len() as f64 / BENCH_DURATION.as_secs_f64();
    (cpp_ops, rust_ops)
}

fn bench_profiler(seqs: &[Vec<u8>]) -> (f64, f64) {
    let cpp_count = count_ops(BENCH_DURATION, || {
        let km = maxbin_rs::original_ffi::OriginalKmerMap::new(4, true);
        for seq in seqs {
            let s = std::str::from_utf8(seq).unwrap();
            let prof = maxbin_rs::original_ffi::OriginalProfiler::new(4, s, &km);
            std::hint::black_box(prof.get_profile(136));
        }
    });

    let rust_count = count_ops(BENCH_DURATION, || {
        let km = maxbin_rs::kmer_map::KmerMap::new(4, true);
        for seq in seqs {
            let prof = maxbin_rs::profiler::Profiler::new(4, seq, &km);
            std::hint::black_box(prof.get_profile());
        }
    });

    let cpp_ops = cpp_count as f64 * seqs.len() as f64 / BENCH_DURATION.as_secs_f64();
    let rust_ops = rust_count as f64 * seqs.len() as f64 / BENCH_DURATION.as_secs_f64();
    (cpp_ops, rust_ops)
}

fn bench_euc_dist(profiles: &[(Vec<f64>, Vec<f64>)]) -> (f64, f64) {
    let cpp_count = count_ops(BENCH_DURATION, || {
        let cpp = maxbin_rs::original_ffi::OriginalEucDist::new(4);
        for (p1, p2) in profiles {
            std::hint::black_box(cpp.get_dist_profile(p1, p2));
        }
    });

    let rust_count = count_ops(BENCH_DURATION, || {
        for (p1, p2) in profiles {
            std::hint::black_box(maxbin_rs::distance::euc_dist_profiles(p1, p2));
        }
    });

    let cpp_ops = cpp_count as f64 * profiles.len() as f64 / BENCH_DURATION.as_secs_f64();
    let rust_ops = rust_count as f64 * profiles.len() as f64 / BENCH_DURATION.as_secs_f64();
    (cpp_ops, rust_ops)
}

fn bench_spearman_dist(profiles: &[(Vec<f64>, Vec<f64>)]) -> (f64, f64) {
    let cpp_count = count_ops(BENCH_DURATION, || {
        let cpp = maxbin_rs::original_ffi::OriginalSpearmanDist::new(4);
        for (p1, p2) in profiles {
            std::hint::black_box(cpp.get_dist_profile(p1, p2));
        }
    });

    let rust_count = count_ops(BENCH_DURATION, || {
        let ctx = maxbin_rs::distance::DistanceContext::new(4);
        for (p1, p2) in profiles {
            std::hint::black_box(maxbin_rs::distance::spearman_dist_profiles(&ctx, p1, p2));
        }
    });

    let cpp_ops = cpp_count as f64 * profiles.len() as f64 / BENCH_DURATION.as_secs_f64();
    let rust_ops = rust_count as f64 * profiles.len() as f64 / BENCH_DURATION.as_secs_f64();
    (cpp_ops, rust_ops)
}

fn bench_fasta_parse(path: &std::path::Path, raw_bytes: &[u8]) -> (f64, f64) {
    let cpp_count = count_ops(BENCH_DURATION, || {
        let reader = maxbin_rs::original_ffi::OriginalFastaReader::new(path);
        std::hint::black_box(reader.num_records());
    });

    let rust_count = count_ops(BENCH_DURATION, || {
        let records = maxbin_rs::fasta::parse(Cursor::new(raw_bytes));
        std::hint::black_box(records.len());
    });

    let cpp_ops = cpp_count as f64 / BENCH_DURATION.as_secs_f64();
    let rust_ops = rust_count as f64 / BENCH_DURATION.as_secs_f64();
    (cpp_ops, rust_ops)
}

fn bench_abundance_parse(path: &std::path::Path, raw_bytes: &[u8]) -> (f64, f64) {
    let cpp_count = count_ops(BENCH_DURATION, || {
        let loader = maxbin_rs::original_ffi::OriginalAbundanceLoader::new(path);
        std::hint::black_box(loader.num_records());
    });

    let rust_count = count_ops(BENCH_DURATION, || {
        let records = maxbin_rs::abundance::parse(Cursor::new(raw_bytes)).unwrap();
        std::hint::black_box(records.len());
    });

    let cpp_ops = cpp_count as f64 / BENCH_DURATION.as_secs_f64();
    let rust_ops = rust_count as f64 / BENCH_DURATION.as_secs_f64();
    (cpp_ops, rust_ops)
}

// --- Main test ---

#[test]
#[ignore] // Run explicitly: cargo nextest run bench_components -- --ignored
fn bench_components() {
    eprintln!();
    eprintln!(
        "=== Component Performance: Rust vs C++ ({:.0}s per component) ===",
        BENCH_DURATION.as_secs_f64()
    );

    print_header();

    // Pure-compute benchmarks (no test data needed)
    let (c, r) = bench_normal_distribution();
    print_row("NormalDist::prob (×1000)", c, r);

    let (c, r) = bench_quicksort();
    print_row("quicksort (136 elems)", c, r);

    let (c, r) = bench_kmermap_lookup();
    print_row("kmermap lookup (×256)", c, r);

    // Data-dependent benchmarks
    let gz_path = match std::env::var("MAXBIN2_TEST_CONTIGS") {
        Ok(p) => std::path::PathBuf::from(p),
        Err(_) => {
            eprintln!();
            eprintln!("  MAXBIN2_TEST_CONTIGS not set — skipping data-dependent benchmarks.");
            eprintln!("  Run inside the devshell for full results.");
            return;
        }
    };

    // Decompress for C++ (can't read gzip)
    let raw_bytes = {
        use std::io::Read;
        let file = std::fs::File::open(&gz_path).unwrap();
        let mut decoder = flate2::read::GzDecoder::new(file);
        let mut buf = Vec::new();
        decoder.read_to_end(&mut buf).unwrap();
        buf
    };
    let plain_path = std::env::temp_dir().join("maxbin_rs_bench_contigs.fa");
    std::fs::write(&plain_path, &raw_bytes).unwrap();

    let (c, r) = bench_fasta_parse(&plain_path, &raw_bytes);
    print_row("FASTA parse (file)", c, r);

    // Build profiles from real contigs for profiler/distance benchmarks
    let records = maxbin_rs::fasta::parse(Cursor::new(&raw_bytes[..]));
    let seqs: Vec<Vec<u8>> = records.iter().map(|r| r.seq.clone()).collect();
    // Take first 50 contigs for profiler benchmark (keeps iteration time reasonable)
    let bench_seqs: Vec<Vec<u8>> = seqs.iter().take(50).cloned().collect();

    let (c, r) = bench_profiler(&bench_seqs);
    print_row("Profiler (×50 contigs)", c, r);

    // Build abundance file for abundance benchmark
    let abund_content: String = records
        .iter()
        .enumerate()
        .map(|(i, rec)| {
            format!(
                "{}\t{:.6}\n",
                rec.header,
                1.0 + (i as f64 * 0.37).sin().abs() * 10.0
            )
        })
        .collect();
    let abund_path = std::env::temp_dir().join("maxbin_rs_bench_abund.txt");
    std::fs::write(&abund_path, &abund_content).unwrap();

    let (c, r) = bench_abundance_parse(&abund_path, abund_content.as_bytes());
    print_row("Abundance parse (file)", c, r);

    // Build profile pairs for distance benchmarks
    let kmap = maxbin_rs::kmer_map::KmerMap::new(4, true);
    let profiles: Vec<Vec<f64>> = bench_seqs
        .iter()
        .map(|seq| {
            let prof = maxbin_rs::profiler::Profiler::new(4, seq, &kmap);
            prof.get_profile().to_vec()
        })
        .collect();
    let profile_pairs: Vec<(Vec<f64>, Vec<f64>)> = profiles
        .windows(2)
        .map(|w| (w[0].clone(), w[1].clone()))
        .collect();

    let (c, r) = bench_euc_dist(&profile_pairs);
    print_row("Euclidean dist (×49 pairs)", c, r);

    let (c, r) = bench_spearman_dist(&profile_pairs);
    print_row("Spearman dist (×49 pairs)", c, r);

    // Cleanup
    std::fs::remove_file(&plain_path).ok();
    std::fs::remove_file(&abund_path).ok();

    eprintln!();
    eprintln!("  Ratio > 1.0 = Rust is faster");
    eprintln!();
}
