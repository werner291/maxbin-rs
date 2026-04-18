/// Pipeline orchestration — replaces run_MaxBin.pl.
///
/// This module coordinates the full MaxBin2 pipeline:
/// 1. Filter contigs by minimum length
/// 2. Map reads to contigs (if reads provided instead of abundance)
/// 3. Run gene caller (FragGeneScan or Prodigal) on contigs
/// 4. Run HMMER to find marker genes
/// 5. Recursive binning: seed → EM → re-seed output bins → repeat (up to 5 levels)
/// 6. Merge noclass, filter small bins, sort by abundance
/// 7. Write output files
use crate::cli::{
    CppEmArgs, EmArgs, FilterArgs, GeneCaller, PipelineArgs, SamToAbundArgs, SeedsArgs,
};
use crate::fasta;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Filter contigs by minimum length.
/// Returns (filtered_path, tooshort_path) where filtered_path contains only
/// contigs >= min_length.
///
/// Matches run_MaxBin.pl:366 (checkContig() call): splits contigs into
/// "$out_f.contig.tmp" and "$out_f.tooshort" based on minimum length.
/// Unlike the original, we don't write temp files next to the input
/// (KNOWN ISSUE: the original writes next to the contig file, failing on
/// read-only filesystems).
pub fn filter_contigs(
    contig: &Path,
    out_prefix: &str,
    min_length: usize,
) -> Result<(PathBuf, PathBuf), String> {
    use std::io::BufRead;

    let filtered_path = PathBuf::from(format!("{out_prefix}.contig.tmp"));
    let tooshort_path = PathBuf::from(format!("{out_prefix}.tooshort"));

    let mut filtered = std::io::BufWriter::new(
        std::fs::File::create(&filtered_path)
            .map_err(|e| format!("Can't create {}: {e}", filtered_path.display()))?,
    );
    let mut tooshort = std::io::BufWriter::new(
        std::fs::File::create(&tooshort_path)
            .map_err(|e| format!("Can't create {}: {e}", tooshort_path.display()))?,
    );

    // Matches run_MaxBin.pl:1175-1311 (checkContig): reads the raw FASTA line by line,
    // accumulates sequence, checks length, writes to filtered or tooshort with the
    // ORIGINAL header (including description) and 70-char line wrapping.
    // We must NOT use our FASTA parser here because it truncates headers at whitespace.
    let file = std::fs::File::open(contig).map_err(|e| format!("Can't open contig file: {e}"))?;
    let reader: Box<dyn BufRead> = if contig.extension().and_then(|e| e.to_str()) == Some("gz") {
        // Matches run_MaxBin.pl:1190-1196: gunzip if .gz
        Box::new(std::io::BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(std::io::BufReader::new(file))
    };

    let mut current_header = String::new();
    let mut current_seq = Vec::<u8>::new();
    let mut n_filtered = 0usize;
    let mut n_tooshort = 0usize;

    fn write_fasta_record(dest: &mut impl Write, header: &str, seq: &[u8]) -> Result<(), String> {
        writeln!(dest, "{header}").map_err(|e| e.to_string())?;
        for chunk in seq.chunks(70) {
            dest.write_all(chunk).map_err(|e| e.to_string())?;
            dest.write_all(b"\n").map_err(|e| e.to_string())?;
        }
        Ok(())
    }

    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        let line = line.trim_end_matches(['\n', '\r']);
        if line.starts_with('>') {
            if !current_header.is_empty() {
                if current_seq.len() >= min_length {
                    write_fasta_record(&mut filtered, &current_header, &current_seq)?;
                    n_filtered += 1;
                } else {
                    write_fasta_record(&mut tooshort, &current_header, &current_seq)?;
                    n_tooshort += 1;
                }
                let total = n_filtered + n_tooshort;
                if total.is_multiple_of(10000) {
                    eprintln!(
                        "  filtering contigs: {} kept, {} too short...",
                        n_filtered, n_tooshort
                    );
                }
            }
            current_header = line.to_string();
            current_seq.clear();
        } else if !line.is_empty() {
            current_seq.extend_from_slice(line.as_bytes());
        }
    }
    if !current_header.is_empty() {
        if current_seq.len() >= min_length {
            write_fasta_record(&mut filtered, &current_header, &current_seq)?;
            n_filtered += 1;
        } else {
            write_fasta_record(&mut tooshort, &current_header, &current_seq)?;
            n_tooshort += 1;
        }
    }
    eprintln!(
        "  filtering contigs: {} kept, {} too short (done)",
        n_filtered, n_tooshort
    );

    Ok((filtered_path, tooshort_path))
}

/// Compute abundance from a SAM file.
/// Faithfully reimplements _getabund.pl's getsam() function:
/// - Parses CIGAR strings to sum matched bases (M operations only)
/// - Filters reads by mismatch count (XM:i: tag + indels from CIGAR) <= 10
/// - Abundance = sum(matched_bases) / contig_length
/// - Output sorted by contig name (Perl's `sort keys`)
pub fn compute_abundance_from_sam(sam_path: &Path, out_path: &Path) -> Result<(), String> {
    use std::collections::BTreeMap; // BTreeMap for sorted output like Perl's `sort keys`
    use std::io::BufRead;

    const MAX_MISMATCH: u32 = 10;

    let file = std::fs::File::open(sam_path).map_err(|e| format!("Can't open SAM: {e}"))?;
    let reader = std::io::BufReader::new(file);

    // Matches _getabund.pl:21-22
    let mut contig_lengths: BTreeMap<String, u64> = BTreeMap::new();
    let mut contig_sums: BTreeMap<String, u64> = BTreeMap::new();

    eprintln!("Reading SAM file to estimate abundance values...");

    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        let fields: Vec<&str> = line.split('\t').collect();

        // Matches _getabund.pl:30-38 — @SQ header parsing
        if fields[0].starts_with('@') {
            if fields[0] == "@SQ" && fields.len() >= 3 {
                let name = fields[1].strip_prefix("SN:").unwrap_or("");
                let len: u64 = fields[2]
                    .strip_prefix("LN:")
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
                if !name.is_empty() && len > 0 {
                    contig_lengths.insert(name.to_string(), len);
                    contig_sums.insert(name.to_string(), 0);
                }
            }
            continue;
        }

        if fields.len() < 6 {
            continue;
        }

        let rname = fields[2];
        if rname == "*" {
            continue;
        }

        let cigar = fields[5];

        // Matches _getabund.pl:43-51 — find XM:i: mismatch tag
        let mut mismatch: u32 = MAX_MISMATCH + 1;
        for field in fields.iter().skip(10) {
            if let Some(rest) = field.strip_prefix("XM:i:") {
                mismatch = rest.parse().unwrap_or(MAX_MISMATCH + 1);
                break;
            }
        }

        // Matches _getabund.pl:52 — add indel mismatches from CIGAR
        mismatch = mismatch.saturating_add(get_mismatch_from_cigar(cigar));

        // Matches _getabund.pl:53-56 — filter by max mismatch, sum CIGAR match length
        if mismatch <= MAX_MISMATCH {
            let matched = get_cigar_match_len(cigar);
            *contig_sums.entry(rname.to_string()).or_insert(0) += matched;
        }
    }

    // Matches _getabund.pl:61-67 — output abundance as sum/length, sorted by name
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(out_path).map_err(|e| format!("Can't create abundance file: {e}"))?,
    );
    for (name, sum) in &contig_sums {
        let len = contig_lengths.get(name).copied().unwrap_or(1);
        let abundance = *sum as f64 / len as f64;
        writeln!(out, "{name}\t{abundance}").map_err(|e| e.to_string())?;
    }

    Ok(())
}

/// Sum of M (match) operation lengths in a CIGAR string.
/// Matches _getabund.pl:70-84 (getCIGARLen)
fn get_cigar_match_len(cigar: &str) -> u64 {
    let mut sum: u64 = 0;
    let mut num_str = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num_str.push(ch);
        } else {
            if ch == 'M' {
                sum += num_str.parse::<u64>().unwrap_or(0);
            }
            num_str.clear();
        }
    }
    sum
}

/// Count indel-based mismatches from CIGAR (I and D operations).
/// Matches _getabund.pl:86-117 (getMismatchFromCIAGR)
fn get_mismatch_from_cigar(cigar: &str) -> u32 {
    let mut ret: u32 = 0;
    let mut num_str = String::new();
    let mut ops: Vec<(u32, char)> = Vec::new();

    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num_str.push(ch);
        } else {
            let n: u32 = num_str.parse().unwrap_or(0);
            ops.push((n, ch));
            num_str.clear();
        }
    }

    for &(n, op) in &ops {
        if op == 'I' || op == 'D' {
            ret += n;
        }
    }
    // The Perl also checks for start/end soft clips but the commented-out
    // code ($ret = 10) is never actually executed, so we skip it.
    ret
}

/// Run the `filter` subcommand: filter contigs by minimum length.
pub fn run_filter(args: &FilterArgs) -> Result<(), String> {
    eprintln!(
        "Filtering contigs shorter than {} bp...",
        args.min_contig_length
    );
    let t = std::time::Instant::now();
    let (filtered, tooshort) = filter_contigs(&args.contig, &args.out, args.min_contig_length)?;
    eprintln!("  done in {:.1}s", t.elapsed().as_secs_f64());
    eprintln!("  filtered: {}", filtered.display());
    eprintln!("  tooshort: {}", tooshort.display());
    Ok(())
}

/// Run the `seeds` subcommand: generate seed file from HMMER marker gene hits.
pub fn run_seeds(args: &SeedsArgs) -> Result<(), String> {
    let seed_file = PathBuf::from(format!("{}.seed", args.out));
    eprintln!("Generating seed file from marker genes...");
    let t = std::time::Instant::now();
    let seed_result = generate_seeds(
        &args.hmmout,
        &args.contig,
        args.min_contig_length,
        &seed_file,
        false,
    );
    if seed_result.is_err() && !args.no_max_effort {
        eprintln!("  try harder to dig out marker genes from contigs.");
        generate_seeds(
            &args.hmmout,
            &args.contig,
            args.min_contig_length,
            &seed_file,
            true,
        )?;
    } else if seed_result.is_err() {
        return seed_result;
    }
    let seed_count = std::fs::read_to_string(&seed_file)
        .map(|s| s.lines().filter(|l| !l.trim().is_empty()).count())
        .unwrap_or(0);
    eprintln!(
        "  done in {:.1}s ({} seeds)",
        t.elapsed().as_secs_f64(),
        seed_count
    );
    eprintln!("  seed file: {}", seed_file.display());
    Ok(())
}

/// Run the `em` subcommand: run the Rust EM algorithm on pre-computed inputs.
pub fn run_em(args: &EmArgs) -> Result<(), String> {
    let abund_files = args.all_abund_files().map_err(|e| e.to_string())?;
    if abund_files.is_empty() {
        return Err("No abundance files provided.".into());
    }

    let seed_content =
        std::fs::read_to_string(&args.seed).map_err(|e| format!("Can't read seed file: {e}"))?;
    let seed_names: Vec<String> = seed_content
        .lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();

    if seed_names.len() <= 1 {
        return Err("Need at least 2 seed contigs to run EM.".into());
    }

    let abund_path_refs: Vec<&Path> = abund_files.iter().map(|p| p.as_path()).collect();

    let params = crate::emanager::EmParams {
        min_seq_length: args.min_contig_length,
        max_em: args.max_iteration,
        min_prob_threshold: args.prob_threshold(),
        ..Default::default()
    };

    eprintln!("Loading data...");
    let t_load = std::time::Instant::now();
    let mut state = crate::emanager::init_em(&args.contig, &abund_path_refs, &seed_names, &params);
    eprintln!(
        "  data loaded in {:.1}s ({} contigs, {} seeds, {} abundance files)",
        t_load.elapsed().as_secs_f64(),
        state.records.len(),
        seed_names.len(),
        abund_files.len()
    );

    eprintln!("Running Rust EM...");
    let t_em = std::time::Instant::now();
    crate::emanager::run_em(&mut state, &params);
    let em_secs = t_em.elapsed().as_secs_f64();
    eprintln!("  EM completed in {em_secs:.1}s");

    eprintln!("Classifying contigs...");
    crate::emanager::classify(&mut state, &params);

    eprintln!("Writing results...");
    let t_write = std::time::Instant::now();
    crate::emanager::write_results(&state, &args.out, &params);
    eprintln!(
        "  results written in {:.1}s",
        t_write.elapsed().as_secs_f64()
    );

    let bin_count = rename_bins_to_3digit(&args.out);
    eprintln!("Done: {bin_count} bins produced (EM: {em_secs:.1}s)");
    Ok(())
}

/// Run the `cpp-em` subcommand: run the original C++ EM via FFI.
/// The C++ EManager handles its own data loading, EM, classification,
/// and file writing internally — we time the whole thing as one block.
pub fn run_cpp_em(args: &CppEmArgs) -> Result<(), String> {
    eprintln!("Running original C++ EM via FFI...");
    eprintln!("  contigs: {}", args.contig.display());
    eprintln!("  abund:   {}", args.abund.display());
    eprintln!("  seed:    {}", args.seed.display());
    eprintln!("  threads: {}", args.thread);

    let t = std::time::Instant::now();
    let em = crate::original_ffi::OriginalEManager::new(&args.contig, &args.abund, &args.out);
    em.set_thread_num(args.thread as i32);
    let result = em.run(&args.seed);
    let secs = t.elapsed().as_secs_f64();

    eprintln!("C++ EM completed in {secs:.1}s (returned: {result})");
    Ok(())
}

/// Run the `sam-to-abund` subcommand: compute abundance from a SAM file.
pub fn run_sam_to_abund(args: &SamToAbundArgs) -> Result<(), String> {
    compute_abundance_from_sam(&args.sam, &args.out)?;
    eprintln!("Abundance written to: {}", args.out.display());
    Ok(())
}

/// Run the full pipeline.
pub fn run_pipeline(cli: &PipelineArgs) -> Result<(), String> {
    eprintln!("maxbin-rs {}", env!("CARGO_PKG_VERSION"));

    cli.validate()?;

    // Step 1: Filter contigs
    eprintln!(
        "Filtering contigs shorter than {} bp...",
        cli.min_contig_length
    );
    let (filtered_contigs, _tooshort_path) =
        filter_contigs(&cli.contig, &cli.out, cli.min_contig_length)?;

    // Step 2: Resolve abundance files
    let mut abund_files = cli.all_abund_files().map_err(|e| e.to_string())?;
    let reads_files = cli.all_reads_files().map_err(|e| e.to_string())?;

    if !reads_files.is_empty() {
        eprintln!("Building Bowtie2 index...");
        let idx_prefix = format!("{}.idx", cli.out);
        crate::external::bowtie2_build(&cli.contig, &idx_prefix)?;

        for (i, reads) in reads_files.iter().enumerate() {
            eprintln!("Running Bowtie2 on reads file [{}]...", reads.display());
            let sam_path = PathBuf::from(format!("{}.sam{i}", cli.out));
            let is_fasta = reads
                .extension()
                .and_then(|e| e.to_str())
                .map(|e| matches!(e, "fa" | "fasta" | "fna"))
                .unwrap_or(false);
            crate::external::bowtie2_map(&idx_prefix, reads, &sam_path, cli.thread, is_fasta)?;

            eprintln!("Computing abundance from SAM...");
            let abund_path = PathBuf::from(format!(
                "{}.reads.abund{}",
                filtered_contigs.display(),
                i + 1
            ));
            compute_abundance_from_sam(&sam_path, &abund_path)?;
            abund_files.push(abund_path);
        }
    }

    if abund_files.is_empty() {
        return Err("No abundance information available.".into());
    }

    // Steps 3-4: Gene calling + HMMER marker gene search
    // When -hmmout is provided, skip gene calling and HMMER entirely —
    // use the pre-computed HMMER output directly.
    let marker_hmm = find_marker_hmm(cli)?;
    let (hmmout, faa_path) = if let Some(ref user_hmmout) = cli.hmmout {
        eprintln!("Using pre-computed HMMER output: {}", user_hmmout.display());
        let hmm_hit_count = std::fs::read_to_string(user_hmmout)
            .map(|s| s.lines().filter(|l| !l.starts_with('#')).count())
            .unwrap_or(0);
        eprintln!("  {} hits", hmm_hit_count);
        (user_hmmout.clone(), None)
    } else {
        // Step 3: Gene calling
        let faa = match cli.gene_caller {
            GeneCaller::Fraggenescan => {
                eprintln!("Running FragGeneScan...");
                let prefix = filtered_contigs.display().to_string();
                crate::external::run_fraggenescan(&filtered_contigs, &prefix, cli.thread)?;
                PathBuf::from(format!("{prefix}.frag.faa"))
            }
            GeneCaller::Prodigal => {
                eprintln!("Running Prodigal...");
                let prefix = filtered_contigs.display().to_string();
                crate::external::run_prodigal(&filtered_contigs, &prefix)?;
                PathBuf::from(format!("{prefix}.faa"))
            }
        };

        // Step 4: HMMER marker gene search
        let hmmout_path = PathBuf::from(format!("{}.hmmout", filtered_contigs.display()));
        eprintln!(
            "Running HMMER hmmsearch (marker HMM: {})...",
            marker_hmm.display()
        );
        crate::external::run_hmmsearch(&faa, &marker_hmm, &hmmout_path, cli.thread, cli.markerset)?;

        let hmm_hit_count = std::fs::read_to_string(&hmmout_path)
            .map(|s| s.lines().filter(|l| !l.starts_with('#')).count())
            .unwrap_or(0);
        eprintln!("  HMMER produced {} hits", hmm_hit_count);
        (hmmout_path, Some(faa))
    };

    // Step 5: Recursive binning
    // Matches run_MaxBin.pl:467-579: BFS loop that re-seeds and re-bins each
    // output bin until no more splits are found or recursion depth is reached.
    //
    // The loop reuses the SAME HMMER output from step 4. It calls generate_seeds
    // on each output bin (filtering the HMMER hits to contigs in that bin) to see
    // if the bin contains multiple organisms, and if so, runs EM again on that bin.
    const RECURSION_MAX: usize = 5;
    const MIN_BIN_SIZE: usize = 100_000; // 100 kbp

    let params = crate::emanager::EmParams {
        min_seq_length: cli.min_contig_length,
        max_em: cli.max_iteration,
        min_prob_threshold: cli.prob_threshold(),
        ..Default::default()
    };
    let abund_path_refs: Vec<&Path> = abund_files.iter().map(|p| p.as_path()).collect();

    // BFS queue: (fasta_path, output_prefix, is_initial)
    struct BinTask {
        fasta: PathBuf,
        out_prefix: String,
        is_initial: bool,
    }

    let mut queue: Vec<BinTask> = vec![BinTask {
        fasta: filtered_contigs.clone(),
        out_prefix: cli.out.clone(),
        is_initial: true,
    }];
    let mut cursor = 0;

    // Track which queue entries were re-classified (their bins were split further)
    let mut reclassified: Vec<bool> = vec![false];
    // Collected output bins: (fasta_path, abundance_value)
    let mut all_bins: Vec<(PathBuf, f64)> = Vec::new();
    // Collected noclass files to merge
    let mut noclass_files: Vec<PathBuf> = Vec::new();

    while cursor < queue.len() {
        let task = &queue[cursor];
        let fasta = task.fasta.clone();
        let out_prefix = task.out_prefix.clone();
        let is_initial = task.is_initial;

        // Count recursion depth by counting ".out" in the output prefix
        // Matches run_MaxBin.pl:493: @tmparr = $currbin =~ /[0-9]{4}.out/g;
        let recursion_depth = out_prefix.matches(".out").count();

        // Generate seeds for this bin from the global HMMER output
        let seed_file = PathBuf::from(format!("{out_prefix}.seed"));
        // Only try harder on the initial contigs (matches run_MaxBin.pl:498-505)

        eprintln!(
            "Generating seeds for {} (depth {recursion_depth})...",
            fasta.file_name().unwrap_or_default().to_string_lossy()
        );
        let seed_result = generate_seeds(&hmmout, &fasta, cli.min_contig_length, &seed_file, false);

        let seeds_found = seed_result.is_ok();
        if !seeds_found && is_initial {
            eprintln!("  no seeds (median ≤ 1), retrying with max_effort...");
            let retry = generate_seeds(&hmmout, &fasta, cli.min_contig_length, &seed_file, true);
            if retry.is_err() {
                return Err(
                    "Marker gene search reveals that the dataset cannot be binned \
                     (the median of marker gene number <= 1). Program stop."
                        .into(),
                );
            }
        } else if !seeds_found && !is_initial {
            // No seeds in this sub-bin — skip it (matches run_MaxBin.pl:510-513)
            eprintln!("  no seeds (median ≤ 1), skipping sub-bin");
            cursor += 1;
            continue;
        }

        // Check recursion limit (matches run_MaxBin.pl:507-510)
        if !is_initial && recursion_depth >= RECURSION_MAX {
            eprintln!("  recursion limit reached, skipping");
            cursor += 1;
            continue;
        }

        // Read seeds
        let seed_content = std::fs::read_to_string(&seed_file)
            .map_err(|e| format!("Can't read seed file: {e}"))?;
        let seed_names: Vec<String> = seed_content
            .lines()
            .map(|l| l.trim().to_string())
            .filter(|l| !l.is_empty())
            .collect();

        eprintln!("  {} seed contigs found", seed_names.len());
        if seed_names.len() <= 1 {
            if is_initial {
                return Err("Need at least 2 seed contigs to run EM.".into());
            }
            cursor += 1;
            continue;
        }

        // Mark this entry as reclassified (its output bins replace it)
        reclassified[cursor] = true;

        // Run EM
        eprintln!("Running EM algorithm (depth {recursion_depth})...");
        crate::emanager::run_pipeline(&fasta, &abund_path_refs, &seed_names, &out_prefix, &params);

        // Collect noclass file
        let noclass_path = PathBuf::from(format!("{out_prefix}.noclass"));
        if noclass_path.exists() {
            noclass_files.push(noclass_path);
        }

        // Read summary to find output bins and their abundances,
        // then push each bin onto the queue for potential re-binning.
        // Matches run_MaxBin.pl:561-577
        let summary_path = format!("{out_prefix}.summary");
        if let Ok(summary) = std::fs::read_to_string(&summary_path) {
            for line in summary.lines() {
                // Format: "Bin [path]\tabundance[\tabundance2...]"
                if let Some(rest) = line.strip_prefix("Bin [")
                    && let Some(bracket_end) = rest.find(']')
                {
                    let bin_path_str = &rest[..bracket_end];
                    let abund_str = rest[bracket_end + 1..].trim_start_matches('\t');
                    let abundance: f64 = abund_str
                        .split('\t')
                        .next()
                        .and_then(|s| s.parse().ok())
                        .unwrap_or(0.0);

                    let bin_path = PathBuf::from(bin_path_str);
                    if bin_path.exists() {
                        // Queue this bin for potential re-binning
                        // Output prefix for sub-bins: strip .fasta, add .out
                        let sub_prefix = bin_path_str
                            .strip_suffix(".fasta")
                            .unwrap_or(bin_path_str)
                            .to_string()
                            + ".out";

                        queue.push(BinTask {
                            fasta: bin_path.clone(),
                            out_prefix: sub_prefix,
                            is_initial: false,
                        });
                        reclassified.push(false);
                        all_bins.push((bin_path, abundance));
                    }
                }
            }
        }

        // Clean up intermediate files from this round
        if !cli.preserve_intermediate {
            let _ = std::fs::remove_file(&seed_file);
        }
        if !is_initial {
            let _ = std::fs::remove_file(format!("{out_prefix}.summary"));
            let _ = std::fs::remove_file(format!("{out_prefix}.log"));
        }

        cursor += 1;
    }

    eprintln!(
        "Recursive binning complete: {} total bins before filtering",
        all_bins.len()
    );

    // Step 6: Post-processing
    // Matches run_MaxBin.pl:600-663

    // Remove reclassified bins (they were split into sub-bins)
    // and bins smaller than MIN_BIN_SIZE
    let mut final_bins: Vec<(PathBuf, f64)> = Vec::new();
    let mut small_bin_files: Vec<PathBuf> = Vec::new();

    // Build set of reclassified bin paths
    let reclassified_paths: std::collections::HashSet<PathBuf> = queue
        .iter()
        .enumerate()
        .filter(|(i, _)| reclassified.get(*i).copied().unwrap_or(false))
        .filter_map(|(_, task)| {
            if !task.is_initial {
                Some(task.fasta.clone())
            } else {
                None
            }
        })
        .collect();

    for (path, abundance) in &all_bins {
        if reclassified_paths.contains(path) {
            // This bin was re-classified (split further) — delete it
            let _ = std::fs::remove_file(path);
            continue;
        }

        // Check bin size (sum of sequence lengths)
        let genome_size = get_bin_genome_size(path);
        if genome_size < MIN_BIN_SIZE {
            small_bin_files.push(path.clone());
            continue;
        }

        final_bins.push((path.clone(), *abundance));
    }

    // Sort by abundance descending (matches run_MaxBin.pl:637-641)
    final_bins.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    eprintln!(
        "  {} bins after filtering ({} removed: {} reclassified, {} < {}bp)",
        final_bins.len(),
        all_bins.len() - final_bins.len(),
        reclassified_paths.len(),
        small_bin_files.len(),
        MIN_BIN_SIZE,
    );

    // Rename bins to final sequential numbering
    // Matches run_MaxBin.pl:648-663
    let bin_count = final_bins.len();
    for (i, (old_path, _)) in final_bins.iter().enumerate() {
        let bin_num = format!("{:03}", i + 1);
        let new_path = format!("{}.{bin_num}.fasta", cli.out);
        // Rename via tmp to avoid collisions
        let tmp_path = format!("{new_path}.tmp");
        let _ = std::fs::rename(old_path, &tmp_path);
        let _ = std::fs::rename(&tmp_path, &new_path);
    }

    // Merge all noclass files + small bins into final noclass
    // Matches run_MaxBin.pl:673-710
    // Write to a .tmp file first since the source noclass files may include
    // the final output path (the initial EM round writes to {out}.noclass).
    {
        let noclass_tmp = format!("{}.tmp.noclass", cli.out);
        let noclass_final = format!("{}.noclass", cli.out);
        let mut noclass_out = std::io::BufWriter::new(
            std::fs::File::create(&noclass_tmp)
                .map_err(|e| format!("Can't create noclass: {e}"))?,
        );
        for nc_file in &noclass_files {
            if let Ok(content) = std::fs::read(nc_file) {
                noclass_out.write_all(&content).map_err(|e| e.to_string())?;
            }
            if !cli.preserve_intermediate {
                let _ = std::fs::remove_file(nc_file);
            }
        }
        // Append small bins to noclass
        for small_bin in &small_bin_files {
            if let Ok(content) = std::fs::read(small_bin) {
                noclass_out.write_all(&content).map_err(|e| e.to_string())?;
            }
            let _ = std::fs::remove_file(small_bin);
        }
        drop(noclass_out); // flush before rename
        let _ = std::fs::rename(&noclass_tmp, &noclass_final);
    }

    // Step 7: Rewrite summary with completeness, genome size, GC content
    // Matches run_MaxBin.pl:831-884
    let total_abund_count = abund_files.len();
    rewrite_summary(&cli.out, bin_count, total_abund_count, &hmmout, &marker_hmm)?;

    // Step 8: Write .marker file (marker gene counts per bin)
    // Matches run_MaxBin.pl:818 (countmarker call)
    write_marker_counts(&hmmout, &cli.out, bin_count, &marker_hmm)?;

    // Step 9: Clean up intermediate files (unless --preserve-intermediate)
    if !cli.preserve_intermediate {
        eprintln!("Deleting intermediate files.");
        let _ = std::fs::remove_file(&filtered_contigs);
        let _ = std::fs::remove_file(&hmmout);
        if let Some(ref faa) = faa_path {
            let _ = std::fs::remove_file(faa);
        }
    }

    // Step 10: Write .log file
    // Matches run_MaxBin.pl:390 (openLOG). The original writes all progress
    // output here. We write a minimal log so downstream pipelines (e.g.
    // nf-core/mag) that glob for *.log don't fail.
    {
        let log_path = format!("{}.log", cli.out);
        let mut log = std::io::BufWriter::new(
            std::fs::File::create(&log_path).map_err(|e| format!("Can't create log file: {e}"))?,
        );
        writeln!(log, "maxbin-rs {}", env!("CARGO_PKG_VERSION")).map_err(|e| e.to_string())?;
        writeln!(log, "Input: {}", cli.contig.display()).map_err(|e| e.to_string())?;
        writeln!(log, "Output: {}", cli.out).map_err(|e| e.to_string())?;
        writeln!(log, "Bins: {bin_count}").map_err(|e| e.to_string())?;
    }

    eprintln!("\n========== Job finished ==========");
    Ok(())
}

/// Get total genome size (sum of sequence lengths) for a bin FASTA file.
/// Matches run_MaxBin.pl getBinInfo (genome size component).
fn get_bin_genome_size(path: &Path) -> usize {
    let records = match fasta::parse_file(path) {
        Ok(r) => r,
        Err(_) => return 0,
    };
    records.iter().map(|r| r.len()).sum()
}

/// Find the marker HMM file. Looks relative to the executable, then falls back
/// to common install locations.
fn find_marker_hmm(cli: &PipelineArgs) -> Result<PathBuf, String> {
    let filename = cli.marker_hmm_filename();
    eprintln!("Looking for marker HMM file: {filename}");

    // Check next to the executable
    if let Ok(exe) = std::env::current_exe()
        && let Some(dir) = exe.parent()
    {
        let candidate = dir.join(filename);
        eprintln!("  checking: {}", candidate.display());
        if candidate.exists() {
            eprintln!("  found: {}", candidate.display());
            return Ok(candidate);
        }
    }

    // Check in the current directory
    let candidate = PathBuf::from(filename);
    eprintln!("  checking: {}", candidate.display());
    if candidate.exists() {
        eprintln!("  found: {}", candidate.display());
        return Ok(candidate);
    }

    Err(format!(
        "Cannot find marker HMM file '{filename}'. \
         Place it next to the maxbin-rs binary or in the current directory."
    ))
}

/// Parse HMMER --domtblout output and generate a seed file.
/// Faithfully reimplements _getmarker.pl:10-243 (gethmmmarker() function):
///
/// 1. Parse domtblout, extract contig names from FragGeneScan gene IDs
/// 2. Filter by alignment coverage >= COV_CUTOFF (0.4) — _getmarker.pl:94-96
/// 3. Apply checkMarker() TIGR ID renames — _getmarker.pl:66
/// 4. Group unique contigs by marker gene — _getmarker.pl:73-116
/// 5. Find the median contig-count across markers — _getmarker.pl:144-165
/// 6. Pick the marker with that median count and shortest query length — _getmarker.pl:225-233
/// 7. Output that marker's contigs as seeds — _getmarker.pl:236-241
fn generate_seeds(
    hmmout: &Path,
    contig_file: &Path,
    min_length: usize,
    seed_file: &Path,
    max_effort: bool,
) -> Result<(), String> {
    use std::collections::{HashMap, HashSet};
    use std::io::BufRead;

    // Matches _getmarker.pl:4: `my $COV_CUTOFF = 0.4;`
    const COV_CUTOFF: f64 = 0.4;

    // Read contig lengths from the FASTA file
    let records =
        fasta::parse_file(contig_file).map_err(|e| format!("Failed to parse contigs: {e}"))?;
    let contig_lengths: HashMap<String, usize> = records
        .iter()
        .map(|r| (r.header.clone(), r.len()))
        .collect();

    // Parse HMM domtblout output
    let file = std::fs::File::open(hmmout).map_err(|e| format!("Can't open HMM output: {e}"))?;
    let reader = std::io::BufReader::new(file);

    // For each marker gene: track unique contigs that hit it, and the query length
    struct MarkerInfo {
        /// Unique contigs in insertion order (matching Perl's @queryseq array).
        /// Order matters: seeds are written in this order, which determines bin numbering.
        contigs: Vec<String>,
        query_len: f64,
    }
    let mut markers: Vec<(String, MarkerInfo)> = Vec::new();
    let mut marker_index: HashMap<String, usize> = HashMap::new();

    // The original iterates and groups by marker name as encountered in order.
    // When a new marker name appears, the previous marker's contigs are finalized.
    let mut current_marker = String::new();
    let mut current_contigs: Vec<String> = Vec::new();
    let mut current_contig_set: HashSet<String> = HashSet::new();
    let mut current_query_len: f64 = 0.0;

    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 17 {
            continue;
        }

        let target = fields[0];
        let marker_raw = fields[3];
        let query_len: f64 = fields[5].parse().unwrap_or(0.0);
        let ali_from: f64 = fields[15].parse().unwrap_or(0.0);
        let ali_to: f64 = fields[16].parse().unwrap_or(0.0);

        let marker = check_marker(marker_raw);

        // Extract contig name from gene prediction ID
        let contig = extract_contig_name(target);

        // Check contig exists and meets min length
        match contig_lengths.get(&contig) {
            Some(&len) if len >= min_length => {}
            _ => continue,
        };

        if marker != current_marker {
            // Flush previous marker
            if !current_marker.is_empty() && !current_contigs.is_empty() {
                let idx = markers.len();
                marker_index.insert(current_marker.clone(), idx);
                markers.push((
                    current_marker.clone(),
                    MarkerInfo {
                        contigs: current_contigs.clone(),
                        query_len: current_query_len,
                    },
                ));
            }
            current_marker = marker.clone();
            current_query_len = query_len;
            current_contigs = Vec::new();
            current_contig_set = HashSet::new();
        }

        // Matches _getmarker.pl:94-96: coverage filter — ($arr[16]-$arr[15]) / $arr[5] >= COV_CUTOFF
        let coverage = (ali_to - ali_from) / query_len;
        if coverage >= COV_CUTOFF {
            // Matches _getmarker.pl:97-100: only add if not already in tmphash
            if current_contig_set.insert(contig.clone()) {
                current_contigs.push(contig);
            }
        }
    }
    // Flush the last marker
    if !current_marker.is_empty() && !current_contigs.is_empty() {
        markers.push((
            current_marker,
            MarkerInfo {
                contigs: current_contigs,
                query_len: current_query_len,
            },
        ));
    }

    if markers.is_empty() {
        return Err("No marker genes found. Cannot generate seeds.".into());
    }

    // Get contig counts per marker and find the median.
    // Matches _getmarker.pl lines 144-165.
    let mut contig_counts: Vec<usize> = markers.iter().map(|(_, m)| m.contigs.len()).collect();
    let median_count;

    if max_effort {
        // max_effort=1: sort descending, take median, if it's 1 look for something bigger
        // Matches _getmarker.pl lines 144-159
        contig_counts.sort_by(|a, b| b.cmp(a));
        let raw_median = contig_counts[contig_counts.len() / 2];
        median_count = if raw_median == 1 {
            contig_counts.iter().copied().find(|&c| c > 1).unwrap_or(1)
        } else {
            raw_median
        };

        if median_count <= 1 {
            // Even with max effort, fall back to longest two contigs
            // Matches _getmarker.pl lines 176-221
            let mut sorted_contigs: Vec<(&String, &usize)> = contig_lengths.iter().collect();
            sorted_contigs.sort_by(|a, b| b.1.cmp(a.1));
            let mut out = std::io::BufWriter::new(
                std::fs::File::create(seed_file)
                    .map_err(|e| format!("Can't create seed file: {e}"))?,
            );

            if let Some((name, _)) = sorted_contigs.first() {
                writeln!(out, "{name}").map_err(|e| e.to_string())?;
            }
            if let Some((name, _)) = sorted_contigs.get(1) {
                writeln!(out, "{name}").map_err(|e| e.to_string())?;
            }
            eprintln!("Falling back to longest 2 contigs as seeds");
            return Ok(());
        }
    } else {
        // Normal mode: sort ascending, take median
        // Matches _getmarker.pl lines 162-165
        contig_counts.sort();
        median_count = contig_counts[contig_counts.len() / 2];

        if median_count <= 1 {
            return Err("Marker gene search reveals dataset cannot be binned (median marker gene count <= 1).".into());
        }
    }

    // Find the marker with median count and shortest query length
    // Matches _getmarker.pl lines 225-233
    let mut best_marker_idx = 0;
    let mut best_query_len = f64::MAX;
    // Log candidates (matches _getmarker.pl debug output)
    for (i, (name, info)) in markers.iter().enumerate() {
        if info.contigs.len() == median_count {
            eprintln!(
                "  candidate: marker={} query_len={} contigs={}",
                name,
                info.query_len,
                info.contigs.len()
            );
            if info.query_len < best_query_len {
                best_query_len = info.query_len;
                best_marker_idx = i;
            }
        }
    }

    let seed_contigs = &markers[best_marker_idx].1.contigs;
    eprintln!(
        "  selected marker={} query_len={}",
        markers[best_marker_idx].0, best_query_len
    );

    // Sort seeds to get deterministic order. The original Perl iterates `keys %tmphash`
    // which has implementation-defined ordering. Sorting alphabetically happens to match
    // the original's output for the test data. If this ever diverges, we'd need to
    // match Perl's exact hash iteration order (version-dependent).
    let mut sorted_seeds = seed_contigs.clone();
    sorted_seeds.sort();

    for s in &sorted_seeds {
        eprintln!("  seed: {s}");
    }

    let mut out = std::io::BufWriter::new(
        std::fs::File::create(seed_file).map_err(|e| format!("Can't create seed file: {e}"))?,
    );

    for contig in &sorted_seeds {
        writeln!(out, "{contig}").map_err(|e| e.to_string())?;
    }

    eprintln!(
        "Found {} seed contigs from marker gene '{}' (median count={}, query_len={})",
        seed_contigs.len(),
        markers[best_marker_idx].0,
        median_count,
        best_query_len as i32,
    );
    Ok(())
}

/// Compute genome size (total ACGT bases) and GC content (%) from a FASTA file.
/// Matches run_MaxBin.pl:1314-1343 (getBinInfo)
fn get_bin_info(fasta_path: &Path) -> (usize, f64) {
    let records = match fasta::parse_file(fasta_path) {
        Ok(r) => r,
        Err(_) => return (0, 0.0),
    };
    let mut total = 0usize;
    let mut gc = 0usize;
    for record in &records {
        for &b in &record.seq {
            match b {
                b'A' | b'a' | b'T' | b't' | b'C' | b'c' | b'G' | b'g' => {
                    total += 1;
                    if b == b'C' || b == b'c' || b == b'G' || b == b'g' {
                        gc += 1;
                    }
                }
                _ => {}
            }
        }
    }
    let gc_pct = if total == 0 {
        0.0
    } else {
        gc as f64 / total as f64 * 100.0
    };
    (total, gc_pct)
}

/// Rewrite the summary file to match the Perl's post-processed format.
/// Matches run_MaxBin.pl:832-884
fn rewrite_summary(
    output_prefix: &str,
    bin_count: usize,
    abund_count: usize,
    hmmout: &Path,
    marker_hmm: &Path,
) -> Result<(), String> {
    // Read the C++ summary to get abundance values
    let summary_path = format!("{output_prefix}.summary");
    let raw_summary =
        std::fs::read_to_string(&summary_path).map_err(|e| format!("Can't read summary: {e}"))?;

    // Parse abundance per bin from the C++ summary format: "Bin [name]\tvalue"
    let mut bin_abundances: Vec<String> = Vec::new();
    for line in raw_summary.lines() {
        if line.starts_with("Bin [") {
            let parts: Vec<&str> = line.splitn(2, '\t').collect();
            if parts.len() >= 2 {
                bin_abundances.push(parts[1].to_string());
            }
        }
    }

    // Compute completeness from marker genes
    let completeness = compute_completeness(hmmout, output_prefix, bin_count, marker_hmm)?;

    // Write final summary
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(&summary_path).map_err(|e| format!("Can't write summary: {e}"))?,
    );

    // Matches run_MaxBin.pl:835-841
    if abund_count == 1 {
        writeln!(
            out,
            "Bin name\tAbundance\tCompleteness\tGenome size\tGC content"
        )
        .map_err(|e| e.to_string())?;
    } else {
        writeln!(out, "Bin name\tCompleteness\tGenome size\tGC content")
            .map_err(|e| e.to_string())?;
    }

    let out_name = Path::new(output_prefix)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or(output_prefix);

    for i in 0..bin_count {
        let bin_num = format!("{:03}", i + 1);
        let bin_fasta = format!("{output_prefix}.{bin_num}.fasta");
        let (genome_size, gc) = get_bin_info(Path::new(&bin_fasta));
        let compl = completeness.get(i).copied().unwrap_or(0.0) * 100.0;
        let abund = bin_abundances.get(i).map(|s| s.as_str()).unwrap_or("0.00");

        // Matches run_MaxBin.pl:859-866
        // The Perl formats abundance with printf "%0.2f"
        let abund_fmt: f64 = abund.parse().unwrap_or(0.0);
        if abund_count == 1 {
            writeln!(
                out,
                "{out_name}.{bin_num}.fasta\t{abund_fmt:.2}\t{compl:.1}%\t{genome_size}\t{gc:.1}"
            )
            .map_err(|e| e.to_string())?;
        } else {
            writeln!(
                out,
                "{out_name}.{bin_num}.fasta\t{compl:.1}%\t{genome_size}\t{gc:.1}"
            )
            .map_err(|e| e.to_string())?;
        }
    }

    Ok(())
}

/// Compute completeness (unique_markers / total_markers) per bin.
/// Matches _getmarker.pl:248-408 (countmarker)
fn compute_completeness(
    hmmout: &Path,
    output_prefix: &str,
    bin_count: usize,
    marker_hmm: &Path,
) -> Result<Vec<f64>, String> {
    use std::collections::{HashMap, HashSet};
    use std::io::BufRead;

    // Read marker gene names from the HMM file
    // Matches _getmarker.pl:306-325
    let hmm_file =
        std::fs::File::open(marker_hmm).map_err(|e| format!("Can't open marker HMM: {e}"))?;
    let mut marker_names: Vec<String> = Vec::new();
    let mut marker_index: HashMap<String, usize> = HashMap::new();
    for line in std::io::BufReader::new(hmm_file).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.starts_with("NAME") {
            let name = line.split_whitespace().nth(1).unwrap_or("").to_string();
            let name = check_marker(&name);
            if !marker_index.contains_key(&name) {
                let idx = marker_names.len();
                marker_index.insert(name.clone(), idx);
                marker_names.push(name);
            }
        }
    }
    let gene_num = marker_names.len();

    // Build contig→bin mapping from output FASTA files
    // Matches _getmarker.pl:276-303
    let mut contig_to_bin: HashMap<String, usize> = HashMap::new();
    for i in 0..bin_count {
        let bin_fasta = format!("{}.{:03}.fasta", output_prefix, i + 1);
        if let Ok(records) = fasta::parse_file(Path::new(&bin_fasta)) {
            for r in &records {
                contig_to_bin.insert(r.header.clone(), i);
            }
        }
    }

    // Parse HMM output, count markers per bin
    // Matches _getmarker.pl:327-365
    let mut unique_markers: Vec<HashSet<usize>> = vec![HashSet::new(); bin_count];
    let file = std::fs::File::open(hmmout).map_err(|e| format!("Can't open HMM output: {e}"))?;
    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 17 {
            continue;
        }
        let target = fields[0];
        let marker_raw = fields[3];
        let query_len: f64 = fields[5].parse().unwrap_or(0.0);
        let ali_from: f64 = fields[15].parse().unwrap_or(0.0);
        let ali_to: f64 = fields[16].parse().unwrap_or(0.0);

        let marker = check_marker(marker_raw);
        let contig = extract_contig_name(target);

        let coverage = (ali_to - ali_from) / query_len;
        if coverage < 0.4 {
            continue;
        }

        if let (Some(&bin_idx), Some(&marker_idx)) =
            (contig_to_bin.get(&contig), marker_index.get(&marker))
        {
            unique_markers[bin_idx].insert(marker_idx);
        }
    }

    // Completeness = unique_markers / gene_num
    // Matches _getmarker.pl:379
    let completeness: Vec<f64> = unique_markers
        .iter()
        .map(|s| s.len() as f64 / gene_num as f64)
        .collect();

    Ok(completeness)
}

/// Write the .marker file (marker gene counts per bin).
/// Matches _getmarker.pl:382-405
fn write_marker_counts(
    hmmout: &Path,
    output_prefix: &str,
    bin_count: usize,
    marker_hmm: &Path,
) -> Result<(), String> {
    use std::collections::HashMap;
    use std::io::BufRead;

    // Read marker names from HMM file
    let hmm_file =
        std::fs::File::open(marker_hmm).map_err(|e| format!("Can't open marker HMM: {e}"))?;
    let mut marker_names: Vec<String> = Vec::new();
    let mut marker_index: HashMap<String, usize> = HashMap::new();
    for line in std::io::BufReader::new(hmm_file).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.starts_with("NAME") {
            let name = line.split_whitespace().nth(1).unwrap_or("").to_string();
            let name = check_marker(&name);
            if !marker_index.contains_key(&name) {
                let idx = marker_names.len();
                marker_index.insert(name.clone(), idx);
                marker_names.push(name);
            }
        }
    }
    let gene_num = marker_names.len();

    // Build contig→bin mapping
    let mut contig_to_bin: HashMap<String, usize> = HashMap::new();
    for i in 0..bin_count {
        let bin_fasta = format!("{}.{:03}.fasta", output_prefix, i + 1);
        if let Ok(records) = fasta::parse_file(Path::new(&bin_fasta)) {
            for r in &records {
                contig_to_bin.insert(r.header.clone(), i);
            }
        }
    }

    // Count markers per bin per gene
    let mut counts: Vec<Vec<u32>> = vec![vec![0; gene_num]; bin_count];
    let file = std::fs::File::open(hmmout).map_err(|e| format!("Can't open HMM output: {e}"))?;
    for line in std::io::BufReader::new(file).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 17 {
            continue;
        }
        let target = fields[0];
        let marker_raw = fields[3];
        let query_len: f64 = fields[5].parse().unwrap_or(0.0);
        let ali_from: f64 = fields[15].parse().unwrap_or(0.0);
        let ali_to: f64 = fields[16].parse().unwrap_or(0.0);

        let marker = check_marker(marker_raw);
        let contig = extract_contig_name(target);

        let coverage = (ali_to - ali_from) / query_len;
        if coverage < 0.4 {
            continue;
        }

        if let (Some(&bin_idx), Some(&marker_idx)) =
            (contig_to_bin.get(&contig), marker_index.get(&marker))
        {
            counts[bin_idx][marker_idx] += 1;
        }
    }

    // Write .marker file
    // Matches _getmarker.pl:382-405
    let marker_path = format!("{output_prefix}.marker");
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(&marker_path)
            .map_err(|e| format!("Can't create marker file: {e}"))?,
    );

    let out_name = Path::new(output_prefix)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or(output_prefix);

    // Header
    write!(out, "\tTotal marker\tUnique marker").map_err(|e| e.to_string())?;
    for name in &marker_names {
        write!(out, "\t{name}").map_err(|e| e.to_string())?;
    }
    writeln!(out).map_err(|e| e.to_string())?;

    // Per-bin rows
    for (i, bin_counts) in counts.iter().enumerate().take(bin_count) {
        let bin_num = format!("{:03}", i + 1);
        let total: u32 = bin_counts.iter().sum();
        let unique: u32 = bin_counts.iter().filter(|&&c| c > 0).count() as u32;
        write!(out, "{out_name}.{bin_num}.fasta\t{total}\t{unique}").map_err(|e| e.to_string())?;
        for &c in bin_counts {
            if c > 0 {
                write!(out, "\t{c}").map_err(|e| e.to_string())?;
            } else {
                write!(out, "\t").map_err(|e| e.to_string())?;
            }
        }
        // KNOWN BUG (reproduced): _getmarker.pl:392 uses `$i <= $genenum` instead of
        // `$i < $genenum`, producing one extra trailing tab per row. The header line
        // has the correct number of columns.
        write!(out, "\t").map_err(|e| e.to_string())?;
        writeln!(out).map_err(|e| e.to_string())?;
    }

    Ok(())
}

/// Rename bin files from 4-digit to 3-digit zero-padding.
/// Matches run_MaxBin.pl:732-740 where the Perl renames .0001.fasta → .001.fasta.
/// Also renames the .summary file references.
fn rename_bins_to_3digit(output_prefix: &str) -> usize {
    // Find all .NNNN.fasta files and rename to .NNN
    let mut count = 0;
    for i in 1..=9999 {
        let old = format!("{}.{:04}.fasta", output_prefix, i);
        let new = format!("{}.{:03}.fasta", output_prefix, i);
        if Path::new(&old).exists() {
            let _ = std::fs::rename(&old, &new);
            count += 1;
        } else {
            break;
        }
    }
    count
}

/// Remap certain TIGR marker gene IDs. Matches _getmarker.pl:410-432 (checkMarker()).
fn check_marker(name: &str) -> String {
    // Matches _getmarker.pl:412-427: four explicit TIGR ID remappings
    match name {
        "TIGR00388" => "TIGR00389".to_string(),
        "TIGR00471" => "TIGR00472".to_string(),
        "TIGR00408" => "TIGR00409".to_string(),
        "TIGR02386" => "TIGR02387".to_string(),
        other => other.to_string(),
    }
}

/// Extract contig name from a FragGeneScan gene prediction ID.
/// Format: contigname_start_end_strand → contigname
/// We strip the last 3 underscore-separated parts.
/// Matches _getmarker.pl:69 (regex): `$arr[0] =~ /([...]+)_[0-9]+_[0-9]+_[+-]$/`
/// — capture group $1 is the contig name (everything before _start_end_strand).
fn extract_contig_name(gene_id: &str) -> String {
    let parts: Vec<&str> = gene_id.rsplitn(4, '_').collect();
    if parts.len() >= 4 {
        parts[3].to_string()
    } else {
        // Fallback: return the whole thing
        gene_id.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_contig_name() {
        assert_eq!(extract_contig_name("k141_0_1_903_+"), "k141_0");
        assert_eq!(extract_contig_name("contig_123_45_678_-"), "contig_123");
        assert_eq!(extract_contig_name("simple"), "simple");
    }
}
