//! External tool invocations — shelling out to HMMER, Bowtie2, and gene callers.
//!
//! These replace the Perl orchestration's subprocess calls in run_MaxBin.pl.
//! Each function wraps a single external tool, handling argument construction
//! and error reporting.
//!
//! In the original, run_MaxBin.pl resolves tool paths via:
//!   1. A `setting` file with manual absolute paths (run_MaxBin.pl:~295-360,
//!      checkSetting/checkProgram subroutines)
//!   2. Fallback to $PATH search
//!
//! Rust finds tools on $PATH only, eliminating the setting file entirely.
//!
//! KNOWN ISSUES preserved from the original:
//! - FragGeneScan is designed for short reads, not assembled contigs (SourceForge
//!   ticket #2). Using it on contigs is a misuse that causes slow initialization
//!   and macOS segfaults. Prodigal is the correct tool for assembled contigs.
//! - The original resolves tool paths via a hand-written `setting` file with
//!   manual absolute paths — a constant source of user pain (Biostars #9473674).
//!   We find tools on $PATH instead.
use std::path::{Path, PathBuf};
use std::process::Command;

/// Run FragGeneScanRs (Rust library) on a contig file to predict genes.
/// Produces {out_prefix}.frag.faa (amino acid FASTA).
///
/// This calls the FragGeneScanRs library directly — no subprocess, no Perl
/// wrapper. The training data directory must contain the standard FGS training
/// files (complete, gene, rgene, noncoding, start, stop, start1, stop1, pwm).
pub fn run_fraggenescan_rs(
    contig: &Path,
    out_prefix: &str,
    threads: usize,
    train_dir: &Path,
) -> Result<(), String> {
    use frag_gene_scan_rs::dna::{Nuc, count_cg_content};
    use frag_gene_scan_rs::hmm;
    use frag_gene_scan_rs::viterbi::viterbi;
    use std::io::Write;

    let (global, locals) =
        hmm::get_train_from_file(train_dir.to_path_buf(), PathBuf::from("complete")).map_err(
            |e| {
                format!(
                    "Failed to load FGSrs training data from {}: {e}",
                    train_dir.display()
                )
            },
        )?;

    let faa_path = format!("{out_prefix}.frag.faa");
    let mut faa_out = std::io::BufWriter::new(
        std::fs::File::create(&faa_path).map_err(|e| format!("Can't create {faa_path}: {e}"))?,
    );

    let file = std::fs::File::open(contig).map_err(|e| format!("Can't open contig file: {e}"))?;
    let reader = seq_io::fasta::Reader::new(file);

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| format!("Failed to create thread pool: {e}"))?;

    // Collect predictions in order (matching original FGS output ordering).
    let records: Vec<_> = reader
        .into_records()
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| format!("FASTA parse error: {e}"))?;

    let predictions: Vec<Vec<u8>> = pool.install(|| {
        use rayon::iter::{IntoParallelIterator, ParallelIterator};
        records
            .into_par_iter()
            .map(|record| {
                let head: Vec<u8> = record
                    .head
                    .iter()
                    .copied()
                    .take_while(u8::is_ascii_graphic)
                    .collect();
                let nseq: Vec<Nuc> = record.seq.iter().copied().map(Nuc::from).collect();
                let mut buf = Vec::new();
                if !nseq.is_empty() {
                    let prediction = viterbi(
                        &global,
                        &locals[count_cg_content(&nseq)],
                        head,
                        nseq,
                        false, // complete=0: short reads / contigs, not whole genome
                    );
                    // Ignore error — matches FGSrs binary behavior
                    let _ = prediction.protein(&mut buf, false);
                }
                buf
            })
            .collect()
    });

    for buf in predictions {
        faa_out
            .write_all(&buf)
            .map_err(|e| format!("Write error: {e}"))?;
    }
    faa_out.flush().map_err(|e| format!("Flush error: {e}"))?;

    if !Path::new(&faa_path).exists() {
        return Err(format!(
            "FragGeneScanRs produced no output ({faa_path} not found)"
        ));
    }
    Ok(())
}

/// Run FragGeneScan (original C) on a contig file to predict genes.
/// Produces {out_prefix}.frag.faa (amino acid FASTA).
///
/// Matches run_MaxBin.pl:~440-460:
///   `run_FragGeneScan.pl -genome=$contig_tmp -out=$contig_tmp.frag -complete=0
///    -train=complete -thread=$thread`
///
/// KNOWN ISSUE: FragGeneScan is optimized for short reads but MaxBin2 feeds it
/// assembled contigs. This is a misuse inherited from the original. Use
/// `run_prodigal` for better results on contigs.
pub fn run_fraggenescan(contig: &Path, out_prefix: &str, threads: usize) -> Result<(), String> {
    let status = Command::new("run_FragGeneScan.pl")
        .arg(format!("-genome={}", contig.display()))
        .arg(format!("-out={out_prefix}.frag"))
        .arg("-complete=0")
        .arg("-train=complete")
        .arg(format!("-thread={threads}"))
        .status()
        .map_err(|e| {
            format!(
                "Failed to run FragGeneScan: {e}\n\
            Hint: use --faa to skip gene calling, or install the full package \
            (nix run github:werner291/maxbin-rs) which bundles all dependencies."
            )
        })?;

    if !status.success() {
        return Err("FragGeneScan failed. Check that run_FragGeneScan.pl is on PATH.".into());
    }

    let faa = format!("{out_prefix}.frag.faa");
    if !Path::new(&faa).exists() {
        return Err(format!("FragGeneScan produced no output ({faa} not found)"));
    }
    Ok(())
}

/// Run HMMER hmmsearch against marker gene HMM profiles.
/// Produces {out_file} (hmmsearch table output).
///
/// Matches run_MaxBin.pl:~1076-1093:
///   markerset=107: `hmmsearch --domtblout $hmmout --cut_tc --cpu $thread $marker_hmm $faa`
///   markerset=40:  `hmmsearch --domtblout $hmmout -E 1e-3 --cpu $thread $marker_hmm $faa`
///
/// The 107-gene marker set (marker.hmm) has TC bit score thresholds; the
/// 40-gene set (bacar_marker.hmm) does not, so the original uses E-value
/// cutoff instead.
pub fn run_hmmsearch(
    faa_file: &Path,
    marker_hmm: &Path,
    out_file: &Path,
    threads: usize,
    markerset: u32,
) -> Result<(), String> {
    let mut cmd = Command::new("hmmsearch");
    cmd.arg("--domtblout").arg(out_file);

    if markerset == 107 {
        cmd.arg("--cut_tc");
    } else {
        cmd.args(["-E", "1e-3"]);
    }

    cmd.arg("--cpu")
        .arg(threads.to_string())
        .arg(marker_hmm)
        .arg(faa_file)
        .stdout(std::process::Stdio::null());

    let output = cmd.output().map_err(|e| {
        format!(
            "Failed to run hmmsearch: {e}\n\
            Hint: use --hmmout to skip HMMER, or install the full package \
            (nix run github:werner291/maxbin-rs) which bundles all dependencies."
        )
    })?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(format!("hmmsearch failed:\n{stderr}"));
    }
    Ok(())
}

/// Build a Bowtie2 index from a contig file.
/// Matches run_MaxBin.pl:~385-395: `bowtie2-build $contig $out_f.idx`
pub fn bowtie2_build(contig: &Path, index_prefix: &str) -> Result<(), String> {
    let output = Command::new("bowtie2-build")
        .arg(contig)
        .arg(index_prefix)
        .stdout(std::process::Stdio::null())
        .output()
        .map_err(|e| format!("Failed to run bowtie2-build: {e}"))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(format!("bowtie2-build failed:\n{stderr}"));
    }
    Ok(())
}

/// Map reads to contigs using Bowtie2, producing a SAM file.
/// Matches run_MaxBin.pl:~400-415: `bowtie2 -p $thread -x $idx -U $reads -S $sam`
///
/// DISCREPANCY: The Perl also checks if the input is FASTA vs FASTQ by file
/// extension and adds `-f` flag for FASTA. Rust reproduces this logic.
pub fn bowtie2_map(
    index_prefix: &str,
    reads: &Path,
    sam_out: &Path,
    threads: usize,
    is_fasta: bool,
) -> Result<(), String> {
    let mut cmd = Command::new("bowtie2");
    cmd.arg("-p")
        .arg(threads.to_string())
        .arg("-x")
        .arg(index_prefix)
        .arg("-U")
        .arg(reads)
        .arg("-S")
        .arg(sam_out);

    if is_fasta {
        cmd.arg("-f");
    }

    cmd.stdout(std::process::Stdio::null());

    let output = cmd
        .output()
        .map_err(|e| format!("Failed to run bowtie2: {e}"))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(format!("bowtie2 mapping failed:\n{stderr}"));
    }
    Ok(())
}
