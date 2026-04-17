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
use std::path::Path;
use std::process::Command;

/// Run FragGeneScan on a contig file to predict genes.
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
        .map_err(|e| format!("Failed to run FragGeneScan: {e}"))?;

    if !status.success() {
        return Err("FragGeneScan failed. Check that run_FragGeneScan.pl is on PATH.".into());
    }

    let faa = format!("{out_prefix}.frag.faa");
    if !Path::new(&faa).exists() {
        return Err(format!("FragGeneScan produced no output ({faa} not found)"));
    }
    Ok(())
}

/// Run Prodigal on a contig file to predict genes.
/// Produces {out_prefix}.faa (amino acid FASTA).
///
/// NOTE: No direct Perl equivalent — the original only supports FragGeneScan.
/// Prodigal support is inspired by mruehlemann/maxbin2_custom fork.
pub fn run_prodigal(contig: &Path, out_prefix: &str) -> Result<(), String> {
    let faa = format!("{out_prefix}.faa");
    let status = Command::new("prodigal")
        .args(["-a", &faa])
        .args(["-i", &contig.display().to_string()])
        .args(["-m", "-o", &format!("{out_prefix}.prodigal.tmp")])
        .args(["-p", "meta", "-q"])
        .status()
        .map_err(|e| format!("Failed to run Prodigal: {e}"))?;

    if !status.success() {
        return Err("Prodigal failed. Check that prodigal is on PATH.".into());
    }

    if !Path::new(&faa).exists() {
        return Err(format!("Prodigal produced no output ({faa} not found)"));
    }
    // Clean up temp file
    let _ = std::fs::remove_file(format!("{out_prefix}.prodigal.tmp"));
    Ok(())
}

/// Run HMMER hmmsearch against marker gene HMM profiles.
/// Produces {out_file} (hmmsearch table output).
///
/// Matches run_MaxBin.pl:~470-490:
///   `hmmsearch --domtblout $hmmout --cut_tc --cpu $thread $marker_hmm $faa`
pub fn run_hmmsearch(
    faa_file: &Path,
    marker_hmm: &Path,
    out_file: &Path,
    threads: usize,
) -> Result<(), String> {
    let status = Command::new("hmmsearch")
        .arg("--domtblout")
        .arg(out_file)
        .arg("--cut_tc")
        .arg("--cpu")
        .arg(threads.to_string())
        .arg(marker_hmm)
        .arg(faa_file)
        .stdout(std::process::Stdio::null())
        .status()
        .map_err(|e| format!("Failed to run hmmsearch: {e}"))?;

    if !status.success() {
        return Err("hmmsearch failed. Check that hmmsearch is on PATH.".into());
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
