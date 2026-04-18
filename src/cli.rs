/// CLI argument parsing.
///
/// v0.3 BREAKING CHANGE: flags now use standard double-dash syntax
/// (`--contig`, `--reads`, etc.) instead of the original's single-dash
/// (`-contig`, `-reads`). Multiple abundance/reads files use repeated
/// `--abund` / `--reads` flags instead of `-abund2`, `-abund3`.
///
/// Subcommands expose individual pipeline stages:
///   maxbin-rs filter       — filter contigs by minimum length
///   maxbin-rs seeds        — generate seed file from HMMER marker gene hits
///   maxbin-rs em           — run the Rust EM algorithm
///   maxbin-rs cpp-em       — run the original C++ EM via FFI (equivalence testing)
///   maxbin-rs sam-to-abund — compute abundance from a SAM file
///   maxbin-rs pipeline     — run the full pipeline (default when no subcommand given)
use std::path::PathBuf;

use clap::{ArgAction, Args, Parser, Subcommand};

const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Parser)]
#[command(
    name = "maxbin-rs",
    version = VERSION,
    about = "MaxBin - a metagenomics binning software.",
    propagate_version = true,
)]
struct ClapCli {
    #[command(subcommand)]
    command: Command,
}

/// Parsed subcommand with its arguments.
#[derive(Debug, Subcommand)]
pub enum Command {
    /// Filter contigs by minimum length.
    Filter(FilterArgs),
    /// Generate seed file from HMMER marker gene hits.
    Seeds(SeedsArgs),
    /// Run the Rust EM algorithm on pre-computed inputs.
    Em(EmArgs),
    /// Run the original C++ EM via FFI (for equivalence testing).
    #[command(name = "cpp-em", alias = "cpp_em")]
    CppEm(CppEmArgs),
    /// Compute abundance from a SAM file.
    #[command(name = "sam-to-abund", alias = "sam_to_abund")]
    SamToAbund(SamToAbundArgs),
    /// Run the full pipeline (default).
    Pipeline(PipelineArgs),
}

#[derive(Debug, Args)]
pub struct FilterArgs {
    #[arg(long)]
    pub contig: PathBuf,
    #[arg(long)]
    pub out: String,
    #[arg(long, default_value_t = 1000)]
    pub min_contig_length: usize,
}

#[derive(Debug, Args)]
pub struct SeedsArgs {
    #[arg(long)]
    pub contig: PathBuf,
    #[arg(long)]
    pub hmmout: PathBuf,
    #[arg(long)]
    pub out: String,
    #[arg(long, default_value_t = 1000)]
    pub min_contig_length: usize,
    #[arg(long, default_value_t = 107)]
    pub markerset: u32,
    #[arg(long, default_value_t = false)]
    pub no_max_effort: bool,
}

#[derive(Debug, Args)]
pub struct EmArgs {
    #[arg(long)]
    pub contig: PathBuf,
    #[arg(long, action = ArgAction::Append)]
    pub abund: Vec<PathBuf>,
    #[arg(long)]
    pub abund_list: Option<PathBuf>,
    #[arg(long)]
    pub seed: PathBuf,
    #[arg(long)]
    pub out: String,
    #[arg(long, default_value_t = 1)]
    pub thread: usize,
    #[arg(long, default_value_t = 50)]
    pub max_iteration: usize,
    /// Probability threshold for EM final classification (default: 0.9).
    /// Note: the original code defaults to 0.5 but the help text says 0.9.
    /// We use 0.9 (the documented value). Use -prob_threshold 0.5 for
    /// bug-for-bug compatibility with the original.
    #[arg(long)]
    pub prob_threshold: Option<f64>,
    #[arg(long, default_value_t = 1000)]
    pub min_contig_length: usize,
}

/// Full pipeline arguments — the original Cli struct, used for both the `pipeline`
/// subcommand and legacy backwards-compatible invocation.
#[derive(Debug, Args)]
pub struct PipelineArgs {
    #[arg(long)]
    pub contig: PathBuf,
    #[arg(long)]
    pub out: String,
    #[arg(long, action = ArgAction::Append)]
    pub abund: Vec<PathBuf>,
    #[arg(long)]
    pub abund_list: Option<PathBuf>,
    #[arg(long, action = ArgAction::Append)]
    pub reads: Vec<PathBuf>,
    #[arg(long)]
    pub reads_list: Option<PathBuf>,
    #[arg(long, default_value_t = 1000)]
    pub min_contig_length: usize,
    #[arg(long, default_value_t = 50)]
    pub max_iteration: usize,
    #[arg(long, default_value_t = 1)]
    pub thread: usize,
    /// Probability threshold for EM final classification (default: 0.9).
    /// Note: the original code defaults to 0.5 but the help text says 0.9.
    /// We use 0.9 (the documented value). Use -prob_threshold 0.5 for
    /// bug-for-bug compatibility with the original.
    #[arg(long)]
    pub prob_threshold: Option<f64>,
    #[arg(long, default_value_t = 107)]
    pub markerset: u32,
    /// Pre-computed protein FASTA (e.g. from Prodigal) — skips gene calling.
    /// HMMER will run on this file instead of calling FragGeneScan.
    #[arg(long)]
    pub faa: Option<PathBuf>,
    #[arg(long)]
    pub plotmarker: bool,
    #[arg(long)]
    pub verbose: bool,
    #[arg(long)]
    pub preserve_intermediate: bool,
    /// Pre-computed HMMER output — skips gene calling and HMMER when provided.
    #[arg(long)]
    pub hmmout: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct CppEmArgs {
    #[arg(long)]
    pub contig: PathBuf,
    #[arg(long)]
    pub abund: PathBuf,
    #[arg(long)]
    pub seed: PathBuf,
    #[arg(long)]
    pub out: String,
    #[arg(long, default_value_t = 1)]
    pub thread: usize,
}

#[derive(Debug, Args)]
pub struct SamToAbundArgs {
    #[arg(long)]
    pub sam: PathBuf,
    #[arg(long)]
    pub out: PathBuf,
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parse CLI arguments. If no subcommand is given, assumes `pipeline`.
pub fn parse() -> Command {
    let args: Vec<String> = std::env::args().collect();
    parse_from(&args[1..])
}

pub fn parse_from(args: &[String]) -> Command {
    // If first arg looks like a flag (not a subcommand), prepend "pipeline"
    let needs_default = args.first().is_some_and(|a| a.starts_with('-'))
        && args[0] != "--version"
        && args[0] != "--help";

    let cli_args: Vec<String> = if needs_default {
        std::iter::once("maxbin-rs".to_string())
            .chain(std::iter::once("pipeline".to_string()))
            .chain(args.iter().cloned())
            .collect()
    } else {
        std::iter::once("maxbin-rs".to_string())
            .chain(args.iter().cloned())
            .collect()
    };

    match ClapCli::try_parse_from(&cli_args) {
        Ok(cli) => cli.command,
        Err(e) => {
            e.print().expect("failed to write to stdout");
            std::process::exit(if e.use_stderr() { 1 } else { 0 });
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Effective prob_threshold: defaults to 0.9 if unset or negative.
///
/// v0.2 BREAKING CHANGE: the original code defaults to 0.5
/// (EManager.cpp:155) but the help text and paper say 0.9. We follow
/// the documented behavior. This changes bin assignments on any run
/// that doesn't explicitly set -prob_threshold.
pub fn effective_prob_threshold(raw: Option<f64>) -> f64 {
    raw.filter(|&v| v >= 0.0).unwrap_or(0.9)
}

/// Merge direct file paths with lines from an optional list file.
fn expand_file_list(direct: &[PathBuf], list: Option<&PathBuf>) -> std::io::Result<Vec<PathBuf>> {
    let mut files = direct.to_vec();
    if let Some(list_path) = list {
        for line in std::fs::read_to_string(list_path)?.lines() {
            let line = line.trim();
            if !line.is_empty() {
                files.push(PathBuf::from(line));
            }
        }
    }
    Ok(files)
}

impl PipelineArgs {
    pub fn prob_threshold(&self) -> f64 {
        effective_prob_threshold(self.prob_threshold)
    }

    pub fn all_abund_files(&self) -> std::io::Result<Vec<PathBuf>> {
        expand_file_list(&self.abund, self.abund_list.as_ref())
    }

    pub fn all_reads_files(&self) -> std::io::Result<Vec<PathBuf>> {
        expand_file_list(&self.reads, self.reads_list.as_ref())
    }

    pub fn validate(&self) -> Result<(), String> {
        if !self.contig.exists() {
            return Err(format!("Contig file not found: {}", self.contig.display()));
        }
        let abund = self.all_abund_files().map_err(|e| e.to_string())?;
        let reads = self.all_reads_files().map_err(|e| e.to_string())?;
        if abund.is_empty() && reads.is_empty() {
            return Err(
                "Please input at least one abundance file or reads file. You may also specify a reads or abundance file list.".into(),
            );
        }
        for f in &abund {
            if !f.exists() {
                return Err(format!("Cannot find abundance file [{}]", f.display()));
            }
        }
        for f in &reads {
            if !f.exists() {
                return Err(format!("Cannot find reads file [{}]", f.display()));
            }
        }
        Ok(())
    }

    pub fn marker_hmm_filename(&self) -> &str {
        if self.markerset == 40 {
            "bacar_marker.hmm"
        } else {
            "marker.hmm"
        }
    }
}

impl EmArgs {
    pub fn prob_threshold(&self) -> f64 {
        effective_prob_threshold(self.prob_threshold)
    }

    pub fn all_abund_files(&self) -> std::io::Result<Vec<PathBuf>> {
        expand_file_list(&self.abund, self.abund_list.as_ref())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn args(strs: &[&str]) -> Vec<String> {
        strs.iter().map(|s| String::from(*s)).collect()
    }

    #[test]
    fn implicit_pipeline_subcommand() {
        let Command::Pipeline(cli) = parse_from(&args(&[
            "--contig", "c.fa", "--reads", "r.fq", "--out", "o",
        ])) else {
            panic!("Expected Pipeline")
        };
        assert_eq!(cli.contig, PathBuf::from("c.fa"));
    }

    #[test]
    fn multiple_abund_flags() {
        let Command::Pipeline(cli) = parse_from(&args(&[
            "--contig", "c.fa", "--abund", "a1.txt", "--abund", "a2.txt", "--out", "o",
        ])) else {
            panic!("Expected Pipeline")
        };
        assert_eq!(cli.abund.len(), 2);
    }

    #[test]
    fn default_prob_threshold() {
        let Command::Pipeline(cli) = parse_from(&args(&[
            "--contig", "c.fa", "--reads", "r.fq", "--out", "o",
        ])) else {
            panic!("Expected Pipeline")
        };
        assert_eq!(cli.prob_threshold(), 0.9);
    }

    #[test]
    fn em_default_prob_threshold() {
        let Command::Em(a) = parse_from(&args(&[
            "em", "--contig", "c.fa", "--abund", "a.txt", "--seed", "s.txt", "--out", "o",
        ])) else {
            panic!("Expected Em")
        };
        assert_eq!(a.prob_threshold(), 0.9);
    }
}
