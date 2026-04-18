/// CLI argument parsing — backwards-compatible with run_MaxBin.pl's flag syntax.
///
/// The original uses single-dash flags (-contig, -reads, -reads2, -reads3, etc.)
/// parsed via a custom loop over @ARGV. We reproduce the exact same interface so
/// that pipelines (nf-core/mag, metaWRAP, atlas, etc.) can substitute maxbin-rs
/// without changing their invocations.
///
/// Subcommands expose individual pipeline stages:
///   maxbin-rs filter     — filter contigs by minimum length
///   maxbin-rs seeds      — generate seed file from HMMER marker gene hits
///   maxbin-rs em         — run the Rust EM algorithm
///   maxbin-rs cpp-em     — run the original C++ EM via FFI (equivalence testing)
///   maxbin-rs sam-to-abund — compute abundance from a SAM file
///   maxbin-rs pipeline   — run the full pipeline (default when no subcommand given)
///
/// When no subcommand is detected, all arguments are treated as the legacy
/// single-dash flag syntax and dispatched to `pipeline` for backwards compatibility.
///
/// ## Preprocessing
///
/// Because the original uses single-dash long flags and numbered flag variants
/// (-reads2, -abund3), raw arguments are preprocessed before clap sees them:
///
/// 1. `-flag` → `--flag` (single-dash long flags → double-dash)
/// 2. `-reads2`, `-abund3` → `--reads`, `--abund` (numbered variants collapsed)
/// 3. `_` → `-` in flag names (underscore/hyphen normalization)
/// 4. `-v` → `--version`
/// 5. If no subcommand detected, `pipeline` is prepended
use std::path::PathBuf;

use clap::{ArgAction, Args, Parser, Subcommand, ValueEnum};

const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Debug, Clone, Copy, PartialEq, ValueEnum)]
pub enum GeneCaller {
    Fraggenescan,
    Prodigal,
}

impl std::fmt::Display for GeneCaller {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Fraggenescan => write!(f, "fraggenescan"),
            Self::Prodigal => write!(f, "prodigal"),
        }
    }
}

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
    #[arg(long, default_value_t = GeneCaller::Fraggenescan)]
    pub gene_caller: GeneCaller,
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
// Preprocessing: massage raw args for clap compatibility
// ---------------------------------------------------------------------------

/// Preprocess CLI arguments for clap compatibility.
///
/// Normalizes single-dash long flags to double-dash, collapses numbered flag
/// variants (-reads2, -abund3) to their base form, and prepends the `pipeline`
/// subcommand when no explicit subcommand is given.
fn preprocess_args(raw: &[String]) -> Vec<String> {
    if raw.is_empty() {
        return vec![];
    }

    let normalized: Vec<String> = raw.iter().map(|a| normalize_arg(a)).collect();
    let first = &normalized[0];

    // If first arg doesn't start with '-', it's a subcommand name — pass through
    if !first.starts_with('-') {
        return normalized;
    }

    // Top-level --version / --help should not get a subcommand prepended
    if first == "--version" || first == "--help" {
        return normalized;
    }

    // Legacy mode: no subcommand given, prepend "pipeline"
    let mut result = Vec::with_capacity(normalized.len() + 1);
    result.push("pipeline".to_string());
    result.extend(normalized);
    result
}

/// Normalize a single CLI argument for clap compatibility.
///
/// - Converts single-dash long flags to double-dash: `-contig` → `--contig`
/// - Normalizes underscores to hyphens: `-min_contig_length` → `--min-contig-length`
/// - Collapses numbered variants: `-reads2` → `--reads`, `-abund3` → `--abund`
/// - Maps `-v` to `--version`
/// - Non-flag arguments pass through unchanged.
fn normalize_arg(arg: &str) -> String {
    if !arg.starts_with('-') {
        return arg.to_string();
    }

    let stripped = arg.trim_start_matches('-');
    if stripped.is_empty() {
        return arg.to_string();
    }

    // -v → --version (matches run_MaxBin.pl's -v flag)
    if stripped == "v" {
        return "--version".to_string();
    }

    // Normalize to underscores for pattern matching, then to hyphens for clap
    let with_underscores = stripped.replace('-', "_");
    let base = collapse_numbered_flag(&with_underscores);
    format!("--{}", base.replace('_', "-"))
}

/// Collapse numbered flag variants to their base form.
///
/// Matches run_MaxBin.pl's regex patterns:
/// - `/^\-reads/` catches -reads, -reads2, -reads3, etc.
/// - `/^\-abund/` catches -abund, -abund2, -abund3, etc.
///
/// Does NOT collapse -reads_list or -abund_list (those are separate flags).
fn collapse_numbered_flag(name: &str) -> String {
    for prefix in &["reads", "abund"] {
        if let Some(suffix) = name.strip_prefix(prefix)
            && !suffix.is_empty()
            && suffix != "_list"
            && suffix.chars().all(|c| c.is_ascii_digit())
        {
            return (*prefix).to_string();
        }
    }
    name.to_string()
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parse CLI arguments, detecting subcommands or falling back to legacy mode.
pub fn parse() -> Command {
    let args: Vec<String> = std::env::args().collect();
    parse_from(&args[1..])
}

pub fn parse_from(args: &[String]) -> Command {
    let preprocessed = preprocess_args(args);
    let cli_args = std::iter::once("maxbin-rs".to_string()).chain(preprocessed);

    match ClapCli::try_parse_from(cli_args) {
        Ok(cli) => cli.command,
        Err(e) => {
            use clap::error::ErrorKind;
            match e.kind() {
                ErrorKind::DisplayHelp | ErrorKind::DisplayVersion => {
                    e.print().expect("failed to write to stdout");
                    std::process::exit(0);
                }
                _ => {
                    // Rewrite clap's "unexpected argument" to "Unrecognized token"
                    // for backwards-compatible error messages.
                    let msg = e.render().to_string();
                    let msg = msg.replace("unexpected argument", "Unrecognized token");
                    eprint!("{msg}");
                    std::process::exit(1);
                }
            }
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
    fn parse_original_syntax() {
        let Command::Pipeline(cli) = parse_from(&args(&[
            "-contig",
            "contigs.fa",
            "-reads",
            "r1.fq",
            "-reads2",
            "r2.fq",
            "-out",
            "output",
            "-thread",
            "4",
        ])) else {
            panic!("Expected Pipeline")
        };
        assert_eq!(cli.contig, PathBuf::from("contigs.fa"));
        assert_eq!(cli.reads.len(), 2);
        assert_eq!(cli.out, "output");
        assert_eq!(cli.thread, 4);
    }

    #[test]
    fn parse_double_dash_syntax() {
        let Command::Pipeline(cli) = parse_from(&args(&[
            "--contig",
            "contigs.fa",
            "--reads",
            "r1.fq",
            "--out",
            "output",
            "--thread",
            "4",
        ])) else {
            panic!("Expected Pipeline")
        };
        assert_eq!(cli.contig, PathBuf::from("contigs.fa"));
        assert_eq!(cli.reads.len(), 1);
        assert_eq!(cli.out, "output");
        assert_eq!(cli.thread, 4);
    }

    #[test]
    fn parse_abund_numbered() {
        let Command::Pipeline(cli) = parse_from(&args(&[
            "-contig", "c.fa", "-abund", "a1.txt", "-abund2", "a2.txt", "-abund3", "a3.txt",
            "-out", "o",
        ])) else {
            panic!("Expected Pipeline")
        };
        assert_eq!(cli.abund.len(), 3);
    }

    #[test]
    fn normalize_args() {
        assert_eq!(normalize_arg("-contig"), "--contig");
        assert_eq!(normalize_arg("--contig"), "--contig");
        assert_eq!(normalize_arg("--min-contig-length"), "--min-contig-length");
        assert_eq!(normalize_arg("-min_contig_length"), "--min-contig-length");
        assert_eq!(
            normalize_arg("-preserve_intermediate"),
            "--preserve-intermediate"
        );
        assert_eq!(normalize_arg("-v"), "--version");
        assert_eq!(normalize_arg("-reads2"), "--reads");
        assert_eq!(normalize_arg("-abund3"), "--abund");
        assert_eq!(normalize_arg("-reads_list"), "--reads-list");
        assert_eq!(normalize_arg("-abund_list"), "--abund-list");
        assert_eq!(normalize_arg("filter"), "filter");
    }

    #[test]
    fn preprocess_inserts_pipeline() {
        assert_eq!(
            preprocess_args(&args(&["-contig", "c.fa", "-out", "o"]))[0],
            "pipeline"
        );
    }

    #[test]
    fn preprocess_preserves_subcommand() {
        assert_eq!(
            preprocess_args(&args(&["filter", "-contig", "c.fa", "-out", "o"]))[0],
            "filter"
        );
    }

    #[test]
    fn preprocess_version_no_pipeline() {
        for flag in ["-version", "--version", "-v"] {
            assert_eq!(
                preprocess_args(&[flag.to_string()]),
                vec!["--version"],
                "flag: {flag}"
            );
        }
    }

    #[test]
    fn default_prob_threshold() {
        let Command::Pipeline(cli) =
            parse_from(&args(&["-contig", "c.fa", "-reads", "r.fq", "-out", "o"]))
        else {
            panic!("Expected Pipeline")
        };
        assert_eq!(cli.prob_threshold(), 0.9);
    }

    #[test]
    fn parse_filter_subcommand() {
        let Command::Filter(a) = parse_from(&args(&[
            "filter",
            "-contig",
            "input.fa.gz",
            "-out",
            "prefix",
            "-min_contig_length",
            "500",
        ])) else {
            panic!("Expected Filter")
        };
        assert_eq!(a.contig, PathBuf::from("input.fa.gz"));
        assert_eq!(a.out, "prefix");
        assert_eq!(a.min_contig_length, 500);
    }

    #[test]
    fn parse_seeds_subcommand() {
        let Command::Seeds(a) = parse_from(&args(&[
            "seeds",
            "-contig",
            "filtered.fa",
            "-hmmout",
            "hits.txt",
            "-out",
            "prefix",
        ])) else {
            panic!("Expected Seeds")
        };
        assert_eq!(a.contig, PathBuf::from("filtered.fa"));
        assert_eq!(a.hmmout, PathBuf::from("hits.txt"));
        assert_eq!(a.out, "prefix");
    }

    #[test]
    fn parse_em_subcommand() {
        let Command::Em(a) = parse_from(&args(&[
            "em",
            "-contig",
            "filtered.fa",
            "-abund",
            "depth.txt",
            "-seed",
            "seeds.txt",
            "-out",
            "prefix",
            "-thread",
            "4",
        ])) else {
            panic!("Expected Em")
        };
        assert_eq!(a.contig, PathBuf::from("filtered.fa"));
        assert_eq!(a.abund.len(), 1);
        assert_eq!(a.seed, PathBuf::from("seeds.txt"));
        assert_eq!(a.out, "prefix");
        assert_eq!(a.thread, 4);
    }

    #[test]
    fn parse_cpp_em_subcommand() {
        let Command::CppEm(a) = parse_from(&args(&[
            "cpp-em", "-contig", "f.fa", "-abund", "d.txt", "-seed", "s.txt", "-out", "o",
        ])) else {
            panic!("Expected CppEm")
        };
        assert_eq!(a.contig, PathBuf::from("f.fa"));
        assert_eq!(a.abund, PathBuf::from("d.txt"));
        assert_eq!(a.seed, PathBuf::from("s.txt"));
        assert_eq!(a.out, "o");
    }

    #[test]
    fn parse_sam_to_abund_subcommand() {
        let Command::SamToAbund(a) = parse_from(&args(&[
            "sam-to-abund",
            "-sam",
            "input.sam",
            "-out",
            "output.txt",
        ])) else {
            panic!("Expected SamToAbund")
        };
        assert_eq!(a.sam, PathBuf::from("input.sam"));
        assert_eq!(a.out, PathBuf::from("output.txt"));
    }

    #[test]
    fn parse_explicit_pipeline_subcommand() {
        let Command::Pipeline(cli) = parse_from(&args(&[
            "pipeline", "-contig", "c.fa", "-reads", "r.fq", "-out", "o",
        ])) else {
            panic!("Expected Pipeline")
        };
        assert_eq!(cli.contig, PathBuf::from("c.fa"));
        assert_eq!(cli.reads.len(), 1);
        assert_eq!(cli.out, "o");
    }

    #[test]
    fn legacy_flags_become_pipeline() {
        assert!(matches!(
            parse_from(&args(&["-contig", "c.fa", "-abund", "a.txt", "-out", "o"])),
            Command::Pipeline(_)
        ));
    }

    #[test]
    fn em_default_prob_threshold() {
        let Command::Em(a) = parse_from(&args(&[
            "em", "-contig", "c.fa", "-abund", "a.txt", "-seed", "s.txt", "-out", "o",
        ])) else {
            panic!("Expected Em")
        };
        assert_eq!(a.prob_threshold(), 0.9);
    }

    #[test]
    fn parse_underscore_subcommand_aliases() {
        assert!(matches!(
            parse_from(&args(&[
                "cpp_em", "-contig", "f.fa", "-abund", "d.txt", "-seed", "s.txt", "-out", "o",
            ])),
            Command::CppEm(_)
        ));
        assert!(matches!(
            parse_from(&args(&["sam_to_abund", "-sam", "i.sam", "-out", "o.txt"])),
            Command::SamToAbund(_)
        ));
    }
}
