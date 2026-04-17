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
/// Matches run_MaxBin.pl:140-290 (argument parsing loop).
use std::path::PathBuf;

const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeneCaller {
    Fraggenescan,
    Prodigal,
}

/// Parsed subcommand with its arguments.
#[derive(Debug)]
pub enum Command {
    /// Filter contigs by minimum length.
    Filter(FilterArgs),
    /// Generate seed file from HMMER marker gene hits.
    Seeds(SeedsArgs),
    /// Run the Rust EM algorithm on pre-computed inputs.
    Em(EmArgs),
    /// Run the original C++ EM via FFI (for equivalence testing).
    CppEm(CppEmArgs),
    /// Compute abundance from a SAM file.
    SamToAbund(SamToAbundArgs),
    /// Run the full pipeline (default).
    Pipeline(PipelineArgs),
}

#[derive(Debug)]
pub struct FilterArgs {
    pub contig: PathBuf,
    pub out: String,
    pub min_contig_length: usize,
}

#[derive(Debug)]
pub struct SeedsArgs {
    pub contig: PathBuf,
    pub hmmout: PathBuf,
    pub out: String,
    pub min_contig_length: usize,
    pub markerset: u32,
}

#[derive(Debug)]
pub struct EmArgs {
    pub contig: PathBuf,
    pub abund: Vec<PathBuf>,
    pub abund_list: Option<PathBuf>,
    pub seed: PathBuf,
    pub out: String,
    pub thread: usize,
    pub max_iteration: usize,
    pub prob_threshold: f64,
    pub min_contig_length: usize,
}

#[derive(Debug)]
pub struct CppEmArgs {
    pub contig: PathBuf,
    pub abund: PathBuf,
    pub seed: PathBuf,
    pub out: String,
    pub thread: usize,
}

#[derive(Debug)]
pub struct SamToAbundArgs {
    pub sam: PathBuf,
    pub out: PathBuf,
}

/// Full pipeline arguments — the original Cli struct, used for both the `pipeline`
/// subcommand and legacy backwards-compatible invocation.
#[derive(Debug)]
pub struct PipelineArgs {
    pub contig: PathBuf,
    pub out: String,
    pub abund: Vec<PathBuf>,
    pub abund_list: Option<PathBuf>,
    pub reads: Vec<PathBuf>,
    pub reads_list: Option<PathBuf>,
    pub min_contig_length: usize,
    pub max_iteration: usize,
    pub thread: usize,
    /// NOTE: the original help text claims 0.9 but the code uses 0.5.
    /// We reproduce the code behavior (0.5) for equivalence.
    pub prob_threshold: f64,
    pub markerset: u32,
    pub gene_caller: GeneCaller,
    pub plotmarker: bool,
    pub verbose: bool,
    pub preserve_intermediate: bool,
}

const USAGE: &str = r#"MaxBin - a metagenomics binning software.
Usage:
  maxbin-rs [subcommand] [options]

Subcommands:
  filter       Filter contigs by minimum length
  seeds        Generate seed file from HMMER marker gene hits
  em           Run the EM algorithm on pre-computed inputs
  cpp-em       Run the original C++ EM via FFI (equivalence testing)
  sam-to-abund Compute abundance from a SAM file
  pipeline     Run the full pipeline (default)

If no subcommand is given, legacy flag syntax is assumed (pipeline mode).

Pipeline options:
  maxbin-rs
    -contig (contig file)
    -out (output file)

   (Input reads and abundance information)
    [-reads (reads file) -reads2 (readsfile) -reads3 (readsfile) -reads4 ... ]
    [-abund (abundance file) -abund2 (abundfile) -abund3 (abundfile) -abund4 ... ]

   (You can also input lists consisting of reads and abundance files)
    [-reads_list (list of reads files)]
    [-abund_list (list of abundance files)]

   (Other parameters)
    [-min_contig_length (minimum contig length. Default 1000)]
    [-max_iteration (maximum Expectation-Maximization algorithm iteration number. Default 50)]
    [-thread (thread num; default 1)]
    [-prob_threshold (probability threshold for EM final classification. Default 0.9)]
    [-plotmarker]
    [-markerset (marker gene sets, 107 (default) or 40.  See README for more information.)]
    [-gene_caller (gene caller: fraggenescan (default) or prodigal)]

  (for debug purpose)
    [-version] [-v] (print version number)
    [-verbose]
    [-preserve_intermediate]

  Please specify either -reads or -abund information.
  You can input multiple reads and/or abundance files at the same time.
"#;

const FILTER_USAGE: &str = r#"maxbin-rs filter — Filter contigs by minimum length.
Usage:
  maxbin-rs filter -contig <input.fa[.gz]> -out <prefix> [-min_contig_length <N>]

Outputs:
  <prefix>.contig.tmp  — filtered contigs (>= min length)
  <prefix>.tooshort    — contigs shorter than min length
"#;

const SEEDS_USAGE: &str = r#"maxbin-rs seeds — Generate seed file from HMMER marker gene hits.
Usage:
  maxbin-rs seeds -contig <filtered.fa> -hmmout <hits.txt> -out <prefix> [-markerset <N>] [-min_contig_length <N>]

Outputs:
  <prefix>.seed — one seed contig name per line
"#;

const EM_USAGE: &str = r#"maxbin-rs em — Run the EM algorithm on pre-computed inputs.
Usage:
  maxbin-rs em -contig <filtered.fa> -abund <depth.txt> [-abund2 ...] -seed <seeds.txt> -out <prefix>
    [-thread <N>] [-max_iteration <N>] [-prob_threshold <F>] [-min_contig_length <N>]

Outputs:
  <prefix>.NNN.fasta — one FASTA per bin
  <prefix>.noclass   — unclassified contigs
  <prefix>.summary   — bin summary
"#;

const CPP_EM_USAGE: &str = r#"maxbin-rs cpp-em — Run the original C++ EM via FFI (equivalence testing).
Usage:
  maxbin-rs cpp-em -contig <filtered.fa> -abund <depth.txt> -seed <seeds.txt> -out <prefix> [-thread <N>]
"#;

const SAM_TO_ABUND_USAGE: &str = r#"maxbin-rs sam-to-abund — Compute abundance from a SAM file.
Usage:
  maxbin-rs sam-to-abund -sam <input.sam> -out <output.txt>
"#;

/// Known subcommand names (used to distinguish subcommands from legacy flags).
const SUBCOMMANDS: &[&str] = &[
    "filter",
    "seeds",
    "em",
    "cpp-em",
    "cpp_em",
    "sam-to-abund",
    "sam_to_abund",
    "pipeline",
];

/// Parse CLI arguments, detecting subcommands or falling back to legacy mode.
pub fn parse() -> Command {
    let args: Vec<String> = std::env::args().collect();
    parse_from(&args[1..])
}

pub fn parse_from(args: &[String]) -> Command {
    if args.is_empty() {
        eprintln!("{USAGE}");
        std::process::exit(1);
    }

    // Check for -version / -v before anything else
    if args.len() == 1 {
        let norm = normalize_flag(&args[0]);
        if norm == "-version" || norm == "-v" {
            eprintln!("MaxBin {VERSION}");
            std::process::exit(0);
        }
    }

    // Check if the first argument is a known subcommand
    let first = &args[0];
    let subcmd = normalize_subcmd(first);
    match subcmd.as_deref() {
        Some("filter") => Command::Filter(parse_filter_args(&args[1..])),
        Some("seeds") => Command::Seeds(parse_seeds_args(&args[1..])),
        Some("em") => Command::Em(parse_em_args(&args[1..])),
        Some("cpp-em") => Command::CppEm(parse_cpp_em_args(&args[1..])),
        Some("sam-to-abund") => Command::SamToAbund(parse_sam_to_abund_args(&args[1..])),
        Some("pipeline") => Command::Pipeline(parse_pipeline_args(&args[1..])),
        _ => {
            // No subcommand detected — treat everything as legacy pipeline flags
            Command::Pipeline(parse_pipeline_args(args))
        }
    }
}

/// Check if an argument is a subcommand name (not a flag).
fn normalize_subcmd(arg: &str) -> Option<String> {
    // Subcommands don't start with a dash
    if arg.starts_with('-') {
        return None;
    }
    let normalized = arg.replace('_', "-");
    if SUBCOMMANDS
        .iter()
        .any(|&s| s.replace('_', "-") == normalized)
    {
        Some(normalized)
    } else {
        None
    }
}

fn parse_filter_args(args: &[String]) -> FilterArgs {
    let mut contig = PathBuf::new();
    let mut out = String::new();
    let mut min_contig_length = 1000usize;

    let mut i = 0;
    while i < args.len() {
        let arg = normalize_flag(&args[i]);
        match arg.as_str() {
            "-contig" => {
                i += 1;
                contig = PathBuf::from(&args[i]);
            }
            "-out" => {
                i += 1;
                out = args[i].clone();
            }
            "-min_contig_length" => {
                i += 1;
                min_contig_length = args[i].parse().unwrap_or(1000);
            }
            "-version" | "-v" => {
                eprintln!("MaxBin {VERSION}");
                std::process::exit(0);
            }
            _ => {
                eprintln!("Unrecognized token [{}]", args[i]);
                eprintln!("{FILTER_USAGE}");
                std::process::exit(1);
            }
        }
        i += 1;
    }

    if contig.as_os_str().is_empty() || out.is_empty() {
        eprintln!("{FILTER_USAGE}");
        std::process::exit(1);
    }

    FilterArgs {
        contig,
        out,
        min_contig_length,
    }
}

fn parse_seeds_args(args: &[String]) -> SeedsArgs {
    let mut contig = PathBuf::new();
    let mut hmmout = PathBuf::new();
    let mut out = String::new();
    let mut min_contig_length = 1000usize;
    let mut markerset = 107u32;

    let mut i = 0;
    while i < args.len() {
        let arg = normalize_flag(&args[i]);
        match arg.as_str() {
            "-contig" => {
                i += 1;
                contig = PathBuf::from(&args[i]);
            }
            "-hmmout" => {
                i += 1;
                hmmout = PathBuf::from(&args[i]);
            }
            "-out" => {
                i += 1;
                out = args[i].clone();
            }
            "-min_contig_length" => {
                i += 1;
                min_contig_length = args[i].parse().unwrap_or(1000);
            }
            "-markerset" => {
                i += 1;
                markerset = args[i].parse().unwrap_or(107);
            }
            "-version" | "-v" => {
                eprintln!("MaxBin {VERSION}");
                std::process::exit(0);
            }
            _ => {
                eprintln!("Unrecognized token [{}]", args[i]);
                eprintln!("{SEEDS_USAGE}");
                std::process::exit(1);
            }
        }
        i += 1;
    }

    if contig.as_os_str().is_empty() || hmmout.as_os_str().is_empty() || out.is_empty() {
        eprintln!("{SEEDS_USAGE}");
        std::process::exit(1);
    }

    SeedsArgs {
        contig,
        hmmout,
        out,
        min_contig_length,
        markerset,
    }
}

fn parse_em_args(args: &[String]) -> EmArgs {
    let mut contig = PathBuf::new();
    let mut abund: Vec<PathBuf> = Vec::new();
    let mut abund_list: Option<PathBuf> = None;
    let mut seed = PathBuf::new();
    let mut out = String::new();
    let mut thread = 1usize;
    let mut max_iteration = 50usize;
    let mut prob_threshold = -1.0f64;
    let mut min_contig_length = 1000usize;

    let mut i = 0;
    while i < args.len() {
        let arg = normalize_flag(&args[i]);
        match arg.as_str() {
            "-contig" => {
                i += 1;
                contig = PathBuf::from(&args[i]);
            }
            "-seed" => {
                i += 1;
                seed = PathBuf::from(&args[i]);
            }
            "-out" => {
                i += 1;
                out = args[i].clone();
            }
            "-thread" => {
                i += 1;
                thread = args[i].parse().unwrap_or(1);
            }
            "-max_iteration" => {
                i += 1;
                max_iteration = args[i].parse().unwrap_or(50);
            }
            "-min_contig_length" => {
                i += 1;
                min_contig_length = args[i].parse().unwrap_or(1000);
            }
            "-abund_list" => {
                i += 1;
                abund_list = Some(PathBuf::from(&args[i]));
            }
            "-prob_threshold" => {
                i += 1;
                let v: f64 = args[i].parse().unwrap_or(-1.0);
                if v >= 0.0 {
                    prob_threshold = v;
                }
            }
            "-version" | "-v" => {
                eprintln!("MaxBin {VERSION}");
                std::process::exit(0);
            }
            // Matches run_MaxBin.pl:178 — regex /^\-abund/ catches -abund, -abund2, etc.
            s if s.starts_with("-abund") => {
                i += 1;
                abund.push(PathBuf::from(&args[i]));
            }
            _ => {
                eprintln!("Unrecognized token [{}]", args[i]);
                eprintln!("{EM_USAGE}");
                std::process::exit(1);
            }
        }
        i += 1;
    }

    // Default prob_threshold to 0.5 if not set (matches EManager.cpp:155)
    if prob_threshold < 0.0 {
        prob_threshold = 0.5;
    }

    if contig.as_os_str().is_empty() || seed.as_os_str().is_empty() || out.is_empty() {
        eprintln!("{EM_USAGE}");
        std::process::exit(1);
    }

    EmArgs {
        contig,
        abund,
        abund_list,
        seed,
        out,
        thread,
        max_iteration,
        prob_threshold,
        min_contig_length,
    }
}

fn parse_cpp_em_args(args: &[String]) -> CppEmArgs {
    let mut contig = PathBuf::new();
    let mut abund = PathBuf::new();
    let mut seed = PathBuf::new();
    let mut out = String::new();
    let mut thread = 1usize;

    let mut i = 0;
    while i < args.len() {
        let arg = normalize_flag(&args[i]);
        match arg.as_str() {
            "-contig" => {
                i += 1;
                contig = PathBuf::from(&args[i]);
            }
            "-abund" => {
                i += 1;
                abund = PathBuf::from(&args[i]);
            }
            "-seed" => {
                i += 1;
                seed = PathBuf::from(&args[i]);
            }
            "-out" => {
                i += 1;
                out = args[i].clone();
            }
            "-thread" => {
                i += 1;
                thread = args[i].parse().unwrap_or(1);
            }
            "-version" | "-v" => {
                eprintln!("MaxBin {VERSION}");
                std::process::exit(0);
            }
            _ => {
                eprintln!("Unrecognized token [{}]", args[i]);
                eprintln!("{CPP_EM_USAGE}");
                std::process::exit(1);
            }
        }
        i += 1;
    }

    if contig.as_os_str().is_empty()
        || abund.as_os_str().is_empty()
        || seed.as_os_str().is_empty()
        || out.is_empty()
    {
        eprintln!("{CPP_EM_USAGE}");
        std::process::exit(1);
    }

    CppEmArgs {
        contig,
        abund,
        seed,
        out,
        thread,
    }
}

fn parse_sam_to_abund_args(args: &[String]) -> SamToAbundArgs {
    let mut sam = PathBuf::new();
    let mut out = PathBuf::new();

    let mut i = 0;
    while i < args.len() {
        let arg = normalize_flag(&args[i]);
        match arg.as_str() {
            "-sam" => {
                i += 1;
                sam = PathBuf::from(&args[i]);
            }
            "-out" => {
                i += 1;
                out = PathBuf::from(&args[i]);
            }
            "-version" | "-v" => {
                eprintln!("MaxBin {VERSION}");
                std::process::exit(0);
            }
            _ => {
                eprintln!("Unrecognized token [{}]", args[i]);
                eprintln!("{SAM_TO_ABUND_USAGE}");
                std::process::exit(1);
            }
        }
        i += 1;
    }

    if sam.as_os_str().is_empty() || out.as_os_str().is_empty() {
        eprintln!("{SAM_TO_ABUND_USAGE}");
        std::process::exit(1);
    }

    SamToAbundArgs { sam, out }
}

fn parse_pipeline_args(args: &[String]) -> PipelineArgs {
    let mut cli = PipelineArgs {
        contig: PathBuf::new(),
        out: String::new(),
        abund: Vec::new(),
        abund_list: None,
        reads: Vec::new(),
        reads_list: None,
        min_contig_length: 1000,
        max_iteration: 50,
        thread: 1,
        prob_threshold: -1.0, // sentinel: not set by user
        markerset: 107,
        gene_caller: GeneCaller::Fraggenescan,
        plotmarker: false,
        verbose: false,
        preserve_intermediate: false,
    };

    let mut i = 0;
    while i < args.len() {
        let arg = normalize_flag(&args[i]);
        match arg.as_str() {
            "-contig" => {
                i += 1;
                cli.contig = PathBuf::from(&args[i]);
            }
            "-out" => {
                i += 1;
                cli.out = args[i].clone();
            }
            "-reads_list" => {
                i += 1;
                cli.reads_list = Some(PathBuf::from(&args[i]));
            }
            "-abund_list" => {
                i += 1;
                cli.abund_list = Some(PathBuf::from(&args[i]));
            }
            "-min_contig_length" => {
                i += 1;
                cli.min_contig_length = args[i].parse().unwrap_or(1000);
            }
            "-max_iteration" => {
                i += 1;
                cli.max_iteration = args[i].parse().unwrap_or(50);
            }
            "-thread" => {
                i += 1;
                cli.thread = args[i].parse().unwrap_or(1);
            }
            "-prob_threshold" => {
                i += 1;
                let v: f64 = args[i].parse().unwrap_or(-1.0);
                if v >= 0.0 {
                    cli.prob_threshold = v;
                }
            }
            "-markerset" => {
                i += 1;
                cli.markerset = args[i].parse().unwrap_or(107);
            }
            "-gene_caller" => {
                i += 1;
                cli.gene_caller = match args[i].as_str() {
                    "prodigal" => GeneCaller::Prodigal,
                    _ => GeneCaller::Fraggenescan,
                };
            }
            "-plotmarker" => cli.plotmarker = true,
            "-verbose" => cli.verbose = true,
            "-preserve_intermediate" => cli.preserve_intermediate = true,
            "-version" | "-v" => {
                eprintln!("MaxBin {VERSION}");
                std::process::exit(0);
            }
            // Matches run_MaxBin.pl:160 — regex /^\-reads/ catches -reads, -reads2, etc.
            s if s.starts_with("-reads") => {
                i += 1;
                cli.reads.push(PathBuf::from(&args[i]));
            }
            // Matches run_MaxBin.pl:178 — regex /^\-abund/ catches -abund, -abund2, etc.
            s if s.starts_with("-abund") => {
                i += 1;
                cli.abund.push(PathBuf::from(&args[i]));
            }
            _ => {
                eprintln!("Unrecognized token [{}]", args[i]);
                eprintln!("{USAGE}");
                std::process::exit(1);
            }
        }
        i += 1;
    }

    // Default prob_threshold to 0.5 if not set (matches EManager.cpp:155,
    // NOT the help text which says 0.9 — known bug)
    if cli.prob_threshold < 0.0 {
        cli.prob_threshold = 0.5;
    }

    if cli.contig.as_os_str().is_empty() {
        eprintln!("No Contig file. Please specify contig file by -contig");
        eprintln!("{USAGE}");
        std::process::exit(1);
    }
    if cli.out.is_empty() {
        eprintln!("Please specify output file by -out.");
        eprintln!("{USAGE}");
        std::process::exit(1);
    }

    cli
}

impl PipelineArgs {
    /// Collect all abundance file paths from -abund flags and -abund_list.
    pub fn all_abund_files(&self) -> std::io::Result<Vec<PathBuf>> {
        let mut files = self.abund.clone();
        if let Some(ref list_path) = self.abund_list {
            let content = std::fs::read_to_string(list_path)?;
            for line in content.lines() {
                let line = line.trim();
                if !line.is_empty() {
                    files.push(PathBuf::from(line));
                }
            }
        }
        Ok(files)
    }

    /// Collect all reads file paths from -reads flags and -reads_list.
    pub fn all_reads_files(&self) -> std::io::Result<Vec<PathBuf>> {
        let mut files = self.reads.clone();
        if let Some(ref list_path) = self.reads_list {
            let content = std::fs::read_to_string(list_path)?;
            for line in content.lines() {
                let line = line.trim();
                if !line.is_empty() {
                    files.push(PathBuf::from(line));
                }
            }
        }
        Ok(files)
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
    /// Collect all abundance file paths from -abund flags and -abund_list.
    pub fn all_abund_files(&self) -> std::io::Result<Vec<PathBuf>> {
        let mut files = self.abund.clone();
        if let Some(ref list_path) = self.abund_list {
            let content = std::fs::read_to_string(list_path)?;
            for line in content.lines() {
                let line = line.trim();
                if !line.is_empty() {
                    files.push(PathBuf::from(line));
                }
            }
        }
        Ok(files)
    }
}

/// Normalize flag: strip leading dashes down to single dash.
/// Accepts both -flag and --flag for forwards compatibility.
/// Also converts --flag-with-hyphens to -flag_with_underscores.
fn normalize_flag(arg: &str) -> String {
    let stripped = arg.trim_start_matches('-');
    let normalized = stripped.replace('-', "_");
    format!("-{normalized}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_original_syntax() {
        let args: Vec<String> = vec![
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
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let cmd = parse_from(&args);
        match cmd {
            Command::Pipeline(cli) => {
                assert_eq!(cli.contig, PathBuf::from("contigs.fa"));
                assert_eq!(cli.reads.len(), 2);
                assert_eq!(cli.out, "output");
                assert_eq!(cli.thread, 4);
            }
            _ => panic!("Expected Pipeline command"),
        }
    }

    #[test]
    fn parse_double_dash_syntax() {
        let args: Vec<String> = vec![
            "--contig",
            "contigs.fa",
            "--reads",
            "r1.fq",
            "--out",
            "output",
            "--thread",
            "4",
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let cmd = parse_from(&args);
        match cmd {
            Command::Pipeline(cli) => {
                assert_eq!(cli.contig, PathBuf::from("contigs.fa"));
                assert_eq!(cli.reads.len(), 1);
                assert_eq!(cli.out, "output");
                assert_eq!(cli.thread, 4);
            }
            _ => panic!("Expected Pipeline command"),
        }
    }

    #[test]
    fn parse_abund_numbered() {
        let args: Vec<String> = vec![
            "-contig", "c.fa", "-abund", "a1.txt", "-abund2", "a2.txt", "-abund3", "a3.txt",
            "-out", "o",
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let cmd = parse_from(&args);
        match cmd {
            Command::Pipeline(cli) => {
                assert_eq!(cli.abund.len(), 3);
            }
            _ => panic!("Expected Pipeline command"),
        }
    }

    #[test]
    fn normalize_flags() {
        assert_eq!(normalize_flag("-contig"), "-contig");
        assert_eq!(normalize_flag("--contig"), "-contig");
        assert_eq!(normalize_flag("--min-contig-length"), "-min_contig_length");
        assert_eq!(normalize_flag("-min_contig_length"), "-min_contig_length");
        assert_eq!(
            normalize_flag("--preserve-intermediate"),
            "-preserve_intermediate"
        );
    }

    #[test]
    fn default_prob_threshold() {
        let args: Vec<String> = vec!["-contig", "c.fa", "-reads", "r.fq", "-out", "o"]
            .into_iter()
            .map(String::from)
            .collect();

        let cmd = parse_from(&args);
        match cmd {
            Command::Pipeline(cli) => {
                // Matches EManager.cpp:155 — default is 0.5, not 0.9
                assert_eq!(cli.prob_threshold, 0.5);
            }
            _ => panic!("Expected Pipeline command"),
        }
    }

    #[test]
    fn parse_filter_subcommand() {
        let args: Vec<String> = vec![
            "filter",
            "-contig",
            "input.fa.gz",
            "-out",
            "prefix",
            "-min_contig_length",
            "500",
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let cmd = parse_from(&args);
        match cmd {
            Command::Filter(a) => {
                assert_eq!(a.contig, PathBuf::from("input.fa.gz"));
                assert_eq!(a.out, "prefix");
                assert_eq!(a.min_contig_length, 500);
            }
            _ => panic!("Expected Filter command"),
        }
    }

    #[test]
    fn parse_seeds_subcommand() {
        let args: Vec<String> = vec![
            "seeds",
            "-contig",
            "filtered.fa",
            "-hmmout",
            "hits.txt",
            "-out",
            "prefix",
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let cmd = parse_from(&args);
        match cmd {
            Command::Seeds(a) => {
                assert_eq!(a.contig, PathBuf::from("filtered.fa"));
                assert_eq!(a.hmmout, PathBuf::from("hits.txt"));
                assert_eq!(a.out, "prefix");
            }
            _ => panic!("Expected Seeds command"),
        }
    }

    #[test]
    fn parse_em_subcommand() {
        let args: Vec<String> = vec![
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
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let cmd = parse_from(&args);
        match cmd {
            Command::Em(a) => {
                assert_eq!(a.contig, PathBuf::from("filtered.fa"));
                assert_eq!(a.abund.len(), 1);
                assert_eq!(a.seed, PathBuf::from("seeds.txt"));
                assert_eq!(a.out, "prefix");
                assert_eq!(a.thread, 4);
            }
            _ => panic!("Expected Em command"),
        }
    }

    #[test]
    fn parse_cpp_em_subcommand() {
        let args: Vec<String> = vec![
            "cpp-em", "-contig", "f.fa", "-abund", "d.txt", "-seed", "s.txt", "-out", "o",
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let cmd = parse_from(&args);
        match cmd {
            Command::CppEm(a) => {
                assert_eq!(a.contig, PathBuf::from("f.fa"));
                assert_eq!(a.abund, PathBuf::from("d.txt"));
                assert_eq!(a.seed, PathBuf::from("s.txt"));
                assert_eq!(a.out, "o");
            }
            _ => panic!("Expected CppEm command"),
        }
    }

    #[test]
    fn parse_sam_to_abund_subcommand() {
        let args: Vec<String> = vec!["sam-to-abund", "-sam", "input.sam", "-out", "output.txt"]
            .into_iter()
            .map(String::from)
            .collect();

        let cmd = parse_from(&args);
        match cmd {
            Command::SamToAbund(a) => {
                assert_eq!(a.sam, PathBuf::from("input.sam"));
                assert_eq!(a.out, PathBuf::from("output.txt"));
            }
            _ => panic!("Expected SamToAbund command"),
        }
    }

    #[test]
    fn parse_explicit_pipeline_subcommand() {
        let args: Vec<String> = vec!["pipeline", "-contig", "c.fa", "-reads", "r.fq", "-out", "o"]
            .into_iter()
            .map(String::from)
            .collect();

        let cmd = parse_from(&args);
        match cmd {
            Command::Pipeline(cli) => {
                assert_eq!(cli.contig, PathBuf::from("c.fa"));
                assert_eq!(cli.reads.len(), 1);
                assert_eq!(cli.out, "o");
            }
            _ => panic!("Expected Pipeline command"),
        }
    }

    #[test]
    fn legacy_flags_become_pipeline() {
        // Flags without a subcommand should parse as Pipeline
        let args: Vec<String> = vec!["-contig", "c.fa", "-abund", "a.txt", "-out", "o"]
            .into_iter()
            .map(String::from)
            .collect();

        let cmd = parse_from(&args);
        assert!(matches!(cmd, Command::Pipeline(_)));
    }

    #[test]
    fn em_default_prob_threshold() {
        let args: Vec<String> = vec![
            "em", "-contig", "c.fa", "-abund", "a.txt", "-seed", "s.txt", "-out", "o",
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let cmd = parse_from(&args);
        match cmd {
            Command::Em(a) => {
                assert_eq!(a.prob_threshold, 0.5);
            }
            _ => panic!("Expected Em command"),
        }
    }
}
