//! Reads the EM stage's input files (the pipeline intermediaries) and computes
//! summary metadata for the "what are the inputs" section.
//!
//! Two abundance formats appear in the wild:
//!   - MaxBin 2-column:  `contig_name <tab> depth`  (one sample per file)
//!   - MetaBAT depth:     a header row, then
//!       `contigName  contigLen  totalAvgDepth  <sampleA>  <sampleA-var>  ...`
//!     (many samples in one file; carries each contig's length directly)
//!
//! TEMPORARY: these parsers live in the frontend so we can put real numbers on
//! screen before `maxbin-core` is extracted and made wasm-clean. Parsing is IO,
//! not the EM algorithm, so the brief duplication is acceptable. Once the core
//! is wasm-able the frontend should call its real parsers instead.

/// Length / GC summary of a contig set.
#[derive(Clone, Debug)]
pub struct ContigStats {
    pub n: usize,
    pub total_bp: u64,
    /// GC fraction in [0, 1], only when sequence is available (None if derived
    /// from a depth file's length column, which carries no sequence).
    pub gc: Option<f64>,
    pub len_min: usize,
    pub len_med: usize,
    pub len_max: usize,
}

/// Per-sample depth summary from an abundance file.
#[derive(Clone, Debug)]
pub struct Sample {
    pub label: String,
    pub mean: f64,
    pub min: f64,
    pub max: f64,
}

/// Everything the inputs section needs to describe one dataset's intermediaries.
#[derive(Clone, Debug)]
pub struct Inputs {
    pub label: &'static str,
    pub blurb: &'static str,
    /// Contig stats, when we have them (from FASTA, or from a depth file's
    /// length column). None when the contig set is too large to ship.
    pub contigs: Option<ContigStats>,
    /// True when we hold the actual contig *sequence* (a FASTA), which the EM
    /// needs for tetranucleotide frequencies. False when contig stats came from
    /// a depth file's length column alone (sequence not bundled).
    pub has_sequence: bool,
    /// Where the contig numbers came from / what is not bundled.
    pub contigs_note: String,
    pub abund_format: &'static str,
    pub abund_rows: usize,
    pub samples: Vec<Sample>,
    /// Number of seed contigs (initial bin centers), when a seed file exists.
    pub seed_n: Option<usize>,
    /// Honest characterization of inputs that don't fit in the browser.
    pub notes: Vec<String>,
}

fn median(sorted: &[usize]) -> usize {
    if sorted.is_empty() {
        0
    } else {
        sorted[sorted.len() / 2]
    }
}

fn finish_lengths(mut lens: Vec<usize>, total_bp: u64, gc: Option<f64>) -> ContigStats {
    lens.sort_unstable();
    ContigStats {
        n: lens.len(),
        total_bp,
        gc,
        len_min: lens.first().copied().unwrap_or(0),
        len_med: median(&lens),
        len_max: lens.last().copied().unwrap_or(0),
    }
}

/// Parse a FASTA into contig stats (length distribution + GC).
pub fn fasta_stats(text: &str) -> ContigStats {
    let mut lens = Vec::new();
    let mut total_bp: u64 = 0;
    let mut gc_bp: u64 = 0;
    let mut cur = 0usize;
    let mut open = false;

    for line in text.lines() {
        if line.starts_with('>') {
            if open {
                lens.push(cur);
            }
            cur = 0;
            open = true;
        } else {
            for b in line.bytes() {
                match b {
                    b'A' | b'a' | b'T' | b't' => {
                        cur += 1;
                        total_bp += 1;
                    }
                    b'G' | b'g' | b'C' | b'c' => {
                        cur += 1;
                        total_bp += 1;
                        gc_bp += 1;
                    }
                    b' ' | b'\t' | b'\r' => {}
                    _ => {
                        cur += 1;
                        total_bp += 1;
                    }
                }
            }
        }
    }
    if open {
        lens.push(cur);
    }

    let gc = if total_bp > 0 {
        Some(gc_bp as f64 / total_bp as f64)
    } else {
        None
    };
    finish_lengths(lens, total_bp, gc)
}

/// Parse a MaxBin 2-column abundance file: `name <whitespace> depth`, one
/// sample. Returns (row count, single-sample summary).
pub fn maxbin_abund(text: &str) -> (usize, Sample) {
    let mut depths = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let mut it = line.split_whitespace();
        let _name = it.next();
        if let Some(v) = it.next().and_then(|s| s.parse::<f64>().ok()) {
            depths.push(v);
        }
    }
    (depths.len(), summarize("sample 1", &depths))
}

/// Parse a MetaBAT depth file. Returns (contig stats from the length column,
/// per-sample summaries, row count). Sample columns are the non-`-var` columns
/// after `contigName contigLen totalAvgDepth`.
pub fn metabat_depth(text: &str) -> (ContigStats, Vec<Sample>, usize) {
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    // Column 0=contigName, 1=contigLen, 2=totalAvgDepth, then sample/var pairs.
    let sample_cols: Vec<(usize, String)> = header
        .iter()
        .enumerate()
        .skip(3)
        .filter(|(_, h)| !h.ends_with("-var"))
        .map(|(i, h)| (i, clean_sample_name(h)))
        .collect();

    let mut lens = Vec::new();
    let mut total_bp: u64 = 0;
    let mut cols: Vec<Vec<f64>> = vec![Vec::new(); sample_cols.len()];
    let mut rows = 0;

    for line in lines {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        rows += 1;
        if let Ok(len) = f[1].parse::<usize>() {
            lens.push(len);
            total_bp += len as u64;
        }
        for (slot, (col, _)) in sample_cols.iter().enumerate() {
            if let Some(v) = f.get(*col).and_then(|s| s.parse::<f64>().ok()) {
                cols[slot].push(v);
            }
        }
    }

    let samples = sample_cols
        .iter()
        .zip(&cols)
        .map(|((_, label), vals)| summarize(label, vals))
        .collect();

    (finish_lengths(lens, total_bp, None), samples, rows)
}

/// MetaBAT sample column headers are the full SAM/BAM path; keep the readable
/// tail (e.g. `RH_S001`) so the table doesn't carry the whole filename.
fn clean_sample_name(h: &str) -> String {
    let stem = h.split('/').next_back().unwrap_or(h);
    // ...GoldStandardAssembly.fasta.gz-RH_S001__insert_270.fq.gz.sam.bam
    if let Some(pos) = stem.find("-RH_") {
        if let Some(tail) = stem[pos + 1..].split("__").next() {
            return tail.to_string();
        }
    }
    stem.chars().take(24).collect()
}

fn summarize(label: &str, vals: &[f64]) -> Sample {
    if vals.is_empty() {
        return Sample {
            label: label.to_string(),
            mean: 0.0,
            min: 0.0,
            max: 0.0,
        };
    }
    let sum: f64 = vals.iter().sum();
    let min = vals.iter().copied().fold(f64::INFINITY, f64::min);
    let max = vals.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    Sample {
        label: label.to_string(),
        mean: sum / vals.len() as f64,
        min,
        max,
    }
}

/// Count seed contigs: one name per non-empty line.
pub fn seed_count(text: &str) -> usize {
    text.lines().filter(|l| !l.trim().is_empty()).count()
}
