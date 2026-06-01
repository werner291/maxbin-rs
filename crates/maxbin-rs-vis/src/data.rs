//! Client-side parsing of the EM subcommand's inputs: contigs (FASTA),
//! abundance (per-contig depth), and seed (contig names).
//!
//! TEMPORARY: these are minimal parsers living in the frontend so we can get
//! real data on screen before `maxbin-core` is extracted and made wasm-clean.
//! Once it is, the frontend should call the core's real parsers (which carry
//! the bug-for-bug abundance-format behavior) instead of these. Parsing is IO,
//! not the algorithm, so the brief duplication is acceptable.

use std::collections::HashSet;

/// One contig with the fields the first plot needs.
#[derive(Clone, Debug)]
pub struct Contig {
    pub name: String,
    pub length: usize,
    /// GC fraction in [0, 1].
    pub gc: f64,
    /// Mean depth from the (first) abundance sample. 0.0 if absent.
    pub abundance: f64,
    /// True if this contig appears in the seed file (an initial bin center).
    pub is_seed: bool,
}

/// A parsed EM-input triple, joined on contig name.
#[derive(Clone, Debug, Default)]
pub struct Dataset {
    pub contigs: Vec<Contig>,
}

/// Parse a FASTA file into (name, length, gc) records. The name is the first
/// whitespace-delimited token of the header (matching how MaxBin keys contigs).
fn parse_fasta(text: &str) -> Vec<(String, usize, f64)> {
    let mut out = Vec::new();
    let mut name: Option<String> = None;
    let mut len = 0usize;
    let mut gc = 0usize;

    let mut flush = |name: &mut Option<String>, len: &mut usize, gc: &mut usize| {
        if let Some(n) = name.take() {
            let frac = if *len > 0 {
                *gc as f64 / *len as f64
            } else {
                0.0
            };
            out.push((n, *len, frac));
        }
        *len = 0;
        *gc = 0;
    };

    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            flush(&mut name, &mut len, &mut gc);
            name = Some(header.split_whitespace().next().unwrap_or("").to_string());
        } else {
            for b in line.bytes() {
                match b {
                    b'A' | b'a' | b'T' | b't' => len += 1,
                    b'G' | b'g' | b'C' | b'c' => {
                        len += 1;
                        gc += 1;
                    }
                    // Skip whitespace; count any other base (N, etc.) toward length only.
                    b' ' | b'\t' | b'\r' => {}
                    _ => len += 1,
                }
            }
        }
    }
    flush(&mut name, &mut len, &mut gc);
    out
}

/// Parse an abundance file into (name -> depth). Splits each non-empty line at
/// the first tab/space/comma/semicolon; the rest is parsed as a float.
fn parse_abund(text: &str) -> Vec<(String, f64)> {
    let mut out = Vec::new();
    for line in text.lines() {
        let line = line.trim_end();
        let line = line.strip_prefix('>').unwrap_or(line);
        if line.is_empty() {
            continue;
        }
        let Some(pos) = line.find(['\t', ' ', ',', ';']) else {
            continue;
        };
        let name = line[..pos].trim_end().to_string();
        let value = line[pos..].trim_start_matches(['\t', ' ', ',', ';']);
        if let Ok(v) = value.trim().parse::<f64>() {
            out.push((name, v));
        }
    }
    out
}

/// Parse a seed file: one contig name per non-empty line.
fn parse_seed(text: &str) -> HashSet<String> {
    text.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(String::from)
        .collect()
}

impl Dataset {
    /// Join the three EM inputs on contig name. Contigs missing from the
    /// abundance file get abundance 0.0; seed membership is set from the seed file.
    pub fn from_inputs(contigs_fa: &str, abund: &str, seed: &str) -> Dataset {
        let depths: std::collections::HashMap<String, f64> =
            parse_abund(abund).into_iter().collect();
        let seeds = parse_seed(seed);

        let contigs = parse_fasta(contigs_fa)
            .into_iter()
            .map(|(name, length, gc)| {
                let abundance = depths.get(&name).copied().unwrap_or(0.0);
                let is_seed = seeds.contains(&name);
                Contig {
                    name,
                    length,
                    gc,
                    abundance,
                    is_seed,
                }
            })
            .collect();

        Dataset { contigs }
    }
}
