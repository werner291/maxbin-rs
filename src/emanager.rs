/// Rust reimplementation of MaxBin2's EManager — the EM algorithm core.
///
/// Decomposed from the original monolithic EManager class into:
/// - EmParams: configuration constants
/// - EmState: mutable algorithm state
/// - Free functions for each phase
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use rayon::prelude::*;

use crate::abundance;
use crate::distance;
use crate::fasta;
use crate::kmer_map::KmerMap;
use crate::normal_distribution::NormalDistribution;
use crate::profiler::Profiler;

/// Configuration constants matching the original EManager::init().
/// Matches EManager.cpp:111-165 (init()): hardcoded defaults.
pub struct EmParams {
    /// Matches EManager.cpp:144: kmer_len = 4
    pub kmer_len: usize,
    /// Matches EManager.cpp:152: min_seq_length = 1000
    pub min_seq_length: usize,
    /// Matches EManager.cpp:153: max_EM = 50
    pub max_em: usize,
    /// Matches EManager.cpp:155: MIN_PROB_THRESHOLD = 0.5
    pub min_prob_threshold: f64,
    /// Matches EManager.cpp:154: VERY_SMALL_DOUBLE = 1e-20
    pub very_small_double: f64,
    /// Matches EManager.cpp:157: FASTA_LINE = 70
    pub fasta_line: usize,
    /// Matches EManager.cpp:158: STABLE_BIN_COUNT = 5
    pub stable_bin_count: usize,
    /// Number of threads for parallel E-step abundance computation.
    /// Matches EManager.cpp's ThreadPool size (set via -thread flag).
    pub thread_num: usize,
}

impl Default for EmParams {
    fn default() -> Self {
        Self {
            kmer_len: 4,
            min_seq_length: 1000,
            max_em: 50,
            min_prob_threshold: 0.5,
            very_small_double: 1e-20,
            fasta_line: 70,
            stable_bin_count: 5,
            thread_num: 1,
        }
    }
}

/// Mutable state for the EM algorithm.
pub struct EmState {
    /// FASTA records (all contigs)
    pub records: Vec<fasta::Record>,
    /// seq_abundance[k][i] = abundance of contig i from abundance file k
    pub seq_abundance: Vec<Vec<f64>>,
    /// Number of abundance files
    pub ab_num: usize,

    /// Tetranucleotide profiles for all contigs
    pub seq_profiles: Vec<Profiler>,
    /// Whether each contig's profile is all-N
    pub is_profile_n: Vec<bool>,

    /// Seed indices (into records)
    pub seed_indices: Vec<usize>,
    /// Seed profiles (separate copies, updated during M-step)
    pub seed_profiles: Vec<Profiler>,
    /// seed_abundance[j][k] = abundance of seed j from abundance file k
    pub seed_abundance: Vec<Vec<f64>>,

    /// seq_prob[i][j] = probability of contig i belonging to bin j
    pub seq_prob: Vec<Vec<f64>>,
    /// seq_bin[i] = assigned bin for contig i (-1 = unclassified)
    pub seq_bin: Vec<i32>,
    /// seed_count[j] = number of contigs assigned to bin j
    pub seed_count: Vec<i32>,
    /// Whether contig i was estimated during run_em
    pub is_estimated: Vec<bool>,
}

/// Hardcoded normal distribution parameters from estimate_normaldistr().
/// These come from simulation tests on 3181 IMG bacterial/archaeal genomes.
/// Matches EManager.cpp:424-456 (estimate_normaldistr()):
/// comment says "Pre-defined mean and standard deviation found from simulation test
/// on 3181 IMG bacterial and archaeal genomes"
pub fn make_intra_normal() -> NormalDistribution {
    // Matches EManager.cpp:435-436: mean=0, std=0.01037897/2
    NormalDistribution::new(0.0, 0.01037897 / 2.0)
}

pub fn make_inter_normal() -> NormalDistribution {
    // Matches EManager.cpp:437-438: mean2=0.0676654, std2=0.03419337
    NormalDistribution::new(0.0676654, 0.03419337)
}

/// get_prob_dist: P(intra) / (P(intra) + P(inter))
/// Matches EManager.cpp:1017-1023 (get_prob_dist()):
/// d_intra / (d_inter + d_intra)
pub fn get_prob_dist(distance: f64, intra: &NormalDistribution, inter: &NormalDistribution) -> f64 {
    let d_intra = intra.prob(distance);
    let d_inter = inter.prob(distance);
    d_intra / (d_inter + d_intra)
}

/// get_prob_abund: Poisson-like probability.
/// Uses log/lgamma/exp, adapted from the original C++ implementation.
/// Note: the original uses `long double` but we use f64 and call the same
/// math functions (log, lgamma, exp).
/// Matches EManager.cpp:1033-1048 (get_prob_abund()):
/// exp(curr_abund * log(lambda) - lgamma(curr_abund+1) - lambda)
pub fn get_prob_abund(curr_abund: f64, lambda: f64) -> f64 {
    // Matches EManager.cpp:1037-1038: return 0 if either is zero
    if lambda == 0.0 || curr_abund == 0.0 {
        0.0
    } else {
        // Matches EManager.cpp:1043-1046:
        // l = log(lambda); ret = exp((curr_abund * l) - lgamma(curr_abund+1) - lambda)
        let l = lambda.ln();
        ((curr_abund * l) - lgamma(curr_abund + 1.0) - lambda).exp()
    }
}

/// C's lgamma function (log of gamma function).
fn lgamma(x: f64) -> f64 {
    // Use libc's lgamma via the extern
    unsafe { libc_lgamma(x) }
}

unsafe extern "C" {
    #[link_name = "lgamma"]
    fn libc_lgamma(x: f64) -> f64;
}

/// Initialize EM state from input files and seed contig names.
/// Matches EManager.cpp:312-422 (init_EM()):
/// parse FASTA, load abundance, build profiles, identify seeds.
pub fn init_em(
    fasta_path: &Path,
    abundance_paths: &[&Path],
    seed_names: &[String],
    params: &EmParams,
) -> EmState {
    let kmap = KmerMap::new(params.kmer_len, true);

    // Matches EManager.cpp:18-19 (constructor): parse FASTA, get seqnum
    let records = fasta::parse_file(fasta_path).expect("failed to parse FASTA");
    let seqnum = records.len();

    // Matches EManager.cpp:214-251 (addAbund()): load abundance per contig
    let ab_num = abundance_paths.len();
    let mut seq_abundance = Vec::with_capacity(ab_num);
    for abund_path in abundance_paths {
        let abund_records =
            abundance::parse_file(abund_path).expect("failed to parse abundance file");
        let abund_map: HashMap<&str, f64> = abund_records
            .iter()
            .map(|r| (r.header.as_str(), r.abundance))
            .collect();

        let mut abund_vec = vec![0.0f64; seqnum];
        for (i, rec) in records.iter().enumerate() {
            abund_vec[i] = *abund_map
                .get(rec.header.as_str())
                .expect(&format!("abundance not found for contig: {}", rec.header));
        }
        seq_abundance.push(abund_vec);
    }

    // Matches EManager.cpp:365-378: build Profiler for each sequence, flag all-N profiles
    let mut seq_profiles = Vec::with_capacity(seqnum);
    let mut is_profile_n = vec![false; seqnum];
    for (i, rec) in records.iter().enumerate() {
        let prof = Profiler::new(params.kmer_len, &rec.seq, &kmap);
        // Matches EManager.cpp:369-372: is_profile_N[i] = true if percent_N == 1
        if prof.get_percent_n() == 1.0 {
            is_profile_n[i] = true;
        }
        seq_profiles.push(prof);
    }

    // Build seed-name -> record-index mapping
    let header_to_idx: HashMap<&str, usize> = records
        .iter()
        .enumerate()
        .map(|(i, r)| (r.header.as_str(), i))
        .collect();

    // Matches EManager.cpp:380-412 (init_EM inner loop for marked seeds):
    // identify seed contigs, copy their profiles and abundances
    let mut seed_indices = Vec::with_capacity(seed_names.len());
    let mut seed_profiles = Vec::with_capacity(seed_names.len());
    let mut seed_abundance = Vec::with_capacity(seed_names.len());

    for name in seed_names {
        let &idx = header_to_idx
            .get(name.as_str())
            .expect(&format!("seed not found in FASTA: {name}"));
        seed_indices.push(idx);

        // Matches EManager.cpp:382-384: seed_profile[j] = new Profiler(...)
        let prof = Profiler::new(params.kmer_len, &records[idx].seq, &kmap);
        seed_profiles.push(prof);

        // Matches EManager.cpp:386-389: copy abundance values for this seed
        let mut abund = vec![0.0f64; ab_num];
        for k in 0..ab_num {
            abund[k] = seq_abundance[k][idx];
        }
        seed_abundance.push(abund);
    }

    let seed_num = seed_indices.len();

    // Matches EManager.cpp:362-363 / 416-421: initialize probability/assignment arrays to zero
    let seq_prob = vec![vec![0.0f64; seed_num]; seqnum];
    let seq_bin = vec![0i32; seqnum];
    let seed_count = vec![0i32; seed_num];
    let is_estimated = vec![false; seqnum];

    // Reproduce the no-op bug at EManager.cpp:404: `seq->isMark(i) == false;`
    // In the original, this is a comparison whose result is discarded (statement with no effect).
    // The only effect is that for non-seed contigs, the mark is NOT cleared.
    // Since we don't use marks (we use seed_indices instead), this is already
    // faithfully reproduced — non-seed contigs are simply not in seed_indices.

    EmState {
        records,
        seq_abundance,
        ab_num,
        seq_profiles,
        is_profile_n,
        seed_indices,
        seed_profiles,
        seed_abundance,
        seq_prob,
        seq_bin,
        seed_count,
        is_estimated,
    }
}

/// Compute abundance probabilities for contig `i` across all abundance files.
/// Matches EManager.cpp:1061-1081 (threadfunc_E()): each abundance file `k`
/// is independent. The C++ parallelizes this over `k` via ThreadPool.
///
/// Returns `abund_prob[j][k]` — normalized probability for seed `j`, abundance file `k`.
fn compute_abund_prob_for_contig(
    i: usize,
    seq_abundance: &[Vec<f64>],
    seed_abundance: &[Vec<f64>],
    very_small_double: f64,
) -> Vec<Vec<f64>> {
    let ab_num = seq_abundance.len();
    let seed_num = seed_abundance.len();
    let mut abund_prob = vec![vec![0.0f64; ab_num]; seed_num];

    for k in 0..ab_num {
        let abund = seq_abundance[k][i];
        let mut asum = 0.0f64;
        for j in 0..seed_num {
            abund_prob[j][k] = get_prob_abund(abund, seed_abundance[j][k]);
            // Matches EManager.cpp:1071-1074: clamp to VERY_SMALL_DOUBLE
            if abund_prob[j][k] < very_small_double {
                abund_prob[j][k] = very_small_double;
            }
            asum += abund_prob[j][k];
        }
        // Matches EManager.cpp:1077-1080: normalize abund_prob
        for j in 0..seed_num {
            abund_prob[j][k] /= asum;
        }
    }
    abund_prob
}

/// Parallel version: compute abundance probabilities across abundance files
/// using Rayon, matching the C++ ThreadPool + threadfunc_E parallelism.
///
/// Each abundance file `k` is processed independently, then results are
/// assembled into the same `abund_prob[j][k]` layout.
fn compute_abund_prob_for_contig_par(
    i: usize,
    seq_abundance: &[Vec<f64>],
    seed_abundance: &[Vec<f64>],
    very_small_double: f64,
) -> Vec<Vec<f64>> {
    let ab_num = seq_abundance.len();
    let seed_num = seed_abundance.len();

    // Each element is a Vec<f64> of length seed_num — the column for abundance file k.
    let columns: Vec<Vec<f64>> = (0..ab_num)
        .into_par_iter()
        .map(|k| {
            let abund = seq_abundance[k][i];
            let mut col = vec![0.0f64; seed_num];
            let mut asum = 0.0f64;
            for j in 0..seed_num {
                col[j] = get_prob_abund(abund, seed_abundance[j][k]);
                if col[j] < very_small_double {
                    col[j] = very_small_double;
                }
                asum += col[j];
            }
            for j in 0..seed_num {
                col[j] /= asum;
            }
            col
        })
        .collect();

    // Transpose columns[k][j] into abund_prob[j][k]
    let mut abund_prob = vec![vec![0.0f64; ab_num]; seed_num];
    for k in 0..ab_num {
        for j in 0..seed_num {
            abund_prob[j][k] = columns[k][j];
        }
    }
    abund_prob
}

/// Run the EM algorithm for up to max_iter iterations.
/// Matches EManager.cpp:459-631 (run_EM())
pub fn run_em(state: &mut EmState, params: &EmParams) {
    // Matches EManager.cpp:451-452: instantiate the two normal distributions
    let intra = make_intra_normal();
    let inter = make_inter_normal();
    let seqnum = state.records.len();
    let seed_num = state.seed_profiles.len();
    let ab_num = state.ab_num;
    let use_parallel = params.thread_num > 1 && ab_num > 1;

    // Build a scoped Rayon thread pool matching the C++ -thread flag.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(params.thread_num)
        .build()
        .expect("failed to build Rayon thread pool");

    // Matches EManager.cpp:492: stable_count = 0
    let mut stable_count = 0;

    let mut dist_prob = vec![0.0f64; seed_num];

    // Matches EManager.cpp:493: for (run = 0; run < run_time; run++)
    for run in 0..params.max_em {
        eprintln!(
            "  EM iteration {}/{} ({} contigs, {} bins)",
            run + 1,
            params.max_em,
            seqnum,
            seed_num
        );
        let mut diff_count = 0;

        // E-step: Matches EManager.cpp:511-592
        for i in 0..seqnum {
            // Matches EManager.cpp:514: if (seq len >= min_seq_length && !is_profile_N)
            if state.records[i].len() >= params.min_seq_length && !state.is_profile_n[i] {
                // Matches EManager.cpp:516-531: compute tetranucleotide dist probabilities
                let mut sum = 0.0f64;
                for j in 0..seed_num {
                    let d = distance::euc_dist_profiles(
                        state.seq_profiles[i].get_profile(),
                        state.seed_profiles[j].get_profile(),
                    );
                    dist_prob[j] = get_prob_dist(d, &intra, &inter);
                    // Matches EManager.cpp:522-525: clamp to VERY_SMALL_DOUBLE
                    if dist_prob[j] < params.very_small_double {
                        dist_prob[j] = params.very_small_double;
                    }
                    sum += dist_prob[j];
                }
                // Matches EManager.cpp:528-531: normalize dist_prob
                for j in 0..seed_num {
                    dist_prob[j] /= sum;
                }

                // Compute abundance probabilities (threadfunc_E equivalent)
                // Matches EManager.cpp:1061-1081 (threadfunc_E()): for each abundance file.
                // The C++ parallelizes this over abundance files via ThreadPool;
                // we do the same with Rayon when thread_num > 1.
                let abund_prob = if use_parallel {
                    pool.install(|| {
                        compute_abund_prob_for_contig_par(
                            i,
                            &state.seq_abundance,
                            &state.seed_abundance,
                            params.very_small_double,
                        )
                    })
                } else {
                    compute_abund_prob_for_contig(
                        i,
                        &state.seq_abundance,
                        &state.seed_abundance,
                        params.very_small_double,
                    )
                };

                // Matches EManager.cpp:555-568: combine dist_prob and abund_prob
                let mut sum = 0.0f64;
                for j in 0..seed_num {
                    // Matches EManager.cpp:562: seq_prob[i][j] = dist_prob[j]
                    state.seq_prob[i][j] = dist_prob[j];
                    // Matches EManager.cpp:564-566: multiply by each abundance probability
                    for k in 0..ab_num {
                        state.seq_prob[i][j] *= abund_prob[j][k];
                    }
                    sum += state.seq_prob[i][j];
                }

                // Matches EManager.cpp:570-589: normalize seq_prob, find best bin
                let mut max = 0.0f64;
                let mut tempbin: i32 = -1;
                for j in 0..seed_num {
                    state.seq_prob[i][j] /= sum;
                    if max < state.seq_prob[i][j] {
                        max = state.seq_prob[i][j];
                        tempbin = j as i32;
                    }
                }

                // Matches EManager.cpp:585-589: update seq_bin, count changes
                if state.seq_bin[i] != tempbin {
                    diff_count += 1;
                    state.seq_bin[i] = tempbin;
                }
                // Matches EManager.cpp:591: is_estimated[i] = true
                state.is_estimated[i] = true;
            }
        }

        // Matches EManager.cpp:594-602: stable_count tracking
        if diff_count == 0 {
            stable_count += 1;
        } else {
            stable_count = 0;
        }

        // M-step: Matches EManager.cpp:1084-1128 (threadfunc_M(int i))
        for si in 0..seed_num {
            // Matches EManager.cpp:1090-1093: reset seed abundance to zero
            for k in 0..ab_num {
                state.seed_abundance[si][k] = 0.0;
            }
            // Matches EManager.cpp:1094: seed_profile[i]->reset()
            state.seed_profiles[si].reset();

            let mut d = 0.0f64;
            for j in 0..seqnum {
                let len = state.records[j].len();
                // Matches EManager.cpp:1105: NaN check via seq_prob[j][i] == seq_prob[j][i]
                // (IEEE 754: NaN != NaN, so this is false iff the value is NaN)
                if len >= params.min_seq_length
                    && state.seq_prob[j][si] == state.seq_prob[j][si]
                    && !state.is_profile_n[j]
                {
                    // Matches EManager.cpp:1108-1111: weight = len * seq_prob; accumulate abundance
                    let weight = len as f64 * state.seq_prob[j][si];
                    for k in 0..ab_num {
                        state.seed_abundance[si][k] += state.seq_abundance[k][j] * weight;
                    }
                    // Matches EManager.cpp:1112: d += len * seq_prob
                    d += weight;
                    // Matches EManager.cpp:1114: seed_profile[i]->addProfile(seq_profile[j], len * prob)
                    state.seed_profiles[si].add_profile(&state.seq_profiles[j], weight);
                }
            }

            // Matches EManager.cpp:1120: seed_profile[i]->calcProfile()
            state.seed_profiles[si].calc_profile(d);

            // Matches EManager.cpp:1124-1126: seed_abundance[i][k] /= d
            for k in 0..ab_num {
                state.seed_abundance[si][k] /= d;
            }
        }

        // Matches EManager.cpp:611-614: break if stable for STABLE_BIN_COUNT iterations
        if stable_count >= params.stable_bin_count {
            eprintln!(
                "\r  EM converged after {} iterations ({} changed last iteration)    ",
                run + 1,
                diff_count
            );
            break;
        }
    }
    eprintln!();
}

/// Classify contigs into bins based on final probabilities.
/// Matches EManager.cpp:633-784 (classify()):
/// re-compute probabilities for unestimated contigs, then assign to best bin or -1.
pub fn classify(state: &mut EmState, params: &EmParams) {
    let intra = make_intra_normal();
    let inter = make_inter_normal();
    let seqnum = state.records.len();
    let seed_num = state.seed_profiles.len();
    let ab_num = state.ab_num;
    let use_parallel = params.thread_num > 1 && ab_num > 1;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(params.thread_num)
        .build()
        .expect("failed to build Rayon thread pool");

    let mut dist_prob = vec![0.0f64; seed_num];

    // Matches EManager.cpp:364: seed_count reset (memset to zero)
    state.seed_count.fill(0);

    // Matches EManager.cpp:665-760: classify each sequence
    for i in 0..seqnum {
        // Matches EManager.cpp:667: if (seq len >= min_seqlen && !is_profile_N)
        if state.records[i].len() >= params.min_seq_length && !state.is_profile_n[i] {
            // Matches EManager.cpp:669: if (is_estimated[i] == false) — re-run E step
            if !state.is_estimated[i] {
                // Matches EManager.cpp:671-680: recompute dist_prob
                let mut sum = 0.0f64;
                for j in 0..seed_num {
                    let d = distance::euc_dist_profiles(
                        state.seq_profiles[i].get_profile(),
                        state.seed_profiles[j].get_profile(),
                    );
                    dist_prob[j] = get_prob_dist(d, &intra, &inter);
                    sum += dist_prob[j];
                }
                for j in 0..seed_num {
                    dist_prob[j] /= sum;
                }

                // Matches EManager.cpp:697-701: recompute abund_prob
                let abund_prob = if use_parallel {
                    pool.install(|| {
                        compute_abund_prob_for_contig_par(
                            i,
                            &state.seq_abundance,
                            &state.seed_abundance,
                            params.very_small_double,
                        )
                    })
                } else {
                    compute_abund_prob_for_contig(
                        i,
                        &state.seq_abundance,
                        &state.seed_abundance,
                        params.very_small_double,
                    )
                };

                // Matches EManager.cpp:703-716: combine and normalize seq_prob
                let mut sum = 0.0f64;
                for j in 0..seed_num {
                    state.seq_prob[i][j] = dist_prob[j];
                    for k in 0..ab_num {
                        state.seq_prob[i][j] *= abund_prob[j][k];
                    }
                    sum += state.seq_prob[i][j];
                }
                for j in 0..seed_num {
                    state.seq_prob[i][j] /= sum;
                    // Matches EManager.cpp:725-728: NaN handling (sum too small)
                    if state.seq_prob[i][j] != state.seq_prob[i][j] {
                        state.seq_prob[i][j] = 0.0;
                    }
                }
            }

            // Matches EManager.cpp:733-749: find best bin (max probability)
            let mut max = 0.0f64;
            for j in 0..seed_num {
                if state.seq_prob[i][j] > max {
                    max = state.seq_prob[i][j];
                    state.seq_bin[i] = j as i32;
                }
            }
            // Matches EManager.cpp:742-749: apply min_prob threshold; increment seed_count
            if max <= params.min_prob_threshold {
                state.seq_bin[i] = -1;
            } else {
                state.seed_count[state.seq_bin[i] as usize] += 1;
            }
        } else {
            // Matches EManager.cpp:757-759: sequences too short or all-N get bin = -1
            state.seq_bin[i] = -1;
        }
    }
}

/// Write output files: .NNN.fasta bins, .noclass, .summary
/// Matches EManager.cpp:867-991 (write_result())
pub fn write_results(state: &EmState, output_prefix: &str, params: &EmParams) {
    let seed_num = state.seed_profiles.len();
    let seqnum = state.records.len();

    // Matches EManager.cpp:872-890: name bins with 4-digit zero-padding, skip empty bins
    // Matches EManager.cpp:883: sprintf(bin_name[i], "%s.%04d.fasta", outputfile, j)
    let mut bin_names: Vec<Option<String>> = vec![None; seed_num];
    let mut j = 1;
    for i in 0..seed_num {
        if state.seed_count[i] > 0 {
            bin_names[i] = Some(format!("{}.{:04}.fasta", output_prefix, j));
            j += 1;
        }
    }

    // Matches EManager.cpp:894-944 (write_result() summary section)
    let summary_path = format!("{}.summary", output_prefix);
    let mut summary = std::io::BufWriter::new(
        std::fs::File::create(&summary_path).expect("failed to create summary"),
    );

    for i in 0..seed_num {
        if state.seed_count[i] > 0 {
            // Matches EManager.cpp:901-907: "Bin [name]\tabundance..."
            write!(summary, "Bin [{}]", bin_names[i].as_ref().unwrap()).unwrap();
            for k in 0..state.ab_num {
                write!(summary, "\t{:.4}", state.seed_abundance[i][k]).unwrap();
            }
            writeln!(summary).unwrap();
            for j in 0..seqnum {
                if state.seq_bin[j] == i as i32 {
                    write!(summary, "\t{}", state.records[j].header).unwrap();
                    for k in 0..state.ab_num {
                        write!(summary, "\t{:.4}", state.seq_abundance[k][j]).unwrap();
                    }
                    writeln!(summary).unwrap();
                }
            }
            writeln!(summary).unwrap();
        }
    }

    // Matches EManager.cpp:926-943: "Bins without any sequences:" section
    write!(summary, "\nBins without any sequences:\n").unwrap();
    let mut empty_j = 1;
    for i in 0..seed_num {
        if state.seed_count[i] == 0 {
            write!(summary, "{}:", empty_j).unwrap();
            for k in 0..state.ab_num {
                write!(summary, " ({:.4})", state.seed_abundance[i][k]).unwrap();
            }
            writeln!(summary).unwrap();
            empty_j += 1;
        }
    }
    drop(summary);

    // Matches EManager.cpp:947-987: open per-bin fasta files and noclass file, write sequences
    let mut bin_files: Vec<Option<std::io::BufWriter<std::fs::File>>> = (0..seed_num)
        .map(|i| {
            bin_names[i].as_ref().map(|name| {
                std::io::BufWriter::new(
                    std::fs::File::create(name).expect("failed to create bin file"),
                )
            })
        })
        .collect();

    // Matches EManager.cpp:956-957: open .noclass file
    let noclass_path = format!("{}.noclass", output_prefix);
    let mut noclass = std::io::BufWriter::new(
        std::fs::File::create(&noclass_path).expect("failed to create noclass"),
    );

    for i in 0..seqnum {
        let header_line = format!(">{}\n", state.records[i].header);
        let seq_str = std::str::from_utf8(&state.records[i].seq).unwrap();

        if state.seq_bin[i] != -1 {
            let bin_idx = state.seq_bin[i] as usize;
            let f = bin_files[bin_idx].as_mut().unwrap();
            f.write_all(header_line.as_bytes()).unwrap();
            write_fasta_seq(f, seq_str, params.fasta_line);
        } else {
            noclass.write_all(header_line.as_bytes()).unwrap();
            write_fasta_seq(&mut noclass, seq_str, params.fasta_line);
        }
    }
}

/// Write a FASTA sequence with line wrapping at `line_width` characters.
/// Matches EManager.cpp:993-1015 (write_fasta()): write FASTA_LINE chars at a time.
fn write_fasta_seq(f: &mut impl Write, seq: &str, line_width: usize) {
    let bytes = seq.as_bytes();
    let len = bytes.len();
    let mut pos = 0;
    while pos < len {
        let end = std::cmp::min(pos + line_width, len);
        f.write_all(&bytes[pos..end]).unwrap();
        f.write_all(b"\n").unwrap();
        pos = end;
    }
}

/// Full pipeline: load data, run EM, classify, write results.
pub fn run_pipeline(
    fasta_path: &Path,
    abundance_paths: &[&Path],
    seed_names: &[String],
    output_prefix: &str,
    params: &EmParams,
) {
    let mut state = init_em(fasta_path, abundance_paths, seed_names, params);
    run_em(&mut state, params);
    classify(&mut state, params);
    write_results(&state, output_prefix, params);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn get_prob_dist_basic() {
        let intra = make_intra_normal();
        let inter = make_inter_normal();

        // At distance 0, intra should dominate (but inter PDF at 0 is nonzero)
        let p = get_prob_dist(0.0, &intra, &inter);
        assert!(p > 0.9, "at d=0, prob should be high, got {p}");

        // At large distance, inter should dominate
        let p = get_prob_dist(0.2, &intra, &inter);
        assert!(p < 0.01, "at d=0.2, prob should be near 0, got {p}");
    }

    #[test]
    fn get_prob_abund_basic() {
        // Zero lambda or zero abundance => 0
        assert_eq!(get_prob_abund(0.0, 5.0), 0.0);
        assert_eq!(get_prob_abund(5.0, 0.0), 0.0);

        // Positive values should return positive
        let p = get_prob_abund(5.0, 5.0);
        assert!(p > 0.0, "expected positive, got {p}");
    }

    // NOTE: get_prob_abund_matches_cpp was removed — it computed the "expected"
    // value using the same lgamma and formula as the function under test, making
    // it tautological. Real coverage comes from proptest_emanager.rs (algebraic
    // properties on the real function) and emanager_equivalence.rs (full pipeline
    // byte-for-byte comparison against C++).
}
