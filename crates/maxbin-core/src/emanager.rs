//! The EM algorithm core: E-step, M-step, convergence, classification.
//!
//! This is the pure heart of maxbin-rs, ported from MaxBin2's EManager class.
//! It depends only on the data handed to it - no file IO, no external tools -
//! so it compiles to wasm. File parsing and result writing live in the
//! maxbin-rs shell; this module takes already-parsed contigs and abundances via
//! [`build_state`] and runs the algorithm over them.
//!
//! Threading: on native targets the E-step parallelizes over contigs with
//! rayon (sized by `thread_num`). On wasm32 there is no thread pool, so the
//! E-step always runs serially.

use std::collections::HashSet;

#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;

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
    /// Ignored on wasm32, where the E-step always runs serially.
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

/// Log of the gamma function.
///
/// Previously linked the system `lgamma` via `extern "C"`, which does not exist
/// on wasm. `libm::lgamma` is a pure-Rust port of musl's libm, so it carries no
/// libc dependency and compiles for wasm32. It is not bit-for-bit identical to
/// the system `lgamma` (up to ~4 ULP apart on ~1% of inputs); see the commit
/// that introduced it for the measured output impact.
fn lgamma(x: f64) -> f64 {
    libm::lgamma(x)
}

/// Build EM state from already-parsed contigs and abundances.
///
/// `records` are the contigs in FASTA order, `seq_abundance[k][i]` is the
/// abundance of contig `i` in abundance file `k`, and `seed_names` are the
/// headers of the seed contigs. This is the pure entry point: the native
/// file-based wrapper (maxbin-rs `init_em`) parses files and calls this, and
/// the wasm frontend builds the same inputs in memory and calls it directly.
///
/// Matches the post-parse portion of EManager.cpp:312-422 (init_EM()).
pub fn build_state(
    records: Vec<fasta::Record>,
    seq_abundance: Vec<Vec<f64>>,
    seed_names: &[String],
    params: &EmParams,
) -> EmState {
    let kmap = KmerMap::new(params.kmer_len, true);
    let seqnum = records.len();
    let ab_num = seq_abundance.len();

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

    // Build set of seed names for O(1) lookup
    let seed_set: HashSet<&str> = seed_names.iter().map(|s| s.as_str()).collect();

    // Matches EManager.cpp:380-412 (init_EM inner loop for marked seeds):
    // identify seed contigs, copy their profiles and abundances.
    // CRITICAL: iterate in FASTA order, not seed file order. The C++ does
    // `for (i = 0; i < seqnum; i++) { if (isMark(i)) { seed[j++] = ... } }`
    // so seed index j corresponds to the j-th marked contig in FASTA order.
    let mut seed_indices = Vec::with_capacity(seed_names.len());
    let mut seed_profiles = Vec::with_capacity(seed_names.len());
    let mut seed_abundance = Vec::with_capacity(seed_names.len());

    for (idx, record) in records.iter().enumerate() {
        if !seed_set.contains(record.header.as_str()) {
            continue;
        }
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
    // faithfully reproduced - non-seed contigs are simply not in seed_indices.

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

/// E-step for a single contig: compute dist_prob, abund_prob, combine into
/// seq_prob, and pick the best bin. Returns None if the contig is skipped
/// (too short or all-N), otherwise returns (seq_prob_row, best_bin).
fn e_step_for_contig(
    i: usize,
    state: &EmState,
    seed_num: usize,
    intra: &NormalDistribution,
    inter: &NormalDistribution,
    params: &EmParams,
) -> Option<(Vec<f64>, i32)> {
    if state.records[i].len() < params.min_seq_length || state.is_profile_n[i] {
        return None;
    }

    // Compute tetranucleotide distance probabilities
    let mut dist_prob = vec![0.0f64; seed_num];
    let mut sum = 0.0f64;
    for (dp, seed_prof) in dist_prob.iter_mut().zip(state.seed_profiles.iter()) {
        let d = distance::euc_dist_profiles(
            state.seq_profiles[i].get_profile(),
            seed_prof.get_profile(),
        );
        *dp = get_prob_dist(d, intra, inter);
        if *dp < params.very_small_double {
            *dp = params.very_small_double;
        }
        sum += *dp;
    }
    for dp in dist_prob.iter_mut() {
        *dp /= sum;
    }

    // Compute abundance probabilities
    let abund_prob = compute_abund_prob_for_contig(
        i,
        &state.seq_abundance,
        &state.seed_abundance,
        params.very_small_double,
    );

    // Combine dist_prob and abund_prob into seq_prob
    let mut probs = vec![0.0f64; seed_num];
    let mut sum = 0.0f64;
    for (j, (&dp, ap_row)) in dist_prob.iter().zip(abund_prob.iter()).enumerate() {
        probs[j] = dp;
        for &ap in ap_row {
            probs[j] *= ap;
        }
        sum += probs[j];
    }

    // Normalize and find best bin
    let mut max = 0.0f64;
    let mut best_bin: i32 = -1;
    for (j, p) in probs.iter_mut().enumerate() {
        *p /= sum;
        if max < *p {
            max = *p;
            best_bin = j as i32;
        }
    }

    Some((probs, best_bin))
}

/// Compute abundance probabilities for contig `i` across all abundance files.
/// Matches EManager.cpp:1061-1081 (threadfunc_E()): each abundance file `k`
/// is independent. The C++ parallelizes this over `k` via ThreadPool.
///
/// Returns `abund_prob[j][k]` - normalized probability for seed `j`, abundance file `k`.
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
        for row in abund_prob.iter_mut() {
            row[k] /= asum;
        }
    }
    abund_prob
}

/// Parallel version: compute abundance probabilities across abundance files
/// using Rayon, matching the C++ ThreadPool + threadfunc_E parallelism.
///
/// Each abundance file `k` is processed independently, then results are
/// assembled into the same `abund_prob[j][k]` layout.
#[cfg(not(target_arch = "wasm32"))]
fn compute_abund_prob_for_contig_par(
    i: usize,
    seq_abundance: &[Vec<f64>],
    seed_abundance: &[Vec<f64>],
    very_small_double: f64,
) -> Vec<Vec<f64>> {
    let ab_num = seq_abundance.len();
    let seed_num = seed_abundance.len();

    // Each element is a Vec<f64> of length seed_num - the column for abundance file k.
    let columns: Vec<Vec<f64>> = (0..ab_num)
        .into_par_iter()
        .map(|k| {
            let abund = seq_abundance[k][i];
            let mut col: Vec<f64> = seed_abundance
                .iter()
                .map(|sa| {
                    let p = get_prob_abund(abund, sa[k]);
                    if p < very_small_double {
                        very_small_double
                    } else {
                        p
                    }
                })
                .collect();
            let asum: f64 = col.iter().sum();
            for c in col.iter_mut() {
                *c /= asum;
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

/// Build the rayon thread pool for parallel abundance computation, or `None`
/// when a single thread / single abundance file makes parallelism pointless.
/// Only present on native targets; wasm32 always runs serially.
#[cfg(not(target_arch = "wasm32"))]
fn maybe_pool(params: &EmParams, ab_num: usize) -> Option<rayon::ThreadPool> {
    if params.thread_num > 1 && ab_num > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(params.thread_num)
                .build()
                .expect("failed to build Rayon thread pool"),
        )
    } else {
        None
    }
}

/// Run the EM algorithm for up to max_iter iterations.
/// Matches EManager.cpp:459-631 (run_EM())
pub fn run_em(state: &mut EmState, params: &EmParams) {
    // Matches EManager.cpp:451-452: instantiate the two normal distributions
    let intra = make_intra_normal();
    let inter = make_inter_normal();
    let seqnum = state.records.len();
    let seed_num = state.seed_profiles.len();

    // Build a scoped Rayon thread pool matching the C++ -thread flag (native only).
    #[cfg(not(target_arch = "wasm32"))]
    let pool = maybe_pool(params, state.ab_num);

    // Matches EManager.cpp:492: stable_count = 0
    let mut stable_count = 0;

    // Matches EManager.cpp:493: for (run = 0; run < run_time; run++)
    for run in 0..params.max_em {
        eprintln!(
            "  EM iteration {}/{} ({} contigs, {} bins)",
            run + 1,
            params.max_em,
            seqnum,
            seed_num
        );
        // E-step: Matches EManager.cpp:511-592
        // Each contig's computation is independent: dist_prob, abund_prob,
        // seq_prob[i], and seq_bin[i] don't depend on other contigs.
        // When thread > 1, parallelize with Rayon. The M-step stays
        // sequential so float accumulation order is deterministic.
        let e_step_result: Vec<Option<(Vec<f64>, i32)>> = {
            #[cfg(not(target_arch = "wasm32"))]
            {
                if let Some(ref pool) = pool {
                    pool.install(|| {
                        (0..seqnum)
                            .into_par_iter()
                            .map(|i| e_step_for_contig(i, state, seed_num, &intra, &inter, params))
                            .collect()
                    })
                } else {
                    (0..seqnum)
                        .map(|i| e_step_for_contig(i, state, seed_num, &intra, &inter, params))
                        .collect()
                }
            }
            #[cfg(target_arch = "wasm32")]
            {
                (0..seqnum)
                    .map(|i| e_step_for_contig(i, state, seed_num, &intra, &inter, params))
                    .collect()
            }
        };

        // Apply E-step results to state (sequential - deterministic)
        let mut diff_count = 0;
        for (i, result) in e_step_result.into_iter().enumerate() {
            if let Some((probs, new_bin)) = result {
                state.seq_prob[i] = probs;
                if state.seq_bin[i] != new_bin {
                    diff_count += 1;
                    state.seq_bin[i] = new_bin;
                }
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
            for v in state.seed_abundance[si].iter_mut() {
                *v = 0.0;
            }
            // Matches EManager.cpp:1094: seed_profile[i]->reset()
            state.seed_profiles[si].reset();

            let mut d = 0.0f64;
            #[allow(clippy::needless_range_loop)] // j indexes 5+ parallel arrays
            for j in 0..seqnum {
                let len = state.records[j].len();
                // Matches EManager.cpp:1105: NaN check via seq_prob[j][i] == seq_prob[j][i]
                // (IEEE 754: NaN != NaN, so this is false iff the value is NaN)
                if len >= params.min_seq_length
                    && !state.seq_prob[j][si].is_nan()
                    && !state.is_profile_n[j]
                {
                    // Matches EManager.cpp:871: seed_abundance[i][k] += seq_abundance[k][j] * len * prob
                    // IMPORTANT: must match C++ evaluation order (abund * len * prob),
                    // NOT (abund * (len * prob)). Float multiplication is non-associative;
                    // different parenthesization causes 1 ULP drift that compounds over
                    // EM iterations and eventually flips threshold decisions.
                    let len_f = len as f64;
                    let prob = state.seq_prob[j][si];
                    for (sa, seq_ab) in state.seed_abundance[si]
                        .iter_mut()
                        .zip(state.seq_abundance.iter())
                    {
                        *sa += seq_ab[j] * len_f * prob;
                    }
                    let weight = len_f * prob;
                    // Matches EManager.cpp:1112: d += len * seq_prob
                    d += weight;
                    // Matches EManager.cpp:1114: seed_profile[i]->addProfile(seq_profile[j], len * prob)
                    state.seed_profiles[si].add_profile(&state.seq_profiles[j], weight);
                }
            }

            // Matches EManager.cpp:1120: seed_profile[i]->calcProfile()
            state.seed_profiles[si].calc_profile(d);

            // Matches EManager.cpp:1124-1126: seed_abundance[i][k] /= d
            for v in state.seed_abundance[si].iter_mut() {
                *v /= d;
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

    #[cfg(not(target_arch = "wasm32"))]
    let pool = maybe_pool(params, state.ab_num);

    let mut dist_prob = vec![0.0f64; seed_num];

    // Matches EManager.cpp:364: seed_count reset (memset to zero)
    state.seed_count.fill(0);

    // Matches EManager.cpp:665-760: classify each sequence
    for i in 0..seqnum {
        // Matches EManager.cpp:667: if (seq len >= min_seqlen && !is_profile_N)
        if state.records[i].len() >= params.min_seq_length && !state.is_profile_n[i] {
            // Matches EManager.cpp:669: if (is_estimated[i] == false) - re-run E step
            if !state.is_estimated[i] {
                // Matches EManager.cpp:671-680: recompute dist_prob
                let mut sum = 0.0f64;
                for (dp, seed_prof) in dist_prob.iter_mut().zip(state.seed_profiles.iter()) {
                    let d = distance::euc_dist_profiles(
                        state.seq_profiles[i].get_profile(),
                        seed_prof.get_profile(),
                    );
                    *dp = get_prob_dist(d, &intra, &inter);
                    sum += *dp;
                }
                for dp in dist_prob.iter_mut() {
                    *dp /= sum;
                }

                // Matches EManager.cpp:697-701: recompute abund_prob
                let abund_prob = {
                    #[cfg(not(target_arch = "wasm32"))]
                    {
                        if let Some(ref pool) = pool {
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
                        }
                    }
                    #[cfg(target_arch = "wasm32")]
                    {
                        compute_abund_prob_for_contig(
                            i,
                            &state.seq_abundance,
                            &state.seed_abundance,
                            params.very_small_double,
                        )
                    }
                };

                // Matches EManager.cpp:703-716: combine and normalize seq_prob
                let mut sum = 0.0f64;
                for (j, (&dp, ap_row)) in dist_prob.iter().zip(abund_prob.iter()).enumerate() {
                    state.seq_prob[i][j] = dp;
                    for &ap in ap_row {
                        state.seq_prob[i][j] *= ap;
                    }
                    sum += state.seq_prob[i][j];
                }
                for sp in state.seq_prob[i].iter_mut() {
                    *sp /= sum;
                    // Matches EManager.cpp:725-728: NaN handling (sum too small)
                    if sp.is_nan() {
                        *sp = 0.0;
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

    // NOTE: get_prob_abund_matches_cpp was removed - it computed the "expected"
    // value using the same lgamma and formula as the function under test, making
    // it tautological. Real coverage lives in the maxbin-rs-equivalence crate:
    // its proptest_emanager.rs (algebraic properties on the real function) and
    // emanager_equivalence.rs (full pipeline byte-for-byte comparison against C++).
}
