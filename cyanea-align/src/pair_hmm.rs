//! Pair Hidden Markov Model for probabilistic sequence alignment.
//!
//! Implements a three-state Pair HMM (Match, InsertX, InsertY) operating
//! entirely in log-space for numerical stability. Provides both the Forward
//! algorithm (total log-likelihood) and the Viterbi algorithm (optimal
//! state path with traceback).
//!
//! # Model
//!
//! The three states are:
//!
//! - **M** — both sequences emit a character (match or mismatch)
//! - **X** — only sequence A emits (gap in B)
//! - **Y** — only sequence B emits (gap in A)
//!
//! All probabilities are stored and computed as natural logarithms.

use cyanea_core::{CyaneaError, Result};

/// Parameters for a three-state Pair HMM.
///
/// The model has three states:
/// - **M** (Match): both sequences emit a character (match or mismatch)
/// - **X** (InsertX): only sequence A emits (gap in B)
/// - **Y** (InsertY): only sequence B emits (gap in A)
#[derive(Debug, Clone)]
pub struct PairHmmParams {
    /// Log-probability of matching bases in state M.
    pub log_match_emission: f64,
    /// Log-probability of mismatching bases in state M.
    pub log_mismatch_emission: f64,
    /// Log-probability of transitioning M -> M.
    pub log_trans_mm: f64,
    /// Log-probability of transitioning M -> X (gap open in B).
    pub log_trans_mx: f64,
    /// Log-probability of transitioning M -> Y (gap open in A).
    pub log_trans_my: f64,
    /// Log-probability of transitioning X -> X (gap extend in B).
    pub log_trans_xx: f64,
    /// Log-probability of transitioning X -> M (gap close from X).
    pub log_trans_xm: f64,
    /// Log-probability of transitioning Y -> Y (gap extend in A).
    pub log_trans_yy: f64,
    /// Log-probability of transitioning Y -> M (gap close from Y).
    pub log_trans_ym: f64,
}

impl Default for PairHmmParams {
    fn default() -> Self {
        // Reasonable defaults for DNA alignment
        Self {
            log_match_emission: (0.98_f64).ln(),
            log_mismatch_emission: (0.02_f64 / 3.0).ln(),
            log_trans_mm: (0.9_f64).ln(),
            log_trans_mx: (0.05_f64).ln(),
            log_trans_my: (0.05_f64).ln(),
            log_trans_xx: (0.3_f64).ln(),
            log_trans_xm: (0.7_f64).ln(),
            log_trans_yy: (0.3_f64).ln(),
            log_trans_ym: (0.7_f64).ln(),
        }
    }
}

/// State in a Pair HMM alignment path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairHmmState {
    Match,
    InsertX,
    InsertY,
}

/// Result of Pair HMM Viterbi alignment.
#[derive(Debug, Clone)]
pub struct PairHmmAlignment {
    /// Log-probability score of the Viterbi path.
    pub score: f64,
    /// Sequence of states in the optimal path.
    pub states: Vec<PairHmmState>,
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Numerically stable log(exp(a) + exp(b)).
fn log_sum_exp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let max = a.max(b);
    max + ((a - max).exp() + (b - max).exp()).ln()
}

/// Three-way log-sum-exp.
fn log_sum_exp3(a: f64, b: f64, c: f64) -> f64 {
    log_sum_exp(log_sum_exp(a, b), c)
}

/// Emission log-probability: match vs mismatch.
fn emit(a: u8, b: u8, params: &PairHmmParams) -> f64 {
    if a == b {
        params.log_match_emission
    } else {
        params.log_mismatch_emission
    }
}

// ---------------------------------------------------------------------------
// Forward algorithm
// ---------------------------------------------------------------------------

/// Forward algorithm returning the total log-likelihood P(seq_a, seq_b | model).
///
/// Computes the sum over all possible state paths (in log-space) using
/// the standard three-matrix DP recurrence.
///
/// # Errors
///
/// Returns an error if either sequence is empty.
pub fn pair_hmm_forward(
    seq_a: &[u8],
    seq_b: &[u8],
    params: &PairHmmParams,
) -> Result<f64> {
    let n = seq_a.len();
    let m = seq_b.len();

    if n == 0 || m == 0 {
        return Err(CyaneaError::InvalidInput(
            "sequences must not be empty".into(),
        ));
    }

    let rows = n + 1;
    let cols = m + 1;

    let mut fm = vec![f64::NEG_INFINITY; rows * cols];
    let mut fx = vec![f64::NEG_INFINITY; rows * cols];
    let mut fy = vec![f64::NEG_INFINITY; rows * cols];

    let idx = |i: usize, j: usize| -> usize { i * cols + j };

    // Base case: start in match state with probability 1 (log(1) = 0).
    fm[idx(0, 0)] = 0.0;

    for i in 0..=n {
        for j in 0..=m {
            // --- M state: requires both i > 0 and j > 0 ---
            if i > 0 && j > 0 {
                let e = emit(seq_a[i - 1], seq_b[j - 1], params);
                fm[idx(i, j)] = e
                    + log_sum_exp3(
                        fm[idx(i - 1, j - 1)] + params.log_trans_mm,
                        fx[idx(i - 1, j - 1)] + params.log_trans_xm,
                        fy[idx(i - 1, j - 1)] + params.log_trans_ym,
                    );
            }

            // --- X state: only sequence A emits (i > 0) ---
            if i > 0 {
                fx[idx(i, j)] = log_sum_exp(
                    fm[idx(i - 1, j)] + params.log_trans_mx,
                    fx[idx(i - 1, j)] + params.log_trans_xx,
                );
            }

            // --- Y state: only sequence B emits (j > 0) ---
            if j > 0 {
                fy[idx(i, j)] = log_sum_exp(
                    fm[idx(i, j - 1)] + params.log_trans_my,
                    fy[idx(i, j - 1)] + params.log_trans_yy,
                );
            }
        }
    }

    Ok(log_sum_exp3(
        fm[idx(n, m)],
        fx[idx(n, m)],
        fy[idx(n, m)],
    ))
}

// ---------------------------------------------------------------------------
// Viterbi algorithm
// ---------------------------------------------------------------------------

/// Backpointer tag used during Viterbi traceback.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BackPtr {
    None,
    FromM,
    FromX,
    FromY,
}

/// Viterbi algorithm returning the optimal state path and its log-probability.
///
/// Uses the same DP structure as [`pair_hmm_forward`] but replaces
/// `log_sum_exp` with `max` and records backpointers for traceback.
///
/// # Errors
///
/// Returns an error if either sequence is empty.
pub fn pair_hmm_viterbi(
    seq_a: &[u8],
    seq_b: &[u8],
    params: &PairHmmParams,
) -> Result<PairHmmAlignment> {
    let n = seq_a.len();
    let m = seq_b.len();

    if n == 0 || m == 0 {
        return Err(CyaneaError::InvalidInput(
            "sequences must not be empty".into(),
        ));
    }

    let rows = n + 1;
    let cols = m + 1;

    let mut vm = vec![f64::NEG_INFINITY; rows * cols];
    let mut vx = vec![f64::NEG_INFINITY; rows * cols];
    let mut vy = vec![f64::NEG_INFINITY; rows * cols];

    let mut bp_m = vec![BackPtr::None; rows * cols];
    let mut bp_x = vec![BackPtr::None; rows * cols];
    let mut bp_y = vec![BackPtr::None; rows * cols];

    let idx = |i: usize, j: usize| -> usize { i * cols + j };

    vm[idx(0, 0)] = 0.0;

    for i in 0..=n {
        for j in 0..=m {
            // --- M state ---
            if i > 0 && j > 0 {
                let e = emit(seq_a[i - 1], seq_b[j - 1], params);
                let from_m = vm[idx(i - 1, j - 1)] + params.log_trans_mm;
                let from_x = vx[idx(i - 1, j - 1)] + params.log_trans_xm;
                let from_y = vy[idx(i - 1, j - 1)] + params.log_trans_ym;

                let (best, bp) = max3(from_m, from_x, from_y);
                vm[idx(i, j)] = e + best;
                bp_m[idx(i, j)] = bp;
            }

            // --- X state ---
            if i > 0 {
                let from_m = vm[idx(i - 1, j)] + params.log_trans_mx;
                let from_x = vx[idx(i - 1, j)] + params.log_trans_xx;

                if from_m >= from_x {
                    vx[idx(i, j)] = from_m;
                    bp_x[idx(i, j)] = BackPtr::FromM;
                } else {
                    vx[idx(i, j)] = from_x;
                    bp_x[idx(i, j)] = BackPtr::FromX;
                }
            }

            // --- Y state ---
            if j > 0 {
                let from_m = vm[idx(i, j - 1)] + params.log_trans_my;
                let from_y = vy[idx(i, j - 1)] + params.log_trans_yy;

                if from_m >= from_y {
                    vy[idx(i, j)] = from_m;
                    bp_y[idx(i, j)] = BackPtr::FromM;
                } else {
                    vy[idx(i, j)] = from_y;
                    bp_y[idx(i, j)] = BackPtr::FromY;
                }
            }
        }
    }

    // Determine which terminal state has the best score.
    let score_m = vm[idx(n, m)];
    let score_x = vx[idx(n, m)];
    let score_y = vy[idx(n, m)];

    let (score, mut current_state) = if score_m >= score_x && score_m >= score_y {
        (score_m, PairHmmState::Match)
    } else if score_x >= score_y {
        (score_x, PairHmmState::InsertX)
    } else {
        (score_y, PairHmmState::InsertY)
    };

    // Traceback
    let mut states = Vec::new();
    let mut i = n;
    let mut j = m;

    while i > 0 || j > 0 {
        states.push(current_state);

        match current_state {
            PairHmmState::Match => {
                let bp = bp_m[idx(i, j)];
                i -= 1;
                j -= 1;
                current_state = backptr_to_state(bp);
            }
            PairHmmState::InsertX => {
                let bp = bp_x[idx(i, j)];
                i -= 1;
                current_state = backptr_to_state(bp);
            }
            PairHmmState::InsertY => {
                let bp = bp_y[idx(i, j)];
                j -= 1;
                current_state = backptr_to_state(bp);
            }
        }
    }

    states.reverse();

    Ok(PairHmmAlignment { score, states })
}

/// Pick the maximum of three values and return the corresponding backpointer.
fn max3(from_m: f64, from_x: f64, from_y: f64) -> (f64, BackPtr) {
    if from_m >= from_x && from_m >= from_y {
        (from_m, BackPtr::FromM)
    } else if from_x >= from_y {
        (from_x, BackPtr::FromX)
    } else {
        (from_y, BackPtr::FromY)
    }
}

/// Convert a backpointer tag to the corresponding HMM state.
fn backptr_to_state(bp: BackPtr) -> PairHmmState {
    match bp {
        BackPtr::FromM | BackPtr::None => PairHmmState::Match,
        BackPtr::FromX => PairHmmState::InsertX,
        BackPtr::FromY => PairHmmState::InsertY,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    fn default_params() -> PairHmmParams {
        PairHmmParams::default()
    }

    // ------------------------------------------------------------------
    // Forward algorithm
    // ------------------------------------------------------------------

    #[test]
    fn forward_identical_sequences() {
        let params = default_params();
        let score = pair_hmm_forward(b"ACGT", b"ACGT", &params).unwrap();
        // Identical sequences should produce a reasonable (finite, negative) log-likelihood.
        assert!(score.is_finite());
        assert!(score < 0.0);
    }

    #[test]
    fn forward_identical_higher_than_mismatch() {
        let params = default_params();
        let identical = pair_hmm_forward(b"ACGT", b"ACGT", &params).unwrap();
        let mismatch = pair_hmm_forward(b"ACGT", b"ACGA", &params).unwrap();
        assert!(
            identical > mismatch,
            "identical ({}) should score higher than mismatch ({})",
            identical,
            mismatch
        );
    }

    #[test]
    fn forward_error_on_empty_a() {
        let params = default_params();
        let result = pair_hmm_forward(b"", b"ACGT", &params);
        assert!(result.is_err());
    }

    #[test]
    fn forward_error_on_empty_b() {
        let params = default_params();
        let result = pair_hmm_forward(b"ACGT", b"", &params);
        assert!(result.is_err());
    }

    #[test]
    fn forward_error_on_both_empty() {
        let params = default_params();
        let result = pair_hmm_forward(b"", b"", &params);
        assert!(result.is_err());
    }

    #[test]
    fn forward_short_length_one() {
        let params = default_params();
        let score = pair_hmm_forward(b"A", b"A", &params).unwrap();
        assert!(score.is_finite());
        assert!(score < 0.0);
    }

    #[test]
    fn forward_short_length_two() {
        let params = default_params();
        let score = pair_hmm_forward(b"AC", b"AC", &params).unwrap();
        assert!(score.is_finite());
    }

    #[test]
    fn forward_different_lengths() {
        let params = default_params();
        // Different lengths require insert states.
        let score = pair_hmm_forward(b"ACGT", b"ACG", &params).unwrap();
        assert!(score.is_finite());
    }

    // ------------------------------------------------------------------
    // Viterbi algorithm
    // ------------------------------------------------------------------

    #[test]
    fn viterbi_identical_all_match() {
        let params = default_params();
        let aln = pair_hmm_viterbi(b"ACGT", b"ACGT", &params).unwrap();
        assert!(aln.score.is_finite());
        assert!(
            aln.states.iter().all(|s| *s == PairHmmState::Match),
            "identical sequences should align entirely through Match states"
        );
    }

    #[test]
    fn viterbi_identical_score_reasonable() {
        let params = default_params();
        let aln = pair_hmm_viterbi(b"ACGT", b"ACGT", &params).unwrap();
        assert!(aln.score < 0.0, "log-probability should be negative");
    }

    #[test]
    fn viterbi_mismatch_lower_score() {
        let params = default_params();
        let identical = pair_hmm_viterbi(b"ACGT", b"ACGT", &params).unwrap();
        let mismatch = pair_hmm_viterbi(b"ACGT", b"ACGA", &params).unwrap();
        assert!(
            identical.score > mismatch.score,
            "identical ({}) should score higher than mismatch ({})",
            identical.score,
            mismatch.score
        );
    }

    #[test]
    fn viterbi_gap_produces_insert_state() {
        // seq_a is longer — a gap in B should produce InsertX states.
        let params = default_params();
        let aln = pair_hmm_viterbi(b"ACGTT", b"ACGT", &params).unwrap();
        assert!(
            aln.states.contains(&PairHmmState::InsertX),
            "longer A should produce InsertX states, got {:?}",
            aln.states
        );
    }

    #[test]
    fn viterbi_gap_in_a_produces_insert_y() {
        // seq_b is longer — a gap in A should produce InsertY states.
        let params = default_params();
        let aln = pair_hmm_viterbi(b"ACGT", b"ACGTT", &params).unwrap();
        assert!(
            aln.states.contains(&PairHmmState::InsertY),
            "longer B should produce InsertY states, got {:?}",
            aln.states
        );
    }

    #[test]
    fn viterbi_error_on_empty() {
        let params = default_params();
        assert!(pair_hmm_viterbi(b"", b"ACGT", &params).is_err());
        assert!(pair_hmm_viterbi(b"ACGT", b"", &params).is_err());
        assert!(pair_hmm_viterbi(b"", b"", &params).is_err());
    }

    #[test]
    fn viterbi_path_length_equals_sequence_consumption() {
        // Each state consumes exactly one position: Match consumes one from
        // each sequence, InsertX consumes one from A, InsertY consumes one
        // from B. So the total consumed must equal (len_a, len_b).
        let params = default_params();

        for (a, b) in &[
            (b"ACGT".as_slice(), b"ACGT".as_slice()),
            (b"ACGTT".as_slice(), b"ACGT".as_slice()),
            (b"ACG".as_slice(), b"ACGTTT".as_slice()),
            (b"A".as_slice(), b"T".as_slice()),
        ] {
            let aln = pair_hmm_viterbi(a, b, &params).unwrap();
            let consumed_a: usize = aln
                .states
                .iter()
                .filter(|s| matches!(s, PairHmmState::Match | PairHmmState::InsertX))
                .count();
            let consumed_b: usize = aln
                .states
                .iter()
                .filter(|s| matches!(s, PairHmmState::Match | PairHmmState::InsertY))
                .count();
            assert_eq!(
                consumed_a,
                a.len(),
                "consumed_a={} != len_a={} for {:?} vs {:?}",
                consumed_a,
                a.len(),
                std::str::from_utf8(a),
                std::str::from_utf8(b),
            );
            assert_eq!(
                consumed_b,
                b.len(),
                "consumed_b={} != len_b={} for {:?} vs {:?}",
                consumed_b,
                b.len(),
                std::str::from_utf8(a),
                std::str::from_utf8(b),
            );
        }
    }

    // ------------------------------------------------------------------
    // Forward >= Viterbi
    // ------------------------------------------------------------------

    #[test]
    fn forward_score_gte_viterbi_score() {
        // The forward algorithm marginalizes over all paths, so its result
        // must be >= the single best path found by Viterbi.
        let params = default_params();

        for (a, b) in &[
            (b"ACGT".as_slice(), b"ACGT".as_slice()),
            (b"ACGT".as_slice(), b"ACGA".as_slice()),
            (b"ACGTT".as_slice(), b"ACGT".as_slice()),
            (b"A".as_slice(), b"C".as_slice()),
        ] {
            let fwd = pair_hmm_forward(a, b, &params).unwrap();
            let vit = pair_hmm_viterbi(a, b, &params).unwrap();
            assert!(
                fwd >= vit.score - 1e-10,
                "forward ({}) should be >= viterbi ({}) for {:?} vs {:?}",
                fwd,
                vit.score,
                std::str::from_utf8(a),
                std::str::from_utf8(b),
            );
        }
    }

    // ------------------------------------------------------------------
    // Custom parameters
    // ------------------------------------------------------------------

    #[test]
    fn custom_params_high_gap_penalty() {
        // With very low gap-open probability, identical-length sequences
        // should strongly prefer the all-Match path.
        let params = PairHmmParams {
            log_match_emission: (0.99_f64).ln(),
            log_mismatch_emission: (0.01_f64 / 3.0).ln(),
            log_trans_mm: (0.99_f64).ln(),
            log_trans_mx: (0.005_f64).ln(),
            log_trans_my: (0.005_f64).ln(),
            log_trans_xx: (0.1_f64).ln(),
            log_trans_xm: (0.9_f64).ln(),
            log_trans_yy: (0.1_f64).ln(),
            log_trans_ym: (0.9_f64).ln(),
        };

        let aln = pair_hmm_viterbi(b"ACGTACGT", b"ACGTACGT", &params).unwrap();
        assert!(aln.states.iter().all(|s| *s == PairHmmState::Match));
    }

    #[test]
    fn custom_params_high_gap_extend() {
        // With high gap-extend probability, a gap that opens should tend
        // to stay open for multiple positions.
        let params = PairHmmParams {
            log_match_emission: (0.98_f64).ln(),
            log_mismatch_emission: (0.02_f64 / 3.0).ln(),
            log_trans_mm: (0.85_f64).ln(),
            log_trans_mx: (0.10_f64).ln(),
            log_trans_my: (0.05_f64).ln(),
            log_trans_xx: (0.8_f64).ln(),
            log_trans_xm: (0.2_f64).ln(),
            log_trans_yy: (0.8_f64).ln(),
            log_trans_ym: (0.2_f64).ln(),
        };

        // A is much longer — should see several consecutive InsertX states.
        let aln = pair_hmm_viterbi(b"ACGTTTTT", b"ACGT", &params).unwrap();
        let insert_x_count = aln
            .states
            .iter()
            .filter(|s| **s == PairHmmState::InsertX)
            .count();
        assert!(
            insert_x_count >= 3,
            "expected at least 3 InsertX states, got {}",
            insert_x_count
        );
    }

    // ------------------------------------------------------------------
    // Short sequences
    // ------------------------------------------------------------------

    #[test]
    fn viterbi_single_base_match() {
        let params = default_params();
        let aln = pair_hmm_viterbi(b"A", b"A", &params).unwrap();
        assert_eq!(aln.states.len(), 1);
        assert_eq!(aln.states[0], PairHmmState::Match);
    }

    #[test]
    fn viterbi_single_base_mismatch() {
        let params = default_params();
        let aln = pair_hmm_viterbi(b"A", b"T", &params).unwrap();
        assert_eq!(aln.states.len(), 1);
        assert_eq!(aln.states[0], PairHmmState::Match);
    }

    #[test]
    fn viterbi_two_bases() {
        let params = default_params();
        let aln = pair_hmm_viterbi(b"AC", b"AC", &params).unwrap();
        assert_eq!(aln.states.len(), 2);
        assert!(aln.states.iter().all(|s| *s == PairHmmState::Match));
    }

    // ------------------------------------------------------------------
    // log_sum_exp correctness
    // ------------------------------------------------------------------

    #[test]
    fn log_sum_exp_basic() {
        let a = 2.0_f64.ln();
        let b = 3.0_f64.ln();
        let result = log_sum_exp(a, b);
        let expected = 5.0_f64.ln();
        assert!(
            (result - expected).abs() < 1e-12,
            "log_sum_exp({}, {}) = {}, expected {}",
            a,
            b,
            result,
            expected
        );
    }

    #[test]
    fn log_sum_exp_neg_infinity() {
        assert_eq!(log_sum_exp(f64::NEG_INFINITY, 0.0), 0.0);
        assert_eq!(log_sum_exp(0.0, f64::NEG_INFINITY), 0.0);
        assert_eq!(
            log_sum_exp(f64::NEG_INFINITY, f64::NEG_INFINITY),
            f64::NEG_INFINITY
        );
    }

    // ------------------------------------------------------------------
    // Default params sanity
    // ------------------------------------------------------------------

    #[test]
    fn default_params_produce_finite_results() {
        let params = default_params();
        let fwd = pair_hmm_forward(b"ACGTACGT", b"ACGAACGT", &params).unwrap();
        assert!(fwd.is_finite());

        let vit = pair_hmm_viterbi(b"ACGTACGT", b"ACGAACGT", &params).unwrap();
        assert!(vit.score.is_finite());
        assert!(!vit.states.is_empty());
    }
}
