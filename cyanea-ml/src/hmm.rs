//! Hidden Markov Models for biological sequence analysis.
//!
//! Provides a discrete HMM implementation suitable for bioinformatics tasks
//! such as gene finding, CpG island detection, protein secondary structure
//! prediction, and profile HMM-based sequence classification.
//!
//! All internal computations use log-space arithmetic to avoid numerical
//! underflow on long observation sequences.
//!
//! # Quick start
//!
//! ```
//! use cyanea_ml::hmm::HmmModel;
//!
//! // 2-state fair/loaded coin model
//! let initial = vec![0.5, 0.5];
//! let transition = vec![
//!     0.9, 0.1,  // fair  -> fair, loaded
//!     0.2, 0.8,  // loaded -> fair, loaded
//! ];
//! let emission = vec![
//!     0.5, 0.5,  // fair:   P(H), P(T)
//!     0.8, 0.2,  // loaded: P(H), P(T)
//! ];
//!
//! let model = HmmModel::new(2, 2, initial, transition, emission).unwrap();
//! let obs = vec![0, 0, 1, 0, 0]; // H H T H H
//! let (path, score) = model.viterbi(&obs).unwrap();
//! assert_eq!(path.len(), obs.len());
//! ```

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Log-space helpers
// ---------------------------------------------------------------------------

/// Numerically stable computation of `log(exp(a) + exp(b))`.
///
/// Handles the cases where `a` or `b` are negative infinity.
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

/// Log-sum-exp over a slice.
fn log_sum_exp_slice(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return f64::NEG_INFINITY;
    }
    let max = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if max == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    let sum: f64 = xs.iter().map(|&x| (x - max).exp()).sum();
    max + sum.ln()
}

// ---------------------------------------------------------------------------
// Epsilon for avoiding log(0)
// ---------------------------------------------------------------------------

const EPSILON: f64 = 1e-10;

// ---------------------------------------------------------------------------
// HmmModel
// ---------------------------------------------------------------------------

/// A discrete Hidden Markov Model.
///
/// Parameters are stored in probability space but all algorithms operate in
/// log-space internally to avoid underflow on long sequences.
#[derive(Debug, Clone)]
pub struct HmmModel {
    /// Number of hidden states.
    n_states: usize,
    /// Number of observable symbols.
    n_symbols: usize,
    /// Initial state probabilities pi[i] (length `n_states`).
    initial: Vec<f64>,
    /// Transition matrix A[i][j] = P(state_j | state_i), stored row-major
    /// as `Vec<f64>` of size `n_states * n_states`.
    transition: Vec<f64>,
    /// Emission matrix B[i][k] = P(symbol_k | state_i), stored row-major
    /// as `Vec<f64>` of size `n_states * n_symbols`.
    emission: Vec<f64>,
}

impl HmmModel {
    /// Create a new HMM after validating dimensions and probability constraints.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - `n_states` or `n_symbols` is zero
    /// - Vector dimensions do not match the declared sizes
    /// - Any probability row does not sum to approximately 1.0 (tolerance 1e-6)
    pub fn new(
        n_states: usize,
        n_symbols: usize,
        initial: Vec<f64>,
        transition: Vec<f64>,
        emission: Vec<f64>,
    ) -> Result<Self> {
        if n_states == 0 {
            return Err(CyaneaError::InvalidInput("n_states must be > 0".into()));
        }
        if n_symbols == 0 {
            return Err(CyaneaError::InvalidInput("n_symbols must be > 0".into()));
        }
        if initial.len() != n_states {
            return Err(CyaneaError::InvalidInput(format!(
                "initial length {} != n_states {}",
                initial.len(),
                n_states
            )));
        }
        if transition.len() != n_states * n_states {
            return Err(CyaneaError::InvalidInput(format!(
                "transition length {} != n_states*n_states {}",
                transition.len(),
                n_states * n_states
            )));
        }
        if emission.len() != n_states * n_symbols {
            return Err(CyaneaError::InvalidInput(format!(
                "emission length {} != n_states*n_symbols {}",
                emission.len(),
                n_states * n_symbols
            )));
        }

        // Validate probability sums.
        let tol = 1e-6;

        let pi_sum: f64 = initial.iter().sum();
        if (pi_sum - 1.0).abs() > tol {
            return Err(CyaneaError::InvalidInput(format!(
                "initial probabilities sum to {pi_sum}, expected ~1.0"
            )));
        }

        for i in 0..n_states {
            let row_sum: f64 = transition[i * n_states..(i + 1) * n_states].iter().sum();
            if (row_sum - 1.0).abs() > tol {
                return Err(CyaneaError::InvalidInput(format!(
                    "transition row {i} sums to {row_sum}, expected ~1.0"
                )));
            }
        }

        for i in 0..n_states {
            let row_sum: f64 = emission[i * n_symbols..(i + 1) * n_symbols].iter().sum();
            if (row_sum - 1.0).abs() > tol {
                return Err(CyaneaError::InvalidInput(format!(
                    "emission row {i} sums to {row_sum}, expected ~1.0"
                )));
            }
        }

        Ok(Self {
            n_states,
            n_symbols,
            initial,
            transition,
            emission,
        })
    }

    /// Number of hidden states.
    pub fn n_states(&self) -> usize {
        self.n_states
    }

    /// Number of observable symbols.
    pub fn n_symbols(&self) -> usize {
        self.n_symbols
    }

    // -----------------------------------------------------------------------
    // Validation helper
    // -----------------------------------------------------------------------

    /// Validate an observation sequence, returning an error if it is empty or
    /// contains out-of-range symbols.
    fn validate_observations(&self, observations: &[usize]) -> Result<()> {
        if observations.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "observation sequence is empty".into(),
            ));
        }
        for (t, &o) in observations.iter().enumerate() {
            if o >= self.n_symbols {
                return Err(CyaneaError::InvalidInput(format!(
                    "observation[{t}] = {o} out of range (n_symbols = {})",
                    self.n_symbols
                )));
            }
        }
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Forward algorithm
    // -----------------------------------------------------------------------

    /// Run the forward algorithm in log-space.
    ///
    /// Returns `(alpha, log_likelihood)` where `alpha[t][i]` is the log
    /// probability of observing `o_0..o_t` and being in state `i` at time `t`.
    ///
    /// # Errors
    ///
    /// Returns an error for empty or invalid observation sequences.
    pub fn forward(&self, observations: &[usize]) -> Result<(Vec<Vec<f64>>, f64)> {
        self.validate_observations(observations)?;

        let n = self.n_states;
        let t_len = observations.len();
        let mut alpha = vec![vec![f64::NEG_INFINITY; n]; t_len];

        // Initialization: alpha[0][i] = log(pi[i]) + log(B[i][o_0])
        let o0 = observations[0];
        for i in 0..n {
            alpha[0][i] =
                (self.initial[i] + EPSILON).ln() + (self.emission[i * self.n_symbols + o0] + EPSILON).ln();
        }

        // Induction
        for t in 1..t_len {
            let ot = observations[t];
            for j in 0..n {
                let mut acc = f64::NEG_INFINITY;
                for i in 0..n {
                    let v = alpha[t - 1][i] + (self.transition[i * n + j] + EPSILON).ln();
                    acc = log_sum_exp(acc, v);
                }
                alpha[t][j] = acc + (self.emission[j * self.n_symbols + ot] + EPSILON).ln();
            }
        }

        // Termination
        let ll = log_sum_exp_slice(&alpha[t_len - 1]);

        Ok((alpha, ll))
    }

    /// Compute the log-likelihood of an observation sequence.
    ///
    /// Convenience wrapper around [`forward`](Self::forward).
    pub fn log_likelihood(&self, observations: &[usize]) -> Result<f64> {
        let (_, ll) = self.forward(observations)?;
        Ok(ll)
    }

    // -----------------------------------------------------------------------
    // Backward algorithm
    // -----------------------------------------------------------------------

    /// Run the backward algorithm in log-space.
    ///
    /// Returns `beta` where `beta[t][i]` is the log probability of observing
    /// `o_{t+1}..o_{T-1}` given state `i` at time `t`.
    ///
    /// # Errors
    ///
    /// Returns an error for empty or invalid observation sequences.
    pub fn backward(&self, observations: &[usize]) -> Result<Vec<Vec<f64>>> {
        self.validate_observations(observations)?;

        let n = self.n_states;
        let t_len = observations.len();
        let mut beta = vec![vec![f64::NEG_INFINITY; n]; t_len];

        // Initialization: beta[T-1][i] = log(1) = 0
        for i in 0..n {
            beta[t_len - 1][i] = 0.0;
        }

        // Induction (backwards)
        for t in (0..t_len - 1).rev() {
            let ot1 = observations[t + 1];
            for i in 0..n {
                let mut acc = f64::NEG_INFINITY;
                for j in 0..n {
                    let v = (self.transition[i * n + j] + EPSILON).ln()
                        + (self.emission[j * self.n_symbols + ot1] + EPSILON).ln()
                        + beta[t + 1][j];
                    acc = log_sum_exp(acc, v);
                }
                beta[t][i] = acc;
            }
        }

        Ok(beta)
    }

    // -----------------------------------------------------------------------
    // Viterbi algorithm
    // -----------------------------------------------------------------------

    /// Viterbi decoding: find the most likely state sequence.
    ///
    /// Returns `(path, log_probability)` where `path[t]` is the most likely
    /// state at time `t` and `log_probability` is the log probability of
    /// that best path.
    ///
    /// # Errors
    ///
    /// Returns an error for empty or invalid observation sequences.
    pub fn viterbi(&self, observations: &[usize]) -> Result<(Vec<usize>, f64)> {
        self.validate_observations(observations)?;

        let n = self.n_states;
        let t_len = observations.len();

        let mut delta = vec![vec![f64::NEG_INFINITY; n]; t_len];
        let mut psi = vec![vec![0usize; n]; t_len];

        // Initialization
        let o0 = observations[0];
        for i in 0..n {
            delta[0][i] =
                (self.initial[i] + EPSILON).ln() + (self.emission[i * self.n_symbols + o0] + EPSILON).ln();
        }

        // Recursion
        for t in 1..t_len {
            let ot = observations[t];
            for j in 0..n {
                let mut best_val = f64::NEG_INFINITY;
                let mut best_state = 0;
                for i in 0..n {
                    let v = delta[t - 1][i] + (self.transition[i * n + j] + EPSILON).ln();
                    if v > best_val {
                        best_val = v;
                        best_state = i;
                    }
                }
                delta[t][j] = best_val + (self.emission[j * self.n_symbols + ot] + EPSILON).ln();
                psi[t][j] = best_state;
            }
        }

        // Termination: find best final state
        let mut best_final = 0usize;
        let mut best_score = f64::NEG_INFINITY;
        for i in 0..n {
            if delta[t_len - 1][i] > best_score {
                best_score = delta[t_len - 1][i];
                best_final = i;
            }
        }

        // Backtrack
        let mut path = vec![0usize; t_len];
        path[t_len - 1] = best_final;
        for t in (0..t_len - 1).rev() {
            path[t] = psi[t + 1][path[t + 1]];
        }

        Ok((path, best_score))
    }

    // -----------------------------------------------------------------------
    // Baum-Welch (EM) training
    // -----------------------------------------------------------------------

    /// Train the HMM parameters using the Baum-Welch (EM) algorithm.
    ///
    /// Iterates until the log-likelihood improvement is below `tolerance` or
    /// `max_iter` iterations have been performed. Returns the final
    /// log-likelihood.
    ///
    /// A small epsilon (1e-10) is added to avoid log(0) during computation.
    ///
    /// # Errors
    ///
    /// Returns an error for empty or invalid observation sequences.
    pub fn baum_welch(
        &mut self,
        observations: &[usize],
        max_iter: usize,
        tolerance: f64,
    ) -> Result<f64> {
        self.validate_observations(observations)?;

        let n = self.n_states;
        let m = self.n_symbols;
        let t_len = observations.len();
        let mut prev_ll = f64::NEG_INFINITY;

        for _iter in 0..max_iter {
            // E-step: forward and backward
            let (alpha, ll) = self.forward(observations)?;
            let beta = self.backward(observations)?;

            // Check convergence
            if (ll - prev_ll).abs() < tolerance && prev_ll != f64::NEG_INFINITY {
                return Ok(ll);
            }
            prev_ll = ll;

            // Compute gamma[t][i] = P(state_i at time t | observations, model)
            // log_gamma[t][i] = alpha[t][i] + beta[t][i] - ll
            let mut gamma = vec![vec![0.0f64; n]; t_len];
            for t in 0..t_len {
                let mut log_gamma_t = vec![f64::NEG_INFINITY; n];
                for i in 0..n {
                    log_gamma_t[i] = alpha[t][i] + beta[t][i] - ll;
                }
                // Convert from log to probability space for accumulation
                for i in 0..n {
                    gamma[t][i] = log_gamma_t[i].exp();
                }
            }

            // Compute xi[t][i][j] = P(state_i at t, state_j at t+1 | obs, model)
            // Stored flat: xi[t] is n*n row-major
            let mut xi = vec![vec![0.0f64; n * n]; t_len.saturating_sub(1)];
            for t in 0..t_len.saturating_sub(1) {
                let ot1 = observations[t + 1];
                for i in 0..n {
                    for j in 0..n {
                        let log_xi = alpha[t][i]
                            + (self.transition[i * n + j] + EPSILON).ln()
                            + (self.emission[j * m + ot1] + EPSILON).ln()
                            + beta[t + 1][j]
                            - ll;
                        xi[t][i * n + j] = log_xi.exp();
                    }
                }
            }

            // M-step: re-estimate parameters

            // Initial probabilities
            for i in 0..n {
                self.initial[i] = gamma[0][i] + EPSILON;
            }
            let pi_sum: f64 = self.initial.iter().sum();
            for i in 0..n {
                self.initial[i] /= pi_sum;
            }

            // Transition probabilities
            for i in 0..n {
                let gamma_sum: f64 =
                    (0..t_len.saturating_sub(1)).map(|t| gamma[t][i]).sum::<f64>() + EPSILON;
                for j in 0..n {
                    let xi_sum: f64 =
                        (0..t_len.saturating_sub(1)).map(|t| xi[t][i * n + j]).sum::<f64>() + EPSILON;
                    self.transition[i * n + j] = xi_sum / gamma_sum;
                }
                // Normalize row
                let row_sum: f64 = self.transition[i * n..(i + 1) * n].iter().sum();
                if row_sum > 0.0 {
                    for j in 0..n {
                        self.transition[i * n + j] /= row_sum;
                    }
                }
            }

            // Emission probabilities
            for i in 0..n {
                let gamma_sum: f64 = (0..t_len).map(|t| gamma[t][i]).sum::<f64>() + EPSILON;
                for k in 0..m {
                    let numer: f64 = (0..t_len)
                        .filter(|&t| observations[t] == k)
                        .map(|t| gamma[t][i])
                        .sum::<f64>()
                        + EPSILON;
                    self.emission[i * m + k] = numer / gamma_sum;
                }
                // Normalize row
                let row_sum: f64 = self.emission[i * m..(i + 1) * m].iter().sum();
                if row_sum > 0.0 {
                    for k in 0..m {
                        self.emission[i * m + k] /= row_sum;
                    }
                }
            }
        }

        // Return final log-likelihood after all iterations
        let (_, ll) = self.forward(observations)?;
        Ok(ll)
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: build a 2-state fair/loaded coin HMM.
    fn coin_hmm() -> HmmModel {
        let initial = vec![0.5, 0.5];
        let transition = vec![
            0.9, 0.1, // fair -> fair, loaded
            0.2, 0.8, // loaded -> fair, loaded
        ];
        let emission = vec![
            0.5, 0.5, // fair:   P(H), P(T)
            0.8, 0.2, // loaded: P(H), P(T)
        ];
        HmmModel::new(2, 2, initial, transition, emission).unwrap()
    }

    /// Helper: build a 2-state "occasionally dishonest casino" with a 6-sided die.
    fn casino_hmm() -> HmmModel {
        let initial = vec![0.5, 0.5];
        let transition = vec![
            0.95, 0.05, // fair -> fair, loaded
            0.10, 0.90, // loaded -> fair, loaded
        ];
        // Fair die: uniform 1/6 each
        // Loaded die: heavily biased toward 6 (symbol index 5)
        let mut emission = vec![0.0; 2 * 6];
        for k in 0..6 {
            emission[0 * 6 + k] = 1.0 / 6.0; // fair
        }
        emission[1 * 6 + 0] = 0.1;
        emission[1 * 6 + 1] = 0.1;
        emission[1 * 6 + 2] = 0.1;
        emission[1 * 6 + 3] = 0.1;
        emission[1 * 6 + 4] = 0.1;
        emission[1 * 6 + 5] = 0.5; // loaded favors 6
        HmmModel::new(2, 6, initial, transition, emission).unwrap()
    }

    // -----------------------------------------------------------------------
    // Forward algorithm tests
    // -----------------------------------------------------------------------

    #[test]
    fn forward_log_likelihood_is_finite_and_negative() {
        let model = coin_hmm();
        let obs = vec![0, 0, 1, 0, 1, 1, 0];
        let (alpha, ll) = model.forward(&obs).unwrap();
        assert!(ll.is_finite(), "log-likelihood should be finite");
        assert!(ll < 0.0, "log-likelihood should be negative");
        assert_eq!(alpha.len(), obs.len());
        assert_eq!(alpha[0].len(), 2);
    }

    // -----------------------------------------------------------------------
    // Viterbi tests
    // -----------------------------------------------------------------------

    #[test]
    fn viterbi_path_length_equals_observations() {
        let model = coin_hmm();
        let obs = vec![0, 1, 0, 0, 1, 0, 1, 1, 0, 0];
        let (path, _score) = model.viterbi(&obs).unwrap();
        assert_eq!(path.len(), obs.len());
    }

    #[test]
    fn viterbi_path_values_are_valid_states() {
        let model = coin_hmm();
        let obs = vec![0, 1, 0, 0, 1];
        let (path, _score) = model.viterbi(&obs).unwrap();
        for &s in &path {
            assert!(s < model.n_states(), "state {s} out of range");
        }
    }

    // -----------------------------------------------------------------------
    // Forward-backward consistency
    // -----------------------------------------------------------------------

    #[test]
    fn forward_backward_log_likelihood_match() {
        let model = coin_hmm();
        let obs = vec![0, 1, 0, 0, 1, 0];

        let (alpha, ll_fwd) = model.forward(&obs).unwrap();
        let beta = model.backward(&obs).unwrap();

        // Compute log-likelihood from alpha[0] + beta[0] (at t=0)
        let mut terms = vec![f64::NEG_INFINITY; model.n_states()];
        for i in 0..model.n_states() {
            terms[i] = alpha[0][i] + beta[0][i];
        }
        let ll_from_ab = log_sum_exp_slice(&terms);

        assert!(
            (ll_fwd - ll_from_ab).abs() < 1e-6,
            "forward LL ({ll_fwd}) and alpha+beta LL ({ll_from_ab}) should match"
        );
    }

    // -----------------------------------------------------------------------
    // Baum-Welch tests
    // -----------------------------------------------------------------------

    #[test]
    fn baum_welch_improves_likelihood_monotonically() {
        // Start with a somewhat off model and verify each iteration improves
        let initial = vec![0.6, 0.4];
        let transition = vec![0.7, 0.3, 0.4, 0.6];
        let emission = vec![0.6, 0.4, 0.3, 0.7];
        let mut model = HmmModel::new(2, 2, initial, transition, emission).unwrap();

        let obs = vec![0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0];

        let mut prev_ll = model.log_likelihood(&obs).unwrap();
        for _ in 0..10 {
            let ll = model.baum_welch(&obs, 1, 0.0).unwrap();
            assert!(
                ll >= prev_ll - 1e-10,
                "LL should not decrease: prev={prev_ll}, curr={ll}"
            );
            prev_ll = ll;
        }
    }

    #[test]
    fn baum_welch_converges() {
        let mut model = coin_hmm();
        let obs = vec![0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0];

        let ll_before = model.log_likelihood(&obs).unwrap();
        let ll_after = model.baum_welch(&obs, 100, 1e-8).unwrap();

        assert!(
            ll_after >= ll_before - 1e-10,
            "Baum-Welch should not decrease LL: before={ll_before}, after={ll_after}"
        );
        assert!(ll_after.is_finite());
    }

    // -----------------------------------------------------------------------
    // Casino HMM (known example)
    // -----------------------------------------------------------------------

    #[test]
    fn casino_viterbi_detects_loaded_runs() {
        let model = casino_hmm();
        // Sequence with a run of 6s (symbol index 5) in the middle
        // Fair region, then loaded region, then fair again
        let obs = vec![0, 2, 3, 1, 4, 0, 5, 5, 5, 5, 5, 5, 2, 1, 0, 3];
        let (path, score) = model.viterbi(&obs).unwrap();

        assert_eq!(path.len(), obs.len());
        assert!(score.is_finite());

        // The run of 5s (loaded die showing 6) should mostly be decoded as
        // the loaded state (state 1)
        let loaded_count: usize = (6..12).filter(|&t| path[t] == 1).count();
        assert!(
            loaded_count >= 4,
            "expected at least 4 of 6 positions in the 6-run to be decoded as loaded, got {loaded_count}"
        );
    }

    // -----------------------------------------------------------------------
    // Error handling tests
    // -----------------------------------------------------------------------

    #[test]
    fn error_on_empty_observations() {
        let model = coin_hmm();
        assert!(model.forward(&[]).is_err());
        assert!(model.viterbi(&[]).is_err());
        assert!(model.backward(&[]).is_err());
        assert!(model.log_likelihood(&[]).is_err());
    }

    #[test]
    fn error_on_invalid_symbol() {
        let model = coin_hmm();
        // Symbol 5 is out of range for n_symbols=2
        assert!(model.forward(&[0, 1, 5]).is_err());
        assert!(model.viterbi(&[0, 1, 5]).is_err());
    }

    #[test]
    fn error_on_dimension_mismatch() {
        // Initial wrong length
        assert!(HmmModel::new(2, 2, vec![1.0], vec![0.5; 4], vec![0.5; 4]).is_err());
        // Transition wrong length
        assert!(HmmModel::new(2, 2, vec![0.5, 0.5], vec![0.5; 3], vec![0.5; 4]).is_err());
        // Emission wrong length
        assert!(HmmModel::new(2, 2, vec![0.5, 0.5], vec![0.5; 4], vec![0.5; 3]).is_err());
        // n_states = 0
        assert!(HmmModel::new(0, 2, vec![], vec![], vec![]).is_err());
        // n_symbols = 0
        assert!(HmmModel::new(2, 0, vec![0.5, 0.5], vec![0.25; 4], vec![]).is_err());
        // Probabilities don't sum to 1
        assert!(HmmModel::new(2, 2, vec![0.3, 0.3], vec![0.5; 4], vec![0.5; 4]).is_err());
    }

    // -----------------------------------------------------------------------
    // log_sum_exp numerical stability
    // -----------------------------------------------------------------------

    #[test]
    fn log_sum_exp_numerical_stability() {
        // Basic property: log(exp(a) + exp(b)) >= max(a, b)
        let result = log_sum_exp(-1000.0, -1001.0);
        assert!(result.is_finite());
        assert!(result >= -1000.0);
        assert!(result < -999.0);

        // log(exp(0) + exp(0)) = log(2)
        let r2 = log_sum_exp(0.0, 0.0);
        assert!((r2 - 2.0_f64.ln()).abs() < 1e-12);

        // NEG_INFINITY identity
        assert_eq!(log_sum_exp(f64::NEG_INFINITY, 5.0), 5.0);
        assert_eq!(log_sum_exp(5.0, f64::NEG_INFINITY), 5.0);

        // Very large values should not overflow
        let big = log_sum_exp(700.0, 700.0);
        assert!(big.is_finite());
        assert!((big - (700.0 + 2.0_f64.ln())).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // Round-trip: generate observations from model, then decode
    // -----------------------------------------------------------------------

    #[test]
    fn round_trip_generate_decode() {
        // Use a very biased model so the Viterbi path is predictable.
        // State 0 strongly emits symbol 0; state 1 strongly emits symbol 1.
        // Transitions strongly self-loop.
        let initial = vec![1.0, 0.0]; // always start in state 0
        let transition = vec![
            0.95, 0.05, // state 0: stay in 0 most of the time
            0.05, 0.95, // state 1: stay in 1 most of the time
        ];
        let emission = vec![
            0.95, 0.05, // state 0 emits symbol 0
            0.05, 0.95, // state 1 emits symbol 1
        ];
        let model = HmmModel::new(2, 2, initial, transition, emission).unwrap();

        // Hand-craft an observation that looks like: state0 region, then state1 region
        let obs: Vec<usize> = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1].to_vec();
        let (path, score) = model.viterbi(&obs).unwrap();

        assert_eq!(path.len(), 10);
        assert!(score.is_finite());

        // First half should mostly be state 0, second half state 1
        let state0_first_half: usize = (0..5).filter(|&t| path[t] == 0).count();
        let state1_second_half: usize = (5..10).filter(|&t| path[t] == 1).count();

        assert!(
            state0_first_half >= 4,
            "expected at least 4 of first 5 decoded as state 0, got {state0_first_half}"
        );
        assert!(
            state1_second_half >= 4,
            "expected at least 4 of last 5 decoded as state 1, got {state1_second_half}"
        );
    }
}
