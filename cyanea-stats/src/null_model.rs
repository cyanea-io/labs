//! Wright-Fisher simulation and null model generators.
//!
//! Provides population genetics simulation via the Wright-Fisher model and
//! resampling-based null model generators for hypothesis testing:
//!
//! - [`wright_fisher`] — discrete-generation drift simulation with binomial sampling
//! - [`permutation_null`] — permutation-based null distribution for any group statistic
//! - [`bootstrap_null`] — bootstrap resampling null distribution for any statistic

use cyanea_core::{CyaneaError, Result};

// ── Xorshift64 PRNG ────────────────────────────────────────────────────────

/// Minimal xorshift64 PRNG for reproducible simulations without external deps.
struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 1 } else { seed },
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state ^= self.state << 13;
        self.state ^= self.state >> 7;
        self.state ^= self.state << 17;
        self.state
    }

    fn next_f64(&mut self) -> f64 {
        self.next_u64() as f64 / u64::MAX as f64
    }
}

// ── Binomial sampling ──────────────────────────────────────────────────────

/// Sample from Binomial(n, p) using direct trials for small n, or a normal
/// approximation (Box-Muller) for large n (> 100).
fn sample_binomial(rng: &mut Xorshift64, n: usize, p: f64) -> usize {
    if p <= 0.0 {
        return 0;
    }
    if p >= 1.0 {
        return n;
    }

    if n <= 100 {
        // Direct sampling: count successes in n Bernoulli trials.
        let mut successes = 0usize;
        for _ in 0..n {
            if rng.next_f64() < p {
                successes += 1;
            }
        }
        successes
    } else {
        // Normal approximation: N(np, sqrt(np(1-p))), rounded and clamped.
        let np = n as f64 * p;
        let std = (np * (1.0 - p)).sqrt();

        // Box-Muller transform for a standard normal variate.
        let u1 = rng.next_f64().max(1e-300); // avoid log(0)
        let u2 = rng.next_f64();
        let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();

        let value = (np + std * z).round();
        // Clamp to [0, n].
        value.max(0.0).min(n as f64) as usize
    }
}

// ── Wright-Fisher result ───────────────────────────────────────────────────

/// Result of a Wright-Fisher drift simulation.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct WrightFisherResult {
    /// Allele frequency at each generation (length = `n_generations + 1`,
    /// including the initial frequency at generation 0).
    pub frequencies: Vec<f64>,
    /// Generation at which fixation occurred (freq == 0.0 or freq == 1.0),
    /// or `None` if the allele is still segregating.
    pub fixation_gen: Option<usize>,
    /// Allele frequency at the final generation.
    pub final_freq: f64,
}

// ── Wright-Fisher simulation ───────────────────────────────────────────────

/// Simulate allele frequency drift under the Wright-Fisher model.
///
/// Each generation, the number of copies of the focal allele is drawn from
/// `Binomial(2 * pop_size, current_freq)`. The simulation runs for
/// `n_generations` discrete generations.
///
/// # Errors
///
/// Returns an error if `pop_size` is zero or `initial_freq` is outside
/// `[0.0, 1.0]`.
///
/// # Example
///
/// ```
/// use cyanea_stats::null_model::wright_fisher;
///
/// let result = wright_fisher(100, 0.5, 200, 42).unwrap();
/// assert_eq!(result.frequencies.len(), 201); // gen 0..=200
/// ```
pub fn wright_fisher(
    pop_size: usize,
    initial_freq: f64,
    n_generations: usize,
    seed: u64,
) -> Result<WrightFisherResult> {
    if pop_size == 0 {
        return Err(CyaneaError::InvalidInput(
            "pop_size must be > 0".to_string(),
        ));
    }
    if !(0.0..=1.0).contains(&initial_freq) {
        return Err(CyaneaError::InvalidInput(
            "initial_freq must be in [0.0, 1.0]".to_string(),
        ));
    }

    let two_n = 2 * pop_size;
    let mut rng = Xorshift64::new(seed);
    let mut frequencies = Vec::with_capacity(n_generations + 1);
    let mut freq = initial_freq;
    let mut fixation_gen: Option<usize> = None;

    frequencies.push(freq);

    for gen in 1..=n_generations {
        if freq == 0.0 || freq == 1.0 {
            // Already fixed — frequency stays constant.
            frequencies.push(freq);
            if fixation_gen.is_none() {
                fixation_gen = Some(gen - 1);
            }
            continue;
        }

        let copies = sample_binomial(&mut rng, two_n, freq);
        freq = copies as f64 / two_n as f64;
        frequencies.push(freq);

        if fixation_gen.is_none() && (freq == 0.0 || freq == 1.0) {
            fixation_gen = Some(gen);
        }
    }

    let final_freq = *frequencies.last().unwrap();

    Ok(WrightFisherResult {
        frequencies,
        fixation_gen,
        final_freq,
    })
}

// ── Permutation null ───────────────────────────────────────────────────────

/// Generate a permutation null distribution for a group statistic.
///
/// On each iteration the pooled `values` are shuffled (Fisher-Yates) and split
/// into groups according to `group_sizes`. The user-supplied `statistic`
/// function is evaluated on each permuted grouping, and the resulting
/// `n_permutations` statistic values are returned.
///
/// # Errors
///
/// Returns an error if `group_sizes` is empty or its sum does not equal
/// `values.len()`.
///
/// # Example
///
/// ```
/// use cyanea_stats::null_model::permutation_null;
///
/// let values = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
/// let sizes = vec![3, 3];
/// let diff_means = |groups: &[&[f64]]| {
///     let m0: f64 = groups[0].iter().sum::<f64>() / groups[0].len() as f64;
///     let m1: f64 = groups[1].iter().sum::<f64>() / groups[1].len() as f64;
///     (m0 - m1).abs()
/// };
/// let null = permutation_null(&values, &sizes, &diff_means, 100, 42).unwrap();
/// assert_eq!(null.len(), 100);
/// ```
pub fn permutation_null(
    values: &[f64],
    group_sizes: &[usize],
    statistic: &dyn Fn(&[&[f64]]) -> f64,
    n_permutations: usize,
    seed: u64,
) -> Result<Vec<f64>> {
    if group_sizes.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "group_sizes must be non-empty".to_string(),
        ));
    }
    let total: usize = group_sizes.iter().sum();
    if total != values.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "sum of group_sizes ({}) != values.len() ({})",
            total,
            values.len()
        )));
    }

    let mut rng = Xorshift64::new(seed);
    let mut buf = values.to_vec();
    let mut results = Vec::with_capacity(n_permutations);

    for _ in 0..n_permutations {
        // Fisher-Yates shuffle.
        let n = buf.len();
        for i in (1..n).rev() {
            let j = (rng.next_u64() as usize) % (i + 1);
            buf.swap(i, j);
        }

        // Split into groups and compute statistic.
        let mut groups: Vec<&[f64]> = Vec::with_capacity(group_sizes.len());
        let mut offset = 0;
        for &sz in group_sizes {
            groups.push(&buf[offset..offset + sz]);
            offset += sz;
        }

        results.push(statistic(&groups));
    }

    Ok(results)
}

// ── Bootstrap null ─────────────────────────────────────────────────────────

/// Generate a bootstrap null distribution for a statistic.
///
/// On each iteration, `data.len()` values are resampled with replacement from
/// `data`, and the user-supplied `statistic` function is evaluated on the
/// resample. The resulting `n_bootstrap` statistic values are returned.
///
/// # Example
///
/// ```
/// use cyanea_stats::null_model::bootstrap_null;
///
/// let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
/// let mean_fn = |d: &[f64]| d.iter().sum::<f64>() / d.len() as f64;
/// let null = bootstrap_null(&data, &mean_fn, 500, 42);
/// assert_eq!(null.len(), 500);
/// ```
pub fn bootstrap_null(
    data: &[f64],
    statistic: &dyn Fn(&[f64]) -> f64,
    n_bootstrap: usize,
    seed: u64,
) -> Vec<f64> {
    let n = data.len();
    if n == 0 {
        return vec![f64::NAN; n_bootstrap];
    }

    let mut rng = Xorshift64::new(seed);
    let mut resample = vec![0.0; n];
    let mut results = Vec::with_capacity(n_bootstrap);

    for _ in 0..n_bootstrap {
        for slot in resample.iter_mut() {
            let idx = (rng.next_u64() as usize) % n;
            *slot = data[idx];
        }
        results.push(statistic(&resample));
    }

    results
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wf_zero_freq_stays_zero() {
        let result = wright_fisher(100, 0.0, 50, 123).unwrap();
        assert_eq!(result.frequencies.len(), 51);
        for &f in &result.frequencies {
            assert_eq!(f, 0.0);
        }
        assert_eq!(result.final_freq, 0.0);
        // Fixation at generation 0 since it starts fixed at loss.
        assert_eq!(result.fixation_gen, Some(0));
    }

    #[test]
    fn wf_one_freq_stays_one() {
        let result = wright_fisher(100, 1.0, 50, 456).unwrap();
        assert_eq!(result.frequencies.len(), 51);
        for &f in &result.frequencies {
            assert_eq!(f, 1.0);
        }
        assert_eq!(result.final_freq, 1.0);
        assert_eq!(result.fixation_gen, Some(0));
    }

    #[test]
    fn permutation_null_returns_correct_count() {
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let sizes = vec![3, 3];
        let diff_means = |groups: &[&[f64]]| {
            let m0: f64 = groups[0].iter().sum::<f64>() / groups[0].len() as f64;
            let m1: f64 = groups[1].iter().sum::<f64>() / groups[1].len() as f64;
            (m0 - m1).abs()
        };
        let null = permutation_null(&values, &sizes, &diff_means, 200, 789).unwrap();
        assert_eq!(null.len(), 200);
        // All values should be finite.
        for &v in &null {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn bootstrap_null_returns_correct_count() {
        let data = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        let mean_fn = |d: &[f64]| d.iter().sum::<f64>() / d.len() as f64;
        let null = bootstrap_null(&data, &mean_fn, 300, 101);
        assert_eq!(null.len(), 300);
        // Bootstrap means must be within the range of the data.
        for &v in &null {
            assert!(v >= 2.0 && v <= 10.0, "bootstrap mean {} out of range", v);
        }
    }
}
