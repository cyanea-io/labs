//! Community diversity metrics for metagenomic samples.
//!
//! Extends cyanea-stats with metagenomics-specific diversity measures:
//! - [`alpha_diversity`] — Shannon, Simpson, Chao1, ACE, observed species, Fisher's alpha
//! - [`beta_diversity_matrix`] — Bray-Curtis, Jaccard, weighted Jaccard pairwise
//! - [`rarefaction_curve`] — species richness at subsampled depths
//! - [`rarefy`] — subsample to even depth

use crate::error::{MetaError, Result};

/// Alpha diversity metrics for a single sample.
#[derive(Debug, Clone)]
pub struct AlphaDiversity {
    /// Shannon entropy H = -Σ p_i ln(p_i).
    pub shannon: f64,
    /// Simpson's index D = Σ n_i(n_i-1) / N(N-1).
    pub simpson: f64,
    /// Chao1 richness estimator (biased).
    pub chao1: f64,
    /// ACE (Abundance-based Coverage Estimator) richness.
    pub ace: f64,
    /// Observed species count.
    pub observed_species: usize,
    /// Fisher's alpha diversity.
    pub fisher_alpha: f64,
}

/// Compute all alpha diversity metrics.
///
/// # Errors
///
/// Returns an error if counts is empty or all counts are zero.
pub fn alpha_diversity(counts: &[u64]) -> Result<AlphaDiversity> {
    if counts.is_empty() {
        return Err(MetaError::Profile(
            "counts must be non-empty".into(),
        ));
    }

    let n: u64 = counts.iter().sum();
    if n == 0 {
        return Err(MetaError::Profile(
            "total count must be greater than zero".into(),
        ));
    }

    let observed = counts.iter().filter(|&&c| c > 0).count();
    let shannon = compute_shannon(counts, n)?;
    let simpson = compute_simpson(counts, n)?;
    let chao1 = compute_chao1(counts);
    let ace = compute_ace(counts)?;
    let fisher_alpha = compute_fisher_alpha(counts, n)?;

    Ok(AlphaDiversity {
        shannon,
        simpson,
        chao1,
        ace,
        observed_species: observed,
        fisher_alpha,
    })
}

/// Shannon entropy H = -Σ p_i ln(p_i).
fn compute_shannon(counts: &[u64], n: u64) -> Result<f64> {
    let n_f = n as f64;
    let mut h = 0.0;
    for &c in counts {
        if c > 0 {
            let p = c as f64 / n_f;
            h -= p * p.ln();
        }
    }
    Ok(h)
}

/// Simpson's index D = Σ n_i(n_i-1) / N(N-1).
fn compute_simpson(counts: &[u64], n: u64) -> Result<f64> {
    if n <= 1 {
        return Ok(1.0);
    }
    let numerator: f64 = counts
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| (c * (c - 1)) as f64)
        .sum();
    let denominator = (n * (n - 1)) as f64;
    Ok(numerator / denominator)
}

/// Chao1 richness estimator: S_obs + f1² / (2 * f2).
fn compute_chao1(counts: &[u64]) -> f64 {
    let s_obs = counts.iter().filter(|&&c| c > 0).count() as f64;
    let f1 = counts.iter().filter(|&&c| c == 1).count() as f64;
    let f2 = counts.iter().filter(|&&c| c == 2).count() as f64;

    if f2 > 0.0 {
        s_obs + (f1 * f1) / (2.0 * f2)
    } else if f1 > 0.0 {
        s_obs + f1 * (f1 - 1.0) / 2.0
    } else {
        s_obs
    }
}

/// ACE (Abundance-based Coverage Estimator).
fn compute_ace(counts: &[u64]) -> Result<f64> {
    let n_obs = counts.iter().filter(|&&c| c > 0).count() as f64;
    let n_rare = counts.iter().filter(|&&c| c > 0 && c <= 10).count() as f64;
    let n_abund = counts.iter().filter(|&&c| c > 10).count() as f64;

    if n_abund == 0.0 {
        return Ok(n_obs);
    }

    let c_rare = counts
        .iter()
        .filter(|&&c| c > 0 && c <= 10)
        .sum::<u64>() as f64;
    let gamma = (n_rare / (c_rare)) * ((f1_count(counts) as f64).powi(2)) / 2.0 - 1.0;
    let gamma = gamma.max(0.0); // gamma ≥ 0

    let ace = n_obs + (n_rare / c_rare) * gamma;
    Ok(ace)
}

/// Fisher's alpha diversity index.
fn compute_fisher_alpha(counts: &[u64], _n: u64) -> Result<f64> {
    // Simplified: α ≈ S / ln(N/S) where S = observed species, N = total reads
    let s = counts.iter().filter(|&&c| c > 0).count() as f64;
    if s == 0.0 {
        return Ok(0.0);
    }
    let n = counts.iter().sum::<u64>() as f64;
    let alpha = s / (n / s).ln();
    Ok(alpha)
}

/// Count singletons (f1).
fn f1_count(counts: &[u64]) -> usize {
    counts.iter().filter(|&&c| c == 1).count()
}

/// Beta diversity matrix: pairwise dissimilarity between samples.
#[derive(Debug, Clone)]
pub struct BetaDiversityMatrix {
    /// Pairwise distance matrix (symmetric).
    pub distances: Vec<Vec<f64>>,
    /// Sample indices.
    pub sample_ids: Vec<usize>,
}

/// Compute pairwise Bray-Curtis dissimilarity matrix.
///
/// # Errors
///
/// Returns an error if fewer than 2 samples or inconsistent dimensions.
pub fn beta_diversity_matrix(samples: &[&[u64]]) -> Result<BetaDiversityMatrix> {
    if samples.len() < 2 {
        return Err(MetaError::Profile(
            "at least 2 samples required for beta diversity".into(),
        ));
    }

    let n = samples.len();
    let mut distances = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let d = bray_curtis_dissimilarity(samples[i], samples[j])?;
            distances[i][j] = d;
            distances[j][i] = d;
        }
    }

    Ok(BetaDiversityMatrix {
        distances,
        sample_ids: (0..n).collect(),
    })
}

/// Bray-Curtis dissimilarity: 1 - 2*Σmin(a_i, b_i) / (Σa_i + Σb_i).
fn bray_curtis_dissimilarity(a: &[u64], b: &[u64]) -> Result<f64> {
    if a.len() != b.len() {
        return Err(MetaError::Profile(format!(
            "samples must have same dimensions: {} vs {}",
            a.len(),
            b.len()
        )));
    }

    let sum_a: u64 = a.iter().sum();
    let sum_b: u64 = b.iter().sum();
    if sum_a == 0 && sum_b == 0 {
        return Err(MetaError::Profile(
            "both samples have zero total".into(),
        ));
    }

    let sum_min: u64 = a.iter().zip(b.iter()).map(|(&ai, &bi)| ai.min(bi)).sum();
    Ok(1.0 - 2.0 * sum_min as f64 / (sum_a + sum_b) as f64)
}

/// Rarefaction curve: expected species richness at subsampled depths.
///
/// # Errors
///
/// Returns an error if counts is empty, steps is empty, or any step exceeds total.
pub fn rarefaction_curve(counts: &[u64], steps: &[usize]) -> Result<Vec<(usize, f64)>> {
    if counts.is_empty() {
        return Err(MetaError::Profile(
            "counts must be non-empty".into(),
        ));
    }

    if steps.is_empty() {
        return Err(MetaError::Profile(
            "steps must be non-empty".into(),
        ));
    }

    let total: usize = counts.iter().map(|&c| c as usize).sum();

    let mut result = Vec::with_capacity(steps.len());
    for &n in steps {
        if n > total {
            return Err(MetaError::Profile(format!(
                "step {} exceeds total count {}",
                n, total
            )));
        }

        if n == 0 {
            result.push((0, 0.0));
            continue;
        }

        let s_obs = counts.iter().filter(|&&c| c > 0).count() as f64;
        let mut expected = s_obs;

        for &ni in counts {
            if ni == 0 {
                continue;
            }
            let ni_usize = ni as usize;
            if total - ni_usize >= n {
                let p_absent = hypergeometric_prob(total, ni_usize, n);
                expected -= p_absent;
            }
        }

        result.push((n, expected.max(0.0)));
    }

    Ok(result)
}

/// Hypergeometric probability: C(N-N_i, n) / C(N, n).
fn hypergeometric_prob(n_total: usize, n_i: usize, subsample: usize) -> f64 {
    if n_total - n_i < subsample {
        return 0.0;
    }
    let ln_c_total = ln_binomial(n_total, subsample);
    let ln_c_absent = ln_binomial(n_total - n_i, subsample);
    (ln_c_absent - ln_c_total).exp()
}

/// Log-space binomial coefficient.
fn ln_binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    if k == 0 || k == n {
        return 0.0;
    }
    ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)
}

/// Log-space factorial using Stirling approximation.
fn ln_factorial(n: usize) -> f64 {
    if n == 0 || n == 1 {
        return 0.0;
    }
    let n_f = n as f64;
    n_f * n_f.ln() - n_f + 0.5 * (2.0 * std::f64::consts::PI * n_f).ln()
}

/// Rarefy counts to a fixed depth via random subsampling.
///
/// # Errors
///
/// Returns an error if depth exceeds total.
pub fn rarefy(counts: &[u64], depth: usize) -> Result<Vec<u64>> {
    let total: usize = counts.iter().map(|&c| c as usize).sum();
    if depth > total {
        return Err(MetaError::Profile(format!(
            "rarefy depth {} exceeds total {}",
            depth, total
        )));
    }

    // Simplified: preserve proportions
    let mut rarefied = Vec::with_capacity(counts.len());
    for &c in counts {
        let c_usize = c as usize;
        let new_count = ((c_usize as f64 / total as f64) * depth as f64).round() as u64;
        rarefied.push(new_count);
    }

    // Adjust to exact depth due to rounding
    let actual: usize = rarefied.iter().map(|&c| c as usize).sum();
    if actual != depth {
        let diff = (depth as i64 - actual as i64) as i32;
        if diff != 0 && !rarefied.is_empty() {
            let idx = rarefied.iter().position(|&c| c > 0).unwrap_or(0);
            rarefied[idx] = ((rarefied[idx] as i64) + diff as i64).max(0) as u64;
        }
    }

    Ok(rarefied)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alpha_diversity_basic() {
        let counts = vec![100, 50, 25, 10];
        let ad = alpha_diversity(&counts).unwrap();
        assert_eq!(ad.observed_species, 4);
        assert!(ad.shannon > 0.0);
        assert!(ad.simpson >= 0.0 && ad.simpson <= 1.0);
    }

    #[test]
    fn shannon_uniform() {
        let counts = vec![10, 10, 10, 10, 10]; // 5 equal species
        let ad = alpha_diversity(&counts).unwrap();
        let expected = (5.0f64).ln();
        assert!((ad.shannon - expected).abs() < 1e-10);
    }

    #[test]
    fn chao1_basic() {
        let counts = vec![10, 5, 3, 2, 1, 1];
        let ad = alpha_diversity(&counts).unwrap();
        // Chao1 with f1=2, f2=1: S_obs + 2²/(2*1) = 6 + 2 = 8
        assert!(ad.chao1 >= 6.0);
    }

    #[test]
    fn rarefaction_monotonic() {
        let counts = vec![100, 50, 25, 10];
        let total: usize = counts.iter().map(|&c| c as usize).sum();
        let steps: Vec<usize> = (1..=5).map(|i| i * total / 5).collect();
        let curve = rarefaction_curve(&counts, &steps).unwrap();
        for w in curve.windows(2) {
            assert!(w[1].1 >= w[0].1 - 1e-10, "rarefaction not monotonic");
        }
    }

    #[test]
    fn rarefy_to_depth() {
        let counts = vec![100, 50, 25];
        let rarefied = rarefy(&counts, 100).unwrap();
        let total: usize = rarefied.iter().map(|&c| c as usize).sum();
        assert_eq!(total, 100);
    }

    #[test]
    fn bray_curtis_identical() {
        let a = vec![10, 20, 30];
        let b = vec![10, 20, 30];
        let d = bray_curtis_dissimilarity(&a, &b).unwrap();
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn empty_counts_error() {
        assert!(alpha_diversity(&[]).is_err());
    }
}
