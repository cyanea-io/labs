//! Alpha and beta diversity metrics for microbial ecology.
//!
//! - **Alpha diversity** — Shannon, Simpson, inverse Simpson, Chao1, observed species
//! - **Beta diversity** — Bray-Curtis dissimilarity and distance matrices
//! - **Rarefaction** — expected species richness at subsampled depths

use cyanea_core::{CyaneaError, Result};

/// Alpha diversity summary for a single sample.
#[derive(Debug, Clone)]
pub struct AlphaDiversity {
    /// Shannon entropy H = -Σ p_i ln(p_i).
    pub shannon: f64,
    /// Simpson's index D = Σ n_i(n_i-1) / N(N-1).
    pub simpson: f64,
    /// Inverse Simpson 1/D.
    pub inverse_simpson: f64,
    /// Chao1 richness estimator.
    pub chao1: f64,
    /// Number of observed species (non-zero counts).
    pub observed_species: usize,
}

/// Compute all alpha diversity metrics for a vector of species counts.
///
/// # Errors
///
/// Returns an error if `counts` is empty or all counts are zero.
pub fn alpha_diversity(counts: &[usize]) -> Result<AlphaDiversity> {
    if counts.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "counts must be non-empty".into(),
        ));
    }
    let n: usize = counts.iter().sum();
    if n == 0 {
        return Err(CyaneaError::InvalidInput(
            "total count must be greater than zero".into(),
        ));
    }

    Ok(AlphaDiversity {
        shannon: shannon_index(counts)?,
        simpson: simpson_index(counts)?,
        inverse_simpson: {
            let d = simpson_index(counts)?;
            if d > 0.0 { 1.0 / d } else { f64::INFINITY }
        },
        chao1: chao1(counts)?,
        observed_species: counts.iter().filter(|&&c| c > 0).count(),
    })
}

/// Shannon diversity index H = -Σ p_i ln(p_i).
///
/// Uses natural logarithm.
///
/// # Errors
///
/// Returns an error if `counts` is empty or total is zero.
pub fn shannon_index(counts: &[usize]) -> Result<f64> {
    validate_counts(counts)?;
    let n: f64 = counts.iter().sum::<usize>() as f64;
    let mut h = 0.0;
    for &c in counts {
        if c > 0 {
            let p = c as f64 / n;
            h -= p * p.ln();
        }
    }
    Ok(h)
}

/// Simpson's diversity index D = Σ n_i(n_i-1) / N(N-1).
///
/// Values range from 0 (infinite diversity) to 1 (single species).
///
/// # Errors
///
/// Returns an error if `counts` is empty or total ≤ 1.
pub fn simpson_index(counts: &[usize]) -> Result<f64> {
    validate_counts(counts)?;
    let n: usize = counts.iter().sum();
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

/// Chao1 richness estimator.
///
/// S_chao1 = S_obs + f1² / (2 * f2), where f1 = singletons, f2 = doubletons.
/// If f2 = 0, uses f1(f1-1)/2 instead.
///
/// # Errors
///
/// Returns an error if `counts` is empty or total is zero.
pub fn chao1(counts: &[usize]) -> Result<f64> {
    validate_counts(counts)?;
    let s_obs = counts.iter().filter(|&&c| c > 0).count() as f64;
    let f1 = counts.iter().filter(|&&c| c == 1).count() as f64;
    let f2 = counts.iter().filter(|&&c| c == 2).count() as f64;

    let estimate = if f2 > 0.0 {
        s_obs + (f1 * f1) / (2.0 * f2)
    } else if f1 > 0.0 {
        s_obs + f1 * (f1 - 1.0) / 2.0
    } else {
        s_obs
    };
    Ok(estimate)
}

/// Bray-Curtis dissimilarity between two samples.
///
/// BC = 1 - 2 * Σ min(a_i, b_i) / (Σ a_i + Σ b_i).
/// Ranges from 0 (identical) to 1 (completely different).
///
/// # Errors
///
/// Returns an error if vectors have different lengths or both are all-zero.
pub fn bray_curtis(a: &[usize], b: &[usize]) -> Result<f64> {
    if a.len() != b.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "count vectors must have the same length: {} vs {}",
            a.len(),
            b.len()
        )));
    }
    let sum_a: usize = a.iter().sum();
    let sum_b: usize = b.iter().sum();
    if sum_a == 0 && sum_b == 0 {
        return Err(CyaneaError::InvalidInput(
            "both samples have zero total counts".into(),
        ));
    }
    let sum_min: usize = a.iter().zip(b.iter()).map(|(&ai, &bi)| ai.min(bi)).sum();
    Ok(1.0 - 2.0 * sum_min as f64 / (sum_a + sum_b) as f64)
}

/// Compute a pairwise Bray-Curtis distance matrix for multiple samples.
///
/// Returns a symmetric matrix where `result[i][j]` is the Bray-Curtis
/// dissimilarity between samples `i` and `j`.
///
/// # Errors
///
/// Returns an error if samples have different lengths or fewer than 2 samples.
pub fn bray_curtis_matrix(samples: &[&[usize]]) -> Result<Vec<Vec<f64>>> {
    if samples.len() < 2 {
        return Err(CyaneaError::InvalidInput(
            "at least 2 samples are required for a distance matrix".into(),
        ));
    }
    let n = samples.len();
    let mut matrix = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = bray_curtis(samples[i], samples[j])?;
            matrix[i][j] = d;
            matrix[j][i] = d;
        }
    }
    Ok(matrix)
}

/// Jaccard dissimilarity between two samples (presence/absence).
///
/// `J = 1 - |A∩B| / |A∪B|` where A and B are the sets of species present.
/// Ranges from 0 (identical) to 1 (completely different).
///
/// # Errors
///
/// Returns an error if vectors have different lengths or both are all-zero.
pub fn jaccard(a: &[usize], b: &[usize]) -> Result<f64> {
    if a.len() != b.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "count vectors must have the same length: {} vs {}",
            a.len(),
            b.len()
        )));
    }
    let mut intersection = 0usize;
    let mut union = 0usize;
    for (&ai, &bi) in a.iter().zip(b.iter()) {
        let pa = ai > 0;
        let pb = bi > 0;
        if pa || pb {
            union += 1;
        }
        if pa && pb {
            intersection += 1;
        }
    }
    if union == 0 {
        return Err(CyaneaError::InvalidInput(
            "both samples have zero total counts".into(),
        ));
    }
    Ok(1.0 - intersection as f64 / union as f64)
}

/// Weighted Jaccard dissimilarity (Ruzicka distance) between two samples.
///
/// `WJ = 1 - Σ min(a_i, b_i) / Σ max(a_i, b_i)`.
/// Ranges from 0 (identical) to 1 (completely different).
///
/// # Errors
///
/// Returns an error if vectors have different lengths or both are all-zero.
pub fn weighted_jaccard(a: &[usize], b: &[usize]) -> Result<f64> {
    if a.len() != b.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "count vectors must have the same length: {} vs {}",
            a.len(),
            b.len()
        )));
    }
    let sum_min: usize = a.iter().zip(b.iter()).map(|(&ai, &bi)| ai.min(bi)).sum();
    let sum_max: usize = a.iter().zip(b.iter()).map(|(&ai, &bi)| ai.max(bi)).sum();
    if sum_max == 0 {
        return Err(CyaneaError::InvalidInput(
            "both samples have zero total counts".into(),
        ));
    }
    Ok(1.0 - sum_min as f64 / sum_max as f64)
}

/// Compute a pairwise Jaccard distance matrix for multiple samples.
///
/// Returns a symmetric matrix where `result[i][j]` is the Jaccard
/// dissimilarity between samples `i` and `j`.
///
/// # Errors
///
/// Returns an error if samples have different lengths or fewer than 2 samples.
pub fn jaccard_matrix(samples: &[&[usize]]) -> Result<Vec<Vec<f64>>> {
    if samples.len() < 2 {
        return Err(CyaneaError::InvalidInput(
            "at least 2 samples are required for a distance matrix".into(),
        ));
    }
    let n = samples.len();
    let mut matrix = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = jaccard(samples[i], samples[j])?;
            matrix[i][j] = d;
            matrix[j][i] = d;
        }
    }
    Ok(matrix)
}

/// Hill numbers (effective number of species) for a range of orders.
///
/// For each order `q`:
/// - q = 0: species richness
/// - q → 1: exp(Shannon entropy)
/// - q = 2: inverse Simpson index
/// - General: `(Σ p_i^q)^(1/(1-q))`
///
/// Hill numbers are monotone non-increasing in q.
///
/// # Errors
///
/// Returns an error if `counts` is empty or all counts are zero.
pub fn hill_numbers(counts: &[usize], orders: &[f64]) -> Result<Vec<(f64, f64)>> {
    validate_counts(counts)?;
    let n: f64 = counts.iter().sum::<usize>() as f64;
    let proportions: Vec<f64> = counts.iter().filter(|&&c| c > 0).map(|&c| c as f64 / n).collect();

    let mut result = Vec::with_capacity(orders.len());
    for &q in orders {
        let hill = if (q - 1.0).abs() < 1e-12 {
            // q → 1: exp(Shannon) = exp(-Σ p_i ln(p_i))
            let h: f64 = proportions.iter().map(|&p| -p * p.ln()).sum();
            h.exp()
        } else if q == 0.0 {
            // q = 0: species richness
            proportions.len() as f64
        } else {
            // General: (Σ p_i^q)^(1/(1-q))
            let sum_pq: f64 = proportions.iter().map(|&p| p.powf(q)).sum();
            sum_pq.powf(1.0 / (1.0 - q))
        };
        result.push((q, hill));
    }
    Ok(result)
}

/// Batch alpha rarefaction: compute rarefaction curves for multiple samples.
///
/// For each sample, computes the expected species richness at each depth
/// using [`rarefaction_curve`].
///
/// # Errors
///
/// Returns an error if any sample is invalid or any depth exceeds a sample's
/// total count.
pub fn alpha_rarefaction(
    samples: &[&[usize]],
    depths: &[usize],
) -> Result<Vec<Vec<(usize, f64)>>> {
    let mut result = Vec::with_capacity(samples.len());
    for (i, sample) in samples.iter().enumerate() {
        let total: usize = sample.iter().sum();
        // Filter depths to those not exceeding this sample's total
        let valid_depths: Vec<usize> = depths.iter().copied().filter(|&d| d <= total).collect();
        if valid_depths.is_empty() {
            return Err(CyaneaError::InvalidInput(format!(
                "no valid depths for sample {} (total count {})",
                i, total
            )));
        }
        result.push(rarefaction_curve(sample, &valid_depths)?);
    }
    Ok(result)
}

/// Compute a rarefaction curve: expected species richness at subsampled depths.
///
/// For each depth `n` in `steps`, estimates the expected number of species
/// using the hypergeometric formula:
/// E[S_n] = S - Σ C(N - N_i, n) / C(N, n)
///
/// Uses log-space binomial coefficients for numerical stability with large N.
///
/// # Errors
///
/// Returns an error if `counts` is empty or any step exceeds the total count.
pub fn rarefaction_curve(counts: &[usize], steps: &[usize]) -> Result<Vec<(usize, f64)>> {
    validate_counts(counts)?;
    let total_n: usize = counts.iter().sum();

    let mut result = Vec::with_capacity(steps.len());
    for &n in steps {
        if n > total_n {
            return Err(CyaneaError::InvalidInput(format!(
                "rarefaction step {} exceeds total count {}",
                n, total_n
            )));
        }
        if n == 0 {
            result.push((0, 0.0));
            continue;
        }

        let s_obs = counts.iter().filter(|&&c| c > 0).count();
        let ln_c_total = ln_binomial(total_n, n);

        let mut expected = s_obs as f64;
        for &ni in counts {
            if ni == 0 {
                continue;
            }
            // Probability that species i is absent from subsample of size n:
            // C(N - N_i, n) / C(N, n)
            if total_n - ni >= n {
                let ln_c_absent = ln_binomial(total_n - ni, n);
                expected -= (ln_c_absent - ln_c_total).exp();
            }
            // If total_n - ni < n, species i is always present.
        }

        result.push((n, expected));
    }
    Ok(result)
}

fn validate_counts(counts: &[usize]) -> Result<()> {
    if counts.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "counts must be non-empty".into(),
        ));
    }
    let n: usize = counts.iter().sum();
    if n == 0 {
        return Err(CyaneaError::InvalidInput(
            "total count must be greater than zero".into(),
        ));
    }
    Ok(())
}

/// Log-space binomial coefficient ln(C(n, k)).
fn ln_binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    if k == 0 || k == n {
        return 0.0;
    }
    // Use Stirling via ln_gamma: ln(C(n,k)) = ln_gamma(n+1) - ln_gamma(k+1) - ln_gamma(n-k+1)
    ln_gamma(n as f64 + 1.0) - ln_gamma(k as f64 + 1.0) - ln_gamma((n - k) as f64 + 1.0)
}

/// Lanczos approximation of ln(Γ(x)).
fn ln_gamma(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::INFINITY;
    }
    // Use the standard Lanczos coefficients (g=7, n=9).
    let coeffs = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];
    if x < 0.5 {
        let pi = std::f64::consts::PI;
        (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut sum = coeffs[0];
        for (i, &c) in coeffs[1..].iter().enumerate() {
            sum += c / (x + i as f64 + 1.0);
        }
        let t = x + 7.5;
        0.5 * (2.0 * std::f64::consts::PI).ln() + (t.ln() * (x + 0.5)) - t + sum.ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shannon_uniform_is_ln_s() {
        // Uniform distribution with S species: H = ln(S).
        let counts = vec![10, 10, 10, 10, 10]; // 5 species, equal abundance
        let h = shannon_index(&counts).unwrap();
        let expected = (5.0f64).ln();
        assert!((h - expected).abs() < 1e-10, "H={}, expected={}", h, expected);
    }

    #[test]
    fn simpson_single_species_is_one() {
        // Single species: D = n(n-1)/n(n-1) = 1.
        let counts = vec![100];
        let d = simpson_index(&counts).unwrap();
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn chao1_all_abundant() {
        // No singletons or doubletons: Chao1 = S_obs.
        let counts = vec![10, 20, 30, 40];
        let c = chao1(&counts).unwrap();
        assert!((c - 4.0).abs() < 1e-10);
    }

    #[test]
    fn bray_curtis_identical_is_zero() {
        let a = vec![10, 20, 30];
        let b = vec![10, 20, 30];
        let bc = bray_curtis(&a, &b).unwrap();
        assert!(bc.abs() < 1e-10);
    }

    #[test]
    fn bray_curtis_disjoint_is_one() {
        let a = vec![10, 0, 0];
        let b = vec![0, 20, 30];
        let bc = bray_curtis(&a, &b).unwrap();
        assert!((bc - 1.0).abs() < 1e-10);
    }

    // ── Jaccard tests ──────────────────────────────────────────────

    #[test]
    fn jaccard_identical_is_zero() {
        let a = vec![10, 20, 30];
        let b = vec![5, 15, 25];
        // Same presence/absence pattern → J = 0
        let j = jaccard(&a, &b).unwrap();
        assert!(j.abs() < 1e-10, "J={}", j);
    }

    #[test]
    fn jaccard_disjoint_is_one() {
        let a = vec![10, 0, 0];
        let b = vec![0, 20, 30];
        let j = jaccard(&a, &b).unwrap();
        assert!((j - 1.0).abs() < 1e-10, "J={}", j);
    }

    #[test]
    fn jaccard_symmetric() {
        let a = vec![10, 20, 0, 5];
        let b = vec![0, 15, 30, 0];
        let j1 = jaccard(&a, &b).unwrap();
        let j2 = jaccard(&b, &a).unwrap();
        assert!((j1 - j2).abs() < 1e-10);
    }

    #[test]
    fn weighted_jaccard_identical_is_zero() {
        let a = vec![10, 20, 30];
        let wj = weighted_jaccard(&a, &a).unwrap();
        assert!(wj.abs() < 1e-10, "WJ={}", wj);
    }

    #[test]
    fn weighted_jaccard_disjoint_is_one() {
        let a = vec![10, 0, 0];
        let b = vec![0, 20, 30];
        let wj = weighted_jaccard(&a, &b).unwrap();
        assert!((wj - 1.0).abs() < 1e-10, "WJ={}", wj);
    }

    #[test]
    fn jaccard_matrix_symmetric() {
        let s1 = vec![10, 20, 0];
        let s2 = vec![0, 15, 30];
        let s3 = vec![5, 10, 15];
        let mat = jaccard_matrix(&[&s1, &s2, &s3]).unwrap();
        assert_eq!(mat.len(), 3);
        for i in 0..3 {
            assert!(mat[i][i].abs() < 1e-10);
            for j in 0..3 {
                assert!((mat[i][j] - mat[j][i]).abs() < 1e-10);
            }
        }
    }

    // ── Hill numbers tests ──────────────────────────────────────────

    #[test]
    fn hill_q0_is_richness() {
        let counts = vec![10, 20, 30, 40, 0];
        let result = hill_numbers(&counts, &[0.0]).unwrap();
        assert!((result[0].1 - 4.0).abs() < 1e-10, "q=0: {}", result[0].1);
    }

    #[test]
    fn hill_q1_is_exp_shannon() {
        let counts = vec![10, 20, 30, 40];
        let result = hill_numbers(&counts, &[1.0]).unwrap();
        let h = shannon_index(&counts).unwrap();
        let expected = h.exp();
        assert!(
            (result[0].1 - expected).abs() < 1e-8,
            "q=1: {} expected {}",
            result[0].1,
            expected
        );
    }

    #[test]
    fn hill_q2_is_inv_simpson() {
        let counts = vec![10, 20, 30, 40];
        let result = hill_numbers(&counts, &[2.0]).unwrap();
        // Hill q=2 = 1 / Σ p_i^2 (probability Simpson, not the n(n-1) form)
        let n: f64 = counts.iter().sum::<usize>() as f64;
        let sum_p2: f64 = counts.iter().filter(|&&c| c > 0).map(|&c| (c as f64 / n).powi(2)).sum();
        let expected = 1.0 / sum_p2;
        assert!(
            (result[0].1 - expected).abs() < 1e-8,
            "q=2: {} expected {}",
            result[0].1,
            expected
        );
    }

    #[test]
    fn hill_monotone_non_increasing() {
        let counts = vec![100, 50, 25, 10, 5, 1];
        let orders: Vec<f64> = (0..=10).map(|i| i as f64 * 0.5).collect();
        let result = hill_numbers(&counts, &orders).unwrap();
        for w in result.windows(2) {
            assert!(
                w[0].1 >= w[1].1 - 1e-8,
                "not monotone: q={}: {} > q={}: {}",
                w[0].0, w[0].1, w[1].0, w[1].1
            );
        }
    }

    // ── Alpha rarefaction tests ─────────────────────────────────────

    #[test]
    fn alpha_rarefaction_multiple_samples() {
        let s1 = vec![100, 50, 25, 10, 5, 1];
        let s2 = vec![80, 40, 20, 10];
        let depths = vec![10, 50, 100];
        let result = alpha_rarefaction(&[&s1, &s2], &depths).unwrap();
        assert_eq!(result.len(), 2);
        // s2 total = 150, so depth 100 should be included
        assert!(result[1].len() >= 2);
    }

    #[test]
    fn rarefaction_monotonic() {
        let counts = vec![100, 50, 25, 10, 5, 1, 1];
        let total: usize = counts.iter().sum();
        let steps: Vec<usize> = (1..=5).map(|i| i * total / 5).collect();
        let curve = rarefaction_curve(&counts, &steps).unwrap();
        // Expected species richness should be non-decreasing.
        for window in curve.windows(2) {
            assert!(
                window[1].1 >= window[0].1 - 1e-10,
                "rarefaction curve is not monotonic: {:?}",
                curve
            );
        }
    }
}
