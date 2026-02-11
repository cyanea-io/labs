//! Hypothesis testing.
//!
//! Provides parametric ([`t_test_one_sample`], [`t_test_two_sample`]) and
//! non-parametric ([`mann_whitney_u`]) statistical tests.

use cyanea_core::{CyaneaError, Result, Scored, Summarizable};

use crate::descriptive;
use crate::distribution::{betai, ln_gamma, ChiSquared, FDistribution, Normal, Distribution};
use crate::rank::{rank, RankMethod};

/// Result of a hypothesis test.
#[derive(Debug, Clone)]
pub struct TestResult {
    /// The test statistic (t, U, z, etc.).
    pub statistic: f64,
    /// Two-tailed p-value.
    pub p_value: f64,
    /// Degrees of freedom, if applicable.
    pub degrees_of_freedom: Option<f64>,
    /// Name of the test method.
    pub method: String,
}

impl Scored for TestResult {
    fn score(&self) -> f64 {
        self.p_value
    }
}

impl Summarizable for TestResult {
    fn summary(&self) -> String {
        match self.degrees_of_freedom {
            Some(df) => format!(
                "{}: statistic={:.4}, df={:.1}, p={:.6}",
                self.method, self.statistic, df, self.p_value,
            ),
            None => format!(
                "{}: statistic={:.4}, p={:.6}",
                self.method, self.statistic, self.p_value,
            ),
        }
    }
}

// ── t-distribution helpers ──────────────────────────────────────────────────

/// Two-tailed p-value for the t-distribution.
fn t_two_tailed_p(t: f64, df: f64) -> f64 {
    let x = df / (df + t * t);
    betai(df / 2.0, 0.5, x).unwrap_or(1.0)
}

// ── One-sample t-test ──────────────────────────────────────────────────────

/// One-sample t-test: test whether the population mean equals `mu`.
///
/// Requires at least 2 observations.
pub fn t_test_one_sample(data: &[f64], mu: f64) -> Result<TestResult> {
    let n = data.len();
    if n < 2 {
        return Err(CyaneaError::InvalidInput(
            "t_test_one_sample: need at least 2 observations".into(),
        ));
    }

    let mean = descriptive::mean(data)?;
    let se = descriptive::std_dev(data, 1)? / (n as f64).sqrt();
    let t = (mean - mu) / se;
    let df = (n - 1) as f64;
    let p = t_two_tailed_p(t, df);

    Ok(TestResult {
        statistic: t,
        p_value: p,
        degrees_of_freedom: Some(df),
        method: "One-sample t-test".into(),
    })
}

// ── Two-sample t-test ──────────────────────────────────────────────────────

/// Two-sample t-test: test whether two populations have the same mean.
///
/// When `equal_var` is `true`, uses pooled variance (Student's t-test).
/// When `false`, uses Welch's t-test (unequal variances).
pub fn t_test_two_sample(x: &[f64], y: &[f64], equal_var: bool) -> Result<TestResult> {
    if x.len() < 2 || y.len() < 2 {
        return Err(CyaneaError::InvalidInput(
            "t_test_two_sample: each group needs at least 2 observations".into(),
        ));
    }

    let nx = x.len() as f64;
    let ny = y.len() as f64;
    let mean_x = descriptive::mean(x)?;
    let mean_y = descriptive::mean(y)?;
    let var_x = descriptive::variance(x, 1)?;
    let var_y = descriptive::variance(y, 1)?;

    let (t, df) = if equal_var {
        // Pooled variance
        let sp2 = ((nx - 1.0) * var_x + (ny - 1.0) * var_y) / (nx + ny - 2.0);
        let se = (sp2 * (1.0 / nx + 1.0 / ny)).sqrt();
        let t = (mean_x - mean_y) / se;
        let df = nx + ny - 2.0;
        (t, df)
    } else {
        // Welch's approximation
        let se = (var_x / nx + var_y / ny).sqrt();
        let t = (mean_x - mean_y) / se;
        let vn_x = var_x / nx;
        let vn_y = var_y / ny;
        let num = (vn_x + vn_y).powi(2);
        let denom = vn_x.powi(2) / (nx - 1.0) + vn_y.powi(2) / (ny - 1.0);
        let df = num / denom;
        (t, df)
    };

    let p = t_two_tailed_p(t, df);
    let method = if equal_var {
        "Two-sample t-test (pooled)"
    } else {
        "Welch's t-test"
    };

    Ok(TestResult {
        statistic: t,
        p_value: p,
        degrees_of_freedom: Some(df),
        method: method.into(),
    })
}

// ── Mann-Whitney U test ────────────────────────────────────────────────────

/// Mann-Whitney U test (Wilcoxon rank-sum test).
///
/// Non-parametric test for whether two independent samples come from the
/// same distribution. Uses normal approximation for the p-value.
///
/// Each group needs at least 1 observation, and total n must be >= 2.
pub fn mann_whitney_u(x: &[f64], y: &[f64]) -> Result<TestResult> {
    if x.is_empty() || y.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "mann_whitney_u: each group must be non-empty".into(),
        ));
    }
    let nx = x.len();
    let ny = y.len();
    let n = nx + ny;
    if n < 2 {
        return Err(CyaneaError::InvalidInput(
            "mann_whitney_u: need at least 2 total observations".into(),
        ));
    }

    // Combine, rank, and sum ranks for x.
    let mut combined: Vec<f64> = Vec::with_capacity(n);
    combined.extend_from_slice(x);
    combined.extend_from_slice(y);
    let ranks = rank(&combined, RankMethod::Average);

    let r1: f64 = ranks[..nx].iter().sum();
    let u1 = r1 - (nx * (nx + 1)) as f64 / 2.0;
    let u2 = (nx * ny) as f64 - u1;
    let u = u1.min(u2);

    // Normal approximation
    let mu_u = (nx * ny) as f64 / 2.0;
    let sigma_u = ((nx * ny * (n + 1)) as f64 / 12.0).sqrt();

    let p = if sigma_u > 0.0 {
        let z = (u - mu_u) / sigma_u;
        // Two-tailed p-value via standard normal
        let normal = Normal::standard();
        (2.0 * normal.cdf(z)).min(1.0) // z <= 0 since u = min(u1,u2) <= mu_u
    } else {
        1.0
    };

    Ok(TestResult {
        statistic: u,
        p_value: p,
        degrees_of_freedom: None,
        method: "Mann-Whitney U test".into(),
    })
}

// ── Fisher's exact test (2×2) ─────────────────────────────────────────────

/// Fisher's exact test for a 2×2 contingency table.
///
/// The table is specified as `[[a, b], [c, d]]`:
///
/// ```text
///           Group 1   Group 2
/// Outcome A    a         b
/// Outcome B    c         d
/// ```
///
/// Returns a two-tailed p-value based on the hypergeometric distribution.
pub fn fisher_exact(table: &[[usize; 2]; 2]) -> Result<TestResult> {
    let a = table[0][0];
    let b = table[0][1];
    let c = table[1][0];
    let d = table[1][1];
    let n = a + b + c + d;

    if n == 0 {
        return Err(CyaneaError::InvalidInput("fisher_exact: table is all zeros".into()));
    }

    let p_observed = hypergeometric_pmf(a, a + b, a + c, n);

    // Two-tailed: sum probabilities of tables as or more extreme than observed
    let row1 = a + b;
    let col1 = a + c;
    let min_a = if row1 + col1 > n { row1 + col1 - n } else { 0 };
    let max_a = row1.min(col1);

    let mut p_value = 0.0;
    for k in min_a..=max_a {
        let p_k = hypergeometric_pmf(k, row1, col1, n);
        if p_k <= p_observed + 1e-12 {
            p_value += p_k;
        }
    }

    Ok(TestResult {
        statistic: p_observed,
        p_value: p_value.min(1.0),
        degrees_of_freedom: None,
        method: "Fisher's exact test".into(),
    })
}

/// Hypergeometric PMF: P(X = k) where X ~ Hypergeometric(N, K, n).
///
/// Probability of drawing exactly `k` successes from a population of `total`
/// containing `success_pop` successes, in a sample of size `sample_size`.
fn hypergeometric_pmf(k: usize, sample_size: usize, success_pop: usize, total: usize) -> f64 {
    // P = C(K,k) * C(N-K, n-k) / C(N, n)
    // Compute in log-space to avoid overflow.
    let log_p = ln_choose(success_pop, k)
        + ln_choose(total - success_pop, sample_size - k)
        - ln_choose(total, sample_size);
    log_p.exp()
}

/// Log of binomial coefficient C(n, k) = ln(n!) - ln(k!) - ln((n-k)!).
fn ln_choose(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    ln_gamma(n as f64 + 1.0) - ln_gamma(k as f64 + 1.0) - ln_gamma((n - k) as f64 + 1.0)
}

// ── Chi-squared test of independence ──────────────────────────────────────

/// Chi-squared test of independence for an r×c contingency table.
///
/// `observed` is a row-major slice of `nrows × ncols` observed counts.
///
/// Uses Pearson's chi-squared statistic: χ² = Σ (O - E)² / E
pub fn chi_squared_test(observed: &[f64], nrows: usize, ncols: usize) -> Result<TestResult> {
    if nrows < 2 || ncols < 2 {
        return Err(CyaneaError::InvalidInput(
            "chi_squared_test: need at least 2×2 table".into(),
        ));
    }
    if observed.len() != nrows * ncols {
        return Err(CyaneaError::InvalidInput(
            "chi_squared_test: observed length must equal nrows × ncols".into(),
        ));
    }

    let total: f64 = observed.iter().sum();
    if total == 0.0 {
        return Err(CyaneaError::InvalidInput("chi_squared_test: all counts are zero".into()));
    }

    // Row and column sums
    let mut row_sums = vec![0.0; nrows];
    let mut col_sums = vec![0.0; ncols];
    for i in 0..nrows {
        for j in 0..ncols {
            let val = observed[i * ncols + j];
            row_sums[i] += val;
            col_sums[j] += val;
        }
    }

    // Compute chi-squared statistic
    let mut chi2 = 0.0;
    for i in 0..nrows {
        for j in 0..ncols {
            let expected = row_sums[i] * col_sums[j] / total;
            if expected > 0.0 {
                let diff = observed[i * ncols + j] - expected;
                chi2 += diff * diff / expected;
            }
        }
    }

    let df = ((nrows - 1) * (ncols - 1)) as f64;
    let chi2_dist = ChiSquared::new(df)?;
    let p_value = 1.0 - chi2_dist.cdf(chi2);

    Ok(TestResult {
        statistic: chi2,
        p_value,
        degrees_of_freedom: Some(df),
        method: "Chi-squared test of independence".into(),
    })
}

// ── One-way ANOVA ─────────────────────────────────────────────────────────

/// One-way analysis of variance (ANOVA).
///
/// Tests whether the means of k groups are equal. Each group must have at
/// least 1 observation, and there must be at least 2 groups.
pub fn anova_oneway(groups: &[&[f64]]) -> Result<TestResult> {
    let k = groups.len();
    if k < 2 {
        return Err(CyaneaError::InvalidInput(
            "anova_oneway: need at least 2 groups".into(),
        ));
    }
    for (i, g) in groups.iter().enumerate() {
        if g.is_empty() {
            return Err(CyaneaError::InvalidInput(
                format!("anova_oneway: group {} is empty", i),
            ));
        }
    }

    let n_total: usize = groups.iter().map(|g| g.len()).sum();
    if n_total <= k {
        return Err(CyaneaError::InvalidInput(
            "anova_oneway: total observations must exceed number of groups".into(),
        ));
    }

    // Grand mean
    let grand_sum: f64 = groups.iter().flat_map(|g| g.iter()).sum();
    let grand_mean = grand_sum / n_total as f64;

    // Between-group sum of squares
    let ss_between: f64 = groups
        .iter()
        .map(|g| {
            let group_mean: f64 = g.iter().sum::<f64>() / g.len() as f64;
            g.len() as f64 * (group_mean - grand_mean).powi(2)
        })
        .sum();

    // Within-group sum of squares
    let ss_within: f64 = groups
        .iter()
        .map(|g| {
            let group_mean: f64 = g.iter().sum::<f64>() / g.len() as f64;
            g.iter().map(|&x| (x - group_mean).powi(2)).sum::<f64>()
        })
        .sum();

    let df_between = (k - 1) as f64;
    let df_within = (n_total - k) as f64;

    let ms_between = ss_between / df_between;
    let ms_within = ss_within / df_within;

    let f_stat = if ms_within > 0.0 {
        ms_between / ms_within
    } else {
        f64::INFINITY
    };

    let f_dist = FDistribution::new(df_between, df_within)?;
    let p_value = 1.0 - f_dist.cdf(f_stat);

    Ok(TestResult {
        statistic: f_stat,
        p_value,
        degrees_of_freedom: Some(df_between),
        method: "One-way ANOVA".into(),
    })
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn t_test_one_sample_mean_equals_mu() {
        // Data centered on mu=0 → large p-value
        let data = [-1.0, -0.5, 0.0, 0.5, 1.0];
        let result = t_test_one_sample(&data, 0.0).unwrap();
        assert!(result.p_value > 0.9, "p={}", result.p_value);
    }

    #[test]
    fn t_test_one_sample_mean_far_from_mu() {
        let data = [10.0, 11.0, 12.0, 13.0, 14.0];
        let result = t_test_one_sample(&data, 0.0).unwrap();
        assert!(result.p_value < 0.001, "p={}", result.p_value);
    }

    #[test]
    fn t_test_one_sample_too_few() {
        assert!(t_test_one_sample(&[1.0], 0.0).is_err());
    }

    #[test]
    fn t_test_two_sample_same_distribution() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [1.5, 2.5, 3.5, 4.5, 5.5];
        let result = t_test_two_sample(&x, &y, true).unwrap();
        // Means are close, p should be moderate to large
        assert!(result.p_value > 0.3, "p={}", result.p_value);
    }

    #[test]
    fn t_test_two_sample_different_means() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [100.0, 101.0, 102.0, 103.0, 104.0];
        let result = t_test_two_sample(&x, &y, true).unwrap();
        assert!(result.p_value < 0.001, "p={}", result.p_value);
    }

    #[test]
    fn t_test_welch() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [100.0, 101.0, 102.0, 103.0, 104.0];
        let result = t_test_two_sample(&x, &y, false).unwrap();
        assert!(result.p_value < 0.001, "p={}", result.p_value);
        assert!(result.method.contains("Welch"));
    }

    #[test]
    fn t_test_two_sample_too_few() {
        assert!(t_test_two_sample(&[1.0], &[2.0, 3.0], true).is_err());
    }

    #[test]
    fn mann_whitney_same() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [1.5, 2.5, 3.5, 4.5, 5.5];
        let result = mann_whitney_u(&x, &y).unwrap();
        assert!(result.p_value > 0.3, "p={}", result.p_value);
    }

    #[test]
    fn mann_whitney_different() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [100.0, 101.0, 102.0, 103.0, 104.0];
        let result = mann_whitney_u(&x, &y).unwrap();
        assert!(result.p_value < 0.05, "p={}", result.p_value);
    }

    #[test]
    fn mann_whitney_empty() {
        assert!(mann_whitney_u(&[], &[1.0]).is_err());
        assert!(mann_whitney_u(&[1.0], &[]).is_err());
    }

    #[test]
    fn test_result_scored() {
        let result = t_test_one_sample(&[1.0, 2.0, 3.0], 2.0).unwrap();
        assert!((result.score() - result.p_value).abs() < 1e-15);
    }

    #[test]
    fn test_result_summary() {
        let result = t_test_one_sample(&[1.0, 2.0, 3.0, 4.0, 5.0], 0.0).unwrap();
        let s = result.summary();
        assert!(s.contains("One-sample t-test"));
        assert!(s.contains("statistic="));
        assert!(s.contains("p="));
    }

    // ── Fisher's exact test ────────────────────────────────────────────

    #[test]
    fn fisher_exact_significant() {
        // Classic lady tasting tea: strong association
        let table = [[8, 1], [1, 8]];
        let result = fisher_exact(&table).unwrap();
        assert!(result.p_value < 0.05, "p={}", result.p_value);
    }

    #[test]
    fn fisher_exact_not_significant() {
        // No association
        let table = [[5, 5], [5, 5]];
        let result = fisher_exact(&table).unwrap();
        assert!(result.p_value > 0.5, "p={}", result.p_value);
    }

    #[test]
    fn fisher_exact_extreme() {
        // Perfect association
        let table = [[10, 0], [0, 10]];
        let result = fisher_exact(&table).unwrap();
        assert!(result.p_value < 0.001, "p={}", result.p_value);
    }

    #[test]
    fn fisher_exact_zero_table() {
        let table = [[0, 0], [0, 0]];
        assert!(fisher_exact(&table).is_err());
    }

    // ── Chi-squared test ───────────────────────────────────────────────

    #[test]
    fn chi_squared_test_independent() {
        // Observed ≈ expected → not significant
        #[rustfmt::skip]
        let observed = [
            50.0, 50.0,
            50.0, 50.0,
        ];
        let result = chi_squared_test(&observed, 2, 2).unwrap();
        assert!(result.p_value > 0.9, "p={}", result.p_value);
    }

    #[test]
    fn chi_squared_test_dependent() {
        // Strong deviation from expected
        #[rustfmt::skip]
        let observed = [
            90.0, 10.0,
            10.0, 90.0,
        ];
        let result = chi_squared_test(&observed, 2, 2).unwrap();
        assert!(result.p_value < 0.001, "p={}", result.p_value);
        assert!((result.degrees_of_freedom.unwrap() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn chi_squared_test_3x3() {
        #[rustfmt::skip]
        let observed = [
            10.0, 20.0, 30.0,
            20.0, 30.0, 10.0,
            30.0, 10.0, 20.0,
        ];
        let result = chi_squared_test(&observed, 3, 3).unwrap();
        assert!((result.degrees_of_freedom.unwrap() - 4.0).abs() < 1e-10);
        assert!(result.p_value < 0.05, "p={}", result.p_value);
    }

    #[test]
    fn chi_squared_test_invalid() {
        assert!(chi_squared_test(&[1.0], 1, 1).is_err());
        assert!(chi_squared_test(&[1.0, 2.0], 2, 2).is_err()); // length mismatch
    }

    // ── ANOVA ──────────────────────────────────────────────────────────

    #[test]
    fn anova_same_groups() {
        let g1 = [1.0, 2.0, 3.0, 4.0, 5.0];
        let g2 = [1.5, 2.5, 3.5, 4.5, 5.5];
        let g3 = [1.0, 2.0, 3.0, 4.0, 5.0];
        let result = anova_oneway(&[&g1, &g2, &g3]).unwrap();
        assert!(result.p_value > 0.3, "p={}", result.p_value);
    }

    #[test]
    fn anova_different_groups() {
        let g1 = [1.0, 2.0, 3.0, 4.0, 5.0];
        let g2 = [100.0, 101.0, 102.0, 103.0, 104.0];
        let g3 = [200.0, 201.0, 202.0, 203.0, 204.0];
        let result = anova_oneway(&[&g1, &g2, &g3]).unwrap();
        assert!(result.p_value < 0.001, "p={}", result.p_value);
        assert!(result.method.contains("ANOVA"));
    }

    #[test]
    fn anova_two_groups_matches_t() {
        // With 2 groups, ANOVA F = t² and p-values should agree
        let g1 = [1.0, 2.0, 3.0, 4.0, 5.0];
        let g2 = [3.0, 4.0, 5.0, 6.0, 7.0];
        let anova_result = anova_oneway(&[&g1, &g2]).unwrap();
        let t_result = t_test_two_sample(&g1, &g2, true).unwrap();
        assert!((anova_result.p_value - t_result.p_value).abs() < 0.01);
    }

    #[test]
    fn anova_too_few_groups() {
        assert!(anova_oneway(&[&[1.0, 2.0]]).is_err());
    }

    #[test]
    fn anova_empty_group() {
        let g1: [f64; 0] = [];
        assert!(anova_oneway(&[&g1, &[1.0, 2.0]]).is_err());
    }
}
