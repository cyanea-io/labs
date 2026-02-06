//! Hypothesis testing.
//!
//! Provides parametric ([`t_test_one_sample`], [`t_test_two_sample`]) and
//! non-parametric ([`mann_whitney_u`]) statistical tests.

use cyanea_core::{CyaneaError, Result, Scored, Summarizable};

use crate::descriptive;
use crate::distribution::{betai, Normal, Distribution};
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
}
