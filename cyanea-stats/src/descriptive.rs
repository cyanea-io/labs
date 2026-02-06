//! Descriptive statistics for numeric data.
//!
//! Provides individual functions ([`mean`], [`median`], [`variance`], etc.) and
//! the aggregate [`describe`] function that computes all common statistics in
//! one pass.

use cyanea_core::{CyaneaError, Result, Summarizable};

/// Aggregate descriptive statistics for a numeric sample.
#[derive(Debug, Clone)]
pub struct DescriptiveStats {
    /// Number of observations.
    pub count: usize,
    /// Arithmetic mean.
    pub mean: f64,
    /// Median (50th percentile).
    pub median: f64,
    /// Population variance (ddof=0).
    pub variance: f64,
    /// Sample variance (ddof=1).
    pub sample_variance: f64,
    /// Population standard deviation.
    pub std_dev: f64,
    /// Sample standard deviation.
    pub sample_std_dev: f64,
    /// Minimum value.
    pub min: f64,
    /// Maximum value.
    pub max: f64,
    /// Range (max - min).
    pub range: f64,
    /// First quartile (25th percentile).
    pub q1: f64,
    /// Third quartile (75th percentile).
    pub q3: f64,
    /// Interquartile range (q3 - q1).
    pub iqr: f64,
    /// Skewness (Fisher's definition, population).
    pub skewness: f64,
    /// Excess kurtosis (population).
    pub kurtosis: f64,
}

impl Summarizable for DescriptiveStats {
    fn summary(&self) -> String {
        format!(
            "n={}, mean={:.4}, std={:.4}, min={:.4}, max={:.4}",
            self.count, self.mean, self.std_dev, self.min, self.max,
        )
    }
}

/// Compute all descriptive statistics for `data`.
///
/// Requires at least 1 element. Sample variance/std_dev require at least 2.
pub fn describe(data: &[f64]) -> Result<DescriptiveStats> {
    if data.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "describe: data must not be empty".into(),
        ));
    }

    let n = data.len();
    let n_f = n as f64;

    // First pass: mean, min, max
    let mut sum = 0.0;
    let mut min_val = f64::INFINITY;
    let mut max_val = f64::NEG_INFINITY;
    for &x in data {
        sum += x;
        if x < min_val {
            min_val = x;
        }
        if x > max_val {
            max_val = x;
        }
    }
    let mean_val = sum / n_f;

    // Second pass: variance, skewness, kurtosis
    let mut m2 = 0.0;
    let mut m3 = 0.0;
    let mut m4 = 0.0;
    for &x in data {
        let d = x - mean_val;
        m2 += d * d;
        m3 += d * d * d;
        m4 += d * d * d * d;
    }

    let pop_var = m2 / n_f;
    let sample_var = if n > 1 { m2 / (n_f - 1.0) } else { f64::NAN };
    let pop_std = pop_var.sqrt();
    let sample_std = sample_var.sqrt();

    let skewness = if pop_std > 0.0 {
        (m3 / n_f) / (pop_std * pop_std * pop_std)
    } else {
        0.0
    };

    let kurtosis = if pop_std > 0.0 {
        (m4 / n_f) / (pop_var * pop_var) - 3.0
    } else {
        0.0
    };

    // Sort a clone for quantiles
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));

    let median_val = compute_quantile_sorted(&sorted, 0.5);
    let q1_val = compute_quantile_sorted(&sorted, 0.25);
    let q3_val = compute_quantile_sorted(&sorted, 0.75);

    Ok(DescriptiveStats {
        count: n,
        mean: mean_val,
        median: median_val,
        variance: pop_var,
        sample_variance: sample_var,
        std_dev: pop_std,
        sample_std_dev: sample_std,
        min: min_val,
        max: max_val,
        range: max_val - min_val,
        q1: q1_val,
        q3: q3_val,
        iqr: q3_val - q1_val,
        skewness,
        kurtosis,
    })
}

// ── Individual functions ───────────────────────────────────────────────────

/// Arithmetic mean.
pub fn mean(data: &[f64]) -> Result<f64> {
    if data.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "mean: data must not be empty".into(),
        ));
    }
    Ok(data.iter().sum::<f64>() / data.len() as f64)
}

/// Median (50th percentile).
pub fn median(data: &[f64]) -> Result<f64> {
    if data.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "median: data must not be empty".into(),
        ));
    }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    Ok(compute_quantile_sorted(&sorted, 0.5))
}

/// Variance with given degrees-of-freedom correction.
///
/// - `ddof = 0` → population variance
/// - `ddof = 1` → sample variance (Bessel's correction)
pub fn variance(data: &[f64], ddof: usize) -> Result<f64> {
    let n = data.len();
    if n <= ddof {
        return Err(CyaneaError::InvalidInput(format!(
            "variance: need more than {} observations (got {})",
            ddof, n,
        )));
    }
    let m = mean(data)?;
    let ss: f64 = data.iter().map(|&x| (x - m).powi(2)).sum();
    Ok(ss / (n - ddof) as f64)
}

/// Standard deviation with given degrees-of-freedom correction.
pub fn std_dev(data: &[f64], ddof: usize) -> Result<f64> {
    Ok(variance(data, ddof)?.sqrt())
}

/// Quantile using linear interpolation (exclusive method).
pub fn quantile(data: &[f64], q: f64) -> Result<f64> {
    if data.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "quantile: data must not be empty".into(),
        ));
    }
    if !(0.0..=1.0).contains(&q) {
        return Err(CyaneaError::InvalidInput(
            "quantile: q must be in [0, 1]".into(),
        ));
    }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    Ok(compute_quantile_sorted(&sorted, q))
}

/// Interquartile range (Q3 - Q1).
pub fn iqr(data: &[f64]) -> Result<f64> {
    let q1 = quantile(data, 0.25)?;
    let q3 = quantile(data, 0.75)?;
    Ok(q3 - q1)
}

/// Median absolute deviation (MAD).
pub fn mad(data: &[f64]) -> Result<f64> {
    let med = median(data)?;
    let deviations: Vec<f64> = data.iter().map(|&x| (x - med).abs()).collect();
    median(&deviations)
}

// ── Internal ───────────────────────────────────────────────────────────────

/// Compute a quantile from a pre-sorted slice using linear interpolation.
fn compute_quantile_sorted(sorted: &[f64], q: f64) -> f64 {
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    let pos = q * (n - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = lo + 1;
    let frac = pos - lo as f64;
    if hi >= n {
        sorted[n - 1]
    } else {
        sorted[lo] * (1.0 - frac) + sorted[hi] * frac
    }
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    #[test]
    fn describe_known_data() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let stats = describe(&data).unwrap();
        assert_eq!(stats.count, 5);
        assert!((stats.mean - 3.0).abs() < TOL);
        assert!((stats.median - 3.0).abs() < TOL);
        assert!((stats.min - 1.0).abs() < TOL);
        assert!((stats.max - 5.0).abs() < TOL);
        assert!((stats.range - 4.0).abs() < TOL);
        // Population variance of [1,2,3,4,5] = 2.0
        assert!((stats.variance - 2.0).abs() < TOL);
        // Sample variance = 2.5
        assert!((stats.sample_variance - 2.5).abs() < TOL);
    }

    #[test]
    fn describe_single() {
        let data = [42.0];
        let stats = describe(&data).unwrap();
        assert_eq!(stats.count, 1);
        assert!((stats.mean - 42.0).abs() < TOL);
        assert!((stats.variance - 0.0).abs() < TOL);
        assert!(stats.sample_variance.is_nan());
    }

    #[test]
    fn describe_empty() {
        assert!(describe(&[]).is_err());
    }

    #[test]
    fn mean_basic() {
        assert!((mean(&[2.0, 4.0, 6.0]).unwrap() - 4.0).abs() < TOL);
    }

    #[test]
    fn mean_empty() {
        assert!(mean(&[]).is_err());
    }

    #[test]
    fn median_odd() {
        assert!((median(&[3.0, 1.0, 2.0]).unwrap() - 2.0).abs() < TOL);
    }

    #[test]
    fn median_even() {
        assert!((median(&[4.0, 1.0, 3.0, 2.0]).unwrap() - 2.5).abs() < TOL);
    }

    #[test]
    fn variance_population() {
        let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        assert!((variance(&data, 0).unwrap() - 4.0).abs() < TOL);
    }

    #[test]
    fn variance_sample() {
        let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let expected = 32.0 / 7.0; // ≈ 4.571
        assert!((variance(&data, 1).unwrap() - expected).abs() < TOL);
    }

    #[test]
    fn variance_too_few() {
        assert!(variance(&[1.0], 1).is_err());
        assert!(variance(&[], 0).is_err());
    }

    #[test]
    fn quantile_basic() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((quantile(&data, 0.0).unwrap() - 1.0).abs() < TOL);
        assert!((quantile(&data, 1.0).unwrap() - 5.0).abs() < TOL);
        assert!((quantile(&data, 0.5).unwrap() - 3.0).abs() < TOL);
    }

    #[test]
    fn quantile_invalid_q() {
        assert!(quantile(&[1.0], -0.1).is_err());
        assert!(quantile(&[1.0], 1.1).is_err());
    }

    #[test]
    fn iqr_basic() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let result = iqr(&data).unwrap();
        // Q1 = 2.75, Q3 = 6.25, IQR = 3.5
        assert!((result - 3.5).abs() < TOL);
    }

    #[test]
    fn mad_basic() {
        let data = [1.0, 1.0, 2.0, 2.0, 4.0, 6.0, 9.0];
        // median = 2, deviations = [1,1,0,0,2,4,7], median of deviations = 1
        assert!((mad(&data).unwrap() - 1.0).abs() < TOL);
    }

    #[test]
    fn describe_skewness_symmetric() {
        // Symmetric data should have skewness ≈ 0
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let stats = describe(&data).unwrap();
        assert!((stats.skewness).abs() < TOL);
    }

    #[test]
    fn summarizable_impl() {
        let stats = describe(&[1.0, 2.0, 3.0]).unwrap();
        let s = stats.summary();
        assert!(s.contains("n=3"));
        assert!(s.contains("mean="));
    }
}
