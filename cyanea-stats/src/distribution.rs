//! Probability distributions and numerical helpers.
//!
//! Provides the [`Distribution`] trait and implementations for [`Normal`] and
//! [`Poisson`] distributions, plus low-level functions ([`erf`], [`ln_gamma`],
//! [`betai`]) used throughout the crate for p-value computation.

use core::f64::consts::PI;

use cyanea_core::{CyaneaError, Result};

// ── Numerical helpers ──────────────────────────────────────────────────────

/// Error function via Abramowitz & Stegun 7.1.26 (max error ~1.5e-7).
pub fn erf(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    sign * (1.0 - poly * (-x * x).exp())
}

/// Natural log of the gamma function via the Lanczos approximation (g=7).
pub fn ln_gamma(x: f64) -> f64 {
    const COEFFS: [f64; 8] = [
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
        // Reflection formula: Γ(x) = π / (sin(πx) · Γ(1-x))
        let log_pi_over_sin = (PI / (PI * x).sin()).ln();
        log_pi_over_sin - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut ag = 0.99999999999980993_f64;
        for (i, &c) in COEFFS.iter().enumerate() {
            ag += c / (x + i as f64 + 1.0);
        }
        let t = x + 7.5; // g + 0.5
        0.5 * (2.0 * PI).ln() + (x + 0.5) * t.ln() - t + ag.ln()
    }
}

/// Regularized incomplete beta function I_x(a, b) via continued fraction
/// (Lentz's method, max 200 iterations).
///
/// Used to compute p-values for the t-distribution and other tests.
pub fn betai(a: f64, b: f64, x: f64) -> Result<f64> {
    if x < 0.0 || x > 1.0 {
        return Err(CyaneaError::InvalidInput(
            "betai: x must be in [0, 1]".into(),
        ));
    }
    if x == 0.0 || x == 1.0 {
        return Ok(x);
    }

    // Use symmetry relation for numerical stability.
    if x > (a + 1.0) / (a + b + 2.0) {
        return Ok(1.0 - betai(b, a, 1.0 - x)?);
    }

    let ln_prefactor = ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b)
        + a * x.ln()
        + b * (1.0 - x).ln();
    let prefactor = ln_prefactor.exp();

    // Evaluate continued fraction with modified Lentz's method.
    let tiny = 1e-30_f64;
    let eps = 1e-10_f64;
    let max_iter = 200;

    let mut c = 1.0_f64;
    let mut d = (1.0 - (a + b) * x / (a + 1.0)).recip();
    if d.abs() < tiny {
        d = tiny;
    }
    let mut h = d;

    for m in 1..=max_iter {
        let m_f64 = m as f64;

        // Even step: d_{2m}
        let num_even = m_f64 * (b - m_f64) * x / ((a + 2.0 * m_f64 - 1.0) * (a + 2.0 * m_f64));
        d = 1.0 + num_even * d;
        if d.abs() < tiny {
            d = tiny;
        }
        d = d.recip();
        c = 1.0 + num_even / c;
        if c.abs() < tiny {
            c = tiny;
        }
        h *= d * c;

        // Odd step: d_{2m+1}
        let num_odd = -((a + m_f64) * (a + b + m_f64) * x)
            / ((a + 2.0 * m_f64) * (a + 2.0 * m_f64 + 1.0));
        d = 1.0 + num_odd * d;
        if d.abs() < tiny {
            d = tiny;
        }
        d = d.recip();
        c = 1.0 + num_odd / c;
        if c.abs() < tiny {
            c = tiny;
        }
        let delta = d * c;
        h *= delta;

        if (delta - 1.0).abs() < eps {
            return Ok(prefactor * h / a);
        }
    }

    Ok(prefactor * h / a)
}

// ── Distribution trait ─────────────────────────────────────────────────────

/// A probability distribution with basic statistical properties.
pub trait Distribution {
    /// Probability density (or mass) function at `x`.
    fn pdf(&self, x: f64) -> f64;

    /// Cumulative distribution function at `x`.
    fn cdf(&self, x: f64) -> f64;

    /// Distribution mean.
    fn mean(&self) -> f64;

    /// Distribution variance.
    fn variance(&self) -> f64;

    /// Distribution standard deviation (default: sqrt of variance).
    fn std_dev(&self) -> f64 {
        self.variance().sqrt()
    }
}

// ── Normal distribution ────────────────────────────────────────────────────

/// Normal (Gaussian) distribution with parameters μ and σ.
#[derive(Debug, Clone, Copy)]
pub struct Normal {
    mu: f64,
    sigma: f64,
}

impl Normal {
    /// Create a new Normal distribution. `sigma` must be positive.
    pub fn new(mu: f64, sigma: f64) -> Result<Self> {
        if sigma <= 0.0 {
            return Err(CyaneaError::InvalidInput(
                "Normal: sigma must be positive".into(),
            ));
        }
        Ok(Self { mu, sigma })
    }

    /// Standard normal distribution N(0, 1).
    pub fn standard() -> Self {
        Self {
            mu: 0.0,
            sigma: 1.0,
        }
    }
}

impl Distribution for Normal {
    fn pdf(&self, x: f64) -> f64 {
        let z = (x - self.mu) / self.sigma;
        (-0.5 * z * z).exp() / (self.sigma * (2.0 * PI).sqrt())
    }

    fn cdf(&self, x: f64) -> f64 {
        let z = (x - self.mu) / self.sigma;
        0.5 * (1.0 + erf(z / core::f64::consts::SQRT_2))
    }

    fn mean(&self) -> f64 {
        self.mu
    }

    fn variance(&self) -> f64 {
        self.sigma * self.sigma
    }
}

// ── Poisson distribution ───────────────────────────────────────────────────

/// Poisson distribution with rate parameter λ.
#[derive(Debug, Clone, Copy)]
pub struct Poisson {
    lambda: f64,
}

impl Poisson {
    /// Create a new Poisson distribution. `lambda` must be positive.
    pub fn new(lambda: f64) -> Result<Self> {
        if lambda <= 0.0 {
            return Err(CyaneaError::InvalidInput(
                "Poisson: lambda must be positive".into(),
            ));
        }
        Ok(Self { lambda })
    }
}

impl Distribution for Poisson {
    fn pdf(&self, x: f64) -> f64 {
        let k = x.round() as i64;
        if k < 0 || (x - k as f64).abs() > 1e-9 {
            return 0.0;
        }
        let k = k as f64;
        // Compute in log-space to avoid overflow.
        (k * self.lambda.ln() - self.lambda - ln_gamma(k + 1.0)).exp()
    }

    fn cdf(&self, x: f64) -> f64 {
        let k_max = x.floor() as i64;
        if k_max < 0 {
            return 0.0;
        }
        let mut sum = 0.0;
        for k in 0..=k_max {
            let k = k as f64;
            sum += (k * self.lambda.ln() - self.lambda - ln_gamma(k + 1.0)).exp();
        }
        sum.min(1.0)
    }

    fn mean(&self) -> f64 {
        self.lambda
    }

    fn variance(&self) -> f64 {
        self.lambda
    }
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    #[test]
    fn erf_zero() {
        assert!((erf(0.0)).abs() < TOL);
    }

    #[test]
    fn erf_one() {
        assert!((erf(1.0) - 0.8427007929).abs() < 1e-5);
    }

    #[test]
    fn erf_negative_symmetry() {
        assert!((erf(-0.5) + erf(0.5)).abs() < TOL);
    }

    #[test]
    fn ln_gamma_integers() {
        // Γ(n) = (n-1)! for positive integers
        assert!((ln_gamma(1.0) - 0.0).abs() < TOL); // 0! = 1
        assert!((ln_gamma(2.0) - 0.0).abs() < TOL); // 1! = 1
        assert!((ln_gamma(5.0) - (24.0_f64).ln()).abs() < TOL); // 4! = 24
        assert!((ln_gamma(7.0) - (720.0_f64).ln()).abs() < TOL); // 6! = 720
    }

    #[test]
    fn ln_gamma_half() {
        // Γ(0.5) = √π
        assert!((ln_gamma(0.5) - 0.5 * PI.ln()).abs() < 1e-5);
    }

    #[test]
    fn betai_boundaries() {
        assert_eq!(betai(1.0, 1.0, 0.0).unwrap(), 0.0);
        assert_eq!(betai(1.0, 1.0, 1.0).unwrap(), 1.0);
    }

    #[test]
    fn betai_uniform() {
        // Beta(1,1) is uniform, so I_x(1,1) = x
        assert!((betai(1.0, 1.0, 0.5).unwrap() - 0.5).abs() < TOL);
        assert!((betai(1.0, 1.0, 0.3).unwrap() - 0.3).abs() < TOL);
    }

    #[test]
    fn betai_symmetry() {
        // I_x(a,b) = 1 - I_{1-x}(b,a)
        let a = 2.0;
        let b = 3.0;
        let x = 0.4;
        let lhs = betai(a, b, x).unwrap();
        let rhs = 1.0 - betai(b, a, 1.0 - x).unwrap();
        assert!((lhs - rhs).abs() < TOL);
    }

    #[test]
    fn betai_invalid_x() {
        assert!(betai(1.0, 1.0, -0.1).is_err());
        assert!(betai(1.0, 1.0, 1.1).is_err());
    }

    #[test]
    fn normal_standard_cdf() {
        let n = Normal::standard();
        assert!((n.cdf(0.0) - 0.5).abs() < TOL);
        assert!((n.cdf(1.0) - 0.8413447).abs() < 1e-5);
        assert!((n.cdf(-1.0) - 0.1586553).abs() < 1e-5);
        assert!((n.cdf(2.0) - 0.9772499).abs() < 1e-5);
    }

    #[test]
    fn normal_standard_pdf_at_zero() {
        let n = Normal::standard();
        let expected = 1.0 / (2.0 * PI).sqrt();
        assert!((n.pdf(0.0) - expected).abs() < TOL);
    }

    #[test]
    fn normal_invalid_sigma() {
        assert!(Normal::new(0.0, 0.0).is_err());
        assert!(Normal::new(0.0, -1.0).is_err());
    }

    #[test]
    fn poisson_pmf() {
        let p = Poisson::new(3.0).unwrap();
        // P(X=0) = e^(-3) ≈ 0.04979
        assert!((p.pdf(0.0) - (-3.0_f64).exp()).abs() < TOL);
        // P(X=3) = 3^3 * e^(-3) / 3! ≈ 0.22404
        let expected = 27.0 * (-3.0_f64).exp() / 6.0;
        assert!((p.pdf(3.0) - expected).abs() < TOL);
    }

    #[test]
    fn poisson_cdf() {
        let p = Poisson::new(1.0).unwrap();
        // P(X <= 0) = e^(-1) ≈ 0.36788
        assert!((p.cdf(0.0) - (-1.0_f64).exp()).abs() < TOL);
        // P(X <= 1) = e^(-1) + e^(-1) = 2 * e^(-1) ≈ 0.73576
        assert!((p.cdf(1.0) - 2.0 * (-1.0_f64).exp()).abs() < TOL);
    }

    #[test]
    fn poisson_invalid_lambda() {
        assert!(Poisson::new(0.0).is_err());
        assert!(Poisson::new(-1.0).is_err());
    }
}
