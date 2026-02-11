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

// ── Regularized lower incomplete gamma function ──────────────────────────

/// Regularized lower incomplete gamma function P(a, x) = γ(a, x) / Γ(a).
///
/// Uses the series expansion when x < a + 1 and the continued fraction
/// representation (computing Q = 1 - P) otherwise.
pub fn gammainc(a: f64, x: f64) -> Result<f64> {
    if a <= 0.0 {
        return Err(CyaneaError::InvalidInput("gammainc: a must be positive".into()));
    }
    if x < 0.0 {
        return Err(CyaneaError::InvalidInput("gammainc: x must be non-negative".into()));
    }
    if x == 0.0 {
        return Ok(0.0);
    }

    if x < a + 1.0 {
        // Series expansion
        gammainc_series(a, x)
    } else {
        // Continued fraction for upper gamma, then P = 1 - Q
        let q = gammainc_cf(a, x)?;
        Ok(1.0 - q)
    }
}

/// Series expansion for P(a, x).
fn gammainc_series(a: f64, x: f64) -> Result<f64> {
    let max_iter = 200;
    let eps = 1e-12;
    let ln_prefix = a * x.ln() - x - ln_gamma(a);

    let mut sum = 1.0 / a;
    let mut term = 1.0 / a;

    for n in 1..=max_iter {
        term *= x / (a + n as f64);
        sum += term;
        if term.abs() < sum.abs() * eps {
            return Ok(sum * ln_prefix.exp());
        }
    }

    Ok(sum * ln_prefix.exp())
}

/// Continued fraction for Q(a, x) = 1 - P(a, x) via modified Lentz's method.
fn gammainc_cf(a: f64, x: f64) -> Result<f64> {
    let max_iter = 200;
    let eps = 1e-12;
    let tiny = 1e-30_f64;
    let ln_prefix = a * x.ln() - x - ln_gamma(a);

    let mut b = x + 1.0 - a;
    let mut c = 1.0 / tiny;
    let mut d = 1.0 / b;
    let mut h = d;

    for i in 1..=max_iter {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < tiny {
            d = tiny;
        }
        c = b + an / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        let delta = d * c;
        h *= delta;
        if (delta - 1.0).abs() < eps {
            break;
        }
    }

    Ok(h * ln_prefix.exp())
}

// ── Chi-squared distribution ──────────────────────────────────────────────

/// Chi-squared distribution with k degrees of freedom.
#[derive(Debug, Clone, Copy)]
pub struct ChiSquared {
    k: f64,
}

impl ChiSquared {
    /// Create a chi-squared distribution with `k` degrees of freedom.
    pub fn new(k: f64) -> Result<Self> {
        if k <= 0.0 {
            return Err(CyaneaError::InvalidInput(
                "ChiSquared: k must be positive".into(),
            ));
        }
        Ok(Self { k })
    }

    /// Degrees of freedom.
    pub fn df(&self) -> f64 {
        self.k
    }
}

impl Distribution for ChiSquared {
    fn pdf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        let half_k = self.k / 2.0;
        let ln_pdf = (half_k - 1.0) * x.ln() - x / 2.0 - half_k * 2.0_f64.ln() - ln_gamma(half_k);
        ln_pdf.exp()
    }

    fn cdf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        gammainc(self.k / 2.0, x / 2.0).unwrap_or(0.0)
    }

    fn mean(&self) -> f64 {
        self.k
    }

    fn variance(&self) -> f64 {
        2.0 * self.k
    }
}

// ── F-distribution ────────────────────────────────────────────────────────

/// F-distribution with d1 and d2 degrees of freedom.
#[derive(Debug, Clone, Copy)]
pub struct FDistribution {
    d1: f64,
    d2: f64,
}

impl FDistribution {
    /// Create an F-distribution with `d1` and `d2` degrees of freedom.
    pub fn new(d1: f64, d2: f64) -> Result<Self> {
        if d1 <= 0.0 || d2 <= 0.0 {
            return Err(CyaneaError::InvalidInput(
                "FDistribution: both d1 and d2 must be positive".into(),
            ));
        }
        Ok(Self { d1, d2 })
    }
}

impl Distribution for FDistribution {
    fn pdf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        let d1 = self.d1;
        let d2 = self.d2;
        let ln_pdf = 0.5 * d1 * (d1 * x / (d1 * x + d2)).ln()
            + 0.5 * d2 * (d2 / (d1 * x + d2)).ln()
            - x.ln()
            - ln_gamma(d1 / 2.0)
            - ln_gamma(d2 / 2.0)
            + ln_gamma((d1 + d2) / 2.0);
        ln_pdf.exp()
    }

    fn cdf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        let ix = self.d1 * x / (self.d1 * x + self.d2);
        betai(self.d1 / 2.0, self.d2 / 2.0, ix).unwrap_or(0.0)
    }

    fn mean(&self) -> f64 {
        if self.d2 > 2.0 {
            self.d2 / (self.d2 - 2.0)
        } else {
            f64::INFINITY
        }
    }

    fn variance(&self) -> f64 {
        if self.d2 > 4.0 {
            let d1 = self.d1;
            let d2 = self.d2;
            2.0 * d2 * d2 * (d1 + d2 - 2.0)
                / (d1 * (d2 - 2.0).powi(2) * (d2 - 4.0))
        } else {
            f64::INFINITY
        }
    }
}

// ── Binomial distribution ─────────────────────────────────────────────────

/// Binomial distribution with parameters n (trials) and p (success probability).
#[derive(Debug, Clone, Copy)]
pub struct Binomial {
    n: usize,
    p: f64,
}

impl Binomial {
    /// Create a binomial distribution. `p` must be in [0, 1].
    pub fn new(n: usize, p: f64) -> Result<Self> {
        if !(0.0..=1.0).contains(&p) {
            return Err(CyaneaError::InvalidInput(
                "Binomial: p must be in [0, 1]".into(),
            ));
        }
        Ok(Self { n, p })
    }

    /// Number of trials.
    pub fn trials(&self) -> usize {
        self.n
    }

    /// Success probability.
    pub fn prob(&self) -> f64 {
        self.p
    }

    /// Probability mass function P(X = k).
    pub fn pmf(&self, k: usize) -> f64 {
        if k > self.n {
            return 0.0;
        }
        let ln_binom = ln_gamma(self.n as f64 + 1.0)
            - ln_gamma(k as f64 + 1.0)
            - ln_gamma((self.n - k) as f64 + 1.0);
        let ln_pmf = ln_binom + k as f64 * self.p.ln() + (self.n - k) as f64 * (1.0 - self.p).ln();
        ln_pmf.exp()
    }
}

impl Distribution for Binomial {
    fn pdf(&self, x: f64) -> f64 {
        let k = x.round() as i64;
        if k < 0 || (x - k as f64).abs() > 1e-9 {
            return 0.0;
        }
        self.pmf(k as usize)
    }

    fn cdf(&self, x: f64) -> f64 {
        let k_max = x.floor() as i64;
        if k_max < 0 {
            return 0.0;
        }
        let k_max = k_max as usize;
        if k_max >= self.n {
            return 1.0;
        }
        // Use regularized incomplete beta: P(X <= k) = I_{1-p}(n-k, k+1)
        betai((self.n - k_max) as f64, (k_max + 1) as f64, 1.0 - self.p).unwrap_or(1.0)
    }

    fn mean(&self) -> f64 {
        self.n as f64 * self.p
    }

    fn variance(&self) -> f64 {
        self.n as f64 * self.p * (1.0 - self.p)
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

    // ── gammainc tests ─────────────────────────────────────────────────

    #[test]
    fn gammainc_zero() {
        assert_eq!(gammainc(1.0, 0.0).unwrap(), 0.0);
    }

    #[test]
    fn gammainc_exponential() {
        // P(1, x) = 1 - e^{-x} for exponential distribution
        let x: f64 = 2.0;
        let expected = 1.0 - (-x).exp();
        assert!((gammainc(1.0, x).unwrap() - expected).abs() < 1e-8);
    }

    #[test]
    fn gammainc_half_integer() {
        // P(0.5, x) = erf(sqrt(x)) for a = 0.5
        let x: f64 = 1.0;
        let expected = erf(x.sqrt());
        assert!((gammainc(0.5, x).unwrap() - expected).abs() < 1e-6);
    }

    #[test]
    fn gammainc_large_x() {
        // For large x, P(a, x) → 1
        assert!((gammainc(2.0, 50.0).unwrap() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn gammainc_invalid() {
        assert!(gammainc(-1.0, 1.0).is_err());
        assert!(gammainc(1.0, -1.0).is_err());
    }

    // ── Chi-squared tests ──────────────────────────────────────────────

    #[test]
    fn chi_squared_cdf_known_values() {
        let chi2 = ChiSquared::new(2.0).unwrap();
        // χ²(2) CDF at x: P = 1 - e^{-x/2}
        let x = 5.991; // ≈ p=0.05 for df=2
        let p = chi2.cdf(x);
        assert!((p - 0.95).abs() < 0.01, "p={}", p);
    }

    #[test]
    fn chi_squared_cdf_df1() {
        let chi2 = ChiSquared::new(1.0).unwrap();
        // χ²(1) at 3.841 ≈ p=0.95
        assert!((chi2.cdf(3.841) - 0.95).abs() < 0.01);
    }

    #[test]
    fn chi_squared_mean_variance() {
        let chi2 = ChiSquared::new(5.0).unwrap();
        assert!((chi2.mean() - 5.0).abs() < TOL);
        assert!((chi2.variance() - 10.0).abs() < TOL);
    }

    #[test]
    fn chi_squared_cdf_at_zero() {
        let chi2 = ChiSquared::new(3.0).unwrap();
        assert_eq!(chi2.cdf(0.0), 0.0);
    }

    #[test]
    fn chi_squared_invalid() {
        assert!(ChiSquared::new(0.0).is_err());
        assert!(ChiSquared::new(-1.0).is_err());
    }

    // ── F-distribution tests ───────────────────────────────────────────

    #[test]
    fn f_dist_cdf_known() {
        let f = FDistribution::new(5.0, 10.0).unwrap();
        // F(5,10) at 3.326 ≈ p=0.95
        let p = f.cdf(3.326);
        assert!((p - 0.95).abs() < 0.02, "p={}", p);
    }

    #[test]
    fn f_dist_cdf_at_zero() {
        let f = FDistribution::new(3.0, 5.0).unwrap();
        assert_eq!(f.cdf(0.0), 0.0);
    }

    #[test]
    fn f_dist_mean() {
        let f = FDistribution::new(4.0, 8.0).unwrap();
        // Mean = d2/(d2-2) = 8/6 ≈ 1.333
        assert!((f.mean() - 8.0 / 6.0).abs() < TOL);
    }

    #[test]
    fn f_dist_invalid() {
        assert!(FDistribution::new(0.0, 5.0).is_err());
        assert!(FDistribution::new(5.0, 0.0).is_err());
    }

    // ── Binomial tests ─────────────────────────────────────────────────

    #[test]
    fn binomial_pmf() {
        let b = Binomial::new(10, 0.5).unwrap();
        // P(X=5) = C(10,5) * 0.5^10 = 252/1024 ≈ 0.24609
        assert!((b.pmf(5) - 0.24609375).abs() < 1e-6);
    }

    #[test]
    fn binomial_pmf_sum() {
        let b = Binomial::new(8, 0.3).unwrap();
        let sum: f64 = (0..=8).map(|k| b.pmf(k)).sum();
        assert!((sum - 1.0).abs() < 1e-8);
    }

    #[test]
    fn binomial_cdf() {
        let b = Binomial::new(10, 0.5).unwrap();
        // P(X <= 5) ≈ 0.623047
        assert!((b.cdf(5.0) - 0.623047).abs() < 0.01);
    }

    #[test]
    fn binomial_cdf_boundaries() {
        let b = Binomial::new(5, 0.5).unwrap();
        assert!(b.cdf(-1.0) == 0.0);
        assert!((b.cdf(5.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn binomial_mean_variance() {
        let b = Binomial::new(20, 0.3).unwrap();
        assert!((b.mean() - 6.0).abs() < TOL);
        assert!((b.variance() - 4.2).abs() < TOL);
    }

    #[test]
    fn binomial_invalid() {
        assert!(Binomial::new(10, -0.1).is_err());
        assert!(Binomial::new(10, 1.1).is_err());
    }
}
