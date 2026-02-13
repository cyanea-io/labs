//! Bayesian conjugate prior distributions.
//!
//! Provides four conjugate prior–likelihood pairs commonly used in
//! bioinformatics:
//!
//! - [`Beta`] — conjugate to binomial likelihood
//! - [`Gamma`] — conjugate to Poisson likelihood
//! - [`NormalConjugate`] — normal prior for normal likelihood (known variance)
//! - [`Dirichlet`] — conjugate to multinomial likelihood

use crate::distribution::{betai, gammainc, ln_gamma, Distribution};
use cyanea_core::{CyaneaError, Result};

// ── Beta distribution ────────────────────────────────────────────────────

/// Beta distribution, conjugate prior for binomial likelihood.
///
/// After observing `s` successes in `n` trials, the posterior is
/// `Beta(α + s, β + n − s)`.
#[derive(Debug, Clone, Copy)]
pub struct Beta {
    alpha: f64,
    beta: f64,
}

impl Beta {
    /// Create a Beta distribution with shape parameters `alpha` and `beta`.
    ///
    /// Both must be positive.
    pub fn new(alpha: f64, beta: f64) -> Result<Self> {
        if alpha <= 0.0 || beta <= 0.0 {
            return Err(CyaneaError::InvalidInput(
                "Beta: alpha and beta must be positive".into(),
            ));
        }
        Ok(Self { alpha, beta })
    }

    /// Alpha parameter.
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    /// Beta parameter.
    pub fn beta(&self) -> f64 {
        self.beta
    }

    /// Compute the posterior after observing binomial data.
    pub fn update_binomial(&self, successes: u64, trials: u64) -> Self {
        Self {
            alpha: self.alpha + successes as f64,
            beta: self.beta + (trials - successes) as f64,
        }
    }
}

impl Distribution for Beta {
    fn pdf(&self, x: f64) -> f64 {
        if x <= 0.0 || x >= 1.0 {
            return 0.0;
        }
        let ln_beta_fn = ln_gamma(self.alpha) + ln_gamma(self.beta)
            - ln_gamma(self.alpha + self.beta);
        let ln_pdf = (self.alpha - 1.0) * x.ln()
            + (self.beta - 1.0) * (1.0 - x).ln()
            - ln_beta_fn;
        ln_pdf.exp()
    }

    fn cdf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        if x >= 1.0 {
            return 1.0;
        }
        betai(self.alpha, self.beta, x).unwrap_or(0.0)
    }

    fn mean(&self) -> f64 {
        self.alpha / (self.alpha + self.beta)
    }

    fn variance(&self) -> f64 {
        let ab = self.alpha + self.beta;
        (self.alpha * self.beta) / (ab * ab * (ab + 1.0))
    }
}

// ── Gamma distribution ───────────────────────────────────────────────────

/// Gamma distribution (shape/rate parameterization), conjugate prior for
/// Poisson likelihood.
///
/// After observing a count `c` (one observation), the posterior is
/// `Gamma(shape + c, rate + 1)`.
#[derive(Debug, Clone, Copy)]
pub struct Gamma {
    shape: f64,
    rate: f64,
}

impl Gamma {
    /// Create a Gamma distribution with given `shape` (α) and `rate` (β).
    ///
    /// Both must be positive.
    pub fn new(shape: f64, rate: f64) -> Result<Self> {
        if shape <= 0.0 || rate <= 0.0 {
            return Err(CyaneaError::InvalidInput(
                "Gamma: shape and rate must be positive".into(),
            ));
        }
        Ok(Self { shape, rate })
    }

    /// Shape parameter (α).
    pub fn shape(&self) -> f64 {
        self.shape
    }

    /// Rate parameter (β).
    pub fn rate(&self) -> f64 {
        self.rate
    }

    /// Posterior after observing a single Poisson count.
    pub fn update_poisson(&self, count: u64) -> Self {
        Self {
            shape: self.shape + count as f64,
            rate: self.rate + 1.0,
        }
    }

    /// Posterior after observing multiple Poisson counts.
    pub fn update_poisson_batch(&self, counts: &[u64]) -> Self {
        let total: u64 = counts.iter().sum();
        Self {
            shape: self.shape + total as f64,
            rate: self.rate + counts.len() as f64,
        }
    }
}

impl Distribution for Gamma {
    fn pdf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        let ln_pdf = self.shape * self.rate.ln() - ln_gamma(self.shape)
            + (self.shape - 1.0) * x.ln()
            - self.rate * x;
        ln_pdf.exp()
    }

    fn cdf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        // P(a, rate * x) but gammainc takes shape and x
        // For Gamma(shape, rate): CDF(x) = P(shape, rate * x) = gammainc(shape, rate * x)
        gammainc(self.shape, self.rate * x).unwrap_or(0.0)
    }

    fn mean(&self) -> f64 {
        self.shape / self.rate
    }

    fn variance(&self) -> f64 {
        self.shape / (self.rate * self.rate)
    }
}

// ── Normal conjugate ─────────────────────────────────────────────────────

/// Normal prior for a normal likelihood with known observation variance.
///
/// Uses the precision (inverse variance) formulation for numerically
/// stable Bayesian updates.
#[derive(Debug, Clone, Copy)]
pub struct NormalConjugate {
    prior_mu: f64,
    prior_var: f64,
    obs_var: f64,
}

impl NormalConjugate {
    /// Create a normal conjugate prior.
    ///
    /// - `prior_mu`: prior mean
    /// - `prior_var`: prior variance (must be positive)
    /// - `obs_var`: known observation variance (must be positive)
    pub fn new(prior_mu: f64, prior_var: f64, obs_var: f64) -> Result<Self> {
        if prior_var <= 0.0 {
            return Err(CyaneaError::InvalidInput(
                "NormalConjugate: prior_var must be positive".into(),
            ));
        }
        if obs_var <= 0.0 {
            return Err(CyaneaError::InvalidInput(
                "NormalConjugate: obs_var must be positive".into(),
            ));
        }
        Ok(Self {
            prior_mu,
            prior_var,
            obs_var,
        })
    }

    /// Update with a single observation.
    pub fn update(&self, observation: f64) -> Self {
        let prior_prec = 1.0 / self.prior_var;
        let obs_prec = 1.0 / self.obs_var;
        let post_prec = prior_prec + obs_prec;
        let post_var = 1.0 / post_prec;
        let post_mu = (prior_prec * self.prior_mu + obs_prec * observation) / post_prec;
        Self {
            prior_mu: post_mu,
            prior_var: post_var,
            obs_var: self.obs_var,
        }
    }

    /// Update with a batch of observations.
    pub fn update_batch(&self, observations: &[f64]) -> Self {
        let n = observations.len() as f64;
        if n == 0.0 {
            return *self;
        }
        let obs_mean: f64 = observations.iter().sum::<f64>() / n;
        let prior_prec = 1.0 / self.prior_var;
        let obs_prec = n / self.obs_var;
        let post_prec = prior_prec + obs_prec;
        let post_var = 1.0 / post_prec;
        let post_mu = (prior_prec * self.prior_mu + obs_prec * obs_mean) / post_prec;
        Self {
            prior_mu: post_mu,
            prior_var: post_var,
            obs_var: self.obs_var,
        }
    }

    /// Posterior mean.
    pub fn posterior_mean(&self) -> f64 {
        self.prior_mu
    }

    /// Posterior variance.
    pub fn posterior_variance(&self) -> f64 {
        self.prior_var
    }
}

// ── Dirichlet distribution ───────────────────────────────────────────────

/// Dirichlet distribution, conjugate prior for multinomial likelihood.
///
/// After observing counts `c₁, ..., cₖ`, the posterior is
/// `Dirichlet(α₁ + c₁, ..., αₖ + cₖ)`.
#[derive(Debug, Clone)]
pub struct Dirichlet {
    alpha: Vec<f64>,
}

impl Dirichlet {
    /// Create a Dirichlet distribution with concentration parameters `alpha`.
    ///
    /// All elements must be positive and the vector must have at least 2 elements.
    pub fn new(alpha: Vec<f64>) -> Result<Self> {
        if alpha.len() < 2 {
            return Err(CyaneaError::InvalidInput(
                "Dirichlet: need at least 2 categories".into(),
            ));
        }
        if alpha.iter().any(|&a| a <= 0.0) {
            return Err(CyaneaError::InvalidInput(
                "Dirichlet: all alpha values must be positive".into(),
            ));
        }
        Ok(Self { alpha })
    }

    /// Create a symmetric Dirichlet with `k` categories, each with
    /// concentration `alpha`.
    pub fn symmetric(k: usize, alpha: f64) -> Result<Self> {
        if k < 2 {
            return Err(CyaneaError::InvalidInput(
                "Dirichlet: need at least 2 categories".into(),
            ));
        }
        if alpha <= 0.0 {
            return Err(CyaneaError::InvalidInput(
                "Dirichlet: alpha must be positive".into(),
            ));
        }
        Ok(Self {
            alpha: vec![alpha; k],
        })
    }

    /// Concentration parameters.
    pub fn alpha(&self) -> &[f64] {
        &self.alpha
    }

    /// Posterior after observing multinomial counts.
    ///
    /// # Panics
    ///
    /// Panics if `counts.len() != self.alpha.len()`.
    pub fn update_multinomial(&self, counts: &[u64]) -> Self {
        assert_eq!(
            counts.len(),
            self.alpha.len(),
            "counts length must match alpha length"
        );
        Self {
            alpha: self
                .alpha
                .iter()
                .zip(counts.iter())
                .map(|(&a, &c)| a + c as f64)
                .collect(),
        }
    }

    /// Expected value (mean) of the Dirichlet: `E[Xᵢ] = αᵢ / Σα`.
    pub fn mean(&self) -> Vec<f64> {
        let sum: f64 = self.alpha.iter().sum();
        self.alpha.iter().map(|&a| a / sum).collect()
    }

    /// Variance of each component: `Var[Xᵢ] = αᵢ(Σα − αᵢ) / (Σα² (Σα + 1))`.
    pub fn variance(&self) -> Vec<f64> {
        let sum: f64 = self.alpha.iter().sum();
        let denom = sum * sum * (sum + 1.0);
        self.alpha.iter().map(|&a| a * (sum - a) / denom).collect()
    }

    /// Log-PDF of the Dirichlet at point `x`.
    ///
    /// # Errors
    ///
    /// Returns an error if `x.len() != alpha.len()` or if values don't
    /// sum to approximately 1.
    pub fn ln_pdf(&self, x: &[f64]) -> Result<f64> {
        if x.len() != self.alpha.len() {
            return Err(CyaneaError::InvalidInput(
                "Dirichlet::ln_pdf: x length must match alpha length".into(),
            ));
        }
        let sum: f64 = x.iter().sum();
        if (sum - 1.0).abs() > 1e-6 {
            return Err(CyaneaError::InvalidInput(
                "Dirichlet::ln_pdf: x must sum to 1".into(),
            ));
        }

        let alpha_sum: f64 = self.alpha.iter().sum();
        let mut ln_b = -ln_gamma(alpha_sum);
        for &a in &self.alpha {
            ln_b += ln_gamma(a);
        }

        let mut result = -ln_b;
        for (xi, &ai) in x.iter().zip(self.alpha.iter()) {
            if *xi <= 0.0 {
                return Err(CyaneaError::InvalidInput(
                    "Dirichlet::ln_pdf: all x values must be positive".into(),
                ));
            }
            result += (ai - 1.0) * xi.ln();
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // ── Beta tests ────────────────────────────────────────────────────

    #[test]
    fn beta_uniform_prior() {
        let prior = Beta::new(1.0, 1.0).unwrap();
        assert!((prior.mean() - 0.5).abs() < TOL);
    }

    #[test]
    fn beta_conjugacy() {
        // Beta(1,1) + 3 successes in 10 trials = Beta(4, 8)
        let prior = Beta::new(1.0, 1.0).unwrap();
        let post = prior.update_binomial(3, 10);
        assert!((post.alpha() - 4.0).abs() < TOL);
        assert!((post.beta() - 8.0).abs() < TOL);
        assert!((post.mean() - 4.0 / 12.0).abs() < TOL);
    }

    #[test]
    fn beta_pdf_at_mode() {
        // Beta(2, 5): mode = (2-1)/(2+5-2) = 1/5 = 0.2
        let b = Beta::new(2.0, 5.0).unwrap();
        let pdf_at_mode = b.pdf(0.2);
        // Should be near the maximum
        assert!(pdf_at_mode > b.pdf(0.1));
        assert!(pdf_at_mode > b.pdf(0.5));
    }

    #[test]
    fn beta_cdf_boundaries() {
        let b = Beta::new(2.0, 3.0).unwrap();
        assert_eq!(b.cdf(0.0), 0.0);
        assert!((b.cdf(1.0) - 1.0).abs() < TOL);
    }

    #[test]
    fn beta_cdf_midpoint() {
        // Beta(1,1) is uniform, so CDF(0.5) = 0.5
        let b = Beta::new(1.0, 1.0).unwrap();
        assert!((b.cdf(0.5) - 0.5).abs() < TOL);
    }

    #[test]
    fn beta_invalid() {
        assert!(Beta::new(0.0, 1.0).is_err());
        assert!(Beta::new(1.0, -1.0).is_err());
    }

    // ── Gamma tests ──────────────────────────────────────────────────

    #[test]
    fn gamma_mean_variance() {
        let g = Gamma::new(3.0, 2.0).unwrap();
        assert!((g.mean() - 1.5).abs() < TOL);
        assert!((g.variance() - 0.75).abs() < TOL);
    }

    #[test]
    fn gamma_conjugacy_poisson() {
        // Gamma(2, 1) + observe count=5 → Gamma(7, 2)
        let prior = Gamma::new(2.0, 1.0).unwrap();
        let post = prior.update_poisson(5);
        assert!((post.shape() - 7.0).abs() < TOL);
        assert!((post.rate() - 2.0).abs() < TOL);
    }

    #[test]
    fn gamma_conjugacy_batch() {
        // Gamma(2, 1) + observe [3, 5, 2] → Gamma(2+10, 1+3) = Gamma(12, 4)
        let prior = Gamma::new(2.0, 1.0).unwrap();
        let post = prior.update_poisson_batch(&[3, 5, 2]);
        assert!((post.shape() - 12.0).abs() < TOL);
        assert!((post.rate() - 4.0).abs() < TOL);
    }

    #[test]
    fn gamma_cdf() {
        // Gamma(1, 1) is Exponential(1): CDF(x) = 1 - e^{-x}
        let g = Gamma::new(1.0, 1.0).unwrap();
        let x = 2.0;
        let expected = 1.0 - (-x as f64).exp();
        assert!((g.cdf(x) - expected).abs() < 1e-8);
    }

    #[test]
    fn gamma_invalid() {
        assert!(Gamma::new(0.0, 1.0).is_err());
        assert!(Gamma::new(1.0, 0.0).is_err());
    }

    // ── NormalConjugate tests ────────────────────────────────────────

    #[test]
    fn normal_conjugate_single_update() {
        let prior = NormalConjugate::new(0.0, 1.0, 1.0).unwrap();
        let post = prior.update(2.0);
        // Posterior mean = (0/1 + 2/1) / (1/1 + 1/1) = 2/2 = 1.0
        assert!((post.posterior_mean() - 1.0).abs() < TOL);
        // Posterior variance = 1/(1+1) = 0.5
        assert!((post.posterior_variance() - 0.5).abs() < TOL);
    }

    #[test]
    fn normal_conjugate_batch_update() {
        let prior = NormalConjugate::new(0.0, 1.0, 1.0).unwrap();
        let post = prior.update_batch(&[2.0, 4.0]);
        // obs_mean = 3.0, n = 2
        // prior_prec = 1, obs_prec = 2/1 = 2
        // post_prec = 3, post_var = 1/3
        // post_mu = (0 + 2*3)/3 = 2.0
        assert!((post.posterior_mean() - 2.0).abs() < TOL);
        assert!((post.posterior_variance() - 1.0 / 3.0).abs() < TOL);
    }

    #[test]
    fn normal_conjugate_empty_batch() {
        let prior = NormalConjugate::new(5.0, 2.0, 1.0).unwrap();
        let post = prior.update_batch(&[]);
        assert!((post.posterior_mean() - 5.0).abs() < TOL);
        assert!((post.posterior_variance() - 2.0).abs() < TOL);
    }

    #[test]
    fn normal_conjugate_precision_shrinkage() {
        // With very precise prior, observation barely moves the mean
        let prior = NormalConjugate::new(0.0, 0.01, 100.0).unwrap();
        let post = prior.update(100.0);
        // prior_prec = 100, obs_prec = 0.01
        // post_mu ≈ 0.0 (barely shifts)
        assert!(post.posterior_mean().abs() < 0.02);
    }

    #[test]
    fn normal_conjugate_invalid() {
        assert!(NormalConjugate::new(0.0, 0.0, 1.0).is_err());
        assert!(NormalConjugate::new(0.0, 1.0, 0.0).is_err());
    }

    // ── Dirichlet tests ──────────────────────────────────────────────

    #[test]
    fn dirichlet_symmetric_mean() {
        let d = Dirichlet::symmetric(4, 1.0).unwrap();
        let mean = d.mean();
        assert_eq!(mean.len(), 4);
        for m in &mean {
            assert!((m - 0.25).abs() < TOL);
        }
    }

    #[test]
    fn dirichlet_conjugacy_multinomial() {
        let prior = Dirichlet::symmetric(3, 1.0).unwrap();
        let post = prior.update_multinomial(&[10, 5, 15]);
        let expected = [11.0, 6.0, 16.0];
        for (a, e) in post.alpha().iter().zip(expected.iter()) {
            assert!((a - e).abs() < TOL);
        }
    }

    #[test]
    fn dirichlet_ln_pdf() {
        // Dirichlet(1,1,1) is uniform on the simplex: ln_pdf = ln(Γ(3)/Γ(1)^3) = ln(2)
        let d = Dirichlet::symmetric(3, 1.0).unwrap();
        let ln_pdf = d.ln_pdf(&[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]).unwrap();
        assert!((ln_pdf - 2.0_f64.ln()).abs() < 1e-6);
    }

    #[test]
    fn dirichlet_invalid() {
        assert!(Dirichlet::new(vec![1.0]).is_err()); // too few
        assert!(Dirichlet::new(vec![1.0, -1.0]).is_err());
        assert!(Dirichlet::symmetric(1, 1.0).is_err());
        assert!(Dirichlet::symmetric(3, 0.0).is_err());
    }

    #[test]
    fn dirichlet_ln_pdf_invalid() {
        let d = Dirichlet::symmetric(3, 1.0).unwrap();
        assert!(d.ln_pdf(&[0.5, 0.5]).is_err()); // wrong length
        assert!(d.ln_pdf(&[0.5, 0.3, 0.1]).is_err()); // doesn't sum to 1
    }
}
