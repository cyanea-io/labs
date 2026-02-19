//! Substitution models for maximum likelihood phylogenetics.
//!
//! Provides the JC69 (Jukes-Cantor 1969) nucleotide substitution model,
//! which assumes equal base frequencies and a single substitution rate.

/// Map a nucleotide byte to an index (A=0, C=1, G=2, T=3).
///
/// Accepts both upper and lower case. Returns `None` for non-standard bases.
pub fn nucleotide_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

/// JC69 transition probability matrix for a given branch length `t`.
///
/// Under JC69, the probability of observing nucleotide `j` at time `t`
/// given nucleotide `i` at time 0 is:
///
/// - P(same) = 1/4 + 3/4 * e^{-4t/3}
/// - P(diff) = 1/4 - 1/4 * e^{-4t/3}
///
/// The matrix is 4x4, indexed by `nucleotide_index` (A=0, C=1, G=2, T=3).
pub fn jc69_probability(t: f64) -> [[f64; 4]; 4] {
    let e = (-4.0 * t / 3.0).exp();
    let p_same = 0.25 + 0.75 * e;
    let p_diff = 0.25 - 0.25 * e;

    [
        [p_same, p_diff, p_diff, p_diff],
        [p_diff, p_same, p_diff, p_diff],
        [p_diff, p_diff, p_same, p_diff],
        [p_diff, p_diff, p_diff, p_same],
    ]
}

/// Number of nucleotide states in DNA models.
pub const NUM_STATES: usize = 4;

/// Equilibrium frequency for JC69 (uniform: 1/4 for each base).
pub const JC69_FREQ: f64 = 0.25;

/// GTR (General Time-Reversible) model parameters.
///
/// The GTR model is parameterized by 6 exchangeability parameters and
/// 4 equilibrium base frequencies.
#[derive(Debug, Clone)]
pub struct GtrParams {
    /// Exchangeability rates: \[AC, AG, AT, CG, CT, GT\].
    pub rates: [f64; 6],
    /// Equilibrium frequencies: \[π_A, π_C, π_G, π_T\].
    pub freqs: [f64; 4],
}

impl GtrParams {
    /// Create new GTR parameters with validation.
    ///
    /// Rates must be positive; frequencies must be positive and sum to ~1.
    pub fn new(rates: [f64; 6], freqs: [f64; 4]) -> cyanea_core::Result<Self> {
        for (i, &r) in rates.iter().enumerate() {
            if r <= 0.0 {
                return Err(cyanea_core::CyaneaError::InvalidInput(format!(
                    "rate[{}] = {} must be positive",
                    i, r
                )));
            }
        }
        for (i, &f) in freqs.iter().enumerate() {
            if f <= 0.0 {
                return Err(cyanea_core::CyaneaError::InvalidInput(format!(
                    "freq[{}] = {} must be positive",
                    i, f
                )));
            }
        }
        let sum: f64 = freqs.iter().sum();
        if (sum - 1.0).abs() > 1e-6 {
            return Err(cyanea_core::CyaneaError::InvalidInput(format!(
                "frequencies sum to {} (expected ~1.0)",
                sum
            )));
        }
        Ok(Self { rates, freqs })
    }

    /// Compute the instantaneous rate matrix Q.
    ///
    /// Q\[i\]\[j\] = s_ij * π_j for i ≠ j, where s_ij is the exchangeability.
    /// Diagonal elements set so rows sum to zero.
    /// Normalized so that -Σ π_i Q_ii = 1 (expected substitutions per unit time).
    pub fn rate_matrix(&self) -> [[f64; 4]; 4] {
        let [r_ac, r_ag, r_at, r_cg, r_ct, r_gt] = self.rates;
        let [pi_a, pi_c, pi_g, pi_t] = self.freqs;

        // Exchangeability matrix (symmetric):
        // S = [[-, AC, AG, AT], [AC, -, CG, CT], [AG, CG, -, GT], [AT, CT, GT, -]]
        let s = [
            [0.0, r_ac, r_ag, r_at],
            [r_ac, 0.0, r_cg, r_ct],
            [r_ag, r_cg, 0.0, r_gt],
            [r_at, r_ct, r_gt, 0.0],
        ];
        let pi = [pi_a, pi_c, pi_g, pi_t];

        let mut q = [[0.0f64; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                if i != j {
                    q[i][j] = s[i][j] * pi[j];
                }
            }
            q[i][i] = -q[i].iter().sum::<f64>() + q[i][i]; // fix: -sum(off-diag)
        }
        // Fix diagonal: sum of off-diagonal in row i
        for i in 0..4 {
            let off_diag: f64 = (0..4).filter(|&j| j != i).map(|j| q[i][j]).sum();
            q[i][i] = -off_diag;
        }

        // Normalize so -Σ π_i Q_ii = 1
        let mu: f64 = (0..4).map(|i| -pi[i] * q[i][i]).sum();
        if mu > 0.0 {
            for i in 0..4 {
                for j in 0..4 {
                    q[i][j] /= mu;
                }
            }
        }

        q
    }

    /// Return a closure that computes P(t) = exp(Qt) for any branch length t.
    ///
    /// Uses eigendecomposition of the symmetrized rate matrix.
    pub fn probability_fn(&self) -> impl Fn(f64) -> [[f64; 4]; 4] {
        let q = self.rate_matrix();
        let pi = self.freqs;

        // Symmetrize: B = diag(√π) Q diag(1/√π)
        // B is real symmetric, so eigendecomposition is real.
        let sqrt_pi: [f64; 4] = [pi[0].sqrt(), pi[1].sqrt(), pi[2].sqrt(), pi[3].sqrt()];
        let inv_sqrt_pi: [f64; 4] = [
            1.0 / sqrt_pi[0],
            1.0 / sqrt_pi[1],
            1.0 / sqrt_pi[2],
            1.0 / sqrt_pi[3],
        ];

        let mut b = [[0.0f64; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                b[i][j] = sqrt_pi[i] * q[i][j] * inv_sqrt_pi[j];
            }
        }

        // Eigendecompose 4×4 symmetric matrix with Jacobi iteration.
        let (eigenvalues, eigenvectors) = jacobi_eigen_4x4(b);

        move |t: f64| {
            // P(t) = diag(1/√π) U exp(Λt) U^T diag(√π)
            let mut p = [[0.0f64; 4]; 4];
            for i in 0..4 {
                for j in 0..4 {
                    let mut sum = 0.0;
                    for k in 0..4 {
                        sum += eigenvectors[i][k]
                            * (eigenvalues[k] * t).exp()
                            * eigenvectors[j][k];
                    }
                    p[i][j] = inv_sqrt_pi[i] * sum * sqrt_pi[j];
                }
            }
            // Clamp small negative values from numerical error.
            for i in 0..4 {
                for j in 0..4 {
                    if p[i][j] < 0.0 {
                        p[i][j] = 0.0;
                    }
                }
            }
            p
        }
    }
}

/// Convenience: compute GTR transition probability matrix for branch length t.
pub fn gtr_probability(params: &GtrParams, t: f64) -> [[f64; 4]; 4] {
    params.probability_fn()(t)
}

/// Construct GTR parameters for the HKY85 model.
///
/// HKY85 has a single transition/transversion ratio κ:
/// transitions (A↔G, C↔T) have rate κ, transversions have rate 1.
pub fn hky85_params(kappa: f64, freqs: [f64; 4]) -> cyanea_core::Result<GtrParams> {
    if kappa <= 0.0 {
        return Err(cyanea_core::CyaneaError::InvalidInput(
            "kappa must be positive".into(),
        ));
    }
    // GTR rates: [AC, AG, AT, CG, CT, GT]
    // Transitions: AG (index 1) and CT (index 4)
    // Transversions: AC (0), AT (2), CG (3), GT (5)
    let rates = [1.0, kappa, 1.0, 1.0, kappa, 1.0];
    GtrParams::new(rates, freqs)
}

/// Discrete gamma rate categories (Yang 1994).
#[derive(Debug, Clone)]
pub struct GammaRates {
    /// Shape parameter alpha (smaller = more rate variation).
    pub alpha: f64,
    /// Number of rate categories.
    pub n_categories: usize,
}

impl GammaRates {
    /// Create new gamma rate parameters.
    pub fn new(alpha: f64, n_categories: usize) -> cyanea_core::Result<Self> {
        if alpha <= 0.0 {
            return Err(cyanea_core::CyaneaError::InvalidInput(
                "alpha must be positive".into(),
            ));
        }
        if n_categories == 0 {
            return Err(cyanea_core::CyaneaError::InvalidInput(
                "n_categories must be > 0".into(),
            ));
        }
        Ok(Self {
            alpha,
            n_categories,
        })
    }

    /// Compute discrete gamma category rates (Yang 1994 quantile method).
    ///
    /// Returns `n_categories` rates that average to 1.0.
    pub fn category_rates(&self) -> Vec<f64> {
        let k = self.n_categories;
        let a = self.alpha;

        // Compute quantile boundaries for k equal-probability categories.
        // Category i has probability mass in [i/k, (i+1)/k].
        // Rate for category i = mean of Gamma(a,a) over that interval.
        // Mean = a * [F^{-1}_{a+1}((i+1)/k) - F^{-1}_{a+1}(i/k)] * k
        // where F^{-1}_{a+1} is the inverse CDF of Gamma(a+1, 1).
        //
        // Actually, the mean of the gamma distribution in quantile interval
        // [p_lo, p_hi] is:
        //   (a / (p_hi - p_lo)) * [G_{a+1}(gamma_inv(p_hi, a)) - G_{a+1}(gamma_inv(p_lo, a))]
        //
        // Simpler approach: use quantile midpoints.
        // rate_i = gamma_inv((2i+1)/(2k), a, a)
        // Then normalize so rates average to 1.

        let mut rates = Vec::with_capacity(k);
        for i in 0..k {
            let p = (2 * i + 1) as f64 / (2 * k) as f64;
            let q = gamma_quantile(p, a);
            // Gamma(a, a) quantile = Gamma(a, 1) quantile / a
            // Actually for Gamma(alpha, beta=alpha), the rate = quantile/1.0
            // since we want mean=1. Gamma(a,a) has mean=1.
            // q is quantile of Gamma(a,1) with mean=a, so divide by a.
            rates.push(q / a);
        }

        // Normalize to ensure average = 1.0.
        let mean: f64 = rates.iter().sum::<f64>() / k as f64;
        if mean > 0.0 {
            for r in &mut rates {
                *r /= mean;
            }
        }

        rates
    }
}

/// Jacobi eigendecomposition for a 4×4 real symmetric matrix.
///
/// Returns (eigenvalues, eigenvectors) where eigenvectors[i][k] is the
/// i-th component of the k-th eigenvector.
fn jacobi_eigen_4x4(mut a: [[f64; 4]; 4]) -> ([f64; 4], [[f64; 4]; 4]) {
    let mut v = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];

    for _ in 0..100 {
        // Find largest off-diagonal element.
        let mut max_val = 0.0f64;
        let mut p = 0;
        let mut q = 1;
        for i in 0..4 {
            for j in (i + 1)..4 {
                if a[i][j].abs() > max_val {
                    max_val = a[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }
        if max_val < 1e-15 {
            break;
        }

        // Compute rotation angle.
        let theta = if (a[p][p] - a[q][q]).abs() < 1e-30 {
            std::f64::consts::FRAC_PI_4
        } else {
            0.5 * (2.0 * a[p][q] / (a[p][p] - a[q][q])).atan()
        };
        let c = theta.cos();
        let s = theta.sin();

        // Apply Givens rotation.
        let mut new_a = a;
        for i in 0..4 {
            new_a[i][p] = c * a[i][p] + s * a[i][q];
            new_a[i][q] = -s * a[i][p] + c * a[i][q];
        }
        // Copy columns, then fix rows.
        let tmp = new_a;
        for j in 0..4 {
            new_a[p][j] = c * tmp[p][j] + s * tmp[q][j];
            new_a[q][j] = -s * tmp[p][j] + c * tmp[q][j];
        }
        // Ensure symmetry of rotated elements.
        new_a[p][q] = 0.0;
        new_a[q][p] = 0.0;
        a = new_a;

        // Update eigenvectors.
        let mut new_v = v;
        for i in 0..4 {
            new_v[i][p] = c * v[i][p] + s * v[i][q];
            new_v[i][q] = -s * v[i][p] + c * v[i][q];
        }
        v = new_v;
    }

    let eigenvalues = [a[0][0], a[1][1], a[2][2], a[3][3]];
    (eigenvalues, v)
}

/// Regularized lower incomplete gamma function P(a, x) = γ(a,x) / Γ(a).
pub(crate) fn gamma_regularized(a: f64, x: f64) -> f64 {
    if x < 0.0 {
        return 0.0;
    }
    if x == 0.0 {
        return 0.0;
    }

    if x < a + 1.0 {
        // Series representation.
        gamma_series(a, x)
    } else {
        // Continued fraction representation.
        1.0 - gamma_cf(a, x)
    }
}

/// Series expansion for lower incomplete gamma.
fn gamma_series(a: f64, x: f64) -> f64 {
    let ln_gamma_a = ln_gamma(a);
    let mut sum = 1.0 / a;
    let mut term = 1.0 / a;
    for n in 1..200 {
        term *= x / (a + n as f64);
        sum += term;
        if term.abs() < sum.abs() * 1e-14 {
            break;
        }
    }
    sum * (-x + a * x.ln() - ln_gamma_a).exp()
}

/// Continued fraction for upper incomplete gamma.
fn gamma_cf(a: f64, x: f64) -> f64 {
    let ln_gamma_a = ln_gamma(a);
    let mut c = 1e-30f64;
    let mut d = 1.0 / (x + 1.0 - a);
    let mut f = d;

    for n in 1..200 {
        let an = -(n as f64) * (n as f64 - a);
        let bn = x + 2.0 * n as f64 + 1.0 - a;
        d = bn + an * d;
        if d.abs() < 1e-30 {
            d = 1e-30;
        }
        c = bn + an / c;
        if c.abs() < 1e-30 {
            c = 1e-30;
        }
        d = 1.0 / d;
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).abs() < 1e-14 {
            break;
        }
    }
    f * (-x + a * x.ln() - ln_gamma_a).exp()
}

/// Log-gamma function (Lanczos approximation).
pub(crate) fn ln_gamma(x: f64) -> f64 {
    let coeffs = [
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5,
    ];
    let y = x;
    let tmp = x + 5.5;
    let tmp = tmp - (x + 0.5) * tmp.ln();
    let mut ser = 1.000000000190015;
    for (j, &c) in coeffs.iter().enumerate() {
        ser += c / (y + 1.0 + j as f64);
    }
    -tmp + (2.5066282746310005 * ser / x).ln()
}

/// Inverse gamma CDF (quantile function) via Newton-Raphson.
///
/// Returns x such that P(a, x) ≈ p, for Gamma(a, 1).
fn gamma_quantile(p: f64, a: f64) -> f64 {
    if p <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return f64::MAX;
    }

    // Initial guess using Wilson-Hilferty approximation.
    let mut x = if a > 1.0 {
        let nu = (2.0 / (9.0 * a)).sqrt();
        let t = normal_quantile(p);
        let y = 1.0 - 1.0 / (9.0 * a) + t * nu;
        if y > 0.0 {
            a * y * y * y
        } else {
            a * 0.01
        }
    } else {
        // For small a, use a different initial guess.
        let t = 1.0 - a * (0.253 + a * 0.12);
        if p < t {
            (p / t).powf(1.0 / a)
        } else {
            1.0 - (1.0 - (p - t) / (1.0 - t)).ln()
        }
    };

    if x < 1e-15 {
        x = 1e-15;
    }

    // Newton-Raphson refinement.
    let ln_gamma_a = ln_gamma(a);
    for _ in 0..30 {
        let cdf = gamma_regularized(a, x);
        let pdf = ((a - 1.0) * x.ln() - x - ln_gamma_a).exp();
        if pdf < 1e-30 {
            break;
        }
        let delta = (cdf - p) / pdf;
        x -= delta;
        if x <= 0.0 {
            x = 1e-15;
        }
        if delta.abs() < x * 1e-12 {
            break;
        }
    }
    x
}

/// Approximate inverse normal CDF (Abramowitz & Stegun 26.2.23).
fn normal_quantile(p: f64) -> f64 {
    if p <= 0.0 {
        return -6.0;
    }
    if p >= 1.0 {
        return 6.0;
    }
    if (p - 0.5).abs() < 1e-15 {
        return 0.0;
    }
    let sign;
    let pp;
    if p < 0.5 {
        pp = p;
        sign = -1.0;
    } else {
        pp = 1.0 - p;
        sign = 1.0;
    };
    let t = (-2.0 * pp.ln()).sqrt();
    let c0 = 2.515517;
    let c1 = 0.802853;
    let c2 = 0.010328;
    let d1 = 1.432788;
    let d2 = 0.189269;
    let d3 = 0.001308;
    let x = t - (c0 + c1 * t + c2 * t * t) / (1.0 + d1 * t + d2 * t * t + d3 * t * t * t);
    sign * x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nucleotide_index_standard() {
        assert_eq!(nucleotide_index(b'A'), Some(0));
        assert_eq!(nucleotide_index(b'C'), Some(1));
        assert_eq!(nucleotide_index(b'G'), Some(2));
        assert_eq!(nucleotide_index(b'T'), Some(3));
    }

    #[test]
    fn nucleotide_index_lowercase() {
        assert_eq!(nucleotide_index(b'a'), Some(0));
        assert_eq!(nucleotide_index(b'c'), Some(1));
        assert_eq!(nucleotide_index(b'g'), Some(2));
        assert_eq!(nucleotide_index(b't'), Some(3));
    }

    #[test]
    fn nucleotide_index_uracil() {
        assert_eq!(nucleotide_index(b'U'), Some(3));
        assert_eq!(nucleotide_index(b'u'), Some(3));
    }

    #[test]
    fn nucleotide_index_invalid() {
        assert_eq!(nucleotide_index(b'N'), None);
        assert_eq!(nucleotide_index(b'-'), None);
        assert_eq!(nucleotide_index(b'X'), None);
    }

    #[test]
    fn jc69_t_zero_is_identity() {
        let p = jc69_probability(0.0);
        for i in 0..4 {
            for j in 0..4 {
                if i == j {
                    assert!(
                        (p[i][j] - 1.0).abs() < 1e-12,
                        "P[{}][{}] = {}, expected 1.0",
                        i,
                        j,
                        p[i][j]
                    );
                } else {
                    assert!(
                        p[i][j].abs() < 1e-12,
                        "P[{}][{}] = {}, expected 0.0",
                        i,
                        j,
                        p[i][j]
                    );
                }
            }
        }
    }

    #[test]
    fn jc69_large_t_approaches_uniform() {
        let p = jc69_probability(1000.0);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (p[i][j] - 0.25).abs() < 1e-10,
                    "P[{}][{}] = {}, expected ~0.25",
                    i,
                    j,
                    p[i][j]
                );
            }
        }
    }

    #[test]
    fn jc69_rows_sum_to_one() {
        for &t in &[0.0, 0.01, 0.1, 0.5, 1.0, 5.0, 100.0] {
            let p = jc69_probability(t);
            for i in 0..4 {
                let row_sum: f64 = p[i].iter().sum();
                assert!(
                    (row_sum - 1.0).abs() < 1e-12,
                    "row {} sum = {} at t = {}",
                    i,
                    row_sum,
                    t
                );
            }
        }
    }

    #[test]
    fn jc69_symmetric() {
        let p = jc69_probability(0.3);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (p[i][j] - p[j][i]).abs() < 1e-12,
                    "P[{}][{}] != P[{}][{}]",
                    i,
                    j,
                    j,
                    i
                );
            }
        }
    }

    #[test]
    fn jc69_probabilities_non_negative() {
        for &t in &[0.0, 0.001, 0.01, 0.1, 1.0, 10.0] {
            let p = jc69_probability(t);
            for i in 0..4 {
                for j in 0..4 {
                    assert!(
                        p[i][j] >= 0.0,
                        "P[{}][{}] = {} < 0 at t = {}",
                        i,
                        j,
                        p[i][j],
                        t
                    );
                }
            }
        }
    }

    #[test]
    fn jc69_known_value() {
        // At t = 0.3: e^{-4*0.3/3} = e^{-0.4}
        let e = (-0.4_f64).exp();
        let p_same = 0.25 + 0.75 * e;
        let p_diff = 0.25 - 0.25 * e;
        let p = jc69_probability(0.3);
        assert!((p[0][0] - p_same).abs() < 1e-12);
        assert!((p[0][1] - p_diff).abs() < 1e-12);
    }

    #[test]
    fn gtr_equal_rates_matches_jc69() {
        // GTR with equal rates and uniform frequencies should match JC69.
        let params =
            GtrParams::new([1.0; 6], [0.25; 4]).unwrap();
        let prob_fn = params.probability_fn();
        for &t in &[0.01, 0.1, 0.5, 1.0] {
            let gtr_p = prob_fn(t);
            let jc_p = jc69_probability(t);
            for i in 0..4 {
                for j in 0..4 {
                    assert!(
                        (gtr_p[i][j] - jc_p[i][j]).abs() < 1e-8,
                        "GTR({},{}) = {}, JC69 = {} at t = {}",
                        i,
                        j,
                        gtr_p[i][j],
                        jc_p[i][j],
                        t
                    );
                }
            }
        }
    }

    #[test]
    fn gtr_rows_sum_to_one() {
        let params =
            GtrParams::new([1.0, 2.0, 1.0, 1.0, 2.0, 1.0], [0.3, 0.2, 0.2, 0.3]).unwrap();
        let prob_fn = params.probability_fn();
        for &t in &[0.0, 0.01, 0.1, 0.5, 1.0, 5.0] {
            let p = prob_fn(t);
            for i in 0..4 {
                let row_sum: f64 = p[i].iter().sum();
                assert!(
                    (row_sum - 1.0).abs() < 1e-6,
                    "row {} sum = {} at t = {}",
                    i,
                    row_sum,
                    t
                );
            }
        }
    }

    #[test]
    fn gtr_t_zero_is_identity() {
        let params =
            GtrParams::new([1.0, 3.0, 1.0, 1.0, 3.0, 1.0], [0.3, 0.2, 0.2, 0.3]).unwrap();
        let p = gtr_probability(&params, 0.0);
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (p[i][j] - expected).abs() < 1e-8,
                    "P[{}][{}] = {}, expected {}",
                    i,
                    j,
                    p[i][j],
                    expected
                );
            }
        }
    }

    #[test]
    fn gamma_rates_average_to_one() {
        for &alpha in &[0.1, 0.5, 1.0, 2.0, 5.0, 10.0] {
            let gamma = GammaRates::new(alpha, 4).unwrap();
            let rates = gamma.category_rates();
            assert_eq!(rates.len(), 4);
            let mean: f64 = rates.iter().sum::<f64>() / 4.0;
            assert!(
                (mean - 1.0).abs() < 1e-6,
                "gamma rates mean = {} at alpha = {}",
                mean,
                alpha
            );
            // All rates should be positive.
            for &r in &rates {
                assert!(r > 0.0, "negative rate {} at alpha = {}", r, alpha);
            }
        }
    }

    #[test]
    fn hky85_as_gtr() {
        let params = hky85_params(2.0, [0.3, 0.2, 0.2, 0.3]).unwrap();
        let p = gtr_probability(&params, 0.1);
        // Transitions (A↔G) should be more likely than transversions (A↔C).
        assert!(
            p[0][2] > p[0][1],
            "transition A→G ({}) should exceed transversion A→C ({})",
            p[0][2],
            p[0][1]
        );
    }
}
