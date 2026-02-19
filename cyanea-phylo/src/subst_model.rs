//! Generic substitution model trait for N-state phylogenetic models.
//!
//! Provides a unified interface for nucleotide (4-state) and protein (20-state)
//! substitution models, wrapping existing model implementations.

use crate::models::{self, GtrParams};

/// Trait for substitution models of any state count.
///
/// Implementors provide rate matrices and transition probability matrices
/// for a given number of character states (e.g., 4 for DNA, 20 for protein).
pub trait SubstitutionModel {
    /// Number of character states (4 for DNA, 20 for protein).
    fn n_states(&self) -> usize;

    /// Equilibrium frequencies, one per state.
    fn frequencies(&self) -> &[f64];

    /// Instantaneous rate matrix Q (N x N).
    fn rate_matrix(&self) -> Vec<Vec<f64>>;

    /// Transition probability matrix P(t) = exp(Qt) for branch length t.
    fn transition_probs(&self, t: f64) -> Vec<Vec<f64>>;

    /// Number of free parameters in the model (excludes branch lengths).
    fn n_free_params(&self) -> usize;
}

/// JC69 model wrapper implementing [`SubstitutionModel`].
pub struct Jc69Model {
    freqs: [f64; 4],
}

impl Jc69Model {
    pub fn new() -> Self {
        Self { freqs: [0.25; 4] }
    }
}

impl Default for Jc69Model {
    fn default() -> Self {
        Self::new()
    }
}

impl SubstitutionModel for Jc69Model {
    fn n_states(&self) -> usize {
        4
    }

    fn frequencies(&self) -> &[f64] {
        &self.freqs
    }

    fn rate_matrix(&self) -> Vec<Vec<f64>> {
        // JC69 rate matrix: off-diag = 1/3, diag = -1 (normalized to mean rate 1)
        let mut q = vec![vec![0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                if i != j {
                    q[i][j] = 1.0 / 3.0;
                }
            }
            q[i][i] = -1.0;
        }
        q
    }

    fn transition_probs(&self, t: f64) -> Vec<Vec<f64>> {
        let p44 = models::jc69_probability(t);
        p44.iter().map(|row| row.to_vec()).collect()
    }

    fn n_free_params(&self) -> usize {
        0
    }
}

/// HKY85 model wrapper implementing [`SubstitutionModel`].
pub struct Hky85Model {
    pub kappa: f64,
    pub freqs: [f64; 4],
    gtr: GtrParams,
}

impl Hky85Model {
    pub fn new(kappa: f64, freqs: [f64; 4]) -> cyanea_core::Result<Self> {
        let gtr = models::hky85_params(kappa, freqs)?;
        Ok(Self { kappa, freqs, gtr })
    }
}

impl SubstitutionModel for Hky85Model {
    fn n_states(&self) -> usize {
        4
    }

    fn frequencies(&self) -> &[f64] {
        &self.freqs
    }

    fn rate_matrix(&self) -> Vec<Vec<f64>> {
        let q44 = self.gtr.rate_matrix();
        q44.iter().map(|row| row.to_vec()).collect()
    }

    fn transition_probs(&self, t: f64) -> Vec<Vec<f64>> {
        let p44 = models::gtr_probability(&self.gtr, t);
        p44.iter().map(|row| row.to_vec()).collect()
    }

    fn n_free_params(&self) -> usize {
        // kappa + 3 free frequencies (4 freqs summing to 1)
        4
    }
}

/// GTR model wrapper implementing [`SubstitutionModel`].
pub struct GtrModel {
    pub params: GtrParams,
}

impl GtrModel {
    pub fn new(params: GtrParams) -> Self {
        Self { params }
    }
}

impl SubstitutionModel for GtrModel {
    fn n_states(&self) -> usize {
        4
    }

    fn frequencies(&self) -> &[f64] {
        &self.params.freqs
    }

    fn rate_matrix(&self) -> Vec<Vec<f64>> {
        let q44 = self.params.rate_matrix();
        q44.iter().map(|row| row.to_vec()).collect()
    }

    fn transition_probs(&self, t: f64) -> Vec<Vec<f64>> {
        let p44 = models::gtr_probability(&self.params, t);
        p44.iter().map(|row| row.to_vec()).collect()
    }

    fn n_free_params(&self) -> usize {
        // 5 rate ratios (6 rates, 1 fixed) + 3 free frequencies
        8
    }
}

/// Generalized Jacobi eigendecomposition for N x N real symmetric matrices.
///
/// Returns (eigenvalues, eigenvectors) where eigenvectors\[i\]\[k\] is the
/// i-th component of the k-th eigenvector.
pub(crate) fn eigen_decompose(
    matrix: &[Vec<f64>],
) -> (Vec<f64>, Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = matrix.len();
    let mut a: Vec<Vec<f64>> = matrix.to_vec();

    // Identity matrix for eigenvectors.
    let mut v: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = vec![0.0; n];
            row[i] = 1.0;
            row
        })
        .collect();

    for _ in 0..200 {
        // Find largest off-diagonal element.
        let mut max_val = 0.0f64;
        let mut p = 0;
        let mut q = 1;
        for i in 0..n {
            for j in (i + 1)..n {
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

        // Apply Givens rotation to columns p and q.
        let mut new_a = a.clone();
        for i in 0..n {
            new_a[i][p] = c * a[i][p] + s * a[i][q];
            new_a[i][q] = -s * a[i][p] + c * a[i][q];
        }
        let tmp = new_a.clone();
        for j in 0..n {
            new_a[p][j] = c * tmp[p][j] + s * tmp[q][j];
            new_a[q][j] = -s * tmp[p][j] + c * tmp[q][j];
        }
        new_a[p][q] = 0.0;
        new_a[q][p] = 0.0;
        a = new_a;

        // Update eigenvectors.
        let mut new_v = v.clone();
        for i in 0..n {
            new_v[i][p] = c * v[i][p] + s * v[i][q];
            new_v[i][q] = -s * v[i][p] + c * v[i][q];
        }
        v = new_v;
    }

    let eigenvalues: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    // For symmetric real matrices, left and right eigenvectors are the same (transpose).
    let v_inv = transpose(&v);
    (eigenvalues, v, v_inv)
}

/// Compute transition probabilities P(t) = exp(Qt) for a general rate matrix Q.
///
/// Uses eigendecomposition of the symmetrized rate matrix:
/// B = diag(sqrt(pi)) * Q * diag(1/sqrt(pi)), then
/// P(t) = diag(1/sqrt(pi)) * U * exp(Lambda*t) * U^T * diag(sqrt(pi))
pub(crate) fn transition_probs_eigen(
    q: &[Vec<f64>],
    freqs: &[f64],
    t: f64,
) -> Vec<Vec<f64>> {
    let n = q.len();

    let sqrt_pi: Vec<f64> = freqs.iter().map(|&f| f.sqrt()).collect();
    let inv_sqrt_pi: Vec<f64> = sqrt_pi.iter().map(|&sp| 1.0 / sp).collect();

    // Symmetrize: B = diag(sqrt(pi)) * Q * diag(1/sqrt(pi))
    let mut b = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            b[i][j] = sqrt_pi[i] * q[i][j] * inv_sqrt_pi[j];
        }
    }

    let (eigenvalues, eigenvectors, _) = eigen_decompose(&b);

    // P(t) = diag(1/sqrt(pi)) * U * exp(Lambda*t) * U^T * diag(sqrt(pi))
    let mut p = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut sum = 0.0;
            for k in 0..n {
                sum += eigenvectors[i][k]
                    * (eigenvalues[k] * t).exp()
                    * eigenvectors[j][k];
            }
            p[i][j] = inv_sqrt_pi[i] * sum * sqrt_pi[j];
            if p[i][j] < 0.0 {
                p[i][j] = 0.0;
            }
        }
    }
    p
}

fn transpose(m: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if m.is_empty() {
        return Vec::new();
    }
    let n = m.len();
    let cols = m[0].len();
    let mut t = vec![vec![0.0; n]; cols];
    for i in 0..n {
        for j in 0..cols {
            t[j][i] = m[i][j];
        }
    }
    t
}

/// Build a normalized rate matrix Q from exchangeability matrix S and frequencies pi.
///
/// Q\[i\]\[j\] = S\[i\]\[j\] * pi\[j\] for i != j, rows sum to 0,
/// normalized so -sum(pi_i * Q_ii) = 1.
pub(crate) fn build_rate_matrix(
    exchangeabilities: &[Vec<f64>],
    freqs: &[f64],
) -> Vec<Vec<f64>> {
    let n = freqs.len();
    let mut q = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in 0..n {
            if i != j {
                q[i][j] = exchangeabilities[i][j] * freqs[j];
            }
        }
        let off_diag: f64 = (0..n).filter(|&j| j != i).map(|j| q[i][j]).sum();
        q[i][i] = -off_diag;
    }

    // Normalize so -sum(pi_i * Q_ii) = 1
    let mu: f64 = (0..n).map(|i| -freqs[i] * q[i][i]).sum();
    if mu > 0.0 {
        for i in 0..n {
            for j in 0..n {
                q[i][j] /= mu;
            }
        }
    }

    q
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jc69_trait_matches_original() {
        let model = Jc69Model::new();
        for &t in &[0.01, 0.1, 0.5, 1.0] {
            let trait_p = model.transition_probs(t);
            let orig_p = models::jc69_probability(t);
            for i in 0..4 {
                for j in 0..4 {
                    assert!(
                        (trait_p[i][j] - orig_p[i][j]).abs() < 1e-12,
                        "JC69 mismatch at ({},{}) t={}: {} vs {}",
                        i, j, t, trait_p[i][j], orig_p[i][j]
                    );
                }
            }
        }
    }

    #[test]
    fn gtr_trait_matches_original() {
        let params = GtrParams::new(
            [1.0, 2.0, 1.0, 1.0, 2.0, 1.0],
            [0.3, 0.2, 0.2, 0.3],
        )
        .unwrap();
        let model = GtrModel::new(params.clone());
        for &t in &[0.01, 0.1, 0.5, 1.0] {
            let trait_p = model.transition_probs(t);
            let orig_p = models::gtr_probability(&params, t);
            for i in 0..4 {
                for j in 0..4 {
                    assert!(
                        (trait_p[i][j] - orig_p[i][j]).abs() < 1e-6,
                        "GTR mismatch at ({},{}) t={}: {} vs {}",
                        i, j, t, trait_p[i][j], orig_p[i][j]
                    );
                }
            }
        }
    }

    #[test]
    fn n_states_returns_4_for_nucleotide() {
        assert_eq!(Jc69Model::new().n_states(), 4);
        assert_eq!(
            Hky85Model::new(2.0, [0.25; 4]).unwrap().n_states(),
            4
        );
        let params = GtrParams::new([1.0; 6], [0.25; 4]).unwrap();
        assert_eq!(GtrModel::new(params).n_states(), 4);
    }

    #[test]
    fn frequencies_sum_to_one() {
        let models: Vec<Box<dyn SubstitutionModel>> = vec![
            Box::new(Jc69Model::new()),
            Box::new(Hky85Model::new(2.0, [0.3, 0.2, 0.2, 0.3]).unwrap()),
            Box::new(GtrModel::new(
                GtrParams::new([1.0, 2.0, 1.0, 1.0, 2.0, 1.0], [0.3, 0.2, 0.2, 0.3])
                    .unwrap(),
            )),
        ];
        for m in &models {
            let sum: f64 = m.frequencies().iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-10,
                "frequencies sum to {}",
                sum
            );
        }
    }

    #[test]
    fn p_zero_is_identity() {
        let models: Vec<Box<dyn SubstitutionModel>> = vec![
            Box::new(Jc69Model::new()),
            Box::new(Hky85Model::new(2.0, [0.3, 0.2, 0.2, 0.3]).unwrap()),
        ];
        for m in &models {
            let p = m.transition_probs(0.0);
            let n = m.n_states();
            for i in 0..n {
                for j in 0..n {
                    let expected = if i == j { 1.0 } else { 0.0 };
                    assert!(
                        (p[i][j] - expected).abs() < 1e-8,
                        "P(0)[{}][{}] = {}, expected {}",
                        i, j, p[i][j], expected
                    );
                }
            }
        }
    }

    #[test]
    fn rows_sum_to_one() {
        let model = Hky85Model::new(3.0, [0.3, 0.2, 0.2, 0.3]).unwrap();
        for &t in &[0.01, 0.1, 0.5, 1.0, 5.0] {
            let p = model.transition_probs(t);
            for (i, row) in p.iter().enumerate() {
                let sum: f64 = row.iter().sum();
                assert!(
                    (sum - 1.0).abs() < 1e-6,
                    "row {} sum = {} at t = {}",
                    i, sum, t
                );
            }
        }
    }

    #[test]
    fn eigen_decompose_identity() {
        let m = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 3.0],
        ];
        let (vals, _, _) = eigen_decompose(&m);
        let mut sorted = vals.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((sorted[0] - 1.0).abs() < 1e-10);
        assert!((sorted[1] - 2.0).abs() < 1e-10);
        assert!((sorted[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn n_free_params_correct() {
        assert_eq!(Jc69Model::new().n_free_params(), 0);
        assert_eq!(
            Hky85Model::new(2.0, [0.25; 4]).unwrap().n_free_params(),
            4
        );
        let params = GtrParams::new([1.0; 6], [0.25; 4]).unwrap();
        assert_eq!(GtrModel::new(params).n_free_params(), 8);
    }

    #[test]
    fn build_rate_matrix_rows_sum_to_zero() {
        let s = vec![
            vec![0.0, 1.0, 2.0, 1.0],
            vec![1.0, 0.0, 1.0, 2.0],
            vec![2.0, 1.0, 0.0, 1.0],
            vec![1.0, 2.0, 1.0, 0.0],
        ];
        let freqs = vec![0.3, 0.2, 0.2, 0.3];
        let q = build_rate_matrix(&s, &freqs);
        for (i, row) in q.iter().enumerate() {
            let sum: f64 = row.iter().sum();
            assert!(
                sum.abs() < 1e-10,
                "row {} sums to {}",
                i, sum
            );
        }
    }
}
