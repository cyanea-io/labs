//! Model inference: prediction primitives for bioinformatics ML.
//!
//! Provides simple, dependency-free models suitable for inference on
//! biological data:
//!
//! - **K-Nearest Neighbors (KNN)** — classification and regression via
//!   majority vote or mean of nearest neighbors
//! - **Linear regression** — ordinary least squares for single-target
//!   regression
//!
//! These implementations are designed for small-to-medium datasets typical
//! of bioinformatics workflows (e.g., expression matrices, k-mer profiles).

use cyanea_core::{CyaneaError, Result, Summarizable};

use crate::distance::{compute_distance, DistanceMetric};

// ---------------------------------------------------------------------------
// KNN
// ---------------------------------------------------------------------------

/// Configuration for K-Nearest Neighbors.
#[derive(Debug, Clone)]
pub struct KnnConfig {
    /// Number of neighbors.
    pub k: usize,
    /// Distance metric to use.
    pub metric: DistanceMetric,
}

impl Default for KnnConfig {
    fn default() -> Self {
        Self {
            k: 5,
            metric: DistanceMetric::Euclidean,
        }
    }
}

/// A fitted KNN model that stores training data for lazy evaluation.
#[derive(Debug, Clone)]
pub struct KnnModel {
    /// Training data, row-major: `n_samples × n_features`.
    data: Vec<f64>,
    n_features: usize,
    n_samples: usize,
    config: KnnConfig,
}

impl KnnModel {
    /// Fit a KNN model from flat row-major training data.
    ///
    /// KNN is a lazy learner — this just stores the data and validates dimensions.
    ///
    /// # Errors
    ///
    /// Returns an error if the data is empty, dimensions are inconsistent,
    /// or k is 0.
    pub fn fit(data: &[f64], n_features: usize, config: KnnConfig) -> Result<Self> {
        if data.is_empty() {
            return Err(CyaneaError::InvalidInput("empty training data".into()));
        }
        if n_features == 0 {
            return Err(CyaneaError::InvalidInput("n_features must be > 0".into()));
        }
        if data.len() % n_features != 0 {
            return Err(CyaneaError::InvalidInput(format!(
                "data length {} not divisible by n_features {}",
                data.len(),
                n_features
            )));
        }
        let n_samples = data.len() / n_features;
        if config.k == 0 {
            return Err(CyaneaError::InvalidInput("k must be > 0".into()));
        }
        if config.k > n_samples {
            return Err(CyaneaError::InvalidInput(format!(
                "k ({}) > n_samples ({})",
                config.k, n_samples
            )));
        }
        Ok(Self {
            data: data.to_vec(),
            n_features,
            n_samples,
            config,
        })
    }

    /// Find the k nearest neighbors of a query point.
    ///
    /// Returns a vector of `(index, distance)` pairs sorted by distance.
    pub fn neighbors(&self, query: &[f64]) -> Result<Vec<(usize, f64)>> {
        if query.len() != self.n_features {
            return Err(CyaneaError::InvalidInput(format!(
                "query has {} features, expected {}",
                query.len(),
                self.n_features
            )));
        }

        let mut dists: Vec<(usize, f64)> = (0..self.n_samples)
            .map(|i| {
                let row = &self.data[i * self.n_features..(i + 1) * self.n_features];
                let d = compute_distance(row, query, self.config.metric).unwrap_or(f64::INFINITY);
                (i, d)
            })
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        dists.truncate(self.config.k);
        Ok(dists)
    }

    /// Classify a query point using majority vote among k neighbors.
    ///
    /// `labels` must have one entry per training sample.
    pub fn classify(&self, query: &[f64], labels: &[i32]) -> Result<i32> {
        if labels.len() != self.n_samples {
            return Err(CyaneaError::InvalidInput(format!(
                "labels length {} != n_samples {}",
                labels.len(),
                self.n_samples
            )));
        }
        let neighbors = self.neighbors(query)?;
        // Majority vote
        let mut votes = std::collections::HashMap::new();
        for &(idx, _) in &neighbors {
            *votes.entry(labels[idx]).or_insert(0usize) += 1;
        }
        votes
            .into_iter()
            .max_by_key(|&(_, count)| count)
            .map(|(label, _)| label)
            .ok_or_else(|| CyaneaError::Other("no neighbors found".into()))
    }

    /// Predict a continuous value as the mean of k neighbors' target values.
    ///
    /// `targets` must have one entry per training sample.
    pub fn regress(&self, query: &[f64], targets: &[f64]) -> Result<f64> {
        if targets.len() != self.n_samples {
            return Err(CyaneaError::InvalidInput(format!(
                "targets length {} != n_samples {}",
                targets.len(),
                self.n_samples
            )));
        }
        let neighbors = self.neighbors(query)?;
        let sum: f64 = neighbors.iter().map(|&(idx, _)| targets[idx]).sum();
        Ok(sum / neighbors.len() as f64)
    }

    /// Number of training samples.
    pub fn n_samples(&self) -> usize {
        self.n_samples
    }

    /// Number of features.
    pub fn n_features(&self) -> usize {
        self.n_features
    }
}

impl Summarizable for KnnModel {
    fn summary(&self) -> String {
        format!(
            "KNN: k={}, {} samples, {} features, {:?}",
            self.config.k, self.n_samples, self.n_features, self.config.metric,
        )
    }
}

// ---------------------------------------------------------------------------
// Linear Regression
// ---------------------------------------------------------------------------

/// A fitted linear regression model: y = X * weights + bias.
#[derive(Debug, Clone)]
pub struct LinearRegression {
    /// Feature weights (one per feature).
    pub weights: Vec<f64>,
    /// Bias (intercept) term.
    pub bias: f64,
    /// Number of features.
    pub n_features: usize,
    /// R-squared on training data.
    pub r_squared: f64,
}

impl LinearRegression {
    /// Fit a linear regression model via ordinary least squares.
    ///
    /// Uses the normal equations with a simple Gaussian elimination solver.
    /// Suitable for small-to-medium datasets where n_features is modest.
    ///
    /// `data` is flat row-major `n_samples × n_features`.
    /// `targets` has one value per sample.
    ///
    /// # Errors
    ///
    /// Returns an error if dimensions are inconsistent or the system is singular.
    pub fn fit(data: &[f64], n_features: usize, targets: &[f64]) -> Result<Self> {
        if data.is_empty() {
            return Err(CyaneaError::InvalidInput("empty training data".into()));
        }
        if n_features == 0 {
            return Err(CyaneaError::InvalidInput("n_features must be > 0".into()));
        }
        if data.len() % n_features != 0 {
            return Err(CyaneaError::InvalidInput(format!(
                "data length {} not divisible by n_features {}",
                data.len(),
                n_features
            )));
        }
        let n_samples = data.len() / n_features;
        if targets.len() != n_samples {
            return Err(CyaneaError::InvalidInput(format!(
                "targets length {} != n_samples {}",
                targets.len(),
                n_samples
            )));
        }
        if n_samples <= n_features {
            return Err(CyaneaError::InvalidInput(
                "need more samples than features for OLS".into(),
            ));
        }

        // Augmented system with bias: X_aug = [X | 1], solving (X^T X) w = X^T y
        let p = n_features + 1; // features + bias column

        // Build X^T X (p × p)
        let mut xtx = vec![0.0; p * p];
        let mut xty = vec![0.0; p];

        for row in 0..n_samples {
            let r = &data[row * n_features..(row + 1) * n_features];
            let y = targets[row];

            for i in 0..n_features {
                for j in 0..n_features {
                    xtx[i * p + j] += r[i] * r[j];
                }
                // Bias column
                xtx[i * p + n_features] += r[i];
                xtx[n_features * p + i] += r[i];
                xty[i] += r[i] * y;
            }
            // Bias-bias
            xtx[n_features * p + n_features] += 1.0;
            xty[n_features] += y;
        }

        // Solve via Gaussian elimination with partial pivoting
        let solution = gauss_solve(&mut xtx, &mut xty, p)?;

        let weights: Vec<f64> = solution[..n_features].to_vec();
        let bias = solution[n_features];

        // Compute R-squared
        let y_mean: f64 = targets.iter().sum::<f64>() / n_samples as f64;
        let mut ss_res = 0.0;
        let mut ss_tot = 0.0;
        for row in 0..n_samples {
            let r = &data[row * n_features..(row + 1) * n_features];
            let mut pred = bias;
            for (i, &w) in weights.iter().enumerate() {
                pred += r[i] * w;
            }
            ss_res += (targets[row] - pred).powi(2);
            ss_tot += (targets[row] - y_mean).powi(2);
        }
        let r_squared = if ss_tot > 0.0 {
            1.0 - ss_res / ss_tot
        } else {
            0.0
        };

        Ok(Self {
            weights,
            bias,
            n_features,
            r_squared,
        })
    }

    /// Predict the target for a single query point.
    pub fn predict(&self, query: &[f64]) -> Result<f64> {
        if query.len() != self.n_features {
            return Err(CyaneaError::InvalidInput(format!(
                "query has {} features, expected {}",
                query.len(),
                self.n_features
            )));
        }
        let mut val = self.bias;
        for (i, &w) in self.weights.iter().enumerate() {
            val += query[i] * w;
        }
        Ok(val)
    }

    /// Predict targets for multiple query points.
    ///
    /// `queries` is flat row-major `n_queries × n_features`.
    pub fn predict_batch(&self, queries: &[f64]) -> Result<Vec<f64>> {
        if queries.is_empty() {
            return Ok(vec![]);
        }
        if queries.len() % self.n_features != 0 {
            return Err(CyaneaError::InvalidInput(format!(
                "queries length {} not divisible by n_features {}",
                queries.len(),
                self.n_features
            )));
        }
        let n = queries.len() / self.n_features;
        let mut predictions = Vec::with_capacity(n);
        for row in 0..n {
            let q = &queries[row * self.n_features..(row + 1) * self.n_features];
            predictions.push(self.predict(q)?);
        }
        Ok(predictions)
    }
}

impl Summarizable for LinearRegression {
    fn summary(&self) -> String {
        format!(
            "LinearRegression: {} features, R²={:.4}",
            self.n_features, self.r_squared,
        )
    }
}

/// Solve Ax = b via Gaussian elimination with partial pivoting.
fn gauss_solve(a: &mut [f64], b: &mut [f64], n: usize) -> Result<Vec<f64>> {
    // Forward elimination
    for col in 0..n {
        // Partial pivoting: find row with largest absolute value in this column
        let mut max_val = a[col * n + col].abs();
        let mut max_row = col;
        for row in (col + 1)..n {
            let val = a[row * n + col].abs();
            if val > max_val {
                max_val = val;
                max_row = row;
            }
        }

        if max_val < 1e-14 {
            return Err(CyaneaError::Other("singular matrix in OLS".into()));
        }

        // Swap rows
        if max_row != col {
            for j in 0..n {
                a.swap(col * n + j, max_row * n + j);
            }
            b.swap(col, max_row);
        }

        // Eliminate below
        let pivot = a[col * n + col];
        for row in (col + 1)..n {
            let factor = a[row * n + col] / pivot;
            for j in col..n {
                a[row * n + j] -= factor * a[col * n + j];
            }
            b[row] -= factor * b[col];
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for col in (0..n).rev() {
        let mut sum = b[col];
        for j in (col + 1)..n {
            sum -= a[col * n + j] * x[j];
        }
        x[col] = sum / a[col * n + col];
    }

    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- KNN tests ---

    fn simple_2d_data() -> Vec<f64> {
        // 6 points: cluster A near (0,0), cluster B near (10,10)
        vec![
            0.0, 0.0,
            0.1, 0.0,
            0.0, 0.1,
            10.0, 10.0,
            10.1, 10.0,
            10.0, 10.1,
        ]
    }

    #[test]
    fn knn_classify_basic() {
        let data = simple_2d_data();
        let labels = vec![0, 0, 0, 1, 1, 1];
        let config = KnnConfig {
            k: 3,
            ..Default::default()
        };
        let model = KnnModel::fit(&data, 2, config).unwrap();
        assert_eq!(model.classify(&[0.05, 0.05], &labels).unwrap(), 0);
        assert_eq!(model.classify(&[10.05, 10.05], &labels).unwrap(), 1);
    }

    #[test]
    fn knn_regress_basic() {
        let data = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let targets = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let config = KnnConfig {
            k: 2,
            metric: DistanceMetric::Euclidean,
        };
        let model = KnnModel::fit(&data, 1, config).unwrap();
        // Query at 1.5 → neighbors 1 and 2 → mean = 1.5
        let pred = model.regress(&[1.5], &targets).unwrap();
        assert!((pred - 1.5).abs() < 1e-10);
    }

    #[test]
    fn knn_neighbors() {
        let data = vec![0.0, 1.0, 2.0, 10.0];
        let config = KnnConfig {
            k: 2,
            metric: DistanceMetric::Euclidean,
        };
        let model = KnnModel::fit(&data, 1, config).unwrap();
        let nb = model.neighbors(&[0.5]).unwrap();
        assert_eq!(nb.len(), 2);
        // Closest should be index 0 (distance 0.5)
        assert_eq!(nb[0].0, 0);
    }

    #[test]
    fn knn_fit_errors() {
        let config = KnnConfig::default();
        assert!(KnnModel::fit(&[], 2, config.clone()).is_err());

        let config_k0 = KnnConfig {
            k: 0,
            ..Default::default()
        };
        assert!(KnnModel::fit(&[1.0, 2.0], 2, config_k0).is_err());

        let config_k_big = KnnConfig {
            k: 10,
            ..Default::default()
        };
        assert!(KnnModel::fit(&[1.0, 2.0], 2, config_k_big).is_err());
    }

    #[test]
    fn knn_summary() {
        let data = simple_2d_data();
        let config = KnnConfig {
            k: 3,
            ..Default::default()
        };
        let model = KnnModel::fit(&data, 2, config).unwrap();
        let s = model.summary();
        assert!(s.contains("KNN"));
        assert!(s.contains("k=3"));
    }

    // --- Linear regression tests ---

    #[test]
    fn linreg_perfect_fit() {
        // y = 2x + 1
        let data = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let targets = vec![1.0, 3.0, 5.0, 7.0, 9.0];
        let model = LinearRegression::fit(&data, 1, &targets).unwrap();
        assert!((model.weights[0] - 2.0).abs() < 1e-8);
        assert!((model.bias - 1.0).abs() < 1e-8);
        assert!(model.r_squared > 0.999);
    }

    #[test]
    fn linreg_predict() {
        // y = 3x + 2
        let data = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let targets = vec![2.0, 5.0, 8.0, 11.0, 14.0];
        let model = LinearRegression::fit(&data, 1, &targets).unwrap();
        let pred = model.predict(&[5.0]).unwrap();
        assert!((pred - 17.0).abs() < 1e-6);
    }

    #[test]
    fn linreg_multivariate() {
        // y = 1*x1 + 2*x2 + 3
        let data = vec![
            1.0, 0.0,
            0.0, 1.0,
            1.0, 1.0,
            2.0, 1.0,
            1.0, 2.0,
        ];
        let targets = vec![4.0, 5.0, 6.0, 7.0, 8.0];
        let model = LinearRegression::fit(&data, 2, &targets).unwrap();
        assert!((model.weights[0] - 1.0).abs() < 1e-6);
        assert!((model.weights[1] - 2.0).abs() < 1e-6);
        assert!((model.bias - 3.0).abs() < 1e-6);
    }

    #[test]
    fn linreg_predict_batch() {
        let data = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let targets = vec![0.0, 2.0, 4.0, 6.0, 8.0];
        let model = LinearRegression::fit(&data, 1, &targets).unwrap();
        let preds = model.predict_batch(&[5.0, 10.0]).unwrap();
        assert_eq!(preds.len(), 2);
        assert!((preds[0] - 10.0).abs() < 1e-6);
        assert!((preds[1] - 20.0).abs() < 1e-6);
    }

    #[test]
    fn linreg_fit_errors() {
        assert!(LinearRegression::fit(&[], 1, &[]).is_err());
        assert!(LinearRegression::fit(&[1.0], 1, &[1.0]).is_err()); // n_samples <= n_features
    }

    #[test]
    fn linreg_summary() {
        let data = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let targets = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let model = LinearRegression::fit(&data, 1, &targets).unwrap();
        let s = model.summary();
        assert!(s.contains("LinearRegression"));
        assert!(s.contains("1 features"));
    }
}
