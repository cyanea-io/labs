//! Dense expression matrix for bulk omics data.
//!
//! [`ExpressionMatrix`] stores a row-major dense matrix of `f64` values
//! (n_features × n_samples) with associated feature and sample names.
//! Typical use cases include RNA-seq gene expression, proteomics intensity
//! values, and metabolomics abundances.

use cyanea_core::{CyaneaError, Result, Summarizable};

/// A dense, row-major expression matrix (features × samples).
#[derive(Debug, Clone)]
pub struct ExpressionMatrix {
    data: Vec<f64>,
    n_features: usize,
    n_samples: usize,
    feature_names: Vec<String>,
    sample_names: Vec<String>,
}

impl ExpressionMatrix {
    /// Create a matrix from row-major 2D data.
    ///
    /// Each inner `Vec` is one feature (row) with `n_samples` values.
    pub fn new(
        data: Vec<Vec<f64>>,
        feature_names: Vec<String>,
        sample_names: Vec<String>,
    ) -> Result<Self> {
        let n_features = data.len();
        let n_samples = sample_names.len();

        if feature_names.len() != n_features {
            return Err(CyaneaError::InvalidInput(format!(
                "feature_names length ({}) does not match row count ({n_features})",
                feature_names.len()
            )));
        }

        let mut flat = Vec::with_capacity(n_features * n_samples);
        for (i, row) in data.iter().enumerate() {
            if row.len() != n_samples {
                return Err(CyaneaError::InvalidInput(format!(
                    "row {i} has {} columns, expected {n_samples}",
                    row.len()
                )));
            }
            flat.extend_from_slice(row);
        }

        Ok(Self {
            data: flat,
            n_features,
            n_samples,
            feature_names,
            sample_names,
        })
    }

    /// Create a zero-filled matrix.
    pub fn zeros(
        n_features: usize,
        n_samples: usize,
        feature_names: Vec<String>,
        sample_names: Vec<String>,
    ) -> Result<Self> {
        if feature_names.len() != n_features {
            return Err(CyaneaError::InvalidInput(format!(
                "feature_names length ({}) does not match n_features ({n_features})",
                feature_names.len()
            )));
        }
        if sample_names.len() != n_samples {
            return Err(CyaneaError::InvalidInput(format!(
                "sample_names length ({}) does not match n_samples ({n_samples})",
                sample_names.len()
            )));
        }
        Ok(Self {
            data: vec![0.0; n_features * n_samples],
            n_features,
            n_samples,
            feature_names,
            sample_names,
        })
    }

    /// (n_features, n_samples).
    pub fn shape(&self) -> (usize, usize) {
        (self.n_features, self.n_samples)
    }

    /// Get a single value by feature and sample index.
    pub fn get(&self, feature_idx: usize, sample_idx: usize) -> Option<f64> {
        if feature_idx < self.n_features && sample_idx < self.n_samples {
            Some(self.data[feature_idx * self.n_samples + sample_idx])
        } else {
            None
        }
    }

    /// Set a single value. Returns an error if indices are out of bounds.
    pub fn set(&mut self, feature_idx: usize, sample_idx: usize, value: f64) -> Result<()> {
        if feature_idx >= self.n_features || sample_idx >= self.n_samples {
            return Err(CyaneaError::InvalidInput(format!(
                "index ({feature_idx}, {sample_idx}) out of bounds for ({}, {})",
                self.n_features, self.n_samples
            )));
        }
        self.data[feature_idx * self.n_samples + sample_idx] = value;
        Ok(())
    }

    /// A slice of one feature's expression across all samples.
    pub fn row(&self, feature_idx: usize) -> Option<&[f64]> {
        if feature_idx < self.n_features {
            let start = feature_idx * self.n_samples;
            Some(&self.data[start..start + self.n_samples])
        } else {
            None
        }
    }

    /// All feature values for a single sample (column copy, since data is row-major).
    pub fn column(&self, sample_idx: usize) -> Option<Vec<f64>> {
        if sample_idx >= self.n_samples {
            return None;
        }
        let col: Vec<f64> = (0..self.n_features)
            .map(|r| self.data[r * self.n_samples + sample_idx])
            .collect();
        Some(col)
    }

    /// Mean expression of a feature across all samples.
    pub fn row_mean(&self, feature_idx: usize) -> Option<f64> {
        let row = self.row(feature_idx)?;
        if row.is_empty() {
            return Some(0.0);
        }
        Some(row.iter().sum::<f64>() / row.len() as f64)
    }

    /// Mean expression across all features for a sample.
    pub fn column_mean(&self, sample_idx: usize) -> Option<f64> {
        let col = self.column(sample_idx)?;
        if col.is_empty() {
            return Some(0.0);
        }
        Some(col.iter().sum::<f64>() / col.len() as f64)
    }

    /// Transpose the matrix, swapping features and samples.
    pub fn transpose(&self) -> ExpressionMatrix {
        let mut transposed = vec![0.0; self.data.len()];
        for r in 0..self.n_features {
            for c in 0..self.n_samples {
                transposed[c * self.n_features + r] = self.data[r * self.n_samples + c];
            }
        }
        ExpressionMatrix {
            data: transposed,
            n_features: self.n_samples,
            n_samples: self.n_features,
            feature_names: self.sample_names.clone(),
            sample_names: self.feature_names.clone(),
        }
    }

    /// Subset the matrix to the given feature (row) indices.
    pub fn filter_features(&self, indices: &[usize]) -> Result<ExpressionMatrix> {
        let mut data = Vec::with_capacity(indices.len() * self.n_samples);
        let mut names = Vec::with_capacity(indices.len());

        for &i in indices {
            if i >= self.n_features {
                return Err(CyaneaError::InvalidInput(format!(
                    "feature index {i} out of bounds (n_features={})",
                    self.n_features
                )));
            }
            let start = i * self.n_samples;
            data.extend_from_slice(&self.data[start..start + self.n_samples]);
            names.push(self.feature_names[i].clone());
        }

        Ok(ExpressionMatrix {
            data,
            n_features: indices.len(),
            n_samples: self.n_samples,
            feature_names: names,
            sample_names: self.sample_names.clone(),
        })
    }

    /// Subset the matrix to the given sample (column) indices.
    pub fn filter_samples(&self, indices: &[usize]) -> Result<ExpressionMatrix> {
        for &i in indices {
            if i >= self.n_samples {
                return Err(CyaneaError::InvalidInput(format!(
                    "sample index {i} out of bounds (n_samples={})",
                    self.n_samples
                )));
            }
        }

        let mut data = Vec::with_capacity(self.n_features * indices.len());
        let mut names = Vec::with_capacity(indices.len());

        for &i in indices {
            names.push(self.sample_names[i].clone());
        }

        for r in 0..self.n_features {
            for &c in indices {
                data.push(self.data[r * self.n_samples + c]);
            }
        }

        Ok(ExpressionMatrix {
            data,
            n_features: self.n_features,
            n_samples: indices.len(),
            feature_names: self.feature_names.clone(),
            sample_names: names,
        })
    }

    /// The underlying flat data as a slice (row-major, n_features × n_samples).
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// Feature (gene/protein) names.
    pub fn feature_names(&self) -> &[String] {
        &self.feature_names
    }

    /// Sample names.
    pub fn sample_names(&self) -> &[String] {
        &self.sample_names
    }

    /// Log2-transform all values: `log2(x + pseudocount)`.
    ///
    /// Commonly used to normalize RNA-seq count data.
    pub fn log_transform(&self, pseudocount: f64) -> ExpressionMatrix {
        let data: Vec<f64> = self
            .data
            .iter()
            .map(|&x| (x + pseudocount).log2())
            .collect();
        ExpressionMatrix {
            data,
            n_features: self.n_features,
            n_samples: self.n_samples,
            feature_names: self.feature_names.clone(),
            sample_names: self.sample_names.clone(),
        }
    }
}

impl Summarizable for ExpressionMatrix {
    fn summary(&self) -> String {
        format!(
            "ExpressionMatrix: {} features \u{00d7} {} samples",
            self.n_features, self.n_samples
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_matrix() -> ExpressionMatrix {
        ExpressionMatrix::new(
            vec![
                vec![1.0, 2.0, 3.0],
                vec![4.0, 5.0, 6.0],
            ],
            vec!["gene1".into(), "gene2".into()],
            vec!["s1".into(), "s2".into(), "s3".into()],
        )
        .unwrap()
    }

    #[test]
    fn test_construction() {
        let m = sample_matrix();
        assert_eq!(m.shape(), (2, 3));
    }

    #[test]
    fn test_dimension_mismatch() {
        let result = ExpressionMatrix::new(
            vec![vec![1.0, 2.0]],
            vec!["gene1".into(), "gene2".into()], // 2 names, 1 row
            vec!["s1".into(), "s2".into()],
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_row_length_mismatch() {
        let result = ExpressionMatrix::new(
            vec![vec![1.0, 2.0], vec![3.0]], // second row too short
            vec!["gene1".into(), "gene2".into()],
            vec!["s1".into(), "s2".into()],
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_zeros() {
        let m = ExpressionMatrix::zeros(
            2,
            3,
            vec!["a".into(), "b".into()],
            vec!["x".into(), "y".into(), "z".into()],
        )
        .unwrap();
        assert_eq!(m.get(0, 0), Some(0.0));
        assert_eq!(m.get(1, 2), Some(0.0));
    }

    #[test]
    fn test_get_set() {
        let mut m = sample_matrix();
        assert_eq!(m.get(0, 0), Some(1.0));
        assert_eq!(m.get(1, 2), Some(6.0));
        assert_eq!(m.get(2, 0), None);

        m.set(0, 0, 99.0).unwrap();
        assert_eq!(m.get(0, 0), Some(99.0));
        assert!(m.set(5, 0, 1.0).is_err());
    }

    #[test]
    fn test_row() {
        let m = sample_matrix();
        assert_eq!(m.row(0), Some(&[1.0, 2.0, 3.0][..]));
        assert_eq!(m.row(1), Some(&[4.0, 5.0, 6.0][..]));
        assert_eq!(m.row(2), None);
    }

    #[test]
    fn test_column() {
        let m = sample_matrix();
        assert_eq!(m.column(0), Some(vec![1.0, 4.0]));
        assert_eq!(m.column(2), Some(vec![3.0, 6.0]));
        assert_eq!(m.column(3), None);
    }

    #[test]
    fn test_row_mean() {
        let m = sample_matrix();
        assert_eq!(m.row_mean(0), Some(2.0)); // (1+2+3)/3
        assert_eq!(m.row_mean(1), Some(5.0)); // (4+5+6)/3
    }

    #[test]
    fn test_column_mean() {
        let m = sample_matrix();
        assert_eq!(m.column_mean(0), Some(2.5)); // (1+4)/2
        assert_eq!(m.column_mean(1), Some(3.5)); // (2+5)/2
    }

    #[test]
    fn test_transpose() {
        let m = sample_matrix();
        let t = m.transpose();
        assert_eq!(t.shape(), (3, 2));
        assert_eq!(t.get(0, 0), Some(1.0));
        assert_eq!(t.get(0, 1), Some(4.0));
        assert_eq!(t.get(2, 1), Some(6.0));
    }

    #[test]
    fn test_filter_features() {
        let m = sample_matrix();
        let filtered = m.filter_features(&[1]).unwrap();
        assert_eq!(filtered.shape(), (1, 3));
        assert_eq!(filtered.get(0, 0), Some(4.0));

        assert!(m.filter_features(&[5]).is_err());
    }

    #[test]
    fn test_filter_samples() {
        let m = sample_matrix();
        let filtered = m.filter_samples(&[0, 2]).unwrap();
        assert_eq!(filtered.shape(), (2, 2));
        assert_eq!(filtered.get(0, 0), Some(1.0));
        assert_eq!(filtered.get(0, 1), Some(3.0));
        assert_eq!(filtered.get(1, 0), Some(4.0));

        assert!(m.filter_samples(&[5]).is_err());
    }

    #[test]
    fn test_log_transform() {
        let m = sample_matrix();
        let logged = m.log_transform(1.0);
        // log2(1.0 + 1.0) = 1.0
        assert!((logged.get(0, 0).unwrap() - 1.0).abs() < 1e-10);
        // log2(4.0 + 1.0) = log2(5) ≈ 2.3219
        assert!((logged.get(1, 0).unwrap() - 5.0_f64.log2()).abs() < 1e-10);
    }

    #[test]
    fn test_summary() {
        let m = sample_matrix();
        assert_eq!(m.summary(), "ExpressionMatrix: 2 features \u{00d7} 3 samples");
    }

    #[test]
    fn test_empty_matrix() {
        let m = ExpressionMatrix::new(
            vec![],
            vec![],
            vec!["s1".into()],
        )
        .unwrap();
        assert_eq!(m.shape(), (0, 1));
    }

    #[test]
    fn test_as_slice() {
        let m = sample_matrix();
        assert_eq!(m.as_slice(), &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    }

    #[test]
    fn test_feature_names() {
        let m = sample_matrix();
        assert_eq!(m.feature_names(), &["gene1", "gene2"]);
    }

    #[test]
    fn test_sample_names() {
        let m = sample_matrix();
        assert_eq!(m.sample_names(), &["s1", "s2", "s3"]);
    }
}
