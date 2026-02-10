//! Distance metrics and pairwise distance matrices.

use cyanea_core::{CyaneaError, Result, Summarizable};

/// Supported distance metrics for clustering and evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum DistanceMetric {
    Euclidean,
    Manhattan,
    Cosine,
}

/// Euclidean (L2) distance between two vectors.
pub fn euclidean(a: &[f64], b: &[f64]) -> Result<f64> {
    validate_pair(a, b)?;
    let sum: f64 = a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum();
    Ok(sum.sqrt())
}

/// Manhattan (L1) distance between two vectors.
pub fn manhattan(a: &[f64], b: &[f64]) -> Result<f64> {
    validate_pair(a, b)?;
    let sum: f64 = a.iter().zip(b).map(|(x, y)| (x - y).abs()).sum();
    Ok(sum)
}

/// Cosine similarity between two vectors.
///
/// Returns 0.0 if either vector is the zero vector.
pub fn cosine_similarity(a: &[f64], b: &[f64]) -> Result<f64> {
    validate_pair(a, b)?;
    let mut dot = 0.0_f64;
    let mut norm_a = 0.0_f64;
    let mut norm_b = 0.0_f64;
    for (x, y) in a.iter().zip(b) {
        dot += x * y;
        norm_a += x * x;
        norm_b += y * y;
    }
    let denom = norm_a.sqrt() * norm_b.sqrt();
    if denom == 0.0 {
        return Ok(0.0);
    }
    Ok(dot / denom)
}

/// Cosine distance: `1.0 - cosine_similarity(a, b)`.
pub fn cosine_distance(a: &[f64], b: &[f64]) -> Result<f64> {
    Ok(1.0 - cosine_similarity(a, b)?)
}

/// Hamming distance between two byte slices (e.g. sequences).
///
/// Counts the number of positions where the bytes differ.
pub fn hamming(a: &[u8], b: &[u8]) -> Result<usize> {
    if a.is_empty() {
        return Err(CyaneaError::InvalidInput("empty vectors".into()));
    }
    if a.len() != b.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "length mismatch: {} vs {}",
            a.len(),
            b.len()
        )));
    }
    Ok(a.iter().zip(b).filter(|(x, y)| x != y).count())
}

/// Compute distance between two vectors using the given metric.
pub fn compute_distance(a: &[f64], b: &[f64], metric: DistanceMetric) -> Result<f64> {
    match metric {
        DistanceMetric::Euclidean => euclidean(a, b),
        DistanceMetric::Manhattan => manhattan(a, b),
        DistanceMetric::Cosine => cosine_distance(a, b),
    }
}

/// Symmetric distance matrix stored in condensed upper-triangle form.
///
/// For `n` points the condensed vector has `n*(n-1)/2` elements.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DistanceMatrix {
    condensed: Vec<f64>,
    n: usize,
}

impl DistanceMatrix {
    /// Build a distance matrix from row-vectors of points.
    pub fn from_points(data: &[&[f64]], metric: DistanceMetric) -> Result<Self> {
        let n = data.len();
        if n < 2 {
            return Err(CyaneaError::InvalidInput(
                "need at least 2 points".into(),
            ));
        }
        let dim = data[0].len();
        if dim == 0 {
            return Err(CyaneaError::InvalidInput("empty vectors".into()));
        }
        for (i, row) in data.iter().enumerate() {
            if row.len() != dim {
                return Err(CyaneaError::InvalidInput(format!(
                    "point {} has length {}, expected {}",
                    i,
                    row.len(),
                    dim
                )));
            }
        }
        #[cfg(feature = "parallel")]
        let condensed = {
            use rayon::prelude::*;
            (0..n)
                .into_par_iter()
                .map(|i| {
                    ((i + 1)..n)
                        .map(|j| compute_distance(data[i], data[j], metric))
                        .collect::<Result<Vec<_>>>()
                })
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .flatten()
                .collect::<Vec<f64>>()
        };
        #[cfg(not(feature = "parallel"))]
        let condensed = {
            let size = n * (n - 1) / 2;
            let mut condensed = Vec::with_capacity(size);
            for i in 0..n {
                for j in (i + 1)..n {
                    condensed.push(compute_distance(data[i], data[j], metric)?);
                }
            }
            condensed
        };
        Ok(Self { condensed, n })
    }

    /// Create from a pre-computed condensed distance vector.
    pub fn from_condensed(condensed: Vec<f64>, n: usize) -> Result<Self> {
        let expected = n * (n - 1) / 2;
        if condensed.len() != expected {
            return Err(CyaneaError::InvalidInput(format!(
                "condensed length {} doesn't match n={} (expected {})",
                condensed.len(),
                n,
                expected
            )));
        }
        Ok(Self { condensed, n })
    }

    /// Get the distance between points `i` and `j`.
    ///
    /// Returns 0.0 when `i == j`.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        if i == j {
            return 0.0;
        }
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.condensed[self.index(a, b)]
    }

    /// Number of points.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Access the raw condensed storage.
    pub fn condensed(&self) -> &[f64] {
        &self.condensed
    }

    /// Map (i, j) where i < j to condensed index.
    fn index(&self, i: usize, j: usize) -> usize {
        // row i starts at position: i*n - i*(i+1)/2
        i * self.n - i * (i + 1) / 2 + (j - i - 1)
    }
}

impl Summarizable for DistanceMatrix {
    fn summary(&self) -> String {
        format!("DistanceMatrix: {}x{}", self.n, self.n)
    }
}

/// Convenience alias for [`DistanceMatrix::from_points`].
pub fn pairwise_distances(
    data: &[&[f64]],
    metric: DistanceMetric,
) -> Result<DistanceMatrix> {
    DistanceMatrix::from_points(data, metric)
}

fn validate_pair(a: &[f64], b: &[f64]) -> Result<()> {
    if a.is_empty() {
        return Err(CyaneaError::InvalidInput("empty vectors".into()));
    }
    if a.len() != b.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "length mismatch: {} vs {}",
            a.len(),
            b.len()
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euclidean_known() {
        let d = euclidean(&[0.0, 0.0], &[3.0, 4.0]).unwrap();
        assert!((d - 5.0).abs() < 1e-12);
    }

    #[test]
    fn euclidean_identical() {
        let d = euclidean(&[1.0, 2.0, 3.0], &[1.0, 2.0, 3.0]).unwrap();
        assert!((d - 0.0).abs() < 1e-12);
    }

    #[test]
    fn euclidean_empty_error() {
        assert!(euclidean(&[], &[]).is_err());
    }

    #[test]
    fn euclidean_length_mismatch() {
        assert!(euclidean(&[1.0], &[1.0, 2.0]).is_err());
    }

    #[test]
    fn manhattan_known() {
        let d = manhattan(&[0.0, 0.0], &[3.0, 4.0]).unwrap();
        assert!((d - 7.0).abs() < 1e-12);
    }

    #[test]
    fn cosine_similarity_orthogonal() {
        let s = cosine_similarity(&[1.0, 0.0], &[0.0, 1.0]).unwrap();
        assert!((s - 0.0).abs() < 1e-12);
    }

    #[test]
    fn cosine_similarity_identical() {
        let s = cosine_similarity(&[1.0, 2.0], &[1.0, 2.0]).unwrap();
        assert!((s - 1.0).abs() < 1e-12);
    }

    #[test]
    fn cosine_similarity_opposite() {
        let s = cosine_similarity(&[1.0, 0.0], &[-1.0, 0.0]).unwrap();
        assert!((s - (-1.0)).abs() < 1e-12);
    }

    #[test]
    fn cosine_similarity_zero_vector() {
        let s = cosine_similarity(&[0.0, 0.0], &[1.0, 2.0]).unwrap();
        assert!((s - 0.0).abs() < 1e-12);
    }

    #[test]
    fn cosine_distance_known() {
        let d = cosine_distance(&[1.0, 0.0], &[0.0, 1.0]).unwrap();
        assert!((d - 1.0).abs() < 1e-12);
    }

    #[test]
    fn hamming_known() {
        assert_eq!(hamming(b"ACGT", b"ACGA").unwrap(), 1);
    }

    #[test]
    fn hamming_identical() {
        assert_eq!(hamming(b"ACGT", b"ACGT").unwrap(), 0);
    }

    #[test]
    fn hamming_all_different() {
        assert_eq!(hamming(b"AAAA", b"TTTT").unwrap(), 4);
    }

    #[test]
    fn hamming_empty_error() {
        assert!(hamming(b"", b"").is_err());
    }

    #[test]
    fn distance_matrix_from_points() {
        let pts: Vec<Vec<f64>> = vec![vec![0.0, 0.0], vec![3.0, 0.0], vec![0.0, 4.0]];
        let refs: Vec<&[f64]> = pts.iter().map(|v| v.as_slice()).collect();
        let dm = DistanceMatrix::from_points(&refs, DistanceMetric::Euclidean).unwrap();
        assert_eq!(dm.n(), 3);
        assert!((dm.get(0, 0) - 0.0).abs() < 1e-12);
        assert!((dm.get(0, 1) - 3.0).abs() < 1e-12);
        assert!((dm.get(0, 2) - 4.0).abs() < 1e-12);
        assert!((dm.get(1, 0) - 3.0).abs() < 1e-12); // symmetric
        assert!((dm.get(2, 0) - 4.0).abs() < 1e-12); // symmetric
    }

    #[test]
    fn distance_matrix_summary() {
        let pts: Vec<Vec<f64>> = vec![vec![0.0], vec![1.0], vec![2.0]];
        let refs: Vec<&[f64]> = pts.iter().map(|v| v.as_slice()).collect();
        let dm = DistanceMatrix::from_points(&refs, DistanceMetric::Euclidean).unwrap();
        assert_eq!(dm.summary(), "DistanceMatrix: 3x3");
    }

    #[test]
    fn distance_matrix_too_few_points() {
        let pts: Vec<Vec<f64>> = vec![vec![1.0]];
        let refs: Vec<&[f64]> = pts.iter().map(|v| v.as_slice()).collect();
        assert!(DistanceMatrix::from_points(&refs, DistanceMetric::Euclidean).is_err());
    }

    #[test]
    fn pairwise_distances_alias() {
        let pts: Vec<Vec<f64>> = vec![vec![0.0], vec![5.0]];
        let refs: Vec<&[f64]> = pts.iter().map(|v| v.as_slice()).collect();
        let dm = pairwise_distances(&refs, DistanceMetric::Manhattan).unwrap();
        assert!((dm.get(0, 1) - 5.0).abs() < 1e-12);
    }
}
