//! Cluster evaluation metrics.

use cyanea_core::{CyaneaError, Result};

use crate::distance::euclidean;

/// Compute the silhouette coefficient for each sample.
///
/// Labels of -1 are treated as noise and excluded from the result (their
/// silhouette value is set to 0.0). Requires at least 2 non-noise clusters.
pub fn silhouette_samples(data: &[&[f64]], labels: &[i32]) -> Result<Vec<f64>> {
    let n = data.len();
    if n != labels.len() {
        return Err(CyaneaError::InvalidInput(
            "data and labels length mismatch".into(),
        ));
    }
    if n < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 samples".into(),
        ));
    }

    // Determine unique non-noise clusters
    let mut unique_labels: Vec<i32> = labels.iter().copied().filter(|&l| l >= 0).collect();
    unique_labels.sort_unstable();
    unique_labels.dedup();

    if unique_labels.len() < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 non-noise clusters".into(),
        ));
    }

    #[cfg(feature = "parallel")]
    let scores = {
        use rayon::prelude::*;
        (0..n)
            .into_par_iter()
            .map(|i| {
                if labels[i] < 0 {
                    return Ok(0.0);
                }

                let mut same_sum = 0.0;
                let mut same_count = 0usize;
                for j in 0..n {
                    if j != i && labels[j] == labels[i] {
                        same_sum += euclidean(data[i], data[j])?;
                        same_count += 1;
                    }
                }
                let a = if same_count > 0 {
                    same_sum / same_count as f64
                } else {
                    0.0
                };

                let mut b = f64::INFINITY;
                for &label in &unique_labels {
                    if label == labels[i] {
                        continue;
                    }
                    let mut other_sum = 0.0;
                    let mut other_count = 0usize;
                    for j in 0..n {
                        if labels[j] == label {
                            other_sum += euclidean(data[i], data[j])?;
                            other_count += 1;
                        }
                    }
                    if other_count > 0 {
                        let mean_dist = other_sum / other_count as f64;
                        if mean_dist < b {
                            b = mean_dist;
                        }
                    }
                }

                let max_ab = a.max(b);
                Ok(if max_ab == 0.0 {
                    0.0
                } else {
                    (b - a) / max_ab
                })
            })
            .collect::<Result<Vec<f64>>>()?
    };
    #[cfg(not(feature = "parallel"))]
    let scores = {
        let mut scores = vec![0.0; n];
        for i in 0..n {
            if labels[i] < 0 {
                continue;
            }

            let mut same_sum = 0.0;
            let mut same_count = 0usize;
            for j in 0..n {
                if j != i && labels[j] == labels[i] {
                    same_sum += euclidean(data[i], data[j])?;
                    same_count += 1;
                }
            }
            let a = if same_count > 0 {
                same_sum / same_count as f64
            } else {
                0.0
            };

            let mut b = f64::INFINITY;
            for &label in &unique_labels {
                if label == labels[i] {
                    continue;
                }
                let mut other_sum = 0.0;
                let mut other_count = 0usize;
                for j in 0..n {
                    if labels[j] == label {
                        other_sum += euclidean(data[i], data[j])?;
                        other_count += 1;
                    }
                }
                if other_count > 0 {
                    let mean_dist = other_sum / other_count as f64;
                    if mean_dist < b {
                        b = mean_dist;
                    }
                }
            }

            let max_ab = a.max(b);
            scores[i] = if max_ab == 0.0 {
                0.0
            } else {
                (b - a) / max_ab
            };
        }
        scores
    };

    Ok(scores)
}

/// Mean silhouette score across all non-noise samples.
///
/// See [`silhouette_samples`] for details.
pub fn silhouette_score(data: &[&[f64]], labels: &[i32]) -> Result<f64> {
    let samples = silhouette_samples(data, labels)?;
    let non_noise: Vec<f64> = samples
        .iter()
        .zip(labels)
        .filter(|(_, &l)| l >= 0)
        .map(|(&s, _)| s)
        .collect();
    if non_noise.is_empty() {
        return Err(CyaneaError::InvalidInput("no non-noise samples".into()));
    }
    Ok(non_noise.iter().sum::<f64>() / non_noise.len() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_refs(data: &[Vec<f64>]) -> Vec<&[f64]> {
        data.iter().map(|v| v.as_slice()).collect()
    }

    #[test]
    fn perfect_separation() {
        // Two well-separated clusters
        let data = vec![
            vec![0.0, 0.0],
            vec![0.1, 0.0],
            vec![0.0, 0.1],
            vec![10.0, 10.0],
            vec![10.1, 10.0],
            vec![10.0, 10.1],
        ];
        let refs = make_refs(&data);
        let labels = vec![0, 0, 0, 1, 1, 1];
        let score = silhouette_score(&refs, &labels).unwrap();
        assert!(score > 0.9, "expected high score, got {}", score);
    }

    #[test]
    fn single_cluster_error() {
        let data = vec![vec![0.0], vec![1.0], vec![2.0]];
        let refs = make_refs(&data);
        let labels = vec![0, 0, 0];
        assert!(silhouette_score(&refs, &labels).is_err());
    }

    #[test]
    fn noise_excluded() {
        let data = vec![
            vec![0.0, 0.0],
            vec![0.1, 0.0],
            vec![10.0, 10.0],
            vec![10.1, 10.0],
            vec![5.0, 5.0], // noise
        ];
        let refs = make_refs(&data);
        let labels = vec![0, 0, 1, 1, -1];
        let samples = silhouette_samples(&refs, &labels).unwrap();
        assert!((samples[4] - 0.0).abs() < 1e-12); // noise gets 0
    }

    #[test]
    fn value_range() {
        let data = vec![
            vec![0.0],
            vec![1.0],
            vec![5.0],
            vec![6.0],
        ];
        let refs = make_refs(&data);
        let labels = vec![0, 0, 1, 1];
        let samples = silhouette_samples(&refs, &labels).unwrap();
        for &s in &samples {
            assert!(s >= -1.0 && s <= 1.0, "silhouette {} out of range", s);
        }
    }

    #[test]
    fn length_mismatch_error() {
        let data = vec![vec![0.0], vec![1.0]];
        let refs = make_refs(&data);
        let labels = vec![0];
        assert!(silhouette_samples(&refs, &labels).is_err());
    }

    #[test]
    fn too_few_samples_error() {
        let data = vec![vec![0.0]];
        let refs = make_refs(&data);
        let labels = vec![0];
        assert!(silhouette_samples(&refs, &labels).is_err());
    }
}
