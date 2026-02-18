//! High-level operation primitives that delegate to a [`Backend`].

use cyanea_core::{CyaneaError, Result};

use crate::backend::{Backend, DistanceMetricGpu};
use crate::buffer::Buffer;

/// Reduces a buffer to a single sum.
pub fn reduce_sum(backend: &dyn Backend, buf: &Buffer) -> Result<f64> {
    backend.reduce_sum(buf)
}

/// Reduces a buffer to its minimum value.
pub fn reduce_min(backend: &dyn Backend, buf: &Buffer) -> Result<f64> {
    backend.reduce_min(buf)
}

/// Reduces a buffer to its maximum value.
pub fn reduce_max(backend: &dyn Backend, buf: &Buffer) -> Result<f64> {
    backend.reduce_max(buf)
}

/// Reduces a buffer to its arithmetic mean.
pub fn reduce_mean(backend: &dyn Backend, buf: &Buffer) -> Result<f64> {
    let len = backend.buffer_len(buf);
    if len == 0 {
        return Err(CyaneaError::InvalidInput(
            "reduce_mean: buffer is empty".to_string(),
        ));
    }
    let sum = backend.reduce_sum(buf)?;
    Ok(sum / len as f64)
}

/// Applies `f` element-wise, returning a new output buffer.
pub fn elementwise_map(
    backend: &dyn Backend,
    input: &Buffer,
    f: &dyn Fn(f64) -> f64,
) -> Result<Buffer> {
    let len = backend.buffer_len(input);
    let mut output = backend.buffer_zeros(len)?;
    backend.elementwise_map(input, &mut output, f)?;
    Ok(output)
}

/// Computes an n×n pairwise distance matrix.
///
/// `data` is a flat row-major matrix of shape `(n, dim)`.
pub fn pairwise_distance_matrix(
    backend: &dyn Backend,
    data: &Buffer,
    n: usize,
    dim: usize,
    metric: DistanceMetricGpu,
) -> Result<Buffer> {
    backend.pairwise_distance_matrix(data, n, dim, metric)
}

/// Multiplies two matrices: C = A × B.
///
/// `a` is row-major (m, k), `b` is row-major (k, n).
pub fn matrix_multiply(
    backend: &dyn Backend,
    a: &Buffer,
    b: &Buffer,
    m: usize,
    k: usize,
    n: usize,
) -> Result<Buffer> {
    backend.matrix_multiply(a, b, m, k, n)
}

/// Computes all-pairs results using a custom function.
///
/// `items` is a flat row-major matrix of shape `(n, item_len)`.
pub fn batch_pairwise(
    backend: &dyn Backend,
    items: &Buffer,
    n: usize,
    item_len: usize,
    f: &dyn Fn(&[f64], &[f64]) -> f64,
) -> Result<Buffer> {
    backend.batch_pairwise(items, n, item_len, f)
}

/// Computes an n×n pairwise distance matrix using tiled evaluation.
///
/// Splits the n rows into tiles of at most `tile_size` rows, computes
/// inter-tile and intra-tile distance sub-matrices using the backend's
/// `pairwise_distance_matrix`, and assembles the full n×n result. This
/// enables larger-than-GPU-memory distance computations by processing
/// tiles sequentially.
///
/// `data` is a flat row-major matrix of shape `(n, dim)`.
pub fn tiled_pairwise_distance(
    backend: &dyn Backend,
    data: &[f64],
    n: usize,
    dim: usize,
    metric: DistanceMetricGpu,
    tile_size: usize,
) -> Result<Vec<f64>> {
    if n == 0 {
        return Ok(vec![]);
    }
    if data.len() != n * dim {
        return Err(CyaneaError::InvalidInput(format!(
            "tiled_pairwise_distance: expected {} elements ({}×{}), got {}",
            n * dim,
            n,
            dim,
            data.len()
        )));
    }
    if tile_size == 0 {
        return Err(CyaneaError::InvalidInput(
            "tiled_pairwise_distance: tile_size must be > 0".to_string(),
        ));
    }

    let mut result = vec![0.0_f64; n * n];
    let num_tiles = (n + tile_size - 1) / tile_size;

    for ti in 0..num_tiles {
        let i_start = ti * tile_size;
        let i_end = (i_start + tile_size).min(n);
        let ni = i_end - i_start;

        for tj in ti..num_tiles {
            let j_start = tj * tile_size;
            let j_end = (j_start + tile_size).min(n);
            let nj = j_end - j_start;

            if ti == tj {
                // Intra-tile: compute ni×ni distance matrix directly
                let tile_data: Vec<f64> = data[i_start * dim..i_end * dim].to_vec();
                let buf = backend.buffer_from_slice(&tile_data)?;
                let dist_buf = backend.pairwise_distance_matrix(&buf, ni, dim, metric)?;
                let dist = backend.read_buffer(&dist_buf)?;
                for li in 0..ni {
                    for lj in 0..ni {
                        result[(i_start + li) * n + (i_start + lj)] = dist[li * ni + lj];
                    }
                }
            } else {
                // Inter-tile: combine tile_I ++ tile_J, compute (ni+nj)×(ni+nj),
                // extract cross-tile sub-matrix
                let mut combined = Vec::with_capacity((ni + nj) * dim);
                combined.extend_from_slice(&data[i_start * dim..i_end * dim]);
                combined.extend_from_slice(&data[j_start * dim..j_end * dim]);
                let total = ni + nj;
                let buf = backend.buffer_from_slice(&combined)?;
                let dist_buf = backend.pairwise_distance_matrix(&buf, total, dim, metric)?;
                let dist = backend.read_buffer(&dist_buf)?;
                // Extract cross-tile distances (rows 0..ni, cols ni..total)
                for li in 0..ni {
                    for lj in 0..nj {
                        let d = dist[li * total + (ni + lj)];
                        result[(i_start + li) * n + (j_start + lj)] = d;
                        result[(j_start + lj) * n + (i_start + li)] = d;
                    }
                }
            }
        }
    }

    Ok(result)
}

/// Column-wise z-score normalization.
///
/// `data` is a flat row-major matrix of shape `(n_rows, n_cols)`.
/// Each column is normalized in-place to have mean 0 and standard deviation 1.
/// Constant columns (std = 0) are set to 0.0.
///
/// Returns a new buffer with the normalized data.
pub fn batch_z_score(
    backend: &dyn Backend,
    data: &Buffer,
    n_rows: usize,
    n_cols: usize,
) -> Result<Buffer> {
    if n_cols == 0 {
        return Err(CyaneaError::InvalidInput(
            "batch_z_score: n_cols must be > 0".to_string(),
        ));
    }
    let len = backend.buffer_len(data);
    if len == 0 {
        return Err(CyaneaError::InvalidInput(
            "batch_z_score: data is empty".to_string(),
        ));
    }
    if len != n_rows * n_cols {
        return Err(CyaneaError::InvalidInput(format!(
            "batch_z_score: expected {} elements ({}×{}), got {}",
            n_rows * n_cols,
            n_rows,
            n_cols,
            len
        )));
    }
    let raw = backend.read_buffer(data)?;
    let mut result = raw.clone();

    for col in 0..n_cols {
        // Compute column mean
        let mut sum = 0.0_f64;
        for row in 0..n_rows {
            sum += raw[row * n_cols + col];
        }
        let mean = sum / n_rows as f64;

        // Compute column variance
        let mut var_sum = 0.0_f64;
        for row in 0..n_rows {
            let diff = raw[row * n_cols + col] - mean;
            var_sum += diff * diff;
        }
        let std = (var_sum / n_rows as f64).sqrt();

        // Normalize
        for row in 0..n_rows {
            result[row * n_cols + col] = if std > 0.0 {
                (raw[row * n_cols + col] - mean) / std
            } else {
                0.0
            };
        }
    }

    backend.buffer_from_slice(&result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cpu::CpuBackend;

    fn backend() -> CpuBackend {
        CpuBackend::new()
    }

    #[test]
    fn ops_reduce_sum() {
        let b = backend();
        let buf = b.buffer_from_slice(&[1.0, 2.0, 3.0]).unwrap();
        assert!((reduce_sum(&b, &buf).unwrap() - 6.0).abs() < 1e-12);
    }

    #[test]
    fn ops_reduce_min() {
        let b = backend();
        let buf = b.buffer_from_slice(&[3.0, 1.0, 2.0]).unwrap();
        assert!((reduce_min(&b, &buf).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn ops_reduce_max() {
        let b = backend();
        let buf = b.buffer_from_slice(&[3.0, 1.0, 2.0]).unwrap();
        assert!((reduce_max(&b, &buf).unwrap() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn ops_reduce_mean() {
        let b = backend();
        let buf = b.buffer_from_slice(&[2.0, 4.0, 6.0]).unwrap();
        assert!((reduce_mean(&b, &buf).unwrap() - 4.0).abs() < 1e-12);
    }

    #[test]
    fn ops_reduce_mean_empty() {
        let b = backend();
        let buf = b.buffer_from_slice(&[]).unwrap();
        assert!(reduce_mean(&b, &buf).is_err());
    }

    #[test]
    fn ops_elementwise_map() {
        let b = backend();
        let input = b.buffer_from_slice(&[1.0, 4.0, 9.0]).unwrap();
        let output = elementwise_map(&b, &input, &|x| x.sqrt()).unwrap();
        let data = b.read_buffer(&output).unwrap();
        assert!((data[0] - 1.0).abs() < 1e-12);
        assert!((data[1] - 2.0).abs() < 1e-12);
        assert!((data[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn ops_pairwise_distance_matrix() {
        let b = backend();
        let data = b.buffer_from_slice(&[0.0, 0.0, 1.0, 0.0]).unwrap();
        let result =
            pairwise_distance_matrix(&b, &data, 2, 2, DistanceMetricGpu::Euclidean).unwrap();
        let mat = b.read_buffer(&result).unwrap();
        assert!((mat[0 * 2 + 1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn ops_matrix_multiply() {
        let b = backend();
        let a = b.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();
        let bm = b.buffer_from_slice(&[5.0, 6.0, 7.0, 8.0]).unwrap();
        let result = matrix_multiply(&b, &a, &bm, 2, 2, 2).unwrap();
        assert_eq!(
            b.read_buffer(&result).unwrap(),
            vec![19.0, 22.0, 43.0, 50.0]
        );
    }

    #[test]
    fn ops_batch_pairwise() {
        let b = backend();
        let items = b.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();
        let result = batch_pairwise(&b, &items, 2, 2, &|a, bv| {
            a.iter().zip(bv).map(|(x, y)| (x - y).abs()).sum()
        })
        .unwrap();
        let mat = b.read_buffer(&result).unwrap();
        assert!((mat[0 * 2 + 1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn ops_batch_z_score() {
        let b = backend();
        // 3 rows × 2 cols: col0 = [1,2,3], col1 = [10,20,30]
        let data = b
            .buffer_from_slice(&[1.0, 10.0, 2.0, 20.0, 3.0, 30.0])
            .unwrap();
        let result = batch_z_score(&b, &data, 3, 2).unwrap();
        let out = b.read_buffer(&result).unwrap();
        // Each column should have mean ≈ 0
        let col0_mean = (out[0] + out[2] + out[4]) / 3.0;
        let col1_mean = (out[1] + out[3] + out[5]) / 3.0;
        assert!(col0_mean.abs() < 1e-12);
        assert!(col1_mean.abs() < 1e-12);
        // First element of col0: (1 - 2) / std = -1/std
        // std of [1,2,3] = sqrt(2/3)
        let expected_std = (2.0_f64 / 3.0).sqrt();
        assert!((out[0] - (-1.0 / expected_std)).abs() < 1e-12);
    }

    #[test]
    fn ops_batch_z_score_constant_column() {
        let b = backend();
        // 2 rows × 2 cols: col0 = [5,5] (constant), col1 = [1,3]
        let data = b
            .buffer_from_slice(&[5.0, 1.0, 5.0, 3.0])
            .unwrap();
        let result = batch_z_score(&b, &data, 2, 2).unwrap();
        let out = b.read_buffer(&result).unwrap();
        // Constant column → all zeros
        assert!((out[0]).abs() < 1e-12);
        assert!((out[2]).abs() < 1e-12);
    }

    #[test]
    fn ops_batch_z_score_dimension_mismatch() {
        let b = backend();
        let data = b.buffer_from_slice(&[1.0, 2.0, 3.0]).unwrap();
        assert!(batch_z_score(&b, &data, 2, 2).is_err());
    }

    // ── Tiled pairwise distance tests ──────────────────────────────

    #[test]
    fn tiled_matches_full() {
        let b = backend();
        // 4 points in 2D
        let data = vec![0.0, 0.0, 3.0, 0.0, 0.0, 4.0, 1.0, 1.0];
        let buf = b.buffer_from_slice(&data).unwrap();
        let full_buf = b
            .pairwise_distance_matrix(&buf, 4, 2, DistanceMetricGpu::Euclidean)
            .unwrap();
        let full = b.read_buffer(&full_buf).unwrap();
        let tiled =
            tiled_pairwise_distance(&b, &data, 4, 2, DistanceMetricGpu::Euclidean, 2).unwrap();
        assert_eq!(full.len(), tiled.len());
        for (i, (f, t)) in full.iter().zip(&tiled).enumerate() {
            assert!(
                (f - t).abs() < 1e-10,
                "mismatch at index {i}: full={f}, tiled={t}"
            );
        }
    }

    #[test]
    fn tiled_single_tile() {
        let b = backend();
        let data = vec![0.0, 0.0, 1.0, 1.0];
        let tiled =
            tiled_pairwise_distance(&b, &data, 2, 2, DistanceMetricGpu::Euclidean, 100).unwrap();
        let expected = (2.0_f64).sqrt();
        assert!((tiled[0 * 2 + 1] - expected).abs() < 1e-10);
        assert!((tiled[1 * 2 + 0] - expected).abs() < 1e-10);
        assert!(tiled[0].abs() < 1e-10);
    }

    #[test]
    fn tiled_exact_division() {
        let b = backend();
        // 6 points in 1D, tile_size=3 → exactly 2 tiles
        let data = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let tiled =
            tiled_pairwise_distance(&b, &data, 6, 1, DistanceMetricGpu::Manhattan, 3).unwrap();
        // d(0,5) = 5.0
        assert!((tiled[0 * 6 + 5] - 5.0).abs() < 1e-10);
        // d(2,3) = 1.0
        assert!((tiled[2 * 6 + 3] - 1.0).abs() < 1e-10);
        // symmetric
        assert!((tiled[5 * 6 + 0] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn tiled_remainder() {
        let b = backend();
        // 5 points in 1D, tile_size=2 → 3 tiles (2,2,1)
        let data = vec![0.0, 10.0, 20.0, 30.0, 40.0];
        let tiled =
            tiled_pairwise_distance(&b, &data, 5, 1, DistanceMetricGpu::Manhattan, 2).unwrap();
        // d(0,4) = 40.0
        assert!((tiled[0 * 5 + 4] - 40.0).abs() < 1e-10);
        // d(1,3) = 20.0
        assert!((tiled[1 * 5 + 3] - 20.0).abs() < 1e-10);
        // diagonal = 0
        for i in 0..5 {
            assert!(tiled[i * 5 + i].abs() < 1e-10);
        }
    }
}
