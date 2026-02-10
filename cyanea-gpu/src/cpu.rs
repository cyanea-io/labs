//! CPU backend — full `Backend` implementation using plain iterators.

use cyanea_core::{CyaneaError, Result};

use crate::backend::{Backend, BackendKind, DeviceInfo, DistanceMetricGpu};
use crate::buffer::Buffer;

/// CPU compute backend.
///
/// All operations run on the host using standard iterators. This backend
/// serves as the reference implementation and fallback when no GPU is
/// available.
pub struct CpuBackend {
    parallelism: usize,
}

impl CpuBackend {
    /// Creates a new CPU backend, detecting available parallelism.
    pub fn new() -> Self {
        let parallelism = std::thread::available_parallelism()
            .map(|p| p.get())
            .unwrap_or(1);
        Self { parallelism }
    }
}

impl Default for CpuBackend {
    fn default() -> Self {
        Self::new()
    }
}

impl Backend for CpuBackend {
    fn device_info(&self) -> DeviceInfo {
        DeviceInfo {
            name: "CPU".to_string(),
            kind: BackendKind::Cpu,
            total_memory: 0,
            max_parallelism: self.parallelism,
        }
    }

    fn buffer_from_slice(&self, data: &[f64]) -> Result<Buffer> {
        Ok(Buffer::from_host(data.to_vec(), BackendKind::Cpu))
    }

    fn buffer_zeros(&self, len: usize) -> Result<Buffer> {
        Ok(Buffer::from_host(vec![0.0; len], BackendKind::Cpu))
    }

    fn read_buffer(&self, buf: &Buffer) -> Result<Vec<f64>> {
        buf.as_host_slice()
            .map(|s| s.to_vec())
            .ok_or_else(|| CyaneaError::Other("buffer has no host data".to_string()))
    }

    fn write_buffer(&self, buf: &mut Buffer, data: &[f64]) -> Result<()> {
        if data.len() != buf.len() {
            return Err(CyaneaError::InvalidInput(format!(
                "write_buffer length mismatch: buffer has {} elements, got {}",
                buf.len(),
                data.len()
            )));
        }
        buf.host_data = Some(data.to_vec());
        Ok(())
    }

    fn buffer_len(&self, buf: &Buffer) -> usize {
        buf.len()
    }

    fn elementwise_map(
        &self,
        input: &Buffer,
        output: &mut Buffer,
        f: &dyn Fn(f64) -> f64,
    ) -> Result<()> {
        let src = input
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("input buffer has no host data".to_string()))?;
        if output.len() != input.len() {
            return Err(CyaneaError::InvalidInput(format!(
                "elementwise_map: output length {} != input length {}",
                output.len(),
                input.len()
            )));
        }
        let mapped: Vec<f64> = src.iter().map(|&x| f(x)).collect();
        output.host_data = Some(mapped);
        Ok(())
    }

    fn reduce_sum(&self, buf: &Buffer) -> Result<f64> {
        let data = buf
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("buffer has no host data".to_string()))?;
        if data.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "reduce_sum: buffer is empty".to_string(),
            ));
        }
        Ok(data.iter().sum())
    }

    fn reduce_min(&self, buf: &Buffer) -> Result<f64> {
        let data = buf
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("buffer has no host data".to_string()))?;
        if data.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "reduce_min: buffer is empty".to_string(),
            ));
        }
        Ok(data.iter().copied().fold(f64::INFINITY, f64::min))
    }

    fn reduce_max(&self, buf: &Buffer) -> Result<f64> {
        let data = buf
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("buffer has no host data".to_string()))?;
        if data.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "reduce_max: buffer is empty".to_string(),
            ));
        }
        Ok(data.iter().copied().fold(f64::NEG_INFINITY, f64::max))
    }

    fn pairwise_distance_matrix(
        &self,
        data: &Buffer,
        n: usize,
        dim: usize,
        metric: DistanceMetricGpu,
    ) -> Result<Buffer> {
        let flat = data
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("buffer has no host data".to_string()))?;
        if flat.len() != n * dim {
            return Err(CyaneaError::InvalidInput(format!(
                "pairwise_distance_matrix: expected {} elements ({}×{}), got {}",
                n * dim,
                n,
                dim,
                flat.len()
            )));
        }
        #[cfg(feature = "parallel")]
        let result = {
            use rayon::prelude::*;
            let rows: Vec<Vec<(usize, f64)>> = (0..n)
                .into_par_iter()
                .map(|i| {
                    let row_i = &flat[i * dim..(i + 1) * dim];
                    ((i + 1)..n)
                        .map(|j| {
                            let row_j = &flat[j * dim..(j + 1) * dim];
                            (j, compute_distance(row_i, row_j, metric))
                        })
                        .collect()
                })
                .collect();
            let mut result = vec![0.0_f64; n * n];
            for (i, row) in rows.into_iter().enumerate() {
                for (j, d) in row {
                    result[i * n + j] = d;
                    result[j * n + i] = d;
                }
            }
            result
        };
        #[cfg(not(feature = "parallel"))]
        let result = {
            let mut result = vec![0.0_f64; n * n];
            for i in 0..n {
                let row_i = &flat[i * dim..(i + 1) * dim];
                for j in (i + 1)..n {
                    let row_j = &flat[j * dim..(j + 1) * dim];
                    let d = compute_distance(row_i, row_j, metric);
                    result[i * n + j] = d;
                    result[j * n + i] = d;
                }
            }
            result
        };
        Ok(Buffer::from_host(result, BackendKind::Cpu))
    }

    fn matrix_multiply(
        &self,
        a: &Buffer,
        b: &Buffer,
        m: usize,
        k: usize,
        n: usize,
    ) -> Result<Buffer> {
        let a_data = a
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("buffer a has no host data".to_string()))?;
        let b_data = b
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("buffer b has no host data".to_string()))?;
        if a_data.len() != m * k {
            return Err(CyaneaError::InvalidInput(format!(
                "matrix_multiply: a has {} elements, expected {}×{}={}",
                a_data.len(),
                m,
                k,
                m * k
            )));
        }
        if b_data.len() != k * n {
            return Err(CyaneaError::InvalidInput(format!(
                "matrix_multiply: b has {} elements, expected {}×{}={}",
                b_data.len(),
                k,
                n,
                k * n
            )));
        }
        #[cfg(feature = "parallel")]
        let c = {
            use rayon::prelude::*;
            (0..m)
                .into_par_iter()
                .flat_map(|i| {
                    (0..n)
                        .map(|j| {
                            let mut sum = 0.0;
                            for p in 0..k {
                                sum += a_data[i * k + p] * b_data[p * n + j];
                            }
                            sum
                        })
                        .collect::<Vec<f64>>()
                })
                .collect::<Vec<f64>>()
        };
        #[cfg(not(feature = "parallel"))]
        let c = {
            let mut c = vec![0.0_f64; m * n];
            for i in 0..m {
                for j in 0..n {
                    let mut sum = 0.0;
                    for p in 0..k {
                        sum += a_data[i * k + p] * b_data[p * n + j];
                    }
                    c[i * n + j] = sum;
                }
            }
            c
        };
        Ok(Buffer::from_host(c, BackendKind::Cpu))
    }

    fn batch_pairwise(
        &self,
        items: &Buffer,
        n: usize,
        item_len: usize,
        f: &dyn Fn(&[f64], &[f64]) -> f64,
    ) -> Result<Buffer> {
        let flat = items
            .as_host_slice()
            .ok_or_else(|| CyaneaError::Other("buffer has no host data".to_string()))?;
        if flat.len() != n * item_len {
            return Err(CyaneaError::InvalidInput(format!(
                "batch_pairwise: expected {} elements ({}×{}), got {}",
                n * item_len,
                n,
                item_len,
                flat.len()
            )));
        }
        let mut result = vec![0.0_f64; n * n];
        for i in 0..n {
            let row_i = &flat[i * item_len..(i + 1) * item_len];
            for j in (i + 1)..n {
                let row_j = &flat[j * item_len..(j + 1) * item_len];
                let val = f(row_i, row_j);
                result[i * n + j] = val;
                result[j * n + i] = val;
            }
        }
        Ok(Buffer::from_host(result, BackendKind::Cpu))
    }
}

/// Computes the distance between two vectors for a given metric.
fn compute_distance(a: &[f64], b: &[f64], metric: DistanceMetricGpu) -> f64 {
    match metric {
        DistanceMetricGpu::Euclidean => {
            let sum_sq: f64 = a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum();
            sum_sq.sqrt()
        }
        DistanceMetricGpu::Manhattan => a.iter().zip(b).map(|(x, y)| (x - y).abs()).sum(),
        DistanceMetricGpu::Cosine => {
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
                0.0
            } else {
                1.0 - (dot / denom)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn backend() -> CpuBackend {
        CpuBackend::new()
    }

    #[test]
    fn device_info_kind() {
        let info = backend().device_info();
        assert_eq!(info.kind, BackendKind::Cpu);
        assert!(info.max_parallelism >= 1);
    }

    #[test]
    fn buffer_round_trip() {
        let b = backend();
        let data = vec![1.0, 2.0, 3.0];
        let buf = b.buffer_from_slice(&data).unwrap();
        assert_eq!(b.buffer_len(&buf), 3);
        assert_eq!(b.read_buffer(&buf).unwrap(), data);
    }

    #[test]
    fn buffer_zeros() {
        let b = backend();
        let buf = b.buffer_zeros(4).unwrap();
        assert_eq!(b.read_buffer(&buf).unwrap(), vec![0.0; 4]);
    }

    #[test]
    fn write_buffer_ok() {
        let b = backend();
        let mut buf = b.buffer_zeros(3).unwrap();
        b.write_buffer(&mut buf, &[7.0, 8.0, 9.0]).unwrap();
        assert_eq!(b.read_buffer(&buf).unwrap(), vec![7.0, 8.0, 9.0]);
    }

    #[test]
    fn write_buffer_length_mismatch() {
        let b = backend();
        let mut buf = b.buffer_zeros(3).unwrap();
        assert!(b.write_buffer(&mut buf, &[1.0, 2.0]).is_err());
    }

    #[test]
    fn elementwise_map_double() {
        let b = backend();
        let input = b.buffer_from_slice(&[1.0, 2.0, 3.0]).unwrap();
        let mut output = b.buffer_zeros(3).unwrap();
        b.elementwise_map(&input, &mut output, &|x| x * 2.0)
            .unwrap();
        assert_eq!(b.read_buffer(&output).unwrap(), vec![2.0, 4.0, 6.0]);
    }

    #[test]
    fn reduce_sum_ok() {
        let b = backend();
        let buf = b.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();
        assert!((b.reduce_sum(&buf).unwrap() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn reduce_sum_empty() {
        let b = backend();
        let buf = b.buffer_from_slice(&[]).unwrap();
        assert!(b.reduce_sum(&buf).is_err());
    }

    #[test]
    fn reduce_min_max() {
        let b = backend();
        let buf = b.buffer_from_slice(&[3.0, 1.0, 4.0, 1.5]).unwrap();
        assert!((b.reduce_min(&buf).unwrap() - 1.0).abs() < 1e-12);
        assert!((b.reduce_max(&buf).unwrap() - 4.0).abs() < 1e-12);
    }

    #[test]
    fn reduce_min_empty() {
        let b = backend();
        let buf = b.buffer_from_slice(&[]).unwrap();
        assert!(b.reduce_min(&buf).is_err());
    }

    #[test]
    fn reduce_max_empty() {
        let b = backend();
        let buf = b.buffer_from_slice(&[]).unwrap();
        assert!(b.reduce_max(&buf).is_err());
    }

    #[test]
    fn pairwise_euclidean() {
        let b = backend();
        // 3 points in 2D: (0,0), (3,0), (0,4)
        let data = b
            .buffer_from_slice(&[0.0, 0.0, 3.0, 0.0, 0.0, 4.0])
            .unwrap();
        let result = b
            .pairwise_distance_matrix(&data, 3, 2, DistanceMetricGpu::Euclidean)
            .unwrap();
        let mat = b.read_buffer(&result).unwrap();
        assert!((mat[0 * 3 + 1] - 3.0).abs() < 1e-12); // d(0,1) = 3
        assert!((mat[0 * 3 + 2] - 4.0).abs() < 1e-12); // d(0,2) = 4
        assert!((mat[1 * 3 + 2] - 5.0).abs() < 1e-12); // d(1,2) = 5
        // Diagonal is zero
        assert!((mat[0]).abs() < 1e-12);
    }

    #[test]
    fn pairwise_manhattan() {
        let b = backend();
        let data = b
            .buffer_from_slice(&[0.0, 0.0, 3.0, 4.0])
            .unwrap();
        let result = b
            .pairwise_distance_matrix(&data, 2, 2, DistanceMetricGpu::Manhattan)
            .unwrap();
        let mat = b.read_buffer(&result).unwrap();
        assert!((mat[0 * 2 + 1] - 7.0).abs() < 1e-12); // |3| + |4| = 7
    }

    #[test]
    fn pairwise_cosine() {
        let b = backend();
        // (1,0) and (0,1) → cosine distance = 1.0
        let data = b
            .buffer_from_slice(&[1.0, 0.0, 0.0, 1.0])
            .unwrap();
        let result = b
            .pairwise_distance_matrix(&data, 2, 2, DistanceMetricGpu::Cosine)
            .unwrap();
        let mat = b.read_buffer(&result).unwrap();
        assert!((mat[0 * 2 + 1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn matrix_multiply_identity() {
        let b = backend();
        // 2×2 identity × [1,2; 3,4] = [1,2; 3,4]
        let eye = b.buffer_from_slice(&[1.0, 0.0, 0.0, 1.0]).unwrap();
        let mat = b.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();
        let result = b.matrix_multiply(&eye, &mat, 2, 2, 2).unwrap();
        assert_eq!(b.read_buffer(&result).unwrap(), vec![1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn matrix_multiply_known() {
        let b = backend();
        // [1,2; 3,4] × [5,6; 7,8] = [19,22; 43,50]
        let a = b.buffer_from_slice(&[1.0, 2.0, 3.0, 4.0]).unwrap();
        let mb = b.buffer_from_slice(&[5.0, 6.0, 7.0, 8.0]).unwrap();
        let result = b.matrix_multiply(&a, &mb, 2, 2, 2).unwrap();
        assert_eq!(
            b.read_buffer(&result).unwrap(),
            vec![19.0, 22.0, 43.0, 50.0]
        );
    }

    #[test]
    fn batch_pairwise_abs_diff_sum() {
        let b = backend();
        // 3 items of length 2
        let items = b
            .buffer_from_slice(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
            .unwrap();
        let result = b
            .batch_pairwise(&items, 3, 2, &|a, bv| {
                a.iter().zip(bv).map(|(x, y)| (x - y).abs()).sum()
            })
            .unwrap();
        let mat = b.read_buffer(&result).unwrap();
        // d(0,1) = |1-3| + |2-4| = 4
        assert!((mat[0 * 3 + 1] - 4.0).abs() < 1e-12);
        // d(0,2) = |1-5| + |2-6| = 8
        assert!((mat[0 * 3 + 2] - 8.0).abs() < 1e-12);
        // symmetric
        assert!((mat[1 * 3 + 0] - 4.0).abs() < 1e-12);
    }
}
