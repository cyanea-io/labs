//! Private 3x3 linear algebra for structural superposition.
//!
//! Implements a Jacobi eigenvalue algorithm and SVD decomposition for 3x3
//! matrices without requiring an external linear algebra crate.

use crate::types::Point3D;

/// A 3x3 matrix stored in row-major order.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Matrix3x3 {
    pub data: [[f64; 3]; 3],
}

impl Matrix3x3 {
    /// Zero matrix.
    pub fn zeros() -> Self {
        Self {
            data: [[0.0; 3]; 3],
        }
    }

    /// Identity matrix.
    pub fn identity() -> Self {
        Self {
            data: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    /// Matrix multiplication: self * other.
    pub fn multiply(&self, other: &Matrix3x3) -> Matrix3x3 {
        let mut result = Matrix3x3::zeros();
        for i in 0..3 {
            for j in 0..3 {
                let mut sum = 0.0;
                for k in 0..3 {
                    sum += self.data[i][k] * other.data[k][j];
                }
                result.data[i][j] = sum;
            }
        }
        result
    }

    /// Transpose.
    pub fn transpose(&self) -> Matrix3x3 {
        let mut result = Matrix3x3::zeros();
        for i in 0..3 {
            for j in 0..3 {
                result.data[i][j] = self.data[j][i];
            }
        }
        result
    }

    /// Determinant.
    pub fn determinant(&self) -> f64 {
        let d = &self.data;
        d[0][0] * (d[1][1] * d[2][2] - d[1][2] * d[2][1])
            - d[0][1] * (d[1][0] * d[2][2] - d[1][2] * d[2][0])
            + d[0][2] * (d[1][0] * d[2][1] - d[1][1] * d[2][0])
    }

    /// Apply this matrix as a rotation/transform to a point: M * p.
    pub fn apply(&self, p: &Point3D) -> Point3D {
        Point3D {
            x: self.data[0][0] * p.x + self.data[0][1] * p.y + self.data[0][2] * p.z,
            y: self.data[1][0] * p.x + self.data[1][1] * p.y + self.data[1][2] * p.z,
            z: self.data[2][0] * p.x + self.data[2][1] * p.y + self.data[2][2] * p.z,
        }
    }
}

/// Result of a 3x3 SVD decomposition: A = U * diag(s) * Vt.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub(crate) struct SVD3x3 {
    pub u: Matrix3x3,
    pub s: [f64; 3],
    pub vt: Matrix3x3,
}

/// Compute SVD of a 3x3 matrix using the Jacobi eigenvalue method.
///
/// Algorithm: compute A^T*A, find eigenvalues/eigenvectors via Jacobi iteration,
/// then recover U = A * V * S^{-1}.
pub(crate) fn svd_3x3(matrix: &Matrix3x3) -> SVD3x3 {
    // Compute A^T * A (symmetric, positive semi-definite)
    let ata = matrix.transpose().multiply(matrix);

    // Jacobi eigenvalue decomposition of A^T*A
    let (eigenvectors, eigenvalues) = jacobi_eigenvalue(&ata);

    // Singular values = sqrt(eigenvalues), and V = eigenvectors
    // Sort by descending singular value
    let mut indices = [0usize, 1, 2];
    indices.sort_by(|&a, &b| eigenvalues[b].partial_cmp(&eigenvalues[a]).unwrap());

    let mut s = [0.0f64; 3];
    let mut v = Matrix3x3::zeros();
    for (col, &idx) in indices.iter().enumerate() {
        s[col] = eigenvalues[idx].max(0.0).sqrt();
        for row in 0..3 {
            v.data[row][col] = eigenvectors.data[row][idx];
        }
    }

    // U = A * V * S^{-1}
    let av = matrix.multiply(&v);
    let mut u = Matrix3x3::zeros();
    for col in 0..3 {
        if s[col] > 1e-10 {
            for row in 0..3 {
                u.data[row][col] = av.data[row][col] / s[col];
            }
        } else {
            // For zero singular values, pick an orthogonal vector
            if col == 2 {
                // u2 = u0 x u1
                let u0 = Point3D::new(u.data[0][0], u.data[1][0], u.data[2][0]);
                let u1 = Point3D::new(u.data[0][1], u.data[1][1], u.data[2][1]);
                let u2 = u0.cross(&u1).normalize();
                u.data[0][2] = u2.x;
                u.data[1][2] = u2.y;
                u.data[2][2] = u2.z;
            }
        }
    }

    SVD3x3 {
        u,
        s,
        vt: v.transpose(),
    }
}

/// Jacobi eigenvalue algorithm for a 3x3 symmetric matrix.
///
/// Returns (eigenvectors as columns of a matrix, eigenvalues).
fn jacobi_eigenvalue(matrix: &Matrix3x3) -> (Matrix3x3, [f64; 3]) {
    let mut a = *matrix;
    let mut v = Matrix3x3::identity();

    let max_iter = 100;
    let tol = 1e-15;

    for _ in 0..max_iter {
        // Find the largest off-diagonal element
        let mut max_val = 0.0f64;
        let mut p = 0;
        let mut q = 1;
        for i in 0..3 {
            for j in (i + 1)..3 {
                if a.data[i][j].abs() > max_val {
                    max_val = a.data[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < tol {
            break;
        }

        // Compute Jacobi rotation angle
        let app = a.data[p][p];
        let aqq = a.data[q][q];
        let apq = a.data[p][q];

        let theta = if (app - aqq).abs() < tol {
            core::f64::consts::FRAC_PI_4
        } else {
            0.5 * (2.0 * apq / (app - aqq)).atan()
        };

        let c = theta.cos();
        let s = theta.sin();

        // Apply Givens rotation: A' = G^T * A * G
        let mut new_a = a;

        // Update rows/columns p and q
        for i in 0..3 {
            if i != p && i != q {
                let aip = a.data[i][p];
                let aiq = a.data[i][q];
                new_a.data[i][p] = c * aip + s * aiq;
                new_a.data[p][i] = new_a.data[i][p];
                new_a.data[i][q] = -s * aip + c * aiq;
                new_a.data[q][i] = new_a.data[i][q];
            }
        }

        new_a.data[p][p] = c * c * app + 2.0 * c * s * apq + s * s * aqq;
        new_a.data[q][q] = s * s * app - 2.0 * c * s * apq + c * c * aqq;
        new_a.data[p][q] = 0.0;
        new_a.data[q][p] = 0.0;

        a = new_a;

        // Accumulate rotation into V
        let mut new_v = v;
        for i in 0..3 {
            let vip = v.data[i][p];
            let viq = v.data[i][q];
            new_v.data[i][p] = c * vip + s * viq;
            new_v.data[i][q] = -s * vip + c * viq;
        }
        v = new_v;
    }

    let eigenvalues = [a.data[0][0], a.data[1][1], a.data[2][2]];
    (v, eigenvalues)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiply_identity() {
        let a = Matrix3x3 {
            data: [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
        };
        let result = a.multiply(&Matrix3x3::identity());
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (result.data[i][j] - a.data[i][j]).abs() < 1e-10,
                    "mismatch at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn test_transpose() {
        let a = Matrix3x3 {
            data: [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
        };
        let t = a.transpose();
        assert!((t.data[0][1] - 4.0).abs() < 1e-10);
        assert!((t.data[1][0] - 2.0).abs() < 1e-10);
        assert!((t.data[2][0] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_determinant() {
        let a = Matrix3x3 {
            data: [[1.0, 2.0, 3.0], [0.0, 1.0, 4.0], [5.0, 6.0, 0.0]],
        };
        assert!((a.determinant() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_svd_known_matrix() {
        // Test with a known rotation matrix — its singular values should all be 1
        let angle: f64 = 0.5;
        let c = angle.cos();
        let s = angle.sin();
        let rot = Matrix3x3 {
            data: [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]],
        };

        let svd = svd_3x3(&rot);

        // All singular values should be 1.0
        for &sv in &svd.s {
            assert!(
                (sv - 1.0).abs() < 1e-6,
                "singular value {} should be ~1.0",
                sv
            );
        }

        // U * diag(s) * Vt should reconstruct the original matrix
        let mut reconstructed = Matrix3x3::zeros();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    reconstructed.data[i][j] += svd.u.data[i][k] * svd.s[k] * svd.vt.data[k][j];
                }
            }
        }
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (reconstructed.data[i][j] - rot.data[i][j]).abs() < 1e-6,
                    "reconstruction failed at [{i}][{j}]"
                );
            }
        }
    }
}
