//! Coordinate-format (COO) sparse matrix.
//!
//! [`SparseMatrix`] stores non-zero entries as `(row, col, value)` triplets.
//! This format is efficient for incremental construction and iteration, and is
//! the natural import format for single-cell count matrices and other sparse
//! omics data.

use cyanea_core::{CyaneaError, Result, Summarizable};

/// A sparse matrix in COO (coordinate) format.
#[derive(Debug, Clone)]
pub struct SparseMatrix {
    rows: Vec<usize>,
    cols: Vec<usize>,
    values: Vec<f64>,
    n_rows: usize,
    n_cols: usize,
}

impl SparseMatrix {
    /// Create an empty sparse matrix with the given dimensions.
    pub fn new(n_rows: usize, n_cols: usize) -> Self {
        Self {
            rows: Vec::new(),
            cols: Vec::new(),
            values: Vec::new(),
            n_rows,
            n_cols,
        }
    }

    /// Create a sparse matrix from triplet vectors.
    ///
    /// All three vectors must have the same length, and all indices must be
    /// within bounds.
    pub fn from_triplets(
        rows: Vec<usize>,
        cols: Vec<usize>,
        values: Vec<f64>,
        n_rows: usize,
        n_cols: usize,
    ) -> Result<Self> {
        if rows.len() != cols.len() || cols.len() != values.len() {
            return Err(CyaneaError::InvalidInput(
                "rows, cols, and values must have the same length".into(),
            ));
        }
        for (i, (&r, &c)) in rows.iter().zip(cols.iter()).enumerate() {
            if r >= n_rows || c >= n_cols {
                return Err(CyaneaError::InvalidInput(format!(
                    "triplet {i} index ({r}, {c}) out of bounds for ({n_rows}, {n_cols})"
                )));
            }
        }
        Ok(Self {
            rows,
            cols,
            values,
            n_rows,
            n_cols,
        })
    }

    /// Insert a single entry. Returns an error if indices are out of bounds.
    pub fn insert(&mut self, row: usize, col: usize, value: f64) -> Result<()> {
        if row >= self.n_rows || col >= self.n_cols {
            return Err(CyaneaError::InvalidInput(format!(
                "index ({row}, {col}) out of bounds for ({}, {})",
                self.n_rows, self.n_cols
            )));
        }
        self.rows.push(row);
        self.cols.push(col);
        self.values.push(value);
        Ok(())
    }

    /// Get the value at `(row, col)`. Returns 0.0 if no entry is stored.
    ///
    /// This is an O(nnz) scan.
    pub fn get(&self, row: usize, col: usize) -> f64 {
        for i in 0..self.values.len() {
            if self.rows[i] == row && self.cols[i] == col {
                return self.values[i];
            }
        }
        0.0
    }

    /// Number of stored (non-zero) entries.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Fraction of entries that are stored: `nnz / (n_rows * n_cols)`.
    pub fn density(&self) -> f64 {
        let total = self.n_rows as f64 * self.n_cols as f64;
        if total == 0.0 {
            return 0.0;
        }
        self.values.len() as f64 / total
    }

    /// (n_rows, n_cols).
    pub fn shape(&self) -> (usize, usize) {
        (self.n_rows, self.n_cols)
    }

    /// Convert to a dense 2D vector.
    pub fn to_dense(&self) -> Vec<Vec<f64>> {
        let mut dense = vec![vec![0.0; self.n_cols]; self.n_rows];
        for i in 0..self.values.len() {
            dense[self.rows[i]][self.cols[i]] = self.values[i];
        }
        dense
    }

    /// Create a sparse matrix from dense data, storing only values where `|value| > threshold`.
    pub fn from_dense(data: &[Vec<f64>], threshold: f64) -> Self {
        let n_rows = data.len();
        let n_cols = data.first().map_or(0, |r| r.len());
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut values = Vec::new();

        for (r, row) in data.iter().enumerate() {
            for (c, &val) in row.iter().enumerate() {
                if val.abs() > threshold {
                    rows.push(r);
                    cols.push(c);
                    values.push(val);
                }
            }
        }

        Self {
            rows,
            cols,
            values,
            n_rows,
            n_cols,
        }
    }

    /// Number of stored entries in a given row.
    pub fn row_nnz(&self, row: usize) -> usize {
        self.rows.iter().filter(|&&r| r == row).count()
    }

    /// Number of stored entries in a given column.
    pub fn col_nnz(&self, col: usize) -> usize {
        self.cols.iter().filter(|&&c| c == col).count()
    }

    /// Convert COO to CSR format.
    ///
    /// Returns `(data, indices, indptr)` where:
    /// - `data` contains the non-zero values
    /// - `indices` contains the column index for each value
    /// - `indptr[i]..indptr[i+1]` gives the range of entries for row `i`
    pub fn to_csr(&self) -> (Vec<f64>, Vec<usize>, Vec<usize>) {
        // Build sorted order by (row, col)
        let nnz = self.values.len();
        let mut order: Vec<usize> = (0..nnz).collect();
        order.sort_by_key(|&i| (self.rows[i], self.cols[i]));

        let mut data = Vec::with_capacity(nnz);
        let mut indices = Vec::with_capacity(nnz);
        let mut indptr = vec![0usize; self.n_rows + 1];

        for &i in &order {
            data.push(self.values[i]);
            indices.push(self.cols[i]);
            indptr[self.rows[i] + 1] += 1;
        }

        // Cumulative sum to build indptr
        for i in 1..=self.n_rows {
            indptr[i] += indptr[i - 1];
        }

        (data, indices, indptr)
    }

    /// Create a sparse matrix from CSR format.
    ///
    /// - `data` — non-zero values
    /// - `indices` — column index for each value
    /// - `indptr` — row pointer array (length `n_rows + 1`)
    pub fn from_csr(
        data: Vec<f64>,
        indices: Vec<usize>,
        indptr: Vec<usize>,
        n_rows: usize,
        n_cols: usize,
    ) -> Result<Self> {
        if data.len() != indices.len() {
            return Err(CyaneaError::InvalidInput(
                "CSR data and indices must have the same length".into(),
            ));
        }
        if indptr.len() != n_rows + 1 {
            return Err(CyaneaError::InvalidInput(format!(
                "CSR indptr length ({}) must be n_rows + 1 ({})",
                indptr.len(),
                n_rows + 1
            )));
        }

        let nnz = data.len();
        let mut rows = Vec::with_capacity(nnz);
        let mut cols = Vec::with_capacity(nnz);

        for row in 0..n_rows {
            let start = indptr[row];
            let end = indptr[row + 1];
            for idx in start..end {
                if idx >= nnz {
                    return Err(CyaneaError::InvalidInput(format!(
                        "CSR indptr references index {idx} but nnz is {nnz}"
                    )));
                }
                if indices[idx] >= n_cols {
                    return Err(CyaneaError::InvalidInput(format!(
                        "CSR column index {} out of bounds for n_cols={}",
                        indices[idx], n_cols
                    )));
                }
                rows.push(row);
                cols.push(indices[idx]);
            }
        }

        Ok(Self {
            rows,
            cols,
            values: data,
            n_rows,
            n_cols,
        })
    }

    /// Iterate over stored triplets `(row, col, value)`.
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize, f64)> + '_ {
        self.rows
            .iter()
            .zip(self.cols.iter())
            .zip(self.values.iter())
            .map(|((&r, &c), &v)| (r, c, v))
    }
}

impl Summarizable for SparseMatrix {
    fn summary(&self) -> String {
        format!(
            "SparseMatrix: {}\u{00d7}{}, {} nonzeros ({:.2}% density)",
            self.n_rows,
            self.n_cols,
            self.nnz(),
            self.density() * 100.0
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_empty() {
        let m = SparseMatrix::new(10, 20);
        assert_eq!(m.shape(), (10, 20));
        assert_eq!(m.nnz(), 0);
        assert_eq!(m.density(), 0.0);
    }

    #[test]
    fn test_from_triplets() {
        let m = SparseMatrix::from_triplets(
            vec![0, 1, 2],
            vec![0, 1, 2],
            vec![1.0, 2.0, 3.0],
            3,
            3,
        )
        .unwrap();
        assert_eq!(m.nnz(), 3);
        assert_eq!(m.get(0, 0), 1.0);
        assert_eq!(m.get(1, 1), 2.0);
        assert_eq!(m.get(0, 1), 0.0);
    }

    #[test]
    fn test_from_triplets_bounds_check() {
        let result = SparseMatrix::from_triplets(
            vec![5],
            vec![0],
            vec![1.0],
            3,
            3,
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_from_triplets_length_mismatch() {
        let result = SparseMatrix::from_triplets(
            vec![0, 1],
            vec![0],
            vec![1.0],
            3,
            3,
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_insert() {
        let mut m = SparseMatrix::new(3, 3);
        m.insert(0, 0, 5.0).unwrap();
        assert_eq!(m.get(0, 0), 5.0);
        assert_eq!(m.nnz(), 1);

        assert!(m.insert(10, 0, 1.0).is_err());
    }

    #[test]
    fn test_density() {
        let m = SparseMatrix::from_triplets(
            vec![0, 1],
            vec![0, 1],
            vec![1.0, 2.0],
            10,
            10,
        )
        .unwrap();
        assert!((m.density() - 0.02).abs() < 1e-10);
    }

    #[test]
    fn test_to_dense() {
        let m = SparseMatrix::from_triplets(
            vec![0, 1],
            vec![1, 0],
            vec![3.0, 7.0],
            2,
            2,
        )
        .unwrap();
        let dense = m.to_dense();
        assert_eq!(dense, vec![vec![0.0, 3.0], vec![7.0, 0.0]]);
    }

    #[test]
    fn test_from_dense() {
        let data = vec![vec![0.0, 3.0], vec![7.0, 0.0]];
        let m = SparseMatrix::from_dense(&data, 0.0);
        assert_eq!(m.nnz(), 2);
        assert_eq!(m.get(0, 1), 3.0);
        assert_eq!(m.get(1, 0), 7.0);
    }

    #[test]
    fn test_from_dense_with_threshold() {
        let data = vec![vec![0.1, 3.0], vec![7.0, 0.05]];
        let m = SparseMatrix::from_dense(&data, 0.5);
        assert_eq!(m.nnz(), 2); // only 3.0 and 7.0
    }

    #[test]
    fn test_row_col_nnz() {
        let m = SparseMatrix::from_triplets(
            vec![0, 0, 1],
            vec![0, 1, 0],
            vec![1.0, 2.0, 3.0],
            2,
            2,
        )
        .unwrap();
        assert_eq!(m.row_nnz(0), 2);
        assert_eq!(m.row_nnz(1), 1);
        assert_eq!(m.col_nnz(0), 2);
        assert_eq!(m.col_nnz(1), 1);
    }

    #[test]
    fn test_iter() {
        let m = SparseMatrix::from_triplets(
            vec![0, 1],
            vec![0, 1],
            vec![1.0, 2.0],
            2,
            2,
        )
        .unwrap();
        let triplets: Vec<_> = m.iter().collect();
        assert_eq!(triplets, vec![(0, 0, 1.0), (1, 1, 2.0)]);
    }

    #[test]
    fn test_summary() {
        let m = SparseMatrix::from_triplets(
            vec![0],
            vec![0],
            vec![1.0],
            100,
            50,
        )
        .unwrap();
        assert_eq!(
            m.summary(),
            "SparseMatrix: 100\u{00d7}50, 1 nonzeros (0.02% density)"
        );
    }

    #[test]
    fn test_zero_dimension_density() {
        let m = SparseMatrix::new(0, 0);
        assert_eq!(m.density(), 0.0);
    }

    #[test]
    fn test_csr_roundtrip() {
        let m = SparseMatrix::from_triplets(
            vec![0, 0, 1, 2, 2],
            vec![0, 2, 1, 0, 2],
            vec![1.0, 2.0, 3.0, 4.0, 5.0],
            3,
            3,
        )
        .unwrap();

        let (data, indices, indptr) = m.to_csr();
        let m2 = SparseMatrix::from_csr(data, indices, indptr, 3, 3).unwrap();

        assert_eq!(m2.shape(), (3, 3));
        assert_eq!(m2.nnz(), 5);
        assert_eq!(m2.get(0, 0), 1.0);
        assert_eq!(m2.get(0, 2), 2.0);
        assert_eq!(m2.get(1, 1), 3.0);
        assert_eq!(m2.get(2, 0), 4.0);
        assert_eq!(m2.get(2, 2), 5.0);
        assert_eq!(m2.get(1, 0), 0.0);
    }

    #[test]
    fn test_csr_empty() {
        let m = SparseMatrix::new(3, 4);
        let (data, indices, indptr) = m.to_csr();
        assert!(data.is_empty());
        assert!(indices.is_empty());
        assert_eq!(indptr, vec![0, 0, 0, 0]);

        let m2 = SparseMatrix::from_csr(data, indices, indptr, 3, 4).unwrap();
        assert_eq!(m2.nnz(), 0);
        assert_eq!(m2.shape(), (3, 4));
    }

    #[test]
    fn test_csr_single_row() {
        let m = SparseMatrix::from_triplets(
            vec![0, 0, 0],
            vec![0, 2, 4],
            vec![1.0, 2.0, 3.0],
            1,
            5,
        )
        .unwrap();

        let (data, indices, indptr) = m.to_csr();
        assert_eq!(data, vec![1.0, 2.0, 3.0]);
        assert_eq!(indices, vec![0, 2, 4]);
        assert_eq!(indptr, vec![0, 3]);

        let m2 = SparseMatrix::from_csr(data, indices, indptr, 1, 5).unwrap();
        assert_eq!(m2.nnz(), 3);
        assert_eq!(m2.get(0, 0), 1.0);
        assert_eq!(m2.get(0, 2), 2.0);
        assert_eq!(m2.get(0, 4), 3.0);
    }
}
