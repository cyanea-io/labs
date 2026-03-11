//! 10X Genomics MTX sparse matrix support.
//!
//! Provides a CSC (Compressed Sparse Column) format [`CscMatrix`] along with
//! readers/writers for Matrix Market (`.mtx`) files and the 10X Genomics
//! directory layout (`matrix.mtx` + `barcodes.tsv` + `features.tsv`).
//!
//! # Example
//!
//! ```
//! use cyanea_omics::sc_mtx::CscMatrix;
//!
//! let dense = vec![vec![1.0, 0.0, 2.0], vec![0.0, 3.0, 0.0]];
//! let csc = CscMatrix::from_dense(&dense);
//! assert_eq!(csc.shape(), (2, 3));
//! assert_eq!(csc.nnz(), 3);
//! assert_eq!(csc.get(0, 0), 1.0);
//! assert_eq!(csc.get(1, 1), 3.0);
//! ```

use std::collections::HashMap;

use cyanea_core::{CyaneaError, Result};

use crate::sparse::SparseMatrix;

// ── CscMatrix ─────────────────────────────────────────────────────────────

/// A sparse matrix in CSC (Compressed Sparse Column) format.
///
/// This is the native format for 10X Genomics count matrices. Columns
/// represent cells (barcodes) and rows represent features (genes).
#[derive(Debug, Clone)]
pub struct CscMatrix {
    /// Column pointer array (length `n_cols + 1`).
    /// `indptr[j]..indptr[j+1]` gives the range of entries for column `j`.
    indptr: Vec<usize>,
    /// Row indices for each non-zero entry.
    indices: Vec<usize>,
    /// Non-zero values.
    data: Vec<f64>,
    /// Number of rows (features/genes).
    n_rows: usize,
    /// Number of columns (cells/barcodes).
    n_cols: usize,
    /// Row names (gene/feature names), if available.
    row_names: Vec<String>,
    /// Column names (barcode names), if available.
    col_names: Vec<String>,
}

impl CscMatrix {
    /// Create a new empty CSC matrix with the given dimensions.
    pub fn new(n_rows: usize, n_cols: usize) -> Self {
        Self {
            indptr: vec![0; n_cols + 1],
            indices: Vec::new(),
            data: Vec::new(),
            n_rows,
            n_cols,
            row_names: Vec::new(),
            col_names: Vec::new(),
        }
    }

    /// Create a CSC matrix from raw components.
    ///
    /// # Arguments
    /// - `indptr` — column pointer array (length `n_cols + 1`)
    /// - `indices` — row indices for each non-zero entry
    /// - `data` — non-zero values
    /// - `n_rows` — number of rows
    /// - `n_cols` — number of columns
    pub fn from_raw(
        indptr: Vec<usize>,
        indices: Vec<usize>,
        data: Vec<f64>,
        n_rows: usize,
        n_cols: usize,
    ) -> Result<Self> {
        if data.len() != indices.len() {
            return Err(CyaneaError::InvalidInput(
                "CSC data and indices must have the same length".into(),
            ));
        }
        if indptr.len() != n_cols + 1 {
            return Err(CyaneaError::InvalidInput(format!(
                "CSC indptr length ({}) must be n_cols + 1 ({})",
                indptr.len(),
                n_cols + 1
            )));
        }
        for &idx in &indices {
            if idx >= n_rows {
                return Err(CyaneaError::InvalidInput(format!(
                    "CSC row index {} out of bounds for n_rows={}",
                    idx, n_rows
                )));
            }
        }
        Ok(Self {
            indptr,
            indices,
            data,
            n_rows,
            n_cols,
            row_names: Vec::new(),
            col_names: Vec::new(),
        })
    }

    /// Set row names (gene/feature names).
    pub fn set_row_names(&mut self, names: Vec<String>) {
        self.row_names = names;
    }

    /// Set column names (barcode names).
    pub fn set_col_names(&mut self, names: Vec<String>) {
        self.col_names = names;
    }

    /// Row names (gene/feature names).
    pub fn row_names(&self) -> &[String] {
        &self.row_names
    }

    /// Column names (barcode names).
    pub fn col_names(&self) -> &[String] {
        &self.col_names
    }

    /// (n_rows, n_cols).
    pub fn shape(&self) -> (usize, usize) {
        (self.n_rows, self.n_cols)
    }

    /// Number of stored non-zero entries.
    pub fn nnz(&self) -> usize {
        self.data.len()
    }

    /// Number of rows.
    pub fn n_rows(&self) -> usize {
        self.n_rows
    }

    /// Number of columns.
    pub fn n_cols(&self) -> usize {
        self.n_cols
    }

    /// Get the value at `(row, col)`. Returns 0.0 if the entry is not stored.
    pub fn get(&self, row: usize, col: usize) -> f64 {
        if col >= self.n_cols || row >= self.n_rows {
            return 0.0;
        }
        let start = self.indptr[col];
        let end = self.indptr[col + 1];
        for idx in start..end {
            if self.indices[idx] == row {
                return self.data[idx];
            }
        }
        0.0
    }

    /// Density: fraction of stored entries.
    pub fn density(&self) -> f64 {
        let total = self.n_rows as f64 * self.n_cols as f64;
        if total == 0.0 {
            return 0.0;
        }
        self.data.len() as f64 / total
    }

    /// Convert to dense 2D vector (rows × cols).
    ///
    /// Returns an error if the matrix would exceed `max_elements` entries
    /// (default 100 million).
    pub fn to_dense(&self) -> Result<Vec<Vec<f64>>> {
        self.to_dense_with_limit(100_000_000)
    }

    /// Convert to dense with a custom element limit.
    pub fn to_dense_with_limit(&self, max_elements: usize) -> Result<Vec<Vec<f64>>> {
        let total = self.n_rows * self.n_cols;
        if total > max_elements {
            return Err(CyaneaError::InvalidInput(format!(
                "dense matrix would have {} elements, exceeding limit of {}",
                total, max_elements
            )));
        }
        let mut dense = vec![vec![0.0; self.n_cols]; self.n_rows];
        for col in 0..self.n_cols {
            let start = self.indptr[col];
            let end = self.indptr[col + 1];
            for idx in start..end {
                dense[self.indices[idx]][col] = self.data[idx];
            }
        }
        Ok(dense)
    }

    /// Create a CSC matrix from a dense 2D vector, storing only non-zero values.
    pub fn from_dense(dense: &[Vec<f64>]) -> Self {
        let n_rows = dense.len();
        let n_cols = dense.first().map_or(0, |r| r.len());
        let mut indptr = vec![0usize; n_cols + 1];
        let mut indices = Vec::new();
        let mut data = Vec::new();

        for col in 0..n_cols {
            for row in 0..n_rows {
                let val = dense[row][col];
                if val != 0.0 {
                    indices.push(row);
                    data.push(val);
                }
            }
            indptr[col + 1] = data.len();
        }

        Self {
            indptr,
            indices,
            data,
            n_rows,
            n_cols,
            row_names: Vec::new(),
            col_names: Vec::new(),
        }
    }

    /// Convert to a COO-format [`SparseMatrix`].
    pub fn to_coo(&self) -> SparseMatrix {
        let mut rows = Vec::with_capacity(self.data.len());
        let mut cols = Vec::with_capacity(self.data.len());
        let mut values = Vec::with_capacity(self.data.len());

        for col in 0..self.n_cols {
            let start = self.indptr[col];
            let end = self.indptr[col + 1];
            for idx in start..end {
                rows.push(self.indices[idx]);
                cols.push(col);
                values.push(self.data[idx]);
            }
        }

        SparseMatrix::from_triplets(rows, cols, values, self.n_rows, self.n_cols)
            .expect("CSC to COO conversion should always succeed")
    }

    /// Create a CSC matrix from a COO-format [`SparseMatrix`].
    pub fn from_coo(coo: &SparseMatrix) -> Self {
        let (n_rows, n_cols) = coo.shape();

        // Collect triplets and sort by (col, row)
        let mut triplets: Vec<(usize, usize, f64)> = coo.iter().collect();
        triplets.sort_by_key(|&(r, c, _)| (c, r));

        let mut indptr = vec![0usize; n_cols + 1];
        let mut indices = Vec::with_capacity(triplets.len());
        let mut data = Vec::with_capacity(triplets.len());

        for &(row, col, val) in &triplets {
            indices.push(row);
            data.push(val);
            indptr[col + 1] += 1;
        }

        // Cumulative sum
        for j in 1..=n_cols {
            indptr[j] += indptr[j - 1];
        }

        Self {
            indptr,
            indices,
            data,
            n_rows,
            n_cols,
            row_names: Vec::new(),
            col_names: Vec::new(),
        }
    }

    /// Sum of values in each row.
    pub fn row_sums(&self) -> Vec<f64> {
        let mut sums = vec![0.0; self.n_rows];
        for col in 0..self.n_cols {
            let start = self.indptr[col];
            let end = self.indptr[col + 1];
            for idx in start..end {
                sums[self.indices[idx]] += self.data[idx];
            }
        }
        sums
    }

    /// Sum of values in each column.
    pub fn col_sums(&self) -> Vec<f64> {
        let mut sums = vec![0.0; self.n_cols];
        for col in 0..self.n_cols {
            let start = self.indptr[col];
            let end = self.indptr[col + 1];
            for idx in start..end {
                sums[col] += self.data[idx];
            }
        }
        sums
    }

    /// Number of non-zero entries per row.
    pub fn row_nnz(&self) -> Vec<usize> {
        let mut counts = vec![0usize; self.n_rows];
        for &row in &self.indices {
            counts[row] += 1;
        }
        counts
    }

    /// Number of non-zero entries per column.
    pub fn col_nnz(&self) -> Vec<usize> {
        let mut counts = vec![0usize; self.n_cols];
        for col in 0..self.n_cols {
            counts[col] = self.indptr[col + 1] - self.indptr[col];
        }
        counts
    }

    /// Element-wise addition with another CSC matrix of the same shape.
    pub fn add(&self, other: &CscMatrix) -> Result<CscMatrix> {
        if self.n_rows != other.n_rows || self.n_cols != other.n_cols {
            return Err(CyaneaError::InvalidInput(format!(
                "shape mismatch: ({}, {}) vs ({}, {})",
                self.n_rows, self.n_cols, other.n_rows, other.n_cols
            )));
        }

        let mut indptr = vec![0usize; self.n_cols + 1];
        let mut indices = Vec::new();
        let mut data = Vec::new();

        for col in 0..self.n_cols {
            let mut col_entries: HashMap<usize, f64> = HashMap::new();

            let s1 = self.indptr[col];
            let e1 = self.indptr[col + 1];
            for idx in s1..e1 {
                *col_entries.entry(self.indices[idx]).or_insert(0.0) += self.data[idx];
            }

            let s2 = other.indptr[col];
            let e2 = other.indptr[col + 1];
            for idx in s2..e2 {
                *col_entries.entry(other.indices[idx]).or_insert(0.0) += other.data[idx];
            }

            let mut sorted_entries: Vec<(usize, f64)> = col_entries
                .into_iter()
                .filter(|&(_, v)| v != 0.0)
                .collect();
            sorted_entries.sort_by_key(|&(r, _)| r);

            for (row, val) in sorted_entries {
                indices.push(row);
                data.push(val);
            }
            indptr[col + 1] = data.len();
        }

        CscMatrix::from_raw(indptr, indices, data, self.n_rows, self.n_cols)
    }

    /// Element-wise multiplication (Hadamard product) with another CSC matrix.
    pub fn multiply_elementwise(&self, other: &CscMatrix) -> Result<CscMatrix> {
        if self.n_rows != other.n_rows || self.n_cols != other.n_cols {
            return Err(CyaneaError::InvalidInput(format!(
                "shape mismatch: ({}, {}) vs ({}, {})",
                self.n_rows, self.n_cols, other.n_rows, other.n_cols
            )));
        }

        let mut indptr = vec![0usize; self.n_cols + 1];
        let mut indices = Vec::new();
        let mut data = Vec::new();

        for col in 0..self.n_cols {
            // Build a map from the first matrix's column entries
            let s1 = self.indptr[col];
            let e1 = self.indptr[col + 1];
            let mut a_entries: HashMap<usize, f64> = HashMap::new();
            for idx in s1..e1 {
                a_entries.insert(self.indices[idx], self.data[idx]);
            }

            // Intersect with second matrix
            let s2 = other.indptr[col];
            let e2 = other.indptr[col + 1];
            let mut col_results: Vec<(usize, f64)> = Vec::new();
            for idx in s2..e2 {
                let row = other.indices[idx];
                if let Some(&a_val) = a_entries.get(&row) {
                    let product = a_val * other.data[idx];
                    if product != 0.0 {
                        col_results.push((row, product));
                    }
                }
            }
            col_results.sort_by_key(|&(r, _)| r);

            for (row, val) in col_results {
                indices.push(row);
                data.push(val);
            }
            indptr[col + 1] = data.len();
        }

        CscMatrix::from_raw(indptr, indices, data, self.n_rows, self.n_cols)
    }

    /// Iterate over stored entries as `(row, col, value)`.
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize, f64)> + '_ {
        (0..self.n_cols).flat_map(move |col| {
            let start = self.indptr[col];
            let end = self.indptr[col + 1];
            (start..end).map(move |idx| (self.indices[idx], col, self.data[idx]))
        })
    }
}

// ── MTX I/O ────────────────────────────────────────────────────────────────

/// Parse a Matrix Market (`.mtx`) format string into a [`CscMatrix`].
///
/// Supports the `%%MatrixMarket matrix coordinate real general` format
/// (integer values are accepted and converted to f64).
pub fn read_mtx(content: &str) -> Result<CscMatrix> {
    let mut lines = content.lines().peekable();

    // Parse header
    let header = lines.next().ok_or_else(|| {
        CyaneaError::Parse("empty MTX file".into())
    })?;

    if !header.starts_with("%%MatrixMarket") {
        return Err(CyaneaError::Parse(
            "invalid MTX header: must start with %%MatrixMarket".into(),
        ));
    }

    let header_lower = header.to_lowercase();
    let is_pattern = header_lower.contains("pattern");
    let is_symmetric = header_lower.contains("symmetric");

    // Skip comments
    while let Some(line) = lines.peek() {
        if line.starts_with('%') {
            lines.next();
        } else {
            break;
        }
    }

    // Parse dimensions
    let dim_line = lines.next().ok_or_else(|| {
        CyaneaError::Parse("missing dimension line in MTX file".into())
    })?;
    let dims: Vec<&str> = dim_line.split_whitespace().collect();
    if dims.len() < 3 {
        return Err(CyaneaError::Parse(format!(
            "invalid dimension line: '{}'",
            dim_line
        )));
    }
    let n_rows: usize = dims[0].parse().map_err(|_| {
        CyaneaError::Parse(format!("invalid n_rows: '{}'", dims[0]))
    })?;
    let n_cols: usize = dims[1].parse().map_err(|_| {
        CyaneaError::Parse(format!("invalid n_cols: '{}'", dims[1]))
    })?;
    let _nnz: usize = dims[2].parse().map_err(|_| {
        CyaneaError::Parse(format!("invalid nnz: '{}'", dims[2]))
    })?;

    // Parse triplets (1-based → 0-based)
    let mut coo_rows = Vec::new();
    let mut coo_cols = Vec::new();
    let mut coo_vals = Vec::new();

    for line in lines {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        }

        let row: usize = parts[0].parse().map_err(|_| {
            CyaneaError::Parse(format!("invalid row index: '{}'", parts[0]))
        })?;
        let col: usize = parts[1].parse().map_err(|_| {
            CyaneaError::Parse(format!("invalid col index: '{}'", parts[1]))
        })?;

        let val = if is_pattern {
            1.0
        } else if parts.len() >= 3 {
            parts[2].parse::<f64>().map_err(|_| {
                CyaneaError::Parse(format!("invalid value: '{}'", parts[2]))
            })?
        } else {
            1.0
        };

        // Convert to 0-based
        let r = row.checked_sub(1).ok_or_else(|| {
            CyaneaError::Parse("MTX row index must be >= 1".into())
        })?;
        let c = col.checked_sub(1).ok_or_else(|| {
            CyaneaError::Parse("MTX col index must be >= 1".into())
        })?;

        coo_rows.push(r);
        coo_cols.push(c);
        coo_vals.push(val);

        // Handle symmetric matrices
        if is_symmetric && r != c {
            coo_rows.push(c);
            coo_cols.push(r);
            coo_vals.push(val);
        }
    }

    // Build CSC from COO triplets
    // Sort by (col, row)
    let nnz = coo_vals.len();
    let mut order: Vec<usize> = (0..nnz).collect();
    order.sort_by_key(|&i| (coo_cols[i], coo_rows[i]));

    let mut indptr = vec![0usize; n_cols + 1];
    let mut indices = Vec::with_capacity(nnz);
    let mut data = Vec::with_capacity(nnz);

    for &i in &order {
        indices.push(coo_rows[i]);
        data.push(coo_vals[i]);
        indptr[coo_cols[i] + 1] += 1;
    }

    for j in 1..=n_cols {
        indptr[j] += indptr[j - 1];
    }

    Ok(CscMatrix {
        indptr,
        indices,
        data,
        n_rows,
        n_cols,
        row_names: Vec::new(),
        col_names: Vec::new(),
    })
}

/// Write a [`CscMatrix`] to Matrix Market format string.
pub fn write_mtx(matrix: &CscMatrix) -> String {
    let mut out = String::new();
    out.push_str("%%MatrixMarket matrix coordinate real general\n");
    out.push_str(&format!("{} {} {}\n", matrix.n_rows, matrix.n_cols, matrix.nnz()));

    for col in 0..matrix.n_cols {
        let start = matrix.indptr[col];
        let end = matrix.indptr[col + 1];
        for idx in start..end {
            // 1-based indices
            out.push_str(&format!(
                "{} {} {}\n",
                matrix.indices[idx] + 1,
                col + 1,
                matrix.data[idx]
            ));
        }
    }

    out
}

/// Read a 10X Genomics directory containing `matrix.mtx`, `barcodes.tsv`,
/// and `features.tsv` (or `genes.tsv`).
pub fn read_10x_directory(
    mtx_content: &str,
    barcodes_content: &str,
    features_content: &str,
) -> Result<CscMatrix> {
    let mut matrix = read_mtx(mtx_content)?;

    // Parse barcodes (one per line, may have tab-separated fields)
    let barcodes: Vec<String> = barcodes_content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split('\t').next().unwrap_or("").to_string())
        .collect();

    // Parse features/genes (tab-separated: id, name, type)
    let features: Vec<String> = features_content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let parts: Vec<&str> = l.split('\t').collect();
            if parts.len() >= 2 {
                parts[1].to_string() // gene name
            } else {
                parts[0].to_string() // gene id
            }
        })
        .collect();

    // In 10X format: rows = genes/features, cols = barcodes/cells
    if features.len() == matrix.n_rows {
        matrix.set_row_names(features);
    }
    if barcodes.len() == matrix.n_cols {
        matrix.set_col_names(barcodes);
    }

    Ok(matrix)
}

// ── Filtering ──────────────────────────────────────────────────────────────

/// Configuration for cell filtering.
#[derive(Debug, Clone)]
pub struct CellFilterConfig {
    /// Minimum number of genes expressed per cell.
    pub min_genes: Option<usize>,
    /// Maximum number of genes expressed per cell.
    pub max_genes: Option<usize>,
    /// Maximum fraction of mitochondrial gene counts (0.0–1.0).
    pub max_mito_pct: Option<f64>,
    /// Prefix for mitochondrial gene names (default: "MT-").
    pub mito_prefix: String,
}

impl Default for CellFilterConfig {
    fn default() -> Self {
        Self {
            min_genes: Some(200),
            max_genes: Some(5000),
            max_mito_pct: Some(0.2),
            mito_prefix: "MT-".into(),
        }
    }
}

/// Filter cells from a CSC matrix based on QC criteria.
pub fn filter_cells(matrix: &CscMatrix, config: &CellFilterConfig) -> (CscMatrix, Vec<usize>) {
    let col_nnz = matrix.col_nnz();

    // Identify mitochondrial genes
    let mito_rows: Vec<bool> = if !matrix.row_names.is_empty() {
        matrix
            .row_names
            .iter()
            .map(|name| name.starts_with(&config.mito_prefix))
            .collect()
    } else {
        vec![false; matrix.n_rows]
    };

    // Compute per-cell total counts
    let col_sums = matrix.col_sums();

    let mut kept_cols: Vec<usize> = Vec::new();

    for col in 0..matrix.n_cols {
        let n_genes = col_nnz[col];

        // Min genes filter
        if let Some(min) = config.min_genes {
            if n_genes < min {
                continue;
            }
        }
        // Max genes filter
        if let Some(max) = config.max_genes {
            if n_genes > max {
                continue;
            }
        }
        // Mito fraction filter
        if let Some(max_mito) = config.max_mito_pct {
            let total = col_sums[col];
            if total > 0.0 {
                let start = matrix.indptr[col];
                let end = matrix.indptr[col + 1];
                let mut mito_sum = 0.0;
                for idx in start..end {
                    if mito_rows[matrix.indices[idx]] {
                        mito_sum += matrix.data[idx];
                    }
                }
                if mito_sum / total > max_mito {
                    continue;
                }
            }
        }

        kept_cols.push(col);
    }

    let filtered = subset_cols(matrix, &kept_cols);
    (filtered, kept_cols)
}

/// Configuration for gene filtering.
#[derive(Debug, Clone)]
pub struct GeneFilterConfig {
    /// Minimum number of cells expressing this gene.
    pub min_cells: usize,
}

impl Default for GeneFilterConfig {
    fn default() -> Self {
        Self { min_cells: 3 }
    }
}

/// Filter genes from a CSC matrix based on the number of cells expressing them.
pub fn filter_genes(matrix: &CscMatrix, config: &GeneFilterConfig) -> (CscMatrix, Vec<usize>) {
    let row_nnz = matrix.row_nnz();

    let kept_rows: Vec<usize> = (0..matrix.n_rows)
        .filter(|&r| row_nnz[r] >= config.min_cells)
        .collect();

    let filtered = subset_rows(matrix, &kept_rows);
    (filtered, kept_rows)
}

/// Subset a CSC matrix to keep only the specified columns.
fn subset_cols(matrix: &CscMatrix, cols: &[usize]) -> CscMatrix {
    let new_n_cols = cols.len();
    let mut indptr = vec![0usize; new_n_cols + 1];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut col_names = Vec::new();

    for (new_col, &old_col) in cols.iter().enumerate() {
        let start = matrix.indptr[old_col];
        let end = matrix.indptr[old_col + 1];
        for idx in start..end {
            indices.push(matrix.indices[idx]);
            data.push(matrix.data[idx]);
        }
        indptr[new_col + 1] = data.len();

        if old_col < matrix.col_names.len() {
            col_names.push(matrix.col_names[old_col].clone());
        }
    }

    CscMatrix {
        indptr,
        indices,
        data,
        n_rows: matrix.n_rows,
        n_cols: new_n_cols,
        row_names: matrix.row_names.clone(),
        col_names,
    }
}

/// Subset a CSC matrix to keep only the specified rows.
fn subset_rows(matrix: &CscMatrix, rows: &[usize]) -> CscMatrix {
    let new_n_rows = rows.len();
    let row_map: HashMap<usize, usize> = rows.iter().enumerate().map(|(new, &old)| (old, new)).collect();

    let mut indptr = vec![0usize; matrix.n_cols + 1];
    let mut indices = Vec::new();
    let mut data = Vec::new();

    for col in 0..matrix.n_cols {
        let start = matrix.indptr[col];
        let end = matrix.indptr[col + 1];
        for idx in start..end {
            if let Some(&new_row) = row_map.get(&matrix.indices[idx]) {
                indices.push(new_row);
                data.push(matrix.data[idx]);
            }
        }
        indptr[col + 1] = data.len();
    }

    let row_names: Vec<String> = rows
        .iter()
        .filter_map(|&r| matrix.row_names.get(r).cloned())
        .collect();

    CscMatrix {
        indptr,
        indices,
        data,
        n_rows: new_n_rows,
        n_cols: matrix.n_cols,
        row_names,
        col_names: matrix.col_names.clone(),
    }
}

impl cyanea_core::Summarizable for CscMatrix {
    fn summary(&self) -> String {
        format!(
            "CscMatrix: {}\u{00d7}{}, {} nonzeros ({:.2}% density)",
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
    fn new_empty() {
        let m = CscMatrix::new(3, 4);
        assert_eq!(m.shape(), (3, 4));
        assert_eq!(m.nnz(), 0);
        assert_eq!(m.get(0, 0), 0.0);
    }

    #[test]
    fn from_dense_roundtrip() {
        let dense = vec![
            vec![1.0, 0.0, 2.0],
            vec![0.0, 3.0, 0.0],
            vec![4.0, 0.0, 5.0],
        ];
        let csc = CscMatrix::from_dense(&dense);
        assert_eq!(csc.shape(), (3, 3));
        assert_eq!(csc.nnz(), 5);
        let roundtrip = csc.to_dense().unwrap();
        assert_eq!(roundtrip, dense);
    }

    #[test]
    fn row_col_sums() {
        let dense = vec![
            vec![1.0, 0.0, 2.0],
            vec![0.0, 3.0, 0.0],
        ];
        let csc = CscMatrix::from_dense(&dense);
        assert_eq!(csc.row_sums(), vec![3.0, 3.0]);
        assert_eq!(csc.col_sums(), vec![1.0, 3.0, 2.0]);
    }

    #[test]
    fn read_mtx_basic() {
        let mtx = "\
%%MatrixMarket matrix coordinate real general
% comment
3 3 4
1 1 1.0
1 3 2.0
2 2 3.0
3 1 4.0
";
        let m = read_mtx(mtx).unwrap();
        assert_eq!(m.shape(), (3, 3));
        assert_eq!(m.nnz(), 4);
    }

    #[test]
    fn write_mtx_roundtrip() {
        let dense = vec![
            vec![1.0, 0.0, 2.0],
            vec![0.0, 3.0, 0.0],
        ];
        let csc = CscMatrix::from_dense(&dense);
        let written = write_mtx(&csc);
        let parsed = read_mtx(&written).unwrap();
        assert_eq!(parsed.shape(), (2, 3));
    }

    #[test]
    fn filter_cells_min_genes() {
        let dense = vec![
            vec![1.0, 1.0, 1.0],
            vec![1.0, 0.0, 1.0],
            vec![0.0, 0.0, 1.0],
        ];
        let csc = CscMatrix::from_dense(&dense);
        let config = CellFilterConfig {
            min_genes: Some(2),
            max_genes: None,
            max_mito_pct: None,
            mito_prefix: "MT-".into(),
        };
        let (filtered, kept) = filter_cells(&csc, &config);
        assert_eq!(kept, vec![0, 2]);
        assert_eq!(filtered.n_cols(), 2);
    }

    #[test]
    fn filter_genes_min_cells() {
        let dense = vec![
            vec![1.0, 1.0, 1.0],
            vec![1.0, 0.0, 0.0],
            vec![1.0, 1.0, 0.0],
        ];
        let csc = CscMatrix::from_dense(&dense);
        let config = GeneFilterConfig { min_cells: 2 };
        let (filtered, kept) = filter_genes(&csc, &config);
        assert_eq!(kept, vec![0, 2]);
        assert_eq!(filtered.n_rows(), 2);
    }

    #[test]
    fn density_calculation() {
        let dense = vec![
            vec![1.0, 0.0],
            vec![0.0, 2.0],
        ];
        let csc = CscMatrix::from_dense(&dense);
        assert!((csc.density() - 0.5).abs() < 1e-10);
    }
}
