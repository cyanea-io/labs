//! AnnData-like container for single-cell omics data.
//!
//! Provides an in-memory representation inspired by the Python AnnData format,
//! the standard data structure in the scverse ecosystem (scanpy, scvi-tools).
//!
//! # Structure
//!
//! - `X` — primary data matrix (cells × genes), dense or sparse
//! - `obs` — per-cell metadata (string key-value pairs)
//! - `var` — per-gene metadata
//! - `obsm` / `varm` — multi-dimensional annotations (e.g. PCA, UMAP embeddings)
//! - `layers` — alternative data matrices (e.g. raw counts, normalized)
//!
//! # Example
//!
//! ```
//! use cyanea_omics::single_cell::{AnnData, MatrixData};
//!
//! let x = MatrixData::Dense(vec![vec![1.0, 2.0], vec![3.0, 4.0]]);
//! let adata = AnnData::new(
//!     x,
//!     vec!["cell_1".into(), "cell_2".into()],
//!     vec!["gene_a".into(), "gene_b".into()],
//! ).unwrap();
//! assert_eq!(adata.n_obs(), 2);
//! assert_eq!(adata.n_vars(), 2);
//! ```

use std::collections::HashMap;

use cyanea_core::{CyaneaError, Result, Summarizable};

use crate::sparse::SparseMatrix;

/// The primary data matrix, either dense or sparse.
#[derive(Debug, Clone)]
pub enum MatrixData {
    /// Dense row-major matrix (n_obs × n_vars).
    Dense(Vec<Vec<f64>>),
    /// Sparse COO matrix.
    Sparse(SparseMatrix),
}

impl MatrixData {
    /// (n_obs, n_vars).
    pub fn shape(&self) -> (usize, usize) {
        match self {
            MatrixData::Dense(rows) => {
                let n_obs = rows.len();
                let n_vars = rows.first().map_or(0, |r| r.len());
                (n_obs, n_vars)
            }
            MatrixData::Sparse(s) => s.shape(),
        }
    }

    /// Get a value at (obs_idx, var_idx).
    pub fn get(&self, obs: usize, var: usize) -> f64 {
        match self {
            MatrixData::Dense(rows) => {
                rows.get(obs).and_then(|r| r.get(var)).copied().unwrap_or(0.0)
            }
            MatrixData::Sparse(s) => s.get(obs, var),
        }
    }

    /// Set a value at (obs_idx, var_idx).
    ///
    /// For dense matrices, sets directly. For sparse matrices, inserts a triplet
    /// (does not deduplicate — the last-inserted value wins on `get()`).
    pub fn set(&mut self, obs: usize, var: usize, val: f64) {
        match self {
            MatrixData::Dense(rows) => {
                if let Some(row) = rows.get_mut(obs) {
                    if let Some(cell) = row.get_mut(var) {
                        *cell = val;
                    }
                }
            }
            MatrixData::Sparse(s) => {
                let _ = s.insert(obs, var, val);
            }
        }
    }

    /// Sum of values in each column.
    pub fn column_sums(&self) -> Vec<f64> {
        match self {
            MatrixData::Dense(rows) => {
                let n_vars = rows.first().map_or(0, |r| r.len());
                let mut sums = vec![0.0; n_vars];
                for row in rows {
                    for (j, &v) in row.iter().enumerate() {
                        sums[j] += v;
                    }
                }
                sums
            }
            MatrixData::Sparse(s) => s.column_sums(),
        }
    }

    /// Mean value of each column.
    pub fn column_means(&self) -> Vec<f64> {
        match self {
            MatrixData::Dense(rows) => {
                let n_obs = rows.len();
                if n_obs == 0 {
                    return vec![];
                }
                let sums = self.column_sums();
                let n = n_obs as f64;
                sums.into_iter().map(|s| s / n).collect()
            }
            MatrixData::Sparse(s) => s.column_means(),
        }
    }

    /// Sum of values in each row.
    pub fn row_sums(&self) -> Vec<f64> {
        match self {
            MatrixData::Dense(rows) => rows.iter().map(|r| r.iter().sum()).collect(),
            MatrixData::Sparse(s) => s.row_sums(),
        }
    }

    /// Flatten to a row-major `Vec<f64>` for interop with cyanea-stats/ml.
    pub fn to_flat_row_major(&self) -> Vec<f64> {
        let (n_obs, n_vars) = self.shape();
        match self {
            MatrixData::Dense(rows) => {
                let mut flat = Vec::with_capacity(n_obs * n_vars);
                for row in rows {
                    flat.extend_from_slice(row);
                }
                flat
            }
            MatrixData::Sparse(s) => {
                let mut flat = vec![0.0; n_obs * n_vars];
                for (r, c, v) in s.iter() {
                    flat[r * n_vars + c] = v;
                }
                flat
            }
        }
    }
}

/// A metadata column with typed data.
///
/// Supports string, numeric, and categorical columns as found in `.h5ad` files.
#[derive(Debug, Clone, PartialEq)]
pub enum ColumnData {
    /// Free-text string values.
    Strings(Vec<String>),
    /// Numeric (f64) values.
    Numeric(Vec<f64>),
    /// Categorical data stored as integer codes indexing into a category list.
    Categorical {
        codes: Vec<i32>,
        categories: Vec<String>,
    },
}

impl ColumnData {
    /// Number of elements in this column.
    pub fn len(&self) -> usize {
        match self {
            ColumnData::Strings(v) => v.len(),
            ColumnData::Numeric(v) => v.len(),
            ColumnData::Categorical { codes, .. } => codes.len(),
        }
    }

    /// Whether the column is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Try to get as string slice. Returns `None` if not `Strings` variant.
    pub fn as_strings(&self) -> Option<&Vec<String>> {
        match self {
            ColumnData::Strings(v) => Some(v),
            _ => None,
        }
    }

    /// Try to get as numeric slice. Returns `None` if not `Numeric` variant.
    pub fn as_numeric(&self) -> Option<&Vec<f64>> {
        match self {
            ColumnData::Numeric(v) => Some(v),
            _ => None,
        }
    }

    /// Subset to the given indices.
    fn subset(&self, indices: &[usize]) -> Self {
        match self {
            ColumnData::Strings(v) => {
                ColumnData::Strings(indices.iter().map(|&i| v[i].clone()).collect())
            }
            ColumnData::Numeric(v) => {
                ColumnData::Numeric(indices.iter().map(|&i| v[i]).collect())
            }
            ColumnData::Categorical { codes, categories } => ColumnData::Categorical {
                codes: indices.iter().map(|&i| codes[i]).collect(),
                categories: categories.clone(),
            },
        }
    }
}

/// Per-cell or per-gene quality control metrics.
#[derive(Debug, Clone)]
pub struct QcMetrics {
    /// Total counts per observation.
    pub total_counts: Vec<f64>,
    /// Number of non-zero features per observation.
    pub n_features: Vec<usize>,
}

/// AnnData-like container for single-cell data.
#[derive(Debug, Clone)]
pub struct AnnData {
    /// Primary data matrix (n_obs × n_vars).
    x: MatrixData,
    /// Observation (cell) names.
    obs_names: Vec<String>,
    /// Variable (gene) names.
    var_names: Vec<String>,
    /// Per-cell metadata.
    obs: HashMap<String, ColumnData>,
    /// Per-gene metadata.
    var: HashMap<String, ColumnData>,
    /// Multi-dimensional observation annotations (e.g. PCA embeddings).
    obsm: HashMap<String, Vec<Vec<f64>>>,
    /// Multi-dimensional variable annotations.
    varm: HashMap<String, Vec<Vec<f64>>>,
    /// Alternative data layers (same shape as X).
    layers: HashMap<String, MatrixData>,
    /// Pairwise observation annotations (e.g. kNN graphs, distance matrices).
    obsp: HashMap<String, SparseMatrix>,
    /// Unstructured metadata (free-form key-value pairs).
    uns: HashMap<String, String>,
}

impl AnnData {
    /// Create a new AnnData container.
    ///
    /// # Errors
    ///
    /// Returns an error if the matrix dimensions don't match the name vectors.
    pub fn new(
        x: MatrixData,
        obs_names: Vec<String>,
        var_names: Vec<String>,
    ) -> Result<Self> {
        let (n_obs, n_vars) = x.shape();
        if obs_names.len() != n_obs {
            return Err(CyaneaError::InvalidInput(format!(
                "obs_names length ({}) does not match n_obs ({})",
                obs_names.len(),
                n_obs
            )));
        }
        if var_names.len() != n_vars {
            return Err(CyaneaError::InvalidInput(format!(
                "var_names length ({}) does not match n_vars ({})",
                var_names.len(),
                n_vars
            )));
        }

        Ok(Self {
            x,
            obs_names,
            var_names,
            obs: HashMap::new(),
            var: HashMap::new(),
            obsm: HashMap::new(),
            varm: HashMap::new(),
            layers: HashMap::new(),
            obsp: HashMap::new(),
            uns: HashMap::new(),
        })
    }

    /// Number of observations (cells).
    pub fn n_obs(&self) -> usize {
        self.obs_names.len()
    }

    /// Number of variables (genes).
    pub fn n_vars(&self) -> usize {
        self.var_names.len()
    }

    /// Shape of the primary data matrix.
    pub fn shape(&self) -> (usize, usize) {
        self.x.shape()
    }

    /// Access the primary data matrix.
    pub fn x(&self) -> &MatrixData {
        &self.x
    }

    /// Observation names.
    pub fn obs_names(&self) -> &[String] {
        &self.obs_names
    }

    /// Variable names.
    pub fn var_names(&self) -> &[String] {
        &self.var_names
    }

    /// Add a per-cell string metadata column.
    pub fn add_obs(&mut self, key: &str, values: Vec<String>) -> Result<()> {
        self.add_obs_column(key, ColumnData::Strings(values))
    }

    /// Add a per-cell numeric metadata column.
    pub fn add_obs_numeric(&mut self, key: &str, values: Vec<f64>) -> Result<()> {
        self.add_obs_column(key, ColumnData::Numeric(values))
    }

    /// Add a per-cell metadata column of any type.
    pub fn add_obs_column(&mut self, key: &str, data: ColumnData) -> Result<()> {
        if data.len() != self.n_obs() {
            return Err(CyaneaError::InvalidInput(format!(
                "obs '{}' length ({}) does not match n_obs ({})",
                key,
                data.len(),
                self.n_obs()
            )));
        }
        self.obs.insert(key.to_string(), data);
        Ok(())
    }

    /// Get per-cell metadata column as typed data.
    pub fn get_obs(&self, key: &str) -> Option<&ColumnData> {
        self.obs.get(key)
    }

    /// Get per-cell metadata column as strings (convenience for backward compat).
    pub fn get_obs_strings(&self, key: &str) -> Option<&Vec<String>> {
        self.obs.get(key).and_then(|c| c.as_strings())
    }

    /// All observation metadata columns.
    pub fn obs_columns(&self) -> &HashMap<String, ColumnData> {
        &self.obs
    }

    /// Add a per-gene string metadata column.
    pub fn add_var(&mut self, key: &str, values: Vec<String>) -> Result<()> {
        self.add_var_column(key, ColumnData::Strings(values))
    }

    /// Add a per-gene numeric metadata column.
    pub fn add_var_numeric(&mut self, key: &str, values: Vec<f64>) -> Result<()> {
        self.add_var_column(key, ColumnData::Numeric(values))
    }

    /// Add a per-gene metadata column of any type.
    pub fn add_var_column(&mut self, key: &str, data: ColumnData) -> Result<()> {
        if data.len() != self.n_vars() {
            return Err(CyaneaError::InvalidInput(format!(
                "var '{}' length ({}) does not match n_vars ({})",
                key,
                data.len(),
                self.n_vars()
            )));
        }
        self.var.insert(key.to_string(), data);
        Ok(())
    }

    /// Get per-gene metadata column as typed data.
    pub fn get_var(&self, key: &str) -> Option<&ColumnData> {
        self.var.get(key)
    }

    /// Get per-gene metadata column as strings (convenience for backward compat).
    pub fn get_var_strings(&self, key: &str) -> Option<&Vec<String>> {
        self.var.get(key).and_then(|c| c.as_strings())
    }

    /// All variable metadata columns.
    pub fn var_columns(&self) -> &HashMap<String, ColumnData> {
        &self.var
    }

    /// Add a multi-dimensional observation annotation (e.g. PCA embedding).
    pub fn add_obsm(&mut self, key: &str, data: Vec<Vec<f64>>) -> Result<()> {
        if data.len() != self.n_obs() {
            return Err(CyaneaError::InvalidInput(format!(
                "obsm '{}' length ({}) does not match n_obs ({})",
                key,
                data.len(),
                self.n_obs()
            )));
        }
        self.obsm.insert(key.to_string(), data);
        Ok(())
    }

    /// Get a multi-dimensional observation annotation.
    pub fn get_obsm(&self, key: &str) -> Option<&Vec<Vec<f64>>> {
        self.obsm.get(key)
    }

    /// Add a multi-dimensional variable annotation.
    pub fn add_varm(&mut self, key: &str, data: Vec<Vec<f64>>) -> Result<()> {
        if data.len() != self.n_vars() {
            return Err(CyaneaError::InvalidInput(format!(
                "varm '{}' length ({}) does not match n_vars ({})",
                key,
                data.len(),
                self.n_vars()
            )));
        }
        self.varm.insert(key.to_string(), data);
        Ok(())
    }

    /// Get a multi-dimensional variable annotation.
    pub fn get_varm(&self, key: &str) -> Option<&Vec<Vec<f64>>> {
        self.varm.get(key)
    }

    /// Add an alternative data layer.
    pub fn add_layer(&mut self, key: &str, layer: MatrixData) -> Result<()> {
        let (n_obs, n_vars) = layer.shape();
        if n_obs != self.n_obs() || n_vars != self.n_vars() {
            return Err(CyaneaError::InvalidInput(format!(
                "layer '{}' shape ({}, {}) does not match ({}, {})",
                key,
                n_obs,
                n_vars,
                self.n_obs(),
                self.n_vars()
            )));
        }
        self.layers.insert(key.to_string(), layer);
        Ok(())
    }

    /// Get an alternative data layer.
    pub fn get_layer(&self, key: &str) -> Option<&MatrixData> {
        self.layers.get(key)
    }

    /// All observation multi-dimensional annotations.
    pub fn obsm_keys(&self) -> &HashMap<String, Vec<Vec<f64>>> {
        &self.obsm
    }

    /// All variable multi-dimensional annotations.
    pub fn varm_keys(&self) -> &HashMap<String, Vec<Vec<f64>>> {
        &self.varm
    }

    /// All alternative data layers.
    pub fn layers_keys(&self) -> &HashMap<String, MatrixData> {
        &self.layers
    }

    /// Mutable access to the primary data matrix.
    pub fn x_mut(&mut self) -> &mut MatrixData {
        &mut self.x
    }

    /// Replace the primary data matrix. The new matrix must have the same shape.
    pub fn set_x(&mut self, new_x: MatrixData) -> Result<()> {
        let (n_obs, n_vars) = new_x.shape();
        if n_obs != self.n_obs() || n_vars != self.n_vars() {
            return Err(CyaneaError::InvalidInput(format!(
                "new X shape ({}, {}) does not match ({}, {})",
                n_obs, n_vars, self.n_obs(), self.n_vars()
            )));
        }
        self.x = new_x;
        Ok(())
    }

    /// Subset to the given variable (gene) indices.
    pub fn subset_vars(&self, indices: &[usize]) -> Result<AnnData> {
        for &i in indices {
            if i >= self.n_vars() {
                return Err(CyaneaError::InvalidInput(format!(
                    "var index {} out of bounds (n_vars={})",
                    i, self.n_vars()
                )));
            }
        }

        let x = subset_matrix_cols(&self.x, indices, self.n_obs());
        let var_names: Vec<String> = indices.iter().map(|&i| self.var_names[i].clone()).collect();

        let mut adata = AnnData::new(x, self.obs_names.clone(), var_names)?;

        // Copy obs metadata
        adata.obs = self.obs.clone();
        // Subset var metadata
        for (key, col) in &self.var {
            adata.var.insert(key.clone(), col.subset(indices));
        }
        // Copy obsm/varm
        adata.obsm = self.obsm.clone();
        // Subset layers
        for (key, layer) in &self.layers {
            let sub = subset_matrix_cols(layer, indices, self.n_obs());
            adata.layers.insert(key.clone(), sub);
        }
        // Copy obsp and uns
        adata.obsp = self.obsp.clone();
        adata.uns = self.uns.clone();

        Ok(adata)
    }

    /// Add a pairwise observation annotation (e.g. kNN graph).
    ///
    /// The matrix must be n_obs × n_obs.
    pub fn add_obsp(&mut self, key: &str, matrix: SparseMatrix) -> Result<()> {
        let (r, c) = matrix.shape();
        if r != self.n_obs() || c != self.n_obs() {
            return Err(CyaneaError::InvalidInput(format!(
                "obsp '{}' shape ({}, {}) does not match n_obs ({})",
                key, r, c, self.n_obs()
            )));
        }
        self.obsp.insert(key.to_string(), matrix);
        Ok(())
    }

    /// Get a pairwise observation annotation.
    pub fn get_obsp(&self, key: &str) -> Option<&SparseMatrix> {
        self.obsp.get(key)
    }

    /// Add unstructured metadata.
    pub fn add_uns(&mut self, key: &str, value: String) {
        self.uns.insert(key.to_string(), value);
    }

    /// Get unstructured metadata.
    pub fn get_uns(&self, key: &str) -> Option<&str> {
        self.uns.get(key).map(|s| s.as_str())
    }

    /// All pairwise observation annotation keys.
    pub fn obsp_keys(&self) -> &HashMap<String, SparseMatrix> {
        &self.obsp
    }

    /// All unstructured metadata.
    pub fn uns_keys(&self) -> &HashMap<String, String> {
        &self.uns
    }

    /// Mutable access to a layer.
    pub fn get_layer_mut(&mut self, key: &str) -> Option<&mut MatrixData> {
        self.layers.get_mut(key)
    }

    /// Subset to the given observation indices.
    pub fn subset_obs(&self, indices: &[usize]) -> Result<AnnData> {
        for &i in indices {
            if i >= self.n_obs() {
                return Err(CyaneaError::InvalidInput(format!(
                    "obs index {} out of bounds (n_obs={})",
                    i,
                    self.n_obs()
                )));
            }
        }

        let x = subset_matrix_rows(&self.x, indices, self.n_vars());
        let obs_names: Vec<String> = indices.iter().map(|&i| self.obs_names[i].clone()).collect();

        let mut adata = AnnData::new(x, obs_names, self.var_names.clone())?;

        // Subset obs metadata
        for (key, col) in &self.obs {
            adata.obs.insert(key.clone(), col.subset(indices));
        }
        // Copy var metadata
        adata.var = self.var.clone();

        // Subset obsm
        for (key, data) in &self.obsm {
            let sub: Vec<Vec<f64>> = indices.iter().map(|&i| data[i].clone()).collect();
            adata.obsm.insert(key.clone(), sub);
        }
        adata.varm = self.varm.clone();
        adata.uns = self.uns.clone();

        Ok(adata)
    }

    /// Compute basic QC metrics.
    pub fn qc_metrics(&self) -> QcMetrics {
        let n = self.n_obs();
        let p = self.n_vars();
        let mut total_counts = vec![0.0; n];
        let mut n_features = vec![0usize; n];

        for i in 0..n {
            for j in 0..p {
                let v = self.x.get(i, j);
                total_counts[i] += v;
                if v > 0.0 {
                    n_features[i] += 1;
                }
            }
        }

        QcMetrics {
            total_counts,
            n_features,
        }
    }
}

fn subset_matrix_cols(x: &MatrixData, col_indices: &[usize], n_obs: usize) -> MatrixData {
    match x {
        MatrixData::Dense(rows) => {
            let sub: Vec<Vec<f64>> = rows
                .iter()
                .map(|row| col_indices.iter().map(|&j| row[j]).collect())
                .collect();
            MatrixData::Dense(sub)
        }
        MatrixData::Sparse(s) => {
            let n_new_cols = col_indices.len();
            let mut new_s = SparseMatrix::new(n_obs, n_new_cols);
            // Build reverse map: old_col -> new_col
            let mut col_map = HashMap::new();
            for (new_j, &old_j) in col_indices.iter().enumerate() {
                col_map.insert(old_j, new_j);
            }
            for (r, c, v) in s.iter() {
                if let Some(&new_c) = col_map.get(&c) {
                    let _ = new_s.insert(r, new_c, v);
                }
            }
            MatrixData::Sparse(new_s)
        }
    }
}

fn subset_matrix_rows(x: &MatrixData, indices: &[usize], n_vars: usize) -> MatrixData {
    match x {
        MatrixData::Dense(rows) => {
            let sub: Vec<Vec<f64>> = indices.iter().map(|&i| rows[i].clone()).collect();
            MatrixData::Dense(sub)
        }
        MatrixData::Sparse(s) => {
            let n_new = indices.len();
            let mut new_s = SparseMatrix::new(n_new, n_vars);
            for (new_row, &old_row) in indices.iter().enumerate() {
                for j in 0..n_vars {
                    let v = s.get(old_row, j);
                    if v != 0.0 {
                        let _ = new_s.insert(new_row, j, v);
                    }
                }
            }
            MatrixData::Sparse(new_s)
        }
    }
}

impl Summarizable for AnnData {
    fn summary(&self) -> String {
        format!(
            "AnnData: {} obs \u{00d7} {} vars, {} layers, {} obsm, {} varm, {} obsp, {} uns",
            self.n_obs(),
            self.n_vars(),
            self.layers.len(),
            self.obsm.len(),
            self.varm.len(),
            self.obsp.len(),
            self.uns.len(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_adata() -> AnnData {
        let x = MatrixData::Dense(vec![
            vec![1.0, 2.0, 0.0],
            vec![3.0, 0.0, 4.0],
            vec![0.0, 5.0, 6.0],
        ]);
        AnnData::new(
            x,
            vec!["cell_1".into(), "cell_2".into(), "cell_3".into()],
            vec!["gene_a".into(), "gene_b".into(), "gene_c".into()],
        )
        .unwrap()
    }

    #[test]
    fn basic_construction() {
        let adata = sample_adata();
        assert_eq!(adata.n_obs(), 3);
        assert_eq!(adata.n_vars(), 3);
        assert_eq!(adata.shape(), (3, 3));
    }

    #[test]
    fn dimension_mismatch_error() {
        let x = MatrixData::Dense(vec![vec![1.0, 2.0]]);
        let result = AnnData::new(
            x,
            vec!["cell_1".into(), "cell_2".into()], // 2 names, 1 row
            vec!["gene_a".into(), "gene_b".into()],
        );
        assert!(result.is_err());
    }

    #[test]
    fn obs_metadata() {
        let mut adata = sample_adata();
        adata
            .add_obs(
                "cell_type",
                vec!["T-cell".into(), "B-cell".into(), "NK".into()],
            )
            .unwrap();
        let ct = adata.get_obs_strings("cell_type").unwrap();
        assert_eq!(ct[0], "T-cell");
        assert!(adata.get_obs("missing").is_none());
    }

    #[test]
    fn obs_metadata_length_mismatch() {
        let mut adata = sample_adata();
        let result = adata.add_obs("bad", vec!["a".into()]);
        assert!(result.is_err());
    }

    #[test]
    fn var_metadata() {
        let mut adata = sample_adata();
        adata
            .add_var(
                "gene_type",
                vec!["coding".into(), "coding".into(), "lncRNA".into()],
            )
            .unwrap();
        let gt = adata.get_var_strings("gene_type").unwrap();
        assert_eq!(gt[2], "lncRNA");
    }

    #[test]
    fn obsm_embedding() {
        let mut adata = sample_adata();
        let pca = vec![vec![0.1, 0.2], vec![0.3, 0.4], vec![0.5, 0.6]];
        adata.add_obsm("X_pca", pca).unwrap();
        let emb = adata.get_obsm("X_pca").unwrap();
        assert_eq!(emb.len(), 3);
        assert_eq!(emb[0], vec![0.1, 0.2]);
    }

    #[test]
    fn layers() {
        let mut adata = sample_adata();
        let raw = MatrixData::Dense(vec![
            vec![10.0, 20.0, 0.0],
            vec![30.0, 0.0, 40.0],
            vec![0.0, 50.0, 60.0],
        ]);
        adata.add_layer("raw_counts", raw).unwrap();
        let layer = adata.get_layer("raw_counts").unwrap();
        assert_eq!(layer.get(0, 0), 10.0);
    }

    #[test]
    fn layer_shape_mismatch() {
        let mut adata = sample_adata();
        let bad = MatrixData::Dense(vec![vec![1.0]]);
        assert!(adata.add_layer("bad", bad).is_err());
    }

    #[test]
    fn subset_obs() {
        let mut adata = sample_adata();
        adata
            .add_obs(
                "label",
                vec!["a".into(), "b".into(), "c".into()],
            )
            .unwrap();
        let sub = adata.subset_obs(&[0, 2]).unwrap();
        assert_eq!(sub.n_obs(), 2);
        assert_eq!(sub.n_vars(), 3);
        assert_eq!(sub.obs_names(), &["cell_1", "cell_3"]);
        let labels = sub.get_obs_strings("label").unwrap();
        assert_eq!(labels, &["a", "c"]);
    }

    #[test]
    fn qc_metrics() {
        let adata = sample_adata();
        let qc = adata.qc_metrics();
        assert_eq!(qc.total_counts, vec![3.0, 7.0, 11.0]);
        assert_eq!(qc.n_features, vec![2, 2, 2]);
    }

    #[test]
    fn sparse_x() {
        let s = SparseMatrix::from_triplets(
            vec![0, 1],
            vec![0, 1],
            vec![5.0, 10.0],
            2,
            2,
        )
        .unwrap();
        let x = MatrixData::Sparse(s);
        let adata = AnnData::new(
            x,
            vec!["c1".into(), "c2".into()],
            vec!["g1".into(), "g2".into()],
        )
        .unwrap();
        assert_eq!(adata.x().get(0, 0), 5.0);
        assert_eq!(adata.x().get(0, 1), 0.0);
    }

    #[test]
    fn summary_format() {
        let adata = sample_adata();
        let s = adata.summary();
        assert!(s.contains("3 obs"));
        assert!(s.contains("3 vars"));
        assert!(s.contains("0 obsp"));
        assert!(s.contains("0 uns"));
    }

    #[test]
    fn matrix_data_set_dense() {
        let mut x = MatrixData::Dense(vec![vec![1.0, 2.0], vec![3.0, 4.0]]);
        x.set(0, 1, 99.0);
        assert_eq!(x.get(0, 1), 99.0);
        assert_eq!(x.get(0, 0), 1.0);
    }

    #[test]
    fn matrix_data_set_sparse() {
        let s = SparseMatrix::new(2, 2);
        let mut x = MatrixData::Sparse(s);
        x.set(0, 1, 5.0);
        assert_eq!(x.get(0, 1), 5.0);
    }

    #[test]
    fn matrix_data_column_sums_dense() {
        let x = MatrixData::Dense(vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
        ]);
        assert_eq!(x.column_sums(), vec![5.0, 7.0, 9.0]);
    }

    #[test]
    fn matrix_data_column_means_dense() {
        let x = MatrixData::Dense(vec![
            vec![2.0, 4.0],
            vec![6.0, 8.0],
        ]);
        let means = x.column_means();
        assert!((means[0] - 4.0).abs() < 1e-10);
        assert!((means[1] - 6.0).abs() < 1e-10);
    }

    #[test]
    fn matrix_data_row_sums_dense() {
        let x = MatrixData::Dense(vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
        ]);
        assert_eq!(x.row_sums(), vec![6.0, 15.0]);
    }

    #[test]
    fn matrix_data_to_flat_row_major_dense() {
        let x = MatrixData::Dense(vec![vec![1.0, 2.0], vec![3.0, 4.0]]);
        assert_eq!(x.to_flat_row_major(), vec![1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn matrix_data_to_flat_row_major_sparse() {
        let s = SparseMatrix::from_triplets(
            vec![0, 1],
            vec![1, 0],
            vec![2.0, 3.0],
            2,
            2,
        )
        .unwrap();
        let x = MatrixData::Sparse(s);
        assert_eq!(x.to_flat_row_major(), vec![0.0, 2.0, 3.0, 0.0]);
    }

    #[test]
    fn x_mut_modify() {
        let mut adata = sample_adata();
        adata.x_mut().set(0, 0, 42.0);
        assert_eq!(adata.x().get(0, 0), 42.0);
    }

    #[test]
    fn set_x_valid() {
        let mut adata = sample_adata();
        let new_x = MatrixData::Dense(vec![
            vec![10.0, 20.0, 30.0],
            vec![40.0, 50.0, 60.0],
            vec![70.0, 80.0, 90.0],
        ]);
        adata.set_x(new_x).unwrap();
        assert_eq!(adata.x().get(0, 0), 10.0);
    }

    #[test]
    fn set_x_shape_mismatch() {
        let mut adata = sample_adata();
        let bad = MatrixData::Dense(vec![vec![1.0]]);
        assert!(adata.set_x(bad).is_err());
    }

    #[test]
    fn subset_vars_basic() {
        let mut adata = sample_adata();
        adata.add_var("type", vec!["a".into(), "b".into(), "c".into()]).unwrap();
        let sub = adata.subset_vars(&[0, 2]).unwrap();
        assert_eq!(sub.n_vars(), 2);
        assert_eq!(sub.n_obs(), 3);
        assert_eq!(sub.var_names(), &["gene_a", "gene_c"]);
        assert_eq!(sub.x().get(0, 0), 1.0);
        assert_eq!(sub.x().get(0, 1), 0.0); // was gene_c at col 2
        let types = sub.get_var_strings("type").unwrap();
        assert_eq!(types, &["a", "c"]);
    }

    #[test]
    fn subset_vars_out_of_bounds() {
        let adata = sample_adata();
        assert!(adata.subset_vars(&[0, 10]).is_err());
    }

    #[test]
    fn obsp_add_get() {
        let mut adata = sample_adata();
        let mut m = SparseMatrix::new(3, 3);
        m.insert(0, 1, 0.5).unwrap();
        m.insert(1, 2, 0.3).unwrap();
        adata.add_obsp("connectivities", m).unwrap();
        let conn = adata.get_obsp("connectivities").unwrap();
        assert_eq!(conn.get(0, 1), 0.5);
        assert!(adata.get_obsp("missing").is_none());
    }

    #[test]
    fn obsp_wrong_shape() {
        let mut adata = sample_adata();
        let m = SparseMatrix::new(2, 2); // 3x3 needed
        assert!(adata.add_obsp("bad", m).is_err());
    }

    #[test]
    fn uns_add_get() {
        let mut adata = sample_adata();
        adata.add_uns("method", "leiden".into());
        assert_eq!(adata.get_uns("method"), Some("leiden"));
        assert_eq!(adata.get_uns("missing"), None);
    }

    #[test]
    fn get_layer_mut() {
        let mut adata = sample_adata();
        let raw = MatrixData::Dense(vec![
            vec![10.0, 20.0, 0.0],
            vec![30.0, 0.0, 40.0],
            vec![0.0, 50.0, 60.0],
        ]);
        adata.add_layer("counts", raw).unwrap();
        let layer = adata.get_layer_mut("counts").unwrap();
        layer.set(0, 0, 99.0);
        assert_eq!(adata.get_layer("counts").unwrap().get(0, 0), 99.0);
    }
}
