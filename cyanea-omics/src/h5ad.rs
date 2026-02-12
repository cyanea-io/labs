//! HDF5-backed AnnData reader/writer for `.h5ad` files.
//!
//! The `.h5ad` format is the standard persistence format for the scverse
//! ecosystem (scanpy, scvi-tools). It stores single-cell data in HDF5 with
//! a defined layout: `X` (dense or CSR sparse matrix), `obs`/`var` (columnar
//! metadata), `obsm`/`varm` (embedding matrices), and `layers` (alternative
//! matrices).
//!
//! Requires the `h5ad` feature flag and a system HDF5 library installation
//! (`brew install hdf5` on macOS, `apt install libhdf5-dev` on Linux).

use std::path::Path;
use std::str::FromStr;

use hdf5::types::VarLenUnicode;
use hdf5::File;

use cyanea_core::{CyaneaError, Result};

use crate::single_cell::{AnnData, ColumnData, MatrixData};
use crate::sparse::SparseMatrix;

fn h5err(e: hdf5::Error) -> CyaneaError {
    CyaneaError::InvalidInput(format!("HDF5 error: {e}"))
}

/// Write a VarLenUnicode scalar attribute on any HDF5 location.
/// Uses a macro because Group, Dataset, etc. all expose `new_attr` via Deref
/// chains but don't share a single trait that function generics can bind on.
macro_rules! write_attr_str {
    ($loc:expr, $key:expr, $val:expr) => {{
        let s = VarLenUnicode::from_str($val).expect("valid UTF-8");
        $loc.new_attr::<VarLenUnicode>()
            .create($key)
            .and_then(|attr| attr.write_scalar(&s))
            .map_err(h5err)
    }};
}

/// Read an `.h5ad` file into an [`AnnData`] container.
pub fn read_h5ad<P: AsRef<Path>>(path: P) -> Result<AnnData> {
    let file = File::open(path.as_ref()).map_err(|e| {
        CyaneaError::InvalidInput(format!("cannot open h5ad file: {e}"))
    })?;

    // Read X
    let x = read_x(&file)?;
    let (n_obs, n_vars) = x.shape();

    // Read obs_names and var_names
    let obs_names = read_index(&file, "obs", n_obs)?;
    let var_names = read_index(&file, "var", n_vars)?;

    let mut adata = AnnData::new(x, obs_names, var_names)?;

    // Read obs columns
    if file.link_exists("obs") {
        let obs_group = file.group("obs").map_err(h5err)?;
        for name in obs_group.member_names().map_err(h5err)? {
            if name == "_index" || name == "__categories" {
                continue;
            }
            if let Ok(col) = read_column(&obs_group, &name) {
                adata.add_obs_column(&name, col)?;
            }
        }
    }

    // Read var columns
    if file.link_exists("var") {
        let var_group = file.group("var").map_err(h5err)?;
        for name in var_group.member_names().map_err(h5err)? {
            if name == "_index" || name == "__categories" {
                continue;
            }
            if let Ok(col) = read_column(&var_group, &name) {
                var_add_column_ignore_error(&mut adata, &name, col);
            }
        }
    }

    // Read obsm
    if file.link_exists("obsm") {
        let group = file.group("obsm").map_err(h5err)?;
        for name in group.member_names().map_err(h5err)? {
            if let Ok(data) = read_embedding(&group, &name) {
                adata.add_obsm(&name, data)?;
            }
        }
    }

    // Read varm
    if file.link_exists("varm") {
        let group = file.group("varm").map_err(h5err)?;
        for name in group.member_names().map_err(h5err)? {
            if let Ok(data) = read_embedding(&group, &name) {
                adata.add_varm(&name, data)?;
            }
        }
    }

    // Read layers
    if file.link_exists("layers") {
        let group = file.group("layers").map_err(h5err)?;
        for name in group.member_names().map_err(h5err)? {
            if let Ok(layer) = read_matrix_member(&group, &name) {
                adata.add_layer(&name, layer)?;
            }
        }
    }

    Ok(adata)
}

/// Write an [`AnnData`] container to an `.h5ad` file.
pub fn write_h5ad<P: AsRef<Path>>(adata: &AnnData, path: P) -> Result<()> {
    let file = File::create(path.as_ref()).map_err(|e| {
        CyaneaError::InvalidInput(format!("cannot create h5ad file: {e}"))
    })?;

    // Write X
    write_x(&file, adata.x())?;

    // Write obs
    {
        let obs = file.create_group("obs").map_err(h5err)?;
        write_index(&obs, adata.obs_names())?;
        for (key, col) in adata.obs_columns() {
            write_column(&obs, key, col)?;
        }
    }

    // Write var
    {
        let var = file.create_group("var").map_err(h5err)?;
        write_index(&var, adata.var_names())?;
        for (key, col) in adata.var_columns() {
            write_column(&var, key, col)?;
        }
    }

    // Write obsm
    if !adata.obsm_keys().is_empty() {
        let group = file.create_group("obsm").map_err(h5err)?;
        for (key, data) in adata.obsm_keys() {
            write_embedding(&group, key, data)?;
        }
    }

    // Write varm
    if !adata.varm_keys().is_empty() {
        let group = file.create_group("varm").map_err(h5err)?;
        for (key, data) in adata.varm_keys() {
            write_embedding(&group, key, data)?;
        }
    }

    // Write layers
    if !adata.layers_keys().is_empty() {
        let group = file.create_group("layers").map_err(h5err)?;
        for (key, layer) in adata.layers_keys() {
            write_matrix_member(&group, key, layer)?;
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Read encoding-type attribute from a group (if present).
fn read_encoding_type(loc: &hdf5::Group) -> Option<String> {
    loc.attr("encoding-type")
        .ok()
        .and_then(|a| a.read_scalar::<VarLenUnicode>().ok())
        .map(|s| s.as_str().to_string())
}

/// Read X from file — supports dense arrays and CSR sparse matrices.
fn read_x(file: &File) -> Result<MatrixData> {
    // Check if X is a group (sparse) or dataset (dense)
    if let Ok(ds) = file.dataset("X") {
        // Dense array
        let arr = ds.read_2d::<f64>().map_err(h5err)?;
        let rows: Vec<Vec<f64>> = arr.rows().into_iter().map(|r| r.to_vec()).collect();
        return Ok(MatrixData::Dense(rows));
    }

    if let Ok(group) = file.group("X") {
        // Sparse matrix — read CSR components
        let encoding = read_encoding_type(&group);
        match encoding.as_deref() {
            Some("csr_matrix") | None => read_csr_group(&group),
            Some(other) => Err(CyaneaError::InvalidInput(format!(
                "unsupported X encoding-type: {other}"
            ))),
        }
    } else {
        Err(CyaneaError::InvalidInput("no X dataset or group found".into()))
    }
}

/// Read a CSR group into a SparseMatrix.
fn read_csr_group(group: &hdf5::Group) -> Result<MatrixData> {
    let data = group
        .dataset("data")
        .and_then(|ds| ds.read_1d::<f64>())
        .map_err(h5err)?
        .to_vec();
    let indices = group
        .dataset("indices")
        .and_then(|ds| ds.read_1d::<i64>())
        .map_err(h5err)?
        .iter()
        .map(|&v| v as usize)
        .collect::<Vec<_>>();
    let indptr = group
        .dataset("indptr")
        .and_then(|ds| ds.read_1d::<i64>())
        .map_err(h5err)?
        .iter()
        .map(|&v| v as usize)
        .collect::<Vec<_>>();

    // Determine shape from the shape attribute or indptr/max(indices)
    let n_rows = indptr.len().saturating_sub(1);
    let n_cols = if group.attr("shape").is_ok() {
        let shape_attr = group
            .attr("shape")
            .and_then(|a| a.read_1d::<i64>())
            .map_err(h5err)?;
        shape_attr[1] as usize
    } else {
        indices.iter().copied().max().map_or(0, |m| m + 1)
    };

    let sm = SparseMatrix::from_csr(data, indices, indptr, n_rows, n_cols)?;
    Ok(MatrixData::Sparse(sm))
}

/// Write X to file.
fn write_x(file: &File, x: &MatrixData) -> Result<()> {
    match x {
        MatrixData::Dense(rows) => {
            let n_rows = rows.len();
            let n_cols = rows.first().map_or(0, |r| r.len());
            let flat: Vec<f64> = rows.iter().flat_map(|r| r.iter().copied()).collect();
            let arr = ndarray::Array2::from_shape_vec((n_rows, n_cols), flat)
                .map_err(|e| CyaneaError::InvalidInput(format!("shape error: {e}")))?;

            let ds = file
                .new_dataset_builder()
                .with_data(&arr)
                .create("X")
                .map_err(h5err)?;

            write_attr_str!(ds, "encoding-type", "array")?;
            write_attr_str!(ds, "encoding-version", "0.2.0")?;
        }
        MatrixData::Sparse(sm) => {
            let (data, indices, indptr) = sm.to_csr();
            let (n_rows, n_cols) = sm.shape();

            let group = file.create_group("X").map_err(h5err)?;
            write_attr_str!(group, "encoding-type", "csr_matrix")?;
            write_attr_str!(group, "encoding-version", "0.1.0")?;

            // Write shape attribute
            let shape_arr = ndarray::arr1(&[n_rows as i64, n_cols as i64]);
            group
                .new_attr_builder()
                .with_data(&shape_arr)
                .create("shape")
                .map_err(h5err)?;

            // Write data, indices, indptr
            let data_arr = ndarray::Array1::from(data);
            group
                .new_dataset_builder()
                .with_data(&data_arr)
                .create("data")
                .map_err(h5err)?;

            let indices_arr =
                ndarray::Array1::from(indices.into_iter().map(|v| v as i64).collect::<Vec<_>>());
            group
                .new_dataset_builder()
                .with_data(&indices_arr)
                .create("indices")
                .map_err(h5err)?;

            let indptr_arr =
                ndarray::Array1::from(indptr.into_iter().map(|v| v as i64).collect::<Vec<_>>());
            group
                .new_dataset_builder()
                .with_data(&indptr_arr)
                .create("indptr")
                .map_err(h5err)?;
        }
    }
    Ok(())
}

/// Read the `_index` dataset from an obs/var group to get names.
fn read_index(file: &File, group_name: &str, expected: usize) -> Result<Vec<String>> {
    if !file.link_exists(group_name) {
        // No group — generate default names
        return Ok((0..expected).map(|i| format!("{i}")).collect());
    }
    let group = file.group(group_name).map_err(h5err)?;
    if !group.link_exists("_index") {
        return Ok((0..expected).map(|i| format!("{i}")).collect());
    }
    let ds = group.dataset("_index").map_err(h5err)?;
    let names = ds
        .read_1d::<VarLenUnicode>()
        .map_err(h5err)?
        .iter()
        .map(|s| s.as_str().to_string())
        .collect();
    Ok(names)
}

/// Write index names to a group.
fn write_index(group: &hdf5::Group, names: &[String]) -> Result<()> {
    let unicode: Vec<VarLenUnicode> = names
        .iter()
        .map(|s| VarLenUnicode::from_str(s).expect("valid UTF-8"))
        .collect();
    let arr = ndarray::Array1::from(unicode);
    group
        .new_dataset_builder()
        .with_data(&arr)
        .create("_index")
        .map_err(h5err)?;
    Ok(())
}

/// Read a metadata column from a group dataset.
fn read_column(group: &hdf5::Group, name: &str) -> std::result::Result<ColumnData, hdf5::Error> {
    let ds = group.dataset(name)?;

    // Check if this is a categorical — look for __categories/<name>
    if group.link_exists("__categories") {
        if let Ok(cat_group) = group.group("__categories") {
            if cat_group.link_exists(name) {
                let cats_ds = cat_group.dataset(name)?;
                let categories: Vec<String> = cats_ds
                    .read_1d::<VarLenUnicode>()?
                    .iter()
                    .map(|s| s.as_str().to_string())
                    .collect();
                let codes: Vec<i32> = ds.read_1d::<i32>()?.to_vec();
                return Ok(ColumnData::Categorical { codes, categories });
            }
        }
    }

    // Try reading as string
    if let Ok(arr) = ds.read_1d::<VarLenUnicode>() {
        let strings: Vec<String> = arr.iter().map(|s| s.as_str().to_string()).collect();
        return Ok(ColumnData::Strings(strings));
    }

    // Try reading as f64
    if let Ok(arr) = ds.read_1d::<f64>() {
        return Ok(ColumnData::Numeric(arr.to_vec()));
    }

    // Try reading as i32 and converting to f64
    if let Ok(arr) = ds.read_1d::<i32>() {
        let vals: Vec<f64> = arr.iter().map(|&v| v as f64).collect();
        return Ok(ColumnData::Numeric(vals));
    }

    Err(hdf5::Error::Internal(format!(
        "unsupported column type for '{name}'"
    )))
}

/// Write a metadata column to a group.
fn write_column(group: &hdf5::Group, name: &str, col: &ColumnData) -> Result<()> {
    match col {
        ColumnData::Strings(vals) => {
            let unicode: Vec<VarLenUnicode> = vals
                .iter()
                .map(|s| VarLenUnicode::from_str(s).expect("valid UTF-8"))
                .collect();
            let arr = ndarray::Array1::from(unicode);
            group
                .new_dataset_builder()
                .with_data(&arr)
                .create(name)
                .map_err(h5err)?;
        }
        ColumnData::Numeric(vals) => {
            let arr = ndarray::Array1::from(vals.clone());
            group
                .new_dataset_builder()
                .with_data(&arr)
                .create(name)
                .map_err(h5err)?;
        }
        ColumnData::Categorical { codes, categories } => {
            // Write codes as the dataset
            let codes_arr = ndarray::Array1::from(codes.clone());
            group
                .new_dataset_builder()
                .with_data(&codes_arr)
                .create(name)
                .map_err(h5err)?;

            // Write categories under __categories/<name>
            if !group.link_exists("__categories") {
                group.create_group("__categories").map_err(h5err)?;
            }
            let cat_group = group.group("__categories").map_err(h5err)?;
            let cat_unicode: Vec<VarLenUnicode> = categories
                .iter()
                .map(|s| VarLenUnicode::from_str(s).expect("valid UTF-8"))
                .collect();
            let cat_arr = ndarray::Array1::from(cat_unicode);
            cat_group
                .new_dataset_builder()
                .with_data(&cat_arr)
                .create(name)
                .map_err(h5err)?;
        }
    }
    Ok(())
}

/// Read a 2D embedding from a group.
fn read_embedding(
    group: &hdf5::Group,
    name: &str,
) -> std::result::Result<Vec<Vec<f64>>, hdf5::Error> {
    let ds = group.dataset(name)?;
    let arr = ds.read_2d::<f64>()?;
    Ok(arr.rows().into_iter().map(|r| r.to_vec()).collect())
}

/// Write a 2D embedding to a group.
fn write_embedding(group: &hdf5::Group, name: &str, data: &[Vec<f64>]) -> Result<()> {
    let n_rows = data.len();
    let n_cols = data.first().map_or(0, |r| r.len());
    let flat: Vec<f64> = data.iter().flat_map(|r| r.iter().copied()).collect();
    let arr = ndarray::Array2::from_shape_vec((n_rows, n_cols), flat)
        .map_err(|e| CyaneaError::InvalidInput(format!("embedding shape error: {e}")))?;
    group
        .new_dataset_builder()
        .with_data(&arr)
        .create(name)
        .map_err(h5err)?;
    Ok(())
}

/// Read a matrix member (dense or sparse) from a group.
fn read_matrix_member(
    group: &hdf5::Group,
    name: &str,
) -> std::result::Result<MatrixData, CyaneaError> {
    // Try as dataset first (dense)
    if let Ok(ds) = group.dataset(name) {
        let arr = ds.read_2d::<f64>().map_err(h5err)?;
        let rows: Vec<Vec<f64>> = arr.rows().into_iter().map(|r| r.to_vec()).collect();
        return Ok(MatrixData::Dense(rows));
    }

    // Try as group (sparse)
    if let Ok(sub_group) = group.group(name) {
        return read_csr_group(&sub_group);
    }

    Err(CyaneaError::InvalidInput(format!(
        "cannot read layer '{name}' as dense or sparse"
    )))
}

/// Write a matrix member to a group.
fn write_matrix_member(group: &hdf5::Group, name: &str, matrix: &MatrixData) -> Result<()> {
    match matrix {
        MatrixData::Dense(rows) => {
            let n_rows = rows.len();
            let n_cols = rows.first().map_or(0, |r| r.len());
            let flat: Vec<f64> = rows.iter().flat_map(|r| r.iter().copied()).collect();
            let arr = ndarray::Array2::from_shape_vec((n_rows, n_cols), flat)
                .map_err(|e| CyaneaError::InvalidInput(format!("shape error: {e}")))?;
            let ds = group
                .new_dataset_builder()
                .with_data(&arr)
                .create(name)
                .map_err(h5err)?;
            write_attr_str!(ds, "encoding-type", "array")?;
            write_attr_str!(ds, "encoding-version", "0.2.0")?;
        }
        MatrixData::Sparse(sm) => {
            let (data, indices, indptr) = sm.to_csr();
            let (n_rows, n_cols) = sm.shape();

            let sub = group.create_group(name).map_err(h5err)?;
            write_attr_str!(sub, "encoding-type", "csr_matrix")?;
            write_attr_str!(sub, "encoding-version", "0.1.0")?;

            let shape_arr = ndarray::arr1(&[n_rows as i64, n_cols as i64]);
            sub.new_attr_builder()
                .with_data(&shape_arr)
                .create("shape")
                .map_err(h5err)?;

            let data_arr = ndarray::Array1::from(data);
            sub.new_dataset_builder()
                .with_data(&data_arr)
                .create("data")
                .map_err(h5err)?;

            let indices_arr =
                ndarray::Array1::from(indices.into_iter().map(|v| v as i64).collect::<Vec<_>>());
            sub.new_dataset_builder()
                .with_data(&indices_arr)
                .create("indices")
                .map_err(h5err)?;

            let indptr_arr =
                ndarray::Array1::from(indptr.into_iter().map(|v| v as i64).collect::<Vec<_>>());
            sub.new_dataset_builder()
                .with_data(&indptr_arr)
                .create("indptr")
                .map_err(h5err)?;
        }
    }
    Ok(())
}

/// Helper: add var column, ignoring errors (used during read when column
/// length might not match due to subgroups like __categories).
fn var_add_column_ignore_error(adata: &mut AnnData, name: &str, col: ColumnData) {
    let _ = adata.add_var_column(name, col);
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    fn temp_path() -> (NamedTempFile, std::path::PathBuf) {
        let f = NamedTempFile::new().unwrap();
        let p = f.path().with_extension("h5ad");
        (f, p)
    }

    fn sample_dense_adata() -> AnnData {
        let x = MatrixData::Dense(vec![
            vec![1.0, 2.0, 0.0],
            vec![0.0, 3.0, 4.0],
        ]);
        AnnData::new(
            x,
            vec!["cell_1".into(), "cell_2".into()],
            vec!["gene_a".into(), "gene_b".into(), "gene_c".into()],
        )
        .unwrap()
    }

    fn sample_sparse_adata() -> AnnData {
        let sm = SparseMatrix::from_triplets(
            vec![0, 0, 1, 1],
            vec![0, 2, 1, 2],
            vec![1.0, 2.0, 3.0, 4.0],
            2,
            3,
        )
        .unwrap();
        AnnData::new(
            MatrixData::Sparse(sm),
            vec!["c1".into(), "c2".into()],
            vec!["g1".into(), "g2".into(), "g3".into()],
        )
        .unwrap()
    }

    #[test]
    fn write_read_roundtrip_dense() {
        let (_tmp, path) = temp_path();
        let adata = sample_dense_adata();
        write_h5ad(&adata, &path).unwrap();

        let loaded = read_h5ad(&path).unwrap();
        assert_eq!(loaded.n_obs(), 2);
        assert_eq!(loaded.n_vars(), 3);
        assert_eq!(loaded.obs_names(), &["cell_1", "cell_2"]);
        assert_eq!(loaded.var_names(), &["gene_a", "gene_b", "gene_c"]);
        assert_eq!(loaded.x().get(0, 0), 1.0);
        assert_eq!(loaded.x().get(0, 2), 0.0);
        assert_eq!(loaded.x().get(1, 1), 3.0);
        assert_eq!(loaded.x().get(1, 2), 4.0);
    }

    #[test]
    fn write_read_roundtrip_sparse() {
        let (_tmp, path) = temp_path();
        let adata = sample_sparse_adata();
        write_h5ad(&adata, &path).unwrap();

        let loaded = read_h5ad(&path).unwrap();
        assert_eq!(loaded.n_obs(), 2);
        assert_eq!(loaded.n_vars(), 3);
        assert_eq!(loaded.x().get(0, 0), 1.0);
        assert_eq!(loaded.x().get(0, 1), 0.0);
        assert_eq!(loaded.x().get(0, 2), 2.0);
        assert_eq!(loaded.x().get(1, 1), 3.0);
        assert_eq!(loaded.x().get(1, 2), 4.0);
    }

    #[test]
    fn obs_var_string_columns() {
        let (_tmp, path) = temp_path();
        let mut adata = sample_dense_adata();
        adata
            .add_obs("cell_type", vec!["T".into(), "B".into()])
            .unwrap();
        adata
            .add_var(
                "gene_type",
                vec!["coding".into(), "lncRNA".into(), "coding".into()],
            )
            .unwrap();
        write_h5ad(&adata, &path).unwrap();

        let loaded = read_h5ad(&path).unwrap();
        let ct = loaded.get_obs_strings("cell_type").unwrap();
        assert_eq!(ct, &["T", "B"]);
        let gt = loaded.get_var_strings("gene_type").unwrap();
        assert_eq!(gt, &["coding", "lncRNA", "coding"]);
    }

    #[test]
    fn obs_var_numeric_columns() {
        let (_tmp, path) = temp_path();
        let mut adata = sample_dense_adata();
        adata
            .add_obs_numeric("score", vec![0.5, 0.9])
            .unwrap();
        write_h5ad(&adata, &path).unwrap();

        let loaded = read_h5ad(&path).unwrap();
        let col = loaded.get_obs("score").unwrap();
        match col {
            ColumnData::Numeric(v) => {
                assert_eq!(v.len(), 2);
                assert!((v[0] - 0.5).abs() < 1e-10);
                assert!((v[1] - 0.9).abs() < 1e-10);
            }
            _ => panic!("expected Numeric column"),
        }
    }

    #[test]
    fn obs_var_categorical_columns() {
        let (_tmp, path) = temp_path();
        let mut adata = sample_dense_adata();
        adata
            .add_obs_column(
                "cluster",
                ColumnData::Categorical {
                    codes: vec![0, 1],
                    categories: vec!["A".into(), "B".into(), "C".into()],
                },
            )
            .unwrap();
        write_h5ad(&adata, &path).unwrap();

        let loaded = read_h5ad(&path).unwrap();
        let col = loaded.get_obs("cluster").unwrap();
        match col {
            ColumnData::Categorical { codes, categories } => {
                assert_eq!(codes, &[0, 1]);
                assert_eq!(categories, &["A", "B", "C"]);
            }
            _ => panic!("expected Categorical column"),
        }
    }

    #[test]
    fn obsm_varm_roundtrip() {
        let (_tmp, path) = temp_path();
        let mut adata = sample_dense_adata();
        let pca = vec![vec![0.1, 0.2], vec![0.3, 0.4]];
        adata.add_obsm("X_pca", pca.clone()).unwrap();
        let loadings = vec![vec![0.5, 0.6], vec![0.7, 0.8], vec![0.9, 1.0]];
        adata.add_varm("PCs", loadings.clone()).unwrap();
        write_h5ad(&adata, &path).unwrap();

        let loaded = read_h5ad(&path).unwrap();
        assert_eq!(loaded.get_obsm("X_pca").unwrap(), &pca);
        assert_eq!(loaded.get_varm("PCs").unwrap(), &loadings);
    }

    #[test]
    fn layers_roundtrip() {
        let (_tmp, path) = temp_path();
        let mut adata = sample_dense_adata();
        let raw = MatrixData::Dense(vec![
            vec![10.0, 20.0, 0.0],
            vec![0.0, 30.0, 40.0],
        ]);
        adata.add_layer("raw_counts", raw).unwrap();
        write_h5ad(&adata, &path).unwrap();

        let loaded = read_h5ad(&path).unwrap();
        let layer = loaded.get_layer("raw_counts").unwrap();
        assert_eq!(layer.get(0, 0), 10.0);
        assert_eq!(layer.get(1, 2), 40.0);
    }

    #[test]
    fn empty_adata() {
        let (_tmp, path) = temp_path();
        let x = MatrixData::Dense(vec![]);
        let adata = AnnData::new(x, vec![], vec![]).unwrap();
        write_h5ad(&adata, &path).unwrap();

        let loaded = read_h5ad(&path).unwrap();
        assert_eq!(loaded.n_obs(), 0);
        assert_eq!(loaded.n_vars(), 0);
    }

    #[test]
    fn nonexistent_file_error() {
        let result = read_h5ad("/nonexistent/path/file.h5ad");
        assert!(result.is_err());
    }
}
