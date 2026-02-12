//! Zarr-backed AnnData reader/writer for `.zarr` directories.
//!
//! The Zarr format is a chunked array storage format used in the scverse
//! ecosystem as an alternative to `.h5ad`. This module implements read/write
//! for AnnData containers using the Zarr v3 format via the `zarrs` crate.
//!
//! Requires the `zarr` feature flag.

use std::num::NonZeroU64;
use std::path::Path;
use std::sync::Arc;

use zarrs::array::{Array, ArrayBuilder, DataType, FillValue};
use zarrs::array::chunk_grid::ChunkGrid;
use zarrs::filesystem::FilesystemStore;
use zarrs::group::GroupBuilder;
use zarrs::storage::ReadableWritableListableStorage;

use cyanea_core::{CyaneaError, Result};

use crate::single_cell::{AnnData, ColumnData, MatrixData};
use crate::sparse::SparseMatrix;

fn zarr_err(e: impl std::fmt::Display) -> CyaneaError {
    CyaneaError::InvalidInput(format!("Zarr error: {e}"))
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Read a `.zarr` directory into an [`AnnData`] container.
pub fn read_zarr<P: AsRef<Path>>(path: P) -> Result<AnnData> {
    let path = path.as_ref();
    if !path.exists() {
        return Err(CyaneaError::Io(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("{}: not found", path.display()),
        )));
    }

    let store: ReadableWritableListableStorage =
        Arc::new(FilesystemStore::new(path).map_err(zarr_err)?);

    // Read X
    let x = read_x(&store)?;
    let (n_obs, n_vars) = x.shape();

    // Read obs/var names
    let obs_names = read_string_array(&store, "/obs/_index", n_obs)?;
    let var_names = read_string_array(&store, "/var/_index", n_vars)?;

    let mut adata = AnnData::new(x, obs_names, var_names)?;

    // Read obs columns
    read_columns(path, &store, "/obs", |name, col| {
        let _ = adata.add_obs_column(name, col);
    })?;

    // Read var columns
    read_columns(path, &store, "/var", |name, col| {
        let _ = adata.add_var_column(name, col);
    })?;

    // Read obsm
    read_embeddings(path, &store, "/obsm", |name, data| {
        let _ = adata.add_obsm(name, data);
    })?;

    // Read varm
    read_embeddings(path, &store, "/varm", |name, data| {
        let _ = adata.add_varm(name, data);
    })?;

    // Read layers
    read_layers(path, &store, "/layers", |name, layer| {
        let _ = adata.add_layer(name, layer);
    })?;

    Ok(adata)
}

/// Write an [`AnnData`] container to a `.zarr` directory.
pub fn write_zarr<P: AsRef<Path>>(adata: &AnnData, path: P) -> Result<()> {
    let path = path.as_ref();

    std::fs::create_dir_all(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;

    let store: ReadableWritableListableStorage =
        Arc::new(FilesystemStore::new(path).map_err(zarr_err)?);

    // Root group
    let root = GroupBuilder::new()
        .build(store.clone(), "/")
        .map_err(zarr_err)?;
    root.store_metadata().map_err(zarr_err)?;

    // Write X
    write_x(&store, adata.x())?;

    // Write obs
    write_group(&store, "/obs")?;
    write_string_array(&store, "/obs/_index", adata.obs_names())?;
    for (key, col) in adata.obs_columns() {
        write_column(&store, "/obs", key, col)?;
    }

    // Write var
    write_group(&store, "/var")?;
    write_string_array(&store, "/var/_index", adata.var_names())?;
    for (key, col) in adata.var_columns() {
        write_column(&store, "/var", key, col)?;
    }

    // Write obsm
    if !adata.obsm_keys().is_empty() {
        write_group(&store, "/obsm")?;
        for (key, data) in adata.obsm_keys() {
            write_embedding(&store, "/obsm", key, data)?;
        }
    }

    // Write varm
    if !adata.varm_keys().is_empty() {
        write_group(&store, "/varm")?;
        for (key, data) in adata.varm_keys() {
            write_embedding(&store, "/varm", key, data)?;
        }
    }

    // Write layers
    if !adata.layers_keys().is_empty() {
        write_group(&store, "/layers")?;
        for (key, layer) in adata.layers_keys() {
            write_layer(&store, "/layers", key, layer)?;
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Internal: groups
// ---------------------------------------------------------------------------

fn write_group(store: &ReadableWritableListableStorage, path: &str) -> Result<()> {
    let group = GroupBuilder::new()
        .build(store.clone(), path)
        .map_err(zarr_err)?;
    group.store_metadata().map_err(zarr_err)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Internal: chunk grid helper
// ---------------------------------------------------------------------------

fn make_chunk_grid(shape: &[u64]) -> ChunkGrid {
    let nz: Vec<NonZeroU64> = shape
        .iter()
        .map(|&s| NonZeroU64::new(if s == 0 { 1 } else { s }).unwrap())
        .collect();
    ChunkGrid::from(nz)
}

// ---------------------------------------------------------------------------
// Internal: X matrix
// ---------------------------------------------------------------------------

fn write_x(store: &ReadableWritableListableStorage, x: &MatrixData) -> Result<()> {
    match x {
        MatrixData::Dense(rows) => {
            let n_rows = rows.len();
            let n_cols = rows.first().map_or(0, |r| r.len());
            if n_rows == 0 || n_cols == 0 {
                write_f64_array(store, "/X", &[], &[0, 0])?;
                return Ok(());
            }
            let flat: Vec<f64> = rows.iter().flat_map(|r| r.iter().copied()).collect();
            write_f64_array(store, "/X", &flat, &[n_rows as u64, n_cols as u64])?;
        }
        MatrixData::Sparse(sm) => {
            let (data, indices, indptr) = sm.to_csr();
            let (n_rows, n_cols) = sm.shape();

            write_group(store, "/X")?;
            write_f64_array(store, "/X/data", &data, &[data.len() as u64])?;

            let indices_i64: Vec<i64> = indices.iter().map(|&v| v as i64).collect();
            write_i64_array(store, "/X/indices", &indices_i64, &[indices_i64.len() as u64])?;

            let indptr_i64: Vec<i64> = indptr.iter().map(|&v| v as i64).collect();
            write_i64_array(store, "/X/indptr", &indptr_i64, &[indptr_i64.len() as u64])?;

            write_i64_array(
                store,
                "/X/_shape",
                &[n_rows as i64, n_cols as i64],
                &[2],
            )?;
        }
    }
    Ok(())
}

fn read_x(store: &ReadableWritableListableStorage) -> Result<MatrixData> {
    // Try dense first
    if let Ok(array) = Array::open(store.clone(), "/X") {
        let shape = array.shape().to_vec();
        if shape.len() == 2 {
            let n_rows = shape[0] as usize;
            let n_cols = shape[1] as usize;
            if n_rows == 0 || n_cols == 0 {
                return Ok(MatrixData::Dense(vec![]));
            }
            let flat: Vec<f64> = array
                .retrieve_array_subset_elements(&array.subset_all())
                .map_err(zarr_err)?;
            let rows: Vec<Vec<f64>> = flat.chunks(n_cols).map(|c| c.to_vec()).collect();
            return Ok(MatrixData::Dense(rows));
        }
    }

    // Try CSR sparse
    let data_arr = Array::open(store.clone(), "/X/data").map_err(zarr_err)?;
    let data: Vec<f64> = data_arr
        .retrieve_array_subset_elements(&data_arr.subset_all())
        .map_err(zarr_err)?;

    let indices_arr = Array::open(store.clone(), "/X/indices").map_err(zarr_err)?;
    let indices_i64: Vec<i64> = indices_arr
        .retrieve_array_subset_elements(&indices_arr.subset_all())
        .map_err(zarr_err)?;
    let indices: Vec<usize> = indices_i64.iter().map(|&v| v as usize).collect();

    let indptr_arr = Array::open(store.clone(), "/X/indptr").map_err(zarr_err)?;
    let indptr_i64: Vec<i64> = indptr_arr
        .retrieve_array_subset_elements(&indptr_arr.subset_all())
        .map_err(zarr_err)?;
    let indptr: Vec<usize> = indptr_i64.iter().map(|&v| v as usize).collect();

    let shape_arr = Array::open(store.clone(), "/X/_shape").map_err(zarr_err)?;
    let shape_data: Vec<i64> = shape_arr
        .retrieve_array_subset_elements(&shape_arr.subset_all())
        .map_err(zarr_err)?;
    let n_rows = shape_data[0] as usize;
    let n_cols = shape_data[1] as usize;

    let sm = SparseMatrix::from_csr(data, indices, indptr, n_rows, n_cols)?;
    Ok(MatrixData::Sparse(sm))
}

// ---------------------------------------------------------------------------
// Internal: typed array writers
// ---------------------------------------------------------------------------

fn write_f64_array(
    store: &ReadableWritableListableStorage,
    path: &str,
    data: &[f64],
    shape: &[u64],
) -> Result<()> {
    let total: u64 = shape.iter().product();
    let array = ArrayBuilder::new(
        shape.to_vec(),
        DataType::Float64,
        make_chunk_grid(shape),
        FillValue::from(0.0f64),
    )
    .build(store.clone(), path)
    .map_err(zarr_err)?;
    array.store_metadata().map_err(zarr_err)?;
    if total > 0 {
        array
            .store_array_subset_elements::<f64>(&array.subset_all(), data)
            .map_err(zarr_err)?;
    }
    Ok(())
}

fn write_i64_array(
    store: &ReadableWritableListableStorage,
    path: &str,
    data: &[i64],
    shape: &[u64],
) -> Result<()> {
    let total: u64 = shape.iter().product();
    let array = ArrayBuilder::new(
        shape.to_vec(),
        DataType::Int64,
        make_chunk_grid(shape),
        FillValue::from(0i64),
    )
    .build(store.clone(), path)
    .map_err(zarr_err)?;
    array.store_metadata().map_err(zarr_err)?;
    if total > 0 {
        array
            .store_array_subset_elements::<i64>(&array.subset_all(), data)
            .map_err(zarr_err)?;
    }
    Ok(())
}

fn write_i32_array(
    store: &ReadableWritableListableStorage,
    path: &str,
    data: &[i32],
    shape: &[u64],
) -> Result<()> {
    let total: u64 = shape.iter().product();
    let array = ArrayBuilder::new(
        shape.to_vec(),
        DataType::Int32,
        make_chunk_grid(shape),
        FillValue::from(0i32),
    )
    .build(store.clone(), path)
    .map_err(zarr_err)?;
    array.store_metadata().map_err(zarr_err)?;
    if total > 0 {
        array
            .store_array_subset_elements::<i32>(&array.subset_all(), data)
            .map_err(zarr_err)?;
    }
    Ok(())
}

fn write_string_array(
    store: &ReadableWritableListableStorage,
    path: &str,
    values: &[String],
) -> Result<()> {
    let n = values.len();
    let chunk_size = if n == 0 { 1u64 } else { n as u64 };
    let array = ArrayBuilder::new(
        vec![n as u64],
        DataType::String,
        make_chunk_grid(&[chunk_size]),
        FillValue::from(""),
    )
    .build(store.clone(), path)
    .map_err(zarr_err)?;
    array.store_metadata().map_err(zarr_err)?;
    if n > 0 {
        array
            .store_array_subset_elements::<String>(&array.subset_all(), &values.to_vec())
            .map_err(zarr_err)?;
    }
    Ok(())
}

fn read_string_array(
    store: &ReadableWritableListableStorage,
    path: &str,
    expected: usize,
) -> Result<Vec<String>> {
    if expected == 0 {
        return Ok(vec![]);
    }
    match Array::open(store.clone(), path) {
        Ok(array) => {
            let data: Vec<String> = array
                .retrieve_array_subset_elements(&array.subset_all())
                .map_err(zarr_err)?;
            Ok(data)
        }
        Err(_) => Ok((0..expected).map(|i| format!("{i}")).collect()),
    }
}

// ---------------------------------------------------------------------------
// Internal: columns
// ---------------------------------------------------------------------------

fn write_column(
    store: &ReadableWritableListableStorage,
    group_path: &str,
    name: &str,
    col: &ColumnData,
) -> Result<()> {
    let col_path = format!("{group_path}/{name}");
    match col {
        ColumnData::Strings(vals) => {
            write_string_array(store, &col_path, vals)?;
        }
        ColumnData::Numeric(vals) => {
            write_f64_array(store, &col_path, vals, &[vals.len() as u64])?;
        }
        ColumnData::Categorical { codes, categories } => {
            write_i32_array(store, &col_path, codes, &[codes.len() as u64])?;
            let cat_group = format!("{group_path}/__categories");
            write_group(store, &cat_group)?;
            write_string_array(store, &format!("{cat_group}/{name}"), categories)?;
        }
    }
    Ok(())
}

fn read_columns(
    root: &Path,
    store: &ReadableWritableListableStorage,
    group_path: &str,
    mut add_fn: impl FnMut(&str, ColumnData),
) -> Result<()> {
    let group_dir = root.join(group_path.trim_start_matches('/'));
    let entries = match std::fs::read_dir(&group_dir) {
        Ok(e) => e,
        Err(_) => return Ok(()),
    };

    for entry in entries {
        let entry = match entry {
            Ok(e) => e,
            Err(_) => continue,
        };
        let name = entry.file_name().to_string_lossy().to_string();
        if name == "_index" || name == "__categories" || name == "zarr.json" {
            continue;
        }

        let col_path = format!("{group_path}/{name}");

        // Try categorical first
        let cat_path = format!("{group_path}/__categories/{name}");
        if let Ok(cat_array) = Array::open(store.clone(), &cat_path) {
            if let Ok(categories) =
                cat_array.retrieve_array_subset_elements::<String>(&cat_array.subset_all())
            {
                if let Ok(code_array) = Array::open(store.clone(), &col_path) {
                    if let Ok(codes) =
                        code_array.retrieve_array_subset_elements::<i32>(&code_array.subset_all())
                    {
                        add_fn(&name, ColumnData::Categorical { codes, categories });
                        continue;
                    }
                }
            }
        }

        // Try opening as an array
        if let Ok(array) = Array::open(store.clone(), &col_path) {
            // Try string
            if array.data_type() == &DataType::String {
                if let Ok(data) =
                    array.retrieve_array_subset_elements::<String>(&array.subset_all())
                {
                    add_fn(&name, ColumnData::Strings(data));
                    continue;
                }
            }

            // Try f64
            if let Ok(data) = array.retrieve_array_subset_elements::<f64>(&array.subset_all()) {
                add_fn(&name, ColumnData::Numeric(data));
                continue;
            }

            // Try i32 → f64
            if let Ok(data) = array.retrieve_array_subset_elements::<i32>(&array.subset_all()) {
                let vals: Vec<f64> = data.iter().map(|&v| v as f64).collect();
                add_fn(&name, ColumnData::Numeric(vals));
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Internal: embeddings
// ---------------------------------------------------------------------------

fn write_embedding(
    store: &ReadableWritableListableStorage,
    group_path: &str,
    name: &str,
    data: &[Vec<f64>],
) -> Result<()> {
    let n_rows = data.len();
    let n_cols = data.first().map_or(0, |r| r.len());
    let flat: Vec<f64> = data.iter().flat_map(|r| r.iter().copied()).collect();
    let emb_path = format!("{group_path}/{name}");
    write_f64_array(store, &emb_path, &flat, &[n_rows as u64, n_cols as u64])
}

fn read_embeddings(
    root: &Path,
    store: &ReadableWritableListableStorage,
    group_path: &str,
    mut add_fn: impl FnMut(&str, Vec<Vec<f64>>),
) -> Result<()> {
    let group_dir = root.join(group_path.trim_start_matches('/'));
    let entries = match std::fs::read_dir(&group_dir) {
        Ok(e) => e,
        Err(_) => return Ok(()),
    };

    for entry in entries {
        let entry = match entry {
            Ok(e) => e,
            Err(_) => continue,
        };
        let name = entry.file_name().to_string_lossy().to_string();
        if name == "zarr.json" {
            continue;
        }

        let emb_path = format!("{group_path}/{name}");
        if let Ok(array) = Array::open(store.clone(), &emb_path) {
            let shape = array.shape().to_vec();
            if shape.len() == 2 {
                let n_cols = shape[1] as usize;
                if let Ok(flat) =
                    array.retrieve_array_subset_elements::<f64>(&array.subset_all())
                {
                    let rows: Vec<Vec<f64>> =
                        flat.chunks(n_cols).map(|c| c.to_vec()).collect();
                    add_fn(&name, rows);
                }
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Internal: layers
// ---------------------------------------------------------------------------

fn write_layer(
    store: &ReadableWritableListableStorage,
    group_path: &str,
    name: &str,
    matrix: &MatrixData,
) -> Result<()> {
    let layer_path = format!("{group_path}/{name}");
    match matrix {
        MatrixData::Dense(rows) => {
            let n_rows = rows.len();
            let n_cols = rows.first().map_or(0, |r| r.len());
            let flat: Vec<f64> = rows.iter().flat_map(|r| r.iter().copied()).collect();
            write_f64_array(store, &layer_path, &flat, &[n_rows as u64, n_cols as u64])
        }
        MatrixData::Sparse(sm) => {
            let (data, indices, indptr) = sm.to_csr();
            let (n_rows, n_cols) = sm.shape();

            write_group(store, &layer_path)?;
            write_f64_array(
                store,
                &format!("{layer_path}/data"),
                &data,
                &[data.len() as u64],
            )?;
            let indices_i64: Vec<i64> = indices.iter().map(|&v| v as i64).collect();
            write_i64_array(
                store,
                &format!("{layer_path}/indices"),
                &indices_i64,
                &[indices_i64.len() as u64],
            )?;
            let indptr_i64: Vec<i64> = indptr.iter().map(|&v| v as i64).collect();
            write_i64_array(
                store,
                &format!("{layer_path}/indptr"),
                &indptr_i64,
                &[indptr_i64.len() as u64],
            )?;
            write_i64_array(
                store,
                &format!("{layer_path}/_shape"),
                &[n_rows as i64, n_cols as i64],
                &[2],
            )?;
            Ok(())
        }
    }
}

fn read_layers(
    root: &Path,
    store: &ReadableWritableListableStorage,
    group_path: &str,
    mut add_fn: impl FnMut(&str, MatrixData),
) -> Result<()> {
    let group_dir = root.join(group_path.trim_start_matches('/'));
    let entries = match std::fs::read_dir(&group_dir) {
        Ok(e) => e,
        Err(_) => return Ok(()),
    };

    for entry in entries {
        let entry = match entry {
            Ok(e) => e,
            Err(_) => continue,
        };
        let name = entry.file_name().to_string_lossy().to_string();
        if name == "zarr.json" {
            continue;
        }

        let layer_path = format!("{group_path}/{name}");

        // Try dense
        if let Ok(array) = Array::open(store.clone(), &layer_path) {
            let shape = array.shape().to_vec();
            if shape.len() == 2 {
                let n_cols = shape[1] as usize;
                if let Ok(flat) =
                    array.retrieve_array_subset_elements::<f64>(&array.subset_all())
                {
                    let rows: Vec<Vec<f64>> =
                        flat.chunks(n_cols).map(|c| c.to_vec()).collect();
                    add_fn(&name, MatrixData::Dense(rows));
                    continue;
                }
            }
        }

        // Try sparse
        if let Ok(data_arr) = Array::open(store.clone(), &format!("{layer_path}/data")) {
            if let Ok(data) =
                data_arr.retrieve_array_subset_elements::<f64>(&data_arr.subset_all())
            {
                let indices_arr =
                    Array::open(store.clone(), &format!("{layer_path}/indices"))
                        .map_err(zarr_err)?;
                let indices_i64: Vec<i64> = indices_arr
                    .retrieve_array_subset_elements(&indices_arr.subset_all())
                    .map_err(zarr_err)?;
                let indices: Vec<usize> = indices_i64.iter().map(|&v| v as usize).collect();

                let indptr_arr =
                    Array::open(store.clone(), &format!("{layer_path}/indptr"))
                        .map_err(zarr_err)?;
                let indptr_i64: Vec<i64> = indptr_arr
                    .retrieve_array_subset_elements(&indptr_arr.subset_all())
                    .map_err(zarr_err)?;
                let indptr: Vec<usize> = indptr_i64.iter().map(|&v| v as usize).collect();

                let shape_arr =
                    Array::open(store.clone(), &format!("{layer_path}/_shape"))
                        .map_err(zarr_err)?;
                let shape_data: Vec<i64> = shape_arr
                    .retrieve_array_subset_elements(&shape_arr.subset_all())
                    .map_err(zarr_err)?;
                let n_rows = shape_data[0] as usize;
                let n_cols = shape_data[1] as usize;

                if let Ok(sm) =
                    SparseMatrix::from_csr(data, indices, indptr, n_rows, n_cols)
                {
                    add_fn(&name, MatrixData::Sparse(sm));
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn zarr_dir() -> (TempDir, std::path::PathBuf) {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.zarr");
        (dir, path)
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
        let (_dir, path) = zarr_dir();
        let adata = sample_dense_adata();
        write_zarr(&adata, &path).unwrap();

        let loaded = read_zarr(&path).unwrap();
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
        let (_dir, path) = zarr_dir();
        let adata = sample_sparse_adata();
        write_zarr(&adata, &path).unwrap();

        let loaded = read_zarr(&path).unwrap();
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
        let (_dir, path) = zarr_dir();
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
        write_zarr(&adata, &path).unwrap();

        let loaded = read_zarr(&path).unwrap();
        let ct = loaded.get_obs_strings("cell_type").unwrap();
        assert_eq!(ct, &["T", "B"]);
        let gt = loaded.get_var_strings("gene_type").unwrap();
        assert_eq!(gt, &["coding", "lncRNA", "coding"]);
    }

    #[test]
    fn obs_var_numeric_columns() {
        let (_dir, path) = zarr_dir();
        let mut adata = sample_dense_adata();
        adata.add_obs_numeric("score", vec![0.5, 0.9]).unwrap();
        write_zarr(&adata, &path).unwrap();

        let loaded = read_zarr(&path).unwrap();
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
        let (_dir, path) = zarr_dir();
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
        write_zarr(&adata, &path).unwrap();

        let loaded = read_zarr(&path).unwrap();
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
        let (_dir, path) = zarr_dir();
        let mut adata = sample_dense_adata();
        let pca = vec![vec![0.1, 0.2], vec![0.3, 0.4]];
        adata.add_obsm("X_pca", pca.clone()).unwrap();
        let loadings = vec![vec![0.5, 0.6], vec![0.7, 0.8], vec![0.9, 1.0]];
        adata.add_varm("PCs", loadings.clone()).unwrap();
        write_zarr(&adata, &path).unwrap();

        let loaded = read_zarr(&path).unwrap();
        assert_eq!(loaded.get_obsm("X_pca").unwrap(), &pca);
        assert_eq!(loaded.get_varm("PCs").unwrap(), &loadings);
    }

    #[test]
    fn layers_roundtrip() {
        let (_dir, path) = zarr_dir();
        let mut adata = sample_dense_adata();
        let raw = MatrixData::Dense(vec![
            vec![10.0, 20.0, 0.0],
            vec![0.0, 30.0, 40.0],
        ]);
        adata.add_layer("raw_counts", raw).unwrap();
        write_zarr(&adata, &path).unwrap();

        let loaded = read_zarr(&path).unwrap();
        let layer = loaded.get_layer("raw_counts").unwrap();
        assert_eq!(layer.get(0, 0), 10.0);
        assert_eq!(layer.get(1, 2), 40.0);
    }

    #[test]
    fn empty_adata() {
        let (_dir, path) = zarr_dir();
        let x = MatrixData::Dense(vec![]);
        let adata = AnnData::new(x, vec![], vec![]).unwrap();
        write_zarr(&adata, &path).unwrap();

        let loaded = read_zarr(&path).unwrap();
        assert_eq!(loaded.n_obs(), 0);
        assert_eq!(loaded.n_vars(), 0);
    }

    #[test]
    fn nonexistent_path_error() {
        let result = read_zarr("/nonexistent/path/file.zarr");
        assert!(result.is_err());
    }
}
