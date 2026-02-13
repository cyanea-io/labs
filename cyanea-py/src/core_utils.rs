//! Python bindings for cyanea-core: hashing and compression utilities.

use pyo3::prelude::*;

use crate::error::IntoPyResult;

/// Compute the SHA-256 hash of in-memory data, returned as a hex string.
#[pyfunction]
fn sha256(data: &[u8]) -> String {
    cyanea_core::hash::sha256(data)
}

/// Compute the SHA-256 hash of a file (streaming), returned as a hex string.
#[pyfunction]
fn sha256_file(path: &str) -> PyResult<String> {
    cyanea_core::hash::sha256_file(path).into_pyresult()
}

/// Compress data using zstd at the given level (1–22).
#[pyfunction]
#[pyo3(signature = (data, *, level=3))]
fn zstd_compress(data: &[u8], level: i32) -> PyResult<Vec<u8>> {
    cyanea_core::compress::zstd_compress(data, level).into_pyresult()
}

/// Decompress zstd-compressed data.
#[pyfunction]
fn zstd_decompress(data: &[u8]) -> PyResult<Vec<u8>> {
    cyanea_core::compress::zstd_decompress(data).into_pyresult()
}

/// Compress data using gzip at the given level (0–9).
#[pyfunction]
#[pyo3(signature = (data, *, level=6))]
fn gzip_compress(data: &[u8], level: u32) -> PyResult<Vec<u8>> {
    cyanea_core::compress::gzip_compress(data, level).into_pyresult()
}

/// Decompress gzip-compressed data.
#[pyfunction]
fn gzip_decompress(data: &[u8]) -> PyResult<Vec<u8>> {
    cyanea_core::compress::gzip_decompress(data).into_pyresult()
}

/// Auto-detect compression format and decompress.
#[pyfunction]
fn decompress(data: &[u8]) -> PyResult<Vec<u8>> {
    cyanea_core::compress::decompress(data).into_pyresult()
}

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "core")?;
    m.add_function(wrap_pyfunction!(sha256, &m)?)?;
    m.add_function(wrap_pyfunction!(sha256_file, &m)?)?;
    m.add_function(wrap_pyfunction!(zstd_compress, &m)?)?;
    m.add_function(wrap_pyfunction!(zstd_decompress, &m)?)?;
    m.add_function(wrap_pyfunction!(gzip_compress, &m)?)?;
    m.add_function(wrap_pyfunction!(gzip_decompress, &m)?)?;
    m.add_function(wrap_pyfunction!(decompress, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
