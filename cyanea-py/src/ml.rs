//! Python bindings for cyanea-ml: distance metrics.

use pyo3::prelude::*;

use crate::error::IntoPyResult;

/// Euclidean distance between two vectors.
#[pyfunction]
fn euclidean_distance(a: Vec<f64>, b: Vec<f64>) -> PyResult<f64> {
    cyanea_ml::distance::euclidean(&a, &b).into_pyresult()
}

/// Manhattan distance between two vectors.
#[pyfunction]
fn manhattan_distance(a: Vec<f64>, b: Vec<f64>) -> PyResult<f64> {
    cyanea_ml::distance::manhattan(&a, &b).into_pyresult()
}

/// Hamming distance between two byte sequences.
#[pyfunction]
fn hamming_distance(a: &[u8], b: &[u8]) -> PyResult<usize> {
    cyanea_ml::distance::hamming(a, b).into_pyresult()
}

/// Cosine similarity between two vectors.
#[pyfunction]
fn cosine_similarity(a: Vec<f64>, b: Vec<f64>) -> PyResult<f64> {
    cyanea_ml::distance::cosine_similarity(&a, &b).into_pyresult()
}

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "ml")?;
    m.add_function(wrap_pyfunction!(euclidean_distance, &m)?)?;
    m.add_function(wrap_pyfunction!(manhattan_distance, &m)?)?;
    m.add_function(wrap_pyfunction!(hamming_distance, &m)?)?;
    m.add_function(wrap_pyfunction!(cosine_similarity, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
