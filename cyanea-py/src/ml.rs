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

/// UMAP dimensionality reduction result.
#[pyclass]
#[derive(Debug, Clone)]
struct UmapResult {
    #[pyo3(get)]
    embedding: Vec<f64>,
    #[pyo3(get)]
    n_samples: usize,
    #[pyo3(get)]
    n_components: usize,
    #[pyo3(get)]
    n_epochs: usize,
}

/// UMAP dimensionality reduction.
#[pyfunction]
#[pyo3(signature = (data, n_features, n_components=2, n_neighbors=15, min_dist=0.1, n_epochs=200, metric="euclidean", seed=42))]
fn umap(
    data: Vec<f64>,
    n_features: usize,
    n_components: usize,
    n_neighbors: usize,
    min_dist: f64,
    n_epochs: usize,
    metric: &str,
    seed: u64,
) -> PyResult<UmapResult> {
    let metric = match metric {
        "euclidean" => cyanea_ml::DistanceMetric::Euclidean,
        "manhattan" => cyanea_ml::DistanceMetric::Manhattan,
        "cosine" => cyanea_ml::DistanceMetric::Cosine,
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unknown metric: {metric}"
            )))
        }
    };
    let config = cyanea_ml::UmapConfig {
        n_components,
        n_neighbors,
        min_dist,
        n_epochs,
        metric,
        seed,
        ..Default::default()
    };
    let r = cyanea_ml::umap(&data, n_features, &config).into_pyresult()?;
    Ok(UmapResult {
        embedding: r.embedding,
        n_samples: r.n_samples,
        n_components: r.n_components,
        n_epochs: r.n_epochs,
    })
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
    m.add_function(wrap_pyfunction!(umap, &m)?)?;
    m.add_class::<UmapResult>()?;
    parent.add_submodule(&m)?;
    Ok(())
}
