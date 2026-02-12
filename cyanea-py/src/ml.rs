//! Python bindings for cyanea-ml: ML, clustering, dimensionality reduction.

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

/// Pairwise distance matrix from row-major data.
///
/// Returns a flat symmetric distance matrix (n×n).
#[pyfunction]
#[pyo3(signature = (data, n_features, metric="euclidean"))]
fn pairwise_distances(data: Vec<f64>, n_features: usize, metric: &str) -> PyResult<Vec<f64>> {
    let metric = parse_metric(metric)?;
    let n = data.len() / n_features;
    let rows: Vec<&[f64]> = data.chunks_exact(n_features).collect();
    let dm = cyanea_ml::distance::pairwise_distances(&rows, metric).into_pyresult()?;
    // Expand condensed to full symmetric matrix.
    let mut full = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            full[i * n + j] = dm.get(i, j);
        }
    }
    Ok(full)
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
    let metric = parse_metric(metric)?;
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

/// PCA result.
#[pyclass(frozen, get_all)]
struct PcaResult {
    components: Vec<f64>,
    explained_variance: Vec<f64>,
    explained_variance_ratio: Vec<f64>,
    transformed: Vec<f64>,
    mean: Vec<f64>,
    n_features: usize,
    n_components: usize,
}

/// PCA dimensionality reduction.
#[pyfunction]
#[pyo3(signature = (data, n_features, n_components=2))]
fn pca(data: Vec<f64>, n_features: usize, n_components: usize) -> PyResult<PcaResult> {
    let config = cyanea_ml::PcaConfig {
        n_components,
        ..Default::default()
    };
    let r = cyanea_ml::pca(&data, n_features, &config).into_pyresult()?;
    Ok(PcaResult {
        components: r.components,
        explained_variance: r.explained_variance,
        explained_variance_ratio: r.explained_variance_ratio,
        transformed: r.transformed,
        mean: r.mean,
        n_features: r.n_features,
        n_components: r.n_components,
    })
}

/// t-SNE result.
#[pyclass(frozen, get_all)]
struct TsneResult {
    embedding: Vec<f64>,
    n_samples: usize,
    n_components: usize,
    kl_divergence: f64,
}

/// t-SNE dimensionality reduction.
#[pyfunction]
#[pyo3(signature = (data, n_features, n_components=2, perplexity=30.0, learning_rate=200.0, n_iter=1000, seed=42))]
fn tsne(
    data: Vec<f64>,
    n_features: usize,
    n_components: usize,
    perplexity: f64,
    learning_rate: f64,
    n_iter: usize,
    seed: u64,
) -> PyResult<TsneResult> {
    let config = cyanea_ml::TsneConfig {
        n_components,
        perplexity,
        learning_rate,
        n_iter,
        seed,
    };
    let r = cyanea_ml::tsne(&data, n_features, &config).into_pyresult()?;
    Ok(TsneResult {
        embedding: r.embedding,
        n_samples: r.n_samples,
        n_components: r.n_components,
        kl_divergence: r.kl_divergence,
    })
}

/// K-Means clustering result.
#[pyclass(frozen, get_all)]
struct KMeansResult {
    centroids: Vec<f64>,
    labels: Vec<usize>,
    inertia: f64,
    n_iter: usize,
    n_features: usize,
}

/// K-Means clustering.
#[pyfunction]
#[pyo3(signature = (data, n_features, n_clusters=2, max_iter=300, tolerance=1e-4, seed=42))]
fn kmeans(
    data: Vec<f64>,
    n_features: usize,
    n_clusters: usize,
    max_iter: usize,
    tolerance: f64,
    seed: u64,
) -> PyResult<KMeansResult> {
    let config = cyanea_ml::KMeansConfig {
        n_clusters,
        max_iter,
        tolerance,
        seed,
    };
    let rows: Vec<&[f64]> = data.chunks_exact(n_features).collect();
    let r = cyanea_ml::kmeans(&rows, &config).into_pyresult()?;
    Ok(KMeansResult {
        centroids: r.centroids,
        labels: r.labels,
        inertia: r.inertia,
        n_iter: r.n_iter,
        n_features: r.n_features,
    })
}

// ---------------------------------------------------------------------------
// NumPy interop (feature-gated)
// ---------------------------------------------------------------------------

#[cfg(feature = "numpy")]
mod np {
    use numpy::{PyArray1, PyArray2};
    use pyo3::prelude::*;

    use crate::error::IntoPyResult;

    /// Pairwise distance matrix as a NumPy 2D array.
    #[pyfunction]
    #[pyo3(signature = (data, n_features, metric="euclidean"))]
    pub fn pairwise_distances_np<'py>(
        py: Python<'py>,
        data: Vec<f64>,
        n_features: usize,
        metric: &str,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let metric = super::parse_metric(metric)?;
        let n = data.len() / n_features;
        let rows: Vec<&[f64]> = data.chunks_exact(n_features).collect();
        let dm = cyanea_ml::distance::pairwise_distances(&rows, metric).into_pyresult()?;
        let mut full = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                full[i * n + j] = dm.get(i, j);
            }
        }
        Ok(PyArray1::from_vec(py, full).reshape([n, n]).unwrap())
    }

    /// UMAP embedding as a NumPy 2D array (n_samples × n_components).
    #[pyfunction]
    #[pyo3(signature = (data, n_features, n_components=2, n_neighbors=15, min_dist=0.1, n_epochs=200, metric="euclidean", seed=42))]
    pub fn umap_np<'py>(
        py: Python<'py>,
        data: Vec<f64>,
        n_features: usize,
        n_components: usize,
        n_neighbors: usize,
        min_dist: f64,
        n_epochs: usize,
        metric: &str,
        seed: u64,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let metric = super::parse_metric(metric)?;
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
        let ns = r.n_samples;
        let nc = r.n_components;
        Ok(PyArray1::from_vec(py, r.embedding).reshape([ns, nc]).unwrap())
    }

    /// PCA transformed data as a NumPy 2D array (n_samples × n_components).
    #[pyfunction]
    #[pyo3(signature = (data, n_features, n_components=2))]
    pub fn pca_np<'py>(
        py: Python<'py>,
        data: Vec<f64>,
        n_features: usize,
        n_components: usize,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let config = cyanea_ml::PcaConfig {
            n_components,
            ..Default::default()
        };
        let r = cyanea_ml::pca(&data, n_features, &config).into_pyresult()?;
        let ns = data.len() / n_features;
        let nc = r.n_components;
        Ok(PyArray1::from_vec(py, r.transformed).reshape([ns, nc]).unwrap())
    }

    /// t-SNE embedding as a NumPy 2D array (n_samples × n_components).
    #[pyfunction]
    #[pyo3(signature = (data, n_features, n_components=2, perplexity=30.0, learning_rate=200.0, n_iter=1000, seed=42))]
    pub fn tsne_np<'py>(
        py: Python<'py>,
        data: Vec<f64>,
        n_features: usize,
        n_components: usize,
        perplexity: f64,
        learning_rate: f64,
        n_iter: usize,
        seed: u64,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let config = cyanea_ml::TsneConfig {
            n_components,
            perplexity,
            learning_rate,
            n_iter,
            seed,
        };
        let r = cyanea_ml::tsne(&data, n_features, &config).into_pyresult()?;
        let ns = r.n_samples;
        let nc = r.n_components;
        Ok(PyArray1::from_vec(py, r.embedding).reshape([ns, nc]).unwrap())
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn parse_metric(metric: &str) -> PyResult<cyanea_ml::DistanceMetric> {
    match metric {
        "euclidean" => Ok(cyanea_ml::DistanceMetric::Euclidean),
        "manhattan" => Ok(cyanea_ml::DistanceMetric::Manhattan),
        "cosine" => Ok(cyanea_ml::DistanceMetric::Cosine),
        _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "unknown metric: {metric}"
        ))),
    }
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
    m.add_function(wrap_pyfunction!(pairwise_distances, &m)?)?;
    m.add_function(wrap_pyfunction!(umap, &m)?)?;
    m.add_function(wrap_pyfunction!(pca, &m)?)?;
    m.add_function(wrap_pyfunction!(tsne, &m)?)?;
    m.add_function(wrap_pyfunction!(kmeans, &m)?)?;
    m.add_class::<UmapResult>()?;
    m.add_class::<PcaResult>()?;
    m.add_class::<TsneResult>()?;
    m.add_class::<KMeansResult>()?;

    #[cfg(feature = "numpy")]
    {
        m.add_function(wrap_pyfunction!(np::pairwise_distances_np, &m)?)?;
        m.add_function(wrap_pyfunction!(np::umap_np, &m)?)?;
        m.add_function(wrap_pyfunction!(np::pca_np, &m)?)?;
        m.add_function(wrap_pyfunction!(np::tsne_np, &m)?)?;
    }

    parent.add_submodule(&m)?;
    Ok(())
}
