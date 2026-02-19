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
    use numpy::{PyArray1, PyArray2, PyArrayMethods};
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
// DBSCAN
// ---------------------------------------------------------------------------

/// DBSCAN clustering result.
#[pyclass(frozen, get_all)]
struct DbscanResult {
    labels: Vec<i32>,
    n_clusters: usize,
}

/// DBSCAN density-based clustering.
///
/// Labels of -1 indicate noise points.
#[pyfunction]
#[pyo3(signature = (data, n_features, eps=0.5, min_samples=5, metric="euclidean"))]
fn dbscan(
    data: Vec<f64>,
    n_features: usize,
    eps: f64,
    min_samples: usize,
    metric: &str,
) -> PyResult<DbscanResult> {
    let metric = parse_metric(metric)?;
    let config = cyanea_ml::DbscanConfig {
        eps,
        min_samples,
        metric,
    };
    let rows: Vec<&[f64]> = data.chunks_exact(n_features).collect();
    let r = cyanea_ml::dbscan(&rows, &config).into_pyresult()?;
    Ok(DbscanResult {
        labels: r.labels,
        n_clusters: r.n_clusters,
    })
}

// ---------------------------------------------------------------------------
// Hierarchical clustering
// ---------------------------------------------------------------------------

/// Hierarchical clustering result.
#[pyclass(frozen, get_all)]
struct HierarchicalResult {
    labels: Vec<usize>,
    merge_distances: Vec<f64>,
}

/// Agglomerative hierarchical clustering.
#[pyfunction]
#[pyo3(signature = (data, n_features, n_clusters=2, linkage="average", metric="euclidean"))]
fn hierarchical(
    data: Vec<f64>,
    n_features: usize,
    n_clusters: usize,
    linkage: &str,
    metric: &str,
) -> PyResult<HierarchicalResult> {
    let metric = parse_metric(metric)?;
    let linkage_enum = match linkage {
        "single" => cyanea_ml::Linkage::Single,
        "complete" => cyanea_ml::Linkage::Complete,
        "average" => cyanea_ml::Linkage::Average,
        "ward" => cyanea_ml::Linkage::Ward,
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unknown linkage: {linkage} (expected 'single', 'complete', 'average', or 'ward')"
            )))
        }
    };
    let rows: Vec<&[f64]> = data.chunks_exact(n_features).collect();
    let dm = cyanea_ml::distance::pairwise_distances(&rows, metric).into_pyresult()?;
    let config = cyanea_ml::HierarchicalConfig {
        n_clusters,
        linkage: linkage_enum,
    };
    let r = cyanea_ml::hierarchical(&dm, &config).into_pyresult()?;
    let merge_distances: Vec<f64> = r.merge_history.iter().map(|m| m.distance).collect();
    Ok(HierarchicalResult {
        labels: r.labels,
        merge_distances,
    })
}

// ---------------------------------------------------------------------------
// KNN
// ---------------------------------------------------------------------------

/// K-Nearest Neighbors model.
#[pyclass]
struct KnnModel {
    inner: cyanea_ml::KnnModel,
}

#[pymethods]
impl KnnModel {
    /// Fit a KNN model from row-major data.
    #[new]
    #[pyo3(signature = (data, n_features, k=5, metric="euclidean"))]
    fn new(data: Vec<f64>, n_features: usize, k: usize, metric: &str) -> PyResult<Self> {
        let metric = parse_metric(metric)?;
        let config = cyanea_ml::KnnConfig { k, metric };
        let inner = cyanea_ml::KnnModel::fit(&data, n_features, config).into_pyresult()?;
        Ok(Self { inner })
    }

    /// Find k nearest neighbors. Returns list of (index, distance).
    fn neighbors(&self, query: Vec<f64>) -> PyResult<Vec<(usize, f64)>> {
        self.inner.neighbors(&query).into_pyresult()
    }

    /// Classify a query point by majority vote.
    fn classify(&self, query: Vec<f64>, labels: Vec<i32>) -> PyResult<i32> {
        self.inner.classify(&query, &labels).into_pyresult()
    }

    /// Regress: predict target as mean of k nearest targets.
    fn regress(&self, query: Vec<f64>, targets: Vec<f64>) -> PyResult<f64> {
        self.inner.regress(&query, &targets).into_pyresult()
    }

    fn __repr__(&self) -> String {
        format!(
            "KnnModel(n_samples={}, n_features={})",
            self.inner.n_samples(),
            self.inner.n_features()
        )
    }
}

// ---------------------------------------------------------------------------
// Linear regression
// ---------------------------------------------------------------------------

/// Linear regression model.
#[pyclass]
struct LinearRegression {
    inner: cyanea_ml::LinearRegression,
}

#[pymethods]
impl LinearRegression {
    /// Fit OLS linear regression from row-major data and targets.
    #[new]
    fn new(data: Vec<f64>, n_features: usize, targets: Vec<f64>) -> PyResult<Self> {
        let inner =
            cyanea_ml::LinearRegression::fit(&data, n_features, &targets).into_pyresult()?;
        Ok(Self { inner })
    }

    /// Learned weights (one per feature).
    #[getter]
    fn weights(&self) -> Vec<f64> {
        self.inner.weights.clone()
    }

    /// Learned bias (intercept).
    #[getter]
    fn bias(&self) -> f64 {
        self.inner.bias
    }

    /// R-squared (coefficient of determination).
    #[getter]
    fn r_squared(&self) -> f64 {
        self.inner.r_squared
    }

    /// Predict a single query.
    fn predict(&self, query: Vec<f64>) -> PyResult<f64> {
        self.inner.predict(&query).into_pyresult()
    }

    /// Predict a batch of queries (row-major).
    fn predict_batch(&self, queries: Vec<f64>) -> PyResult<Vec<f64>> {
        self.inner.predict_batch(&queries).into_pyresult()
    }

    fn __repr__(&self) -> String {
        format!(
            "LinearRegression(n_features={}, r_squared={:.4})",
            self.inner.n_features, self.inner.r_squared
        )
    }
}

// ---------------------------------------------------------------------------
// Decision tree
// ---------------------------------------------------------------------------

/// Decision tree classifier.
#[pyclass]
struct DecisionTree {
    inner: cyanea_ml::DecisionTree,
}

#[pymethods]
impl DecisionTree {
    /// Fit a decision tree classifier.
    #[new]
    #[pyo3(signature = (data, n_features, labels, max_depth=10))]
    fn new(data: Vec<f64>, n_features: usize, labels: Vec<usize>, max_depth: usize) -> PyResult<Self> {
        let inner =
            cyanea_ml::DecisionTree::fit(&data, n_features, &labels, max_depth).into_pyresult()?;
        Ok(Self { inner })
    }

    /// Predict class label for a single sample.
    fn predict(&self, sample: Vec<f64>) -> usize {
        self.inner.predict(&sample)
    }

    /// Predict class labels for a batch (row-major).
    fn predict_batch(&self, data: Vec<f64>, n_features: usize) -> Vec<usize> {
        self.inner.predict_batch(&data, n_features)
    }
}

// ---------------------------------------------------------------------------
// Random forest
// ---------------------------------------------------------------------------

/// Random forest classifier.
#[pyclass]
struct RandomForest {
    inner: cyanea_ml::RandomForest,
}

#[pymethods]
impl RandomForest {
    /// Fit a random forest classifier.
    #[new]
    #[pyo3(signature = (data, n_features, labels, n_trees=10, max_depth=5, seed=42))]
    fn new(
        data: Vec<f64>,
        n_features: usize,
        labels: Vec<usize>,
        n_trees: usize,
        max_depth: usize,
        seed: u64,
    ) -> PyResult<Self> {
        let config = cyanea_ml::RandomForestConfig {
            n_trees,
            max_depth,
            max_features: None,
            seed,
        };
        let inner =
            cyanea_ml::RandomForest::fit(&data, n_features, &labels, &config).into_pyresult()?;
        Ok(Self { inner })
    }

    /// Predict class label for a single sample.
    fn predict(&self, sample: Vec<f64>) -> usize {
        self.inner.predict(&sample)
    }

    /// Predict class labels for a batch (row-major).
    fn predict_batch(&self, data: Vec<f64>, n_features: usize) -> Vec<usize> {
        self.inner.predict_batch(&data, n_features)
    }

    /// Feature importance scores.
    fn feature_importances(&self, n_features: usize) -> Vec<f64> {
        self.inner.feature_importance(n_features)
    }
}

// ---------------------------------------------------------------------------
// HMM
// ---------------------------------------------------------------------------

/// Hidden Markov Model.
#[pyclass]
struct HmmModel {
    inner: cyanea_ml::HmmModel,
}

#[pymethods]
impl HmmModel {
    /// Create an HMM with given parameters.
    ///
    /// - initial: flat array of initial probabilities (n_states)
    /// - transition: flat row-major transition matrix (n_states × n_states)
    /// - emission: flat row-major emission matrix (n_states × n_symbols)
    #[new]
    fn new(
        n_states: usize,
        n_symbols: usize,
        initial: Vec<f64>,
        transition: Vec<f64>,
        emission: Vec<f64>,
    ) -> PyResult<Self> {
        let inner =
            cyanea_ml::HmmModel::new(n_states, n_symbols, initial, transition, emission)
                .into_pyresult()?;
        Ok(Self { inner })
    }

    /// Viterbi decoding. Returns (most_likely_path, log_probability).
    fn viterbi(&self, observations: Vec<usize>) -> PyResult<(Vec<usize>, f64)> {
        self.inner.viterbi(&observations).into_pyresult()
    }

    /// Forward algorithm. Returns log-likelihood of observations.
    fn log_likelihood(&self, observations: Vec<usize>) -> PyResult<f64> {
        self.inner.log_likelihood(&observations).into_pyresult()
    }

    fn __repr__(&self) -> String {
        format!(
            "HmmModel(n_states={}, n_symbols={})",
            self.inner.n_states(),
            self.inner.n_symbols()
        )
    }
}

// ---------------------------------------------------------------------------
// Normalization
// ---------------------------------------------------------------------------

/// Min-max normalize data to [0, 1].
#[pyfunction]
fn normalize_min_max(data: Vec<f64>) -> PyResult<Vec<f64>> {
    let mut out = data;
    cyanea_ml::normalize::min_max(&mut out).into_pyresult()?;
    Ok(out)
}

/// Z-score normalize data (zero mean, unit variance).
#[pyfunction]
fn normalize_z_score(data: Vec<f64>) -> PyResult<Vec<f64>> {
    let mut out = data;
    cyanea_ml::normalize::z_score(&mut out).into_pyresult()?;
    Ok(out)
}

/// L2-normalize data to unit norm.
#[pyfunction]
fn normalize_l2(data: Vec<f64>) -> PyResult<Vec<f64>> {
    let mut out = data;
    cyanea_ml::normalize::l2_normalize(&mut out).into_pyresult()?;
    Ok(out)
}

// ---------------------------------------------------------------------------
// Cluster evaluation
// ---------------------------------------------------------------------------

/// Silhouette score (mean across all non-noise samples).
#[pyfunction]
#[pyo3(signature = (data, n_features, labels))]
fn silhouette_score(data: Vec<f64>, n_features: usize, labels: Vec<i32>) -> PyResult<f64> {
    let rows: Vec<&[f64]> = data.chunks_exact(n_features).collect();
    cyanea_ml::silhouette_score(&rows, &labels).into_pyresult()
}

/// Per-sample silhouette coefficients.
#[pyfunction]
#[pyo3(signature = (data, n_features, labels))]
fn silhouette_samples(data: Vec<f64>, n_features: usize, labels: Vec<i32>) -> PyResult<Vec<f64>> {
    let rows: Vec<&[f64]> = data.chunks_exact(n_features).collect();
    cyanea_ml::silhouette_samples(&rows, &labels).into_pyresult()
}

// ---------------------------------------------------------------------------
// Gradient Boosted Trees
// ---------------------------------------------------------------------------

/// Gradient boosted decision tree model.
#[pyclass]
struct GradientBoostedTrees {
    inner: cyanea_ml::GradientBoostedTrees,
}

#[pymethods]
impl GradientBoostedTrees {
    /// Fit a GBDT regression model.
    #[staticmethod]
    #[pyo3(signature = (data, n_features, targets, n_trees=100, max_depth=5, learning_rate=0.1, seed=42))]
    fn fit_regression(
        data: Vec<f64>,
        n_features: usize,
        targets: Vec<f64>,
        n_trees: usize,
        max_depth: usize,
        learning_rate: f64,
        seed: u64,
    ) -> PyResult<Self> {
        let config = cyanea_ml::GbdtConfig {
            n_estimators: n_trees,
            max_depth,
            learning_rate,
            seed,
            ..Default::default()
        };
        let inner = cyanea_ml::GradientBoostedTrees::fit_regression(
            &data, n_features, &targets, &config,
        )
        .into_pyresult()?;
        Ok(Self { inner })
    }

    /// Fit a GBDT classification model.
    #[staticmethod]
    #[pyo3(signature = (data, n_features, labels, n_trees=100, max_depth=5, learning_rate=0.1, seed=42))]
    fn fit_classification(
        data: Vec<f64>,
        n_features: usize,
        labels: Vec<usize>,
        n_trees: usize,
        max_depth: usize,
        learning_rate: f64,
        seed: u64,
    ) -> PyResult<Self> {
        let config = cyanea_ml::GbdtConfig {
            n_estimators: n_trees,
            max_depth,
            learning_rate,
            seed,
            ..Default::default()
        };
        let inner = cyanea_ml::GradientBoostedTrees::fit_classification(
            &data, n_features, &labels, &config,
        )
        .into_pyresult()?;
        Ok(Self { inner })
    }

    /// Predict the raw output for a single sample.
    fn predict(&self, sample: Vec<f64>) -> f64 {
        self.inner.predict(&sample)
    }

    /// Predict raw outputs for multiple samples (row-major).
    fn predict_batch(&self, data: Vec<f64>, n_features: usize) -> Vec<f64> {
        self.inner.predict_batch(&data, n_features)
    }

    /// Predict class label for a single sample.
    fn predict_class(&self, sample: Vec<f64>) -> usize {
        self.inner.predict_class(&sample)
    }

    /// Predict class labels for multiple samples (row-major).
    fn predict_class_batch(&self, data: Vec<f64>, n_features: usize) -> Vec<usize> {
        self.inner.predict_class_batch(&data, n_features)
    }

    /// Impurity-based feature importance scores.
    fn feature_importances(&self) -> Vec<f64> {
        self.inner.feature_importance()
    }
}

// ---------------------------------------------------------------------------
// Classification metrics
// ---------------------------------------------------------------------------

/// Confusion matrix result.
#[pyclass(frozen, get_all)]
struct PyConfusionMatrix {
    matrix: Vec<Vec<usize>>,
    labels: Vec<usize>,
    accuracy: f64,
}

/// Compute a confusion matrix from actual and predicted labels.
#[pyfunction]
fn confusion_matrix(actual: Vec<usize>, predicted: Vec<usize>) -> PyResult<PyConfusionMatrix> {
    let cm =
        cyanea_ml::ConfusionMatrix::from_labels(&actual, &predicted, None).into_pyresult()?;
    let nc = cm.n_classes;
    let mut matrix = Vec::with_capacity(nc);
    for i in 0..nc {
        let row: Vec<usize> = (0..nc).map(|j| cm.get(i, j)).collect();
        matrix.push(row);
    }
    let labels: Vec<usize> = (0..nc).collect();
    let accuracy = cm.accuracy();
    Ok(PyConfusionMatrix {
        matrix,
        labels,
        accuracy,
    })
}

/// Area under the ROC curve.
#[pyfunction]
fn roc_auc(scores: Vec<f64>, labels: Vec<bool>) -> PyResult<f64> {
    cyanea_ml::roc_auc(&scores, &labels).into_pyresult()
}

/// ROC curve as a list of (fpr, tpr, threshold) tuples.
#[pyfunction]
fn roc_curve(scores: Vec<f64>, labels: Vec<bool>) -> PyResult<Vec<(f64, f64, f64)>> {
    let curve = cyanea_ml::roc_curve(&scores, &labels).into_pyresult()?;
    Ok(curve
        .points
        .iter()
        .map(|p| (p.fpr, p.tpr, p.threshold))
        .collect())
}

/// Precision-recall curve as a list of (precision, recall, threshold) tuples.
#[pyfunction]
fn pr_curve(scores: Vec<f64>, labels: Vec<bool>) -> PyResult<Vec<(f64, f64, f64)>> {
    let curve = cyanea_ml::pr_curve(&scores, &labels).into_pyresult()?;
    Ok(curve
        .points
        .iter()
        .map(|p| (p.precision, p.recall, p.threshold))
        .collect())
}

// ---------------------------------------------------------------------------
// Cross-validation
// ---------------------------------------------------------------------------

/// Cross-validation result.
#[pyclass(frozen, get_all)]
struct PyCvResult {
    fold_scores: Vec<f64>,
    mean: f64,
    std_dev: f64,
}

/// K-fold cross-validation using a random forest classifier.
///
/// Trains a random forest on each training fold and evaluates accuracy
/// on the held-out fold.
#[pyfunction]
#[pyo3(signature = (data, n_features, labels, k=5, seed=42))]
fn cross_validate_kfold(
    data: Vec<f64>,
    n_features: usize,
    labels: Vec<usize>,
    k: usize,
    seed: u64,
) -> PyResult<PyCvResult> {
    let n_samples = data.len() / n_features;
    let result = cyanea_ml::cross_validate_kfold(n_samples, k, seed, |train_idx, test_idx| {
        // Build training data
        let train_data: Vec<f64> = train_idx
            .iter()
            .flat_map(|&i| data[i * n_features..(i + 1) * n_features].iter().copied())
            .collect();
        let train_labels: Vec<usize> = train_idx.iter().map(|&i| labels[i]).collect();

        let config = cyanea_ml::RandomForestConfig {
            n_trees: 10,
            max_depth: 5,
            max_features: None,
            seed,
        };
        let forest =
            cyanea_ml::RandomForest::fit(&train_data, n_features, &train_labels, &config)?;

        // Evaluate on test fold
        let correct: usize = test_idx
            .iter()
            .filter(|&&i| {
                let sample = &data[i * n_features..(i + 1) * n_features];
                forest.predict(sample) == labels[i]
            })
            .count();
        Ok(correct as f64 / test_idx.len() as f64)
    })
    .into_pyresult()?;

    let fold_scores: Vec<f64> = result.folds.iter().map(|f| f.score).collect();
    Ok(PyCvResult {
        fold_scores,
        mean: result.mean_score,
        std_dev: result.std_score,
    })
}

// ---------------------------------------------------------------------------
// Feature selection
// ---------------------------------------------------------------------------

/// Select features by variance threshold.
///
/// Returns indices of features whose variance exceeds the threshold.
#[pyfunction]
#[pyo3(signature = (data, n_features, threshold=0.0))]
fn variance_threshold(data: Vec<f64>, n_features: usize, threshold: f64) -> PyResult<Vec<usize>> {
    let result =
        cyanea_ml::variance_threshold(&data, n_features, threshold).into_pyresult()?;
    Ok(result.selected)
}

/// Mutual information between features and discrete labels.
///
/// Returns a vector of MI scores (one per feature, in nats).
#[pyfunction]
#[pyo3(signature = (data, n_features, labels, n_bins=10))]
fn mutual_information(
    data: Vec<f64>,
    n_features: usize,
    labels: Vec<usize>,
    n_bins: usize,
) -> PyResult<Vec<f64>> {
    let result =
        cyanea_ml::mutual_information(&data, n_features, &labels, n_bins).into_pyresult()?;
    Ok(result.scores)
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
    // Distance functions
    m.add_function(wrap_pyfunction!(euclidean_distance, &m)?)?;
    m.add_function(wrap_pyfunction!(manhattan_distance, &m)?)?;
    m.add_function(wrap_pyfunction!(hamming_distance, &m)?)?;
    m.add_function(wrap_pyfunction!(cosine_similarity, &m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_distances, &m)?)?;
    // Dimensionality reduction
    m.add_function(wrap_pyfunction!(umap, &m)?)?;
    m.add_function(wrap_pyfunction!(pca, &m)?)?;
    m.add_function(wrap_pyfunction!(tsne, &m)?)?;
    // Clustering
    m.add_function(wrap_pyfunction!(kmeans, &m)?)?;
    m.add_function(wrap_pyfunction!(dbscan, &m)?)?;
    m.add_function(wrap_pyfunction!(hierarchical, &m)?)?;
    // Normalization
    m.add_function(wrap_pyfunction!(normalize_min_max, &m)?)?;
    m.add_function(wrap_pyfunction!(normalize_z_score, &m)?)?;
    m.add_function(wrap_pyfunction!(normalize_l2, &m)?)?;
    // Evaluation
    m.add_function(wrap_pyfunction!(silhouette_score, &m)?)?;
    m.add_function(wrap_pyfunction!(silhouette_samples, &m)?)?;
    // Classification metrics
    m.add_function(wrap_pyfunction!(confusion_matrix, &m)?)?;
    m.add_function(wrap_pyfunction!(roc_auc, &m)?)?;
    m.add_function(wrap_pyfunction!(roc_curve, &m)?)?;
    m.add_function(wrap_pyfunction!(pr_curve, &m)?)?;
    // Cross-validation
    m.add_function(wrap_pyfunction!(cross_validate_kfold, &m)?)?;
    // Feature selection
    m.add_function(wrap_pyfunction!(variance_threshold, &m)?)?;
    m.add_function(wrap_pyfunction!(mutual_information, &m)?)?;
    // Classes
    m.add_class::<UmapResult>()?;
    m.add_class::<PcaResult>()?;
    m.add_class::<TsneResult>()?;
    m.add_class::<KMeansResult>()?;
    m.add_class::<DbscanResult>()?;
    m.add_class::<HierarchicalResult>()?;
    m.add_class::<KnnModel>()?;
    m.add_class::<LinearRegression>()?;
    m.add_class::<DecisionTree>()?;
    m.add_class::<RandomForest>()?;
    m.add_class::<HmmModel>()?;
    m.add_class::<GradientBoostedTrees>()?;
    m.add_class::<PyConfusionMatrix>()?;
    m.add_class::<PyCvResult>()?;

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
