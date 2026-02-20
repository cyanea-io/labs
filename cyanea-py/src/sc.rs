//! Python bindings for cyanea-omics single-cell analysis pipeline.

use std::collections::HashMap;

use pyo3::prelude::*;

use crate::error::IntoPyResult;

use cyanea_omics::sc_cluster::{ClusterConfig, DistanceMetric, NeighborsConfig};
use cyanea_omics::sc_integrate::{CombatConfig, HarmonyConfig, MnnConfig};
use cyanea_omics::sc_markers::{MarkerConfig, MarkerMethod};
use cyanea_omics::sc_preprocess::{HvgConfig, HvgMethod, NormalizeConfig};
use cyanea_omics::sc_trajectory::{DiffusionConfig, DptConfig};
use cyanea_omics::single_cell::{AnnData, MatrixData};
use cyanea_omics::sparse::SparseMatrix;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn build_adata(data: Vec<f64>, n_features: usize) -> PyResult<AnnData> {
    if data.is_empty() || n_features == 0 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "data and n_features must be non-empty",
        ));
    }
    if data.len() % n_features != 0 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "data length ({}) is not divisible by n_features ({})",
            data.len(),
            n_features
        )));
    }
    let n_obs = data.len() / n_features;
    let rows: Vec<Vec<f64>> = data.chunks(n_features).map(|c| c.to_vec()).collect();
    let obs_names: Vec<String> = (0..n_obs).map(|i| format!("cell_{i}")).collect();
    let var_names: Vec<String> = (0..n_features).map(|i| format!("gene_{i}")).collect();
    AnnData::new(MatrixData::Dense(rows), obs_names, var_names).into_pyresult()
}

// ---------------------------------------------------------------------------
// Result classes
// ---------------------------------------------------------------------------

/// Highly variable genes result.
#[pyclass(frozen, get_all)]
pub struct PyHvgResult {
    pub gene_indices: Vec<usize>,
    pub dispersions: Vec<f64>,
}

/// kNN graph result.
#[pyclass(frozen, get_all)]
pub struct PyNeighborsResult {
    pub distances: Vec<(usize, usize, f64)>,
    pub connectivities: Vec<(usize, usize, f64)>,
}

/// Clustering result.
#[pyclass(frozen, get_all)]
pub struct PyClusterResult {
    pub labels: Vec<usize>,
    pub modularity: f64,
}

/// Diffusion map result.
#[pyclass(frozen, get_all)]
pub struct PyDiffusionResult {
    pub components: Vec<Vec<f64>>,
    pub eigenvalues: Vec<f64>,
}

/// Diffusion pseudotime result.
#[pyclass(frozen, get_all)]
pub struct PyDptResult {
    pub pseudotime: Vec<f64>,
}

/// PAGA graph abstraction result.
#[pyclass(frozen, get_all)]
pub struct PyPagaResult {
    pub connectivities: Vec<Vec<f64>>,
    pub groups: Vec<String>,
    pub cluster_sizes: Vec<usize>,
}

/// A marker gene identified by differential expression.
#[pyclass(frozen, get_all)]
#[derive(Clone)]
pub struct PyMarkerGene {
    pub gene_index: usize,
    pub gene_name: String,
    pub score: f64,
    pub pvalue: f64,
    pub padj: f64,
    pub log2fc: f64,
    pub pct_in: f64,
    pub pct_out: f64,
}

/// Differential expression results grouped by cluster.
#[pyclass(frozen, get_all)]
pub struct PyMarkerResult {
    pub markers: HashMap<String, Vec<PyMarkerGene>>,
    pub n_clusters: usize,
}

/// Batch correction result.
#[pyclass(frozen, get_all)]
pub struct PyCorrectedResult {
    pub corrected: Vec<f64>,
    pub n_obs: usize,
    pub n_vars: usize,
}

// ---------------------------------------------------------------------------
// Preprocessing
// ---------------------------------------------------------------------------

/// Normalize each cell to a target sum and optionally log-transform.
///
/// Returns the corrected expression matrix as a flat row-major list.
#[pyfunction]
#[pyo3(signature = (data, n_features, target_sum=10000.0, log1p=true))]
fn normalize_total(
    data: Vec<f64>,
    n_features: usize,
    target_sum: f64,
    log1p: bool,
) -> PyResult<Vec<f64>> {
    let mut adata = build_adata(data, n_features)?;
    let config = NormalizeConfig {
        target_sum,
        log_transform: log1p,
        save_raw: false,
    };
    cyanea_omics::sc_preprocess::normalize_total(&mut adata, &config).into_pyresult()?;
    Ok(adata.x().to_flat_row_major())
}

/// Identify highly variable genes.
///
/// Returns a PyHvgResult with gene indices and normalized dispersions.
#[pyfunction]
#[pyo3(signature = (data, n_features, method="seurat_v3", n_top_genes=2000))]
fn hvg(data: Vec<f64>, n_features: usize, method: &str, n_top_genes: usize) -> PyResult<PyHvgResult> {
    let mut adata = build_adata(data, n_features)?;
    let hvg_method = match method {
        "cell_ranger" => HvgMethod::CellRanger,
        _ => HvgMethod::SeuratV3,
    };
    let config = HvgConfig {
        n_top_genes,
        method: hvg_method,
        ..HvgConfig::default()
    };
    cyanea_omics::sc_preprocess::highly_variable_genes(&mut adata, &config).into_pyresult()?;
    let hvg_col = adata.get_var("highly_variable");
    let disp_col = adata.get_var("dispersions_norm");
    let mut gene_indices = Vec::new();
    let mut dispersions = Vec::new();
    if let Some(cyanea_omics::ColumnData::Numeric(hvg_vals)) = hvg_col {
        for (i, &v) in hvg_vals.iter().enumerate() {
            if v > 0.0 {
                gene_indices.push(i);
                if let Some(cyanea_omics::ColumnData::Numeric(d)) = disp_col {
                    dispersions.push(d[i]);
                }
            }
        }
    }
    Ok(PyHvgResult {
        gene_indices,
        dispersions,
    })
}

/// Regress out covariates from the expression matrix.
///
/// Returns the corrected expression matrix as a flat row-major list.
#[pyfunction]
fn regress_out(data: Vec<f64>, n_features: usize, covariates: Vec<Vec<f64>>) -> PyResult<Vec<f64>> {
    let mut adata = build_adata(data, n_features)?;
    let mut keys = Vec::new();
    for (i, cov) in covariates.iter().enumerate() {
        let key = format!("covariate_{i}");
        adata.add_obs_numeric(&key, cov.clone()).into_pyresult()?;
        keys.push(key);
    }
    let key_refs: Vec<&str> = keys.iter().map(|s| s.as_str()).collect();
    cyanea_omics::sc_preprocess::regress_out(&mut adata, &key_refs).into_pyresult()?;
    Ok(adata.x().to_flat_row_major())
}

/// Score a set of genes per cell.
///
/// Returns a list of scores (one per cell).
#[pyfunction]
fn score_genes(data: Vec<f64>, n_features: usize, gene_indices: Vec<usize>) -> PyResult<Vec<f64>> {
    let mut adata = build_adata(data, n_features)?;
    cyanea_omics::sc_preprocess::score_genes(&mut adata, &gene_indices, 50, "score")
        .into_pyresult()?;
    match adata.get_obs("score") {
        Some(cyanea_omics::ColumnData::Numeric(scores)) => Ok(scores.clone()),
        _ => Err(pyo3::exceptions::PyRuntimeError::new_err(
            "score_genes did not produce scores",
        )),
    }
}

// ---------------------------------------------------------------------------
// Clustering
// ---------------------------------------------------------------------------

/// Compute a k-nearest neighbors graph.
///
/// Returns a PyNeighborsResult with distance and connectivity sparse matrices.
#[pyfunction]
#[pyo3(signature = (data, n_features, n_neighbors=15, metric="euclidean"))]
fn neighbors(
    data: Vec<f64>,
    n_features: usize,
    n_neighbors: usize,
    metric: &str,
) -> PyResult<PyNeighborsResult> {
    let mut adata = build_adata(data, n_features)?;
    let dist_metric = match metric {
        "cosine" => DistanceMetric::Cosine,
        _ => DistanceMetric::Euclidean,
    };
    let config = NeighborsConfig {
        n_neighbors,
        n_pcs: 0,
        metric: dist_metric,
        seed: 42,
    };
    cyanea_omics::sc_cluster::neighbors(&mut adata, &config).into_pyresult()?;
    let distances = adata
        .get_obsp("distances")
        .map(|s| s.iter().collect::<Vec<_>>())
        .unwrap_or_default();
    let connectivities = adata
        .get_obsp("connectivities")
        .map(|s| s.iter().collect::<Vec<_>>())
        .unwrap_or_default();
    Ok(PyNeighborsResult {
        distances,
        connectivities,
    })
}

/// Leiden clustering on a kNN graph.
///
/// Takes sparse distance and connectivity matrices as lists of (row, col, value) tuples.
/// Returns a PyClusterResult with cluster labels.
#[pyfunction]
#[pyo3(signature = (distances, connectivities, n_obs, resolution=1.0))]
fn leiden(
    distances: Vec<(usize, usize, f64)>,
    connectivities: Vec<(usize, usize, f64)>,
    n_obs: usize,
    resolution: f64,
) -> PyResult<PyClusterResult> {
    cluster_impl(distances, connectivities, n_obs, resolution, true)
}

/// Louvain clustering on a kNN graph.
///
/// Takes sparse distance and connectivity matrices as lists of (row, col, value) tuples.
/// Returns a PyClusterResult with cluster labels.
#[pyfunction]
#[pyo3(signature = (distances, connectivities, n_obs, resolution=1.0))]
fn louvain(
    distances: Vec<(usize, usize, f64)>,
    connectivities: Vec<(usize, usize, f64)>,
    n_obs: usize,
    resolution: f64,
) -> PyResult<PyClusterResult> {
    cluster_impl(distances, connectivities, n_obs, resolution, false)
}

fn cluster_impl(
    distances: Vec<(usize, usize, f64)>,
    connectivities: Vec<(usize, usize, f64)>,
    n_obs: usize,
    resolution: f64,
    use_leiden: bool,
) -> PyResult<PyClusterResult> {
    let dummy = vec![vec![0.0; 1]; n_obs];
    let obs_names: Vec<String> = (0..n_obs).map(|i| format!("cell_{i}")).collect();
    let var_names = vec!["dummy".to_string()];
    let mut adata =
        AnnData::new(MatrixData::Dense(dummy), obs_names, var_names).into_pyresult()?;

    let mut dist_sm = SparseMatrix::new(n_obs, n_obs);
    for (r, c, v) in distances {
        let _ = dist_sm.insert(r, c, v);
    }
    adata.add_obsp("distances", dist_sm).into_pyresult()?;

    let mut conn_sm = SparseMatrix::new(n_obs, n_obs);
    for (r, c, v) in connectivities {
        let _ = conn_sm.insert(r, c, v);
    }
    adata.add_obsp("connectivities", conn_sm).into_pyresult()?;

    let config = ClusterConfig {
        resolution,
        n_iterations: 10,
        seed: 42,
        key_added: "cluster".to_string(),
    };
    if use_leiden {
        cyanea_omics::sc_cluster::leiden(&mut adata, &config).into_pyresult()?;
    } else {
        cyanea_omics::sc_cluster::louvain(&mut adata, &config).into_pyresult()?;
    }

    extract_cluster_labels(&adata)
}

fn extract_cluster_labels(adata: &AnnData) -> PyResult<PyClusterResult> {
    match adata.get_obs("cluster") {
        Some(cyanea_omics::ColumnData::Categorical { codes, .. }) => {
            let labels: Vec<usize> = codes.iter().map(|&c| c as usize).collect();
            let n_clusters = labels.iter().copied().max().map_or(0, |m| m + 1);
            Ok(PyClusterResult {
                labels,
                modularity: n_clusters as f64 / adata.n_obs() as f64,
            })
        }
        Some(cyanea_omics::ColumnData::Strings(labels)) => {
            let mut label_map = HashMap::new();
            let mut next_id = 0usize;
            let int_labels: Vec<usize> = labels
                .iter()
                .map(|l| {
                    *label_map.entry(l.clone()).or_insert_with(|| {
                        let id = next_id;
                        next_id += 1;
                        id
                    })
                })
                .collect();
            Ok(PyClusterResult {
                labels: int_labels,
                modularity: 0.0,
            })
        }
        _ => Err(pyo3::exceptions::PyRuntimeError::new_err(
            "clustering did not produce labels",
        )),
    }
}

// ---------------------------------------------------------------------------
// Trajectory
// ---------------------------------------------------------------------------

/// Compute a diffusion map embedding.
///
/// Returns a PyDiffusionResult with per-cell components and eigenvalues.
#[pyfunction]
#[pyo3(signature = (data, n_features, n_components=10))]
fn diffusion_map(
    data: Vec<f64>,
    n_features: usize,
    n_components: usize,
) -> PyResult<PyDiffusionResult> {
    let mut adata = build_adata(data, n_features)?;
    let config = DiffusionConfig {
        n_components,
        alpha: 1.0,
    };
    let result =
        cyanea_omics::sc_trajectory::diffusion_map(&mut adata, &config).into_pyresult()?;
    Ok(PyDiffusionResult {
        components: result.components,
        eigenvalues: result.eigenvalues,
    })
}

/// Compute diffusion pseudotime from a diffusion map result.
///
/// Returns a PyDptResult with pseudotime values.
#[pyfunction]
fn dpt(components: Vec<Vec<f64>>, eigenvalues: Vec<f64>, root_cell: usize) -> PyResult<PyDptResult> {
    let n_obs = components.len();
    if n_obs == 0 {
        return Err(pyo3::exceptions::PyValueError::new_err("empty components"));
    }
    let dummy = vec![vec![0.0; 1]; n_obs];
    let obs_names: Vec<String> = (0..n_obs).map(|i| format!("cell_{i}")).collect();
    let var_names = vec!["dummy".to_string()];
    let mut adata =
        AnnData::new(MatrixData::Dense(dummy), obs_names, var_names).into_pyresult()?;
    adata.add_obsm("X_diffmap", components).into_pyresult()?;
    let evals_str = eigenvalues
        .iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(",");
    adata.add_uns("diffmap_evals", format!("[{evals_str}]"));
    let config = DptConfig {
        root_cell,
        n_branchings: 0,
    };
    cyanea_omics::sc_trajectory::dpt(&mut adata, &config).into_pyresult()?;
    match adata.get_obs("dpt_pseudotime") {
        Some(cyanea_omics::ColumnData::Numeric(pt)) => Ok(PyDptResult {
            pseudotime: pt.clone(),
        }),
        _ => Err(pyo3::exceptions::PyRuntimeError::new_err(
            "DPT did not produce pseudotime",
        )),
    }
}

/// Compute PAGA graph abstraction.
///
/// Takes sparse connectivities and cluster labels, returns connectivity between clusters.
#[pyfunction]
fn paga(
    connectivities: Vec<(usize, usize, f64)>,
    clusters: Vec<String>,
) -> PyResult<PyPagaResult> {
    let n_obs = clusters.len();
    let dummy = vec![vec![0.0; 1]; n_obs];
    let obs_names: Vec<String> = (0..n_obs).map(|i| format!("cell_{i}")).collect();
    let var_names = vec!["dummy".to_string()];
    let mut adata =
        AnnData::new(MatrixData::Dense(dummy), obs_names, var_names).into_pyresult()?;
    adata.add_obs("leiden", clusters).into_pyresult()?;
    let mut sm = SparseMatrix::new(n_obs, n_obs);
    for (r, c, v) in connectivities {
        let _ = sm.insert(r, c, v);
    }
    adata.add_obsp("connectivities", sm).into_pyresult()?;
    let result = cyanea_omics::sc_trajectory::paga(&adata, "leiden").into_pyresult()?;
    Ok(PyPagaResult {
        connectivities: result.connectivities,
        groups: result.groups,
        cluster_sizes: result.cluster_sizes,
    })
}

// ---------------------------------------------------------------------------
// Markers
// ---------------------------------------------------------------------------

/// Rank genes per cluster by differential expression.
///
/// Returns a PyMarkerResult with per-cluster marker gene lists.
#[pyfunction]
#[pyo3(signature = (data, n_features, clusters, method="t-test"))]
fn rank_genes_groups(
    data: Vec<f64>,
    n_features: usize,
    clusters: Vec<String>,
    method: &str,
) -> PyResult<PyMarkerResult> {
    let mut adata = build_adata(data, n_features)?;
    adata.add_obs("cluster", clusters).into_pyresult()?;
    let marker_method = match method {
        "wilcoxon" => MarkerMethod::Wilcoxon,
        "logistic" => MarkerMethod::LogisticRegression,
        _ => MarkerMethod::TTest,
    };
    let config = MarkerConfig {
        method: marker_method,
        cluster_key: "cluster".to_string(),
        log2fc_threshold: 0.0,
        min_pct: 0.0,
        padj_threshold: 1.0,
        n_genes: None,
    };
    let results =
        cyanea_omics::sc_markers::rank_genes_groups(&adata, &config).into_pyresult()?;
    let markers: HashMap<String, Vec<PyMarkerGene>> = results
        .markers
        .iter()
        .map(|(k, genes)| {
            let py_genes: Vec<PyMarkerGene> = genes
                .iter()
                .map(|g| PyMarkerGene {
                    gene_index: g.gene_index,
                    gene_name: g.gene_name.clone(),
                    score: g.statistic,
                    pvalue: g.p_value,
                    padj: g.p_adjusted,
                    log2fc: g.log2_fold_change,
                    pct_in: g.pct_in,
                    pct_out: g.pct_out,
                })
                .collect();
            (k.clone(), py_genes)
        })
        .collect();
    Ok(PyMarkerResult {
        markers,
        n_clusters: results.n_clusters,
    })
}

/// Filter marker genes by fold-change, pct, and p-value thresholds.
///
/// Returns a filtered PyMarkerResult.
#[pyfunction]
#[pyo3(signature = (markers, log2fc_threshold=1.0, min_pct=0.1, padj_threshold=0.05))]
fn filter_markers(
    markers: &PyMarkerResult,
    log2fc_threshold: f64,
    min_pct: f64,
    padj_threshold: f64,
) -> PyResult<PyMarkerResult> {
    let filtered: HashMap<String, Vec<PyMarkerGene>> = markers
        .markers
        .iter()
        .map(|(k, genes)| {
            let filt: Vec<PyMarkerGene> = genes
                .iter()
                .filter(|g| {
                    g.log2fc.abs() >= log2fc_threshold
                        && g.pct_in >= min_pct
                        && g.padj <= padj_threshold
                })
                .map(|g| PyMarkerGene {
                    gene_index: g.gene_index,
                    gene_name: g.gene_name.clone(),
                    score: g.score,
                    pvalue: g.pvalue,
                    padj: g.padj,
                    log2fc: g.log2fc,
                    pct_in: g.pct_in,
                    pct_out: g.pct_out,
                })
                .collect();
            (k.clone(), filt)
        })
        .collect();
    Ok(PyMarkerResult {
        markers: filtered,
        n_clusters: markers.n_clusters,
    })
}

// ---------------------------------------------------------------------------
// Integration
// ---------------------------------------------------------------------------

/// Run Harmony batch correction.
///
/// Returns a PyCorrectedResult with the corrected expression matrix.
#[pyfunction]
#[pyo3(signature = (data, n_features, batch, n_clusters=None))]
fn harmony(
    data: Vec<f64>,
    n_features: usize,
    batch: Vec<String>,
    n_clusters: Option<usize>,
) -> PyResult<PyCorrectedResult> {
    let mut adata = build_adata(data, n_features)?;
    adata.add_obs("batch", batch).into_pyresult()?;
    let config = HarmonyConfig {
        batch_key: "batch".to_string(),
        n_clusters,
        theta: 2.0,
        sigma: 0.1,
        max_iter: 10,
    };
    cyanea_omics::sc_integrate::harmony(&mut adata, &config).into_pyresult()?;
    let (n_obs, n_vars) = adata.shape();
    Ok(PyCorrectedResult {
        corrected: adata.x().to_flat_row_major(),
        n_obs,
        n_vars,
    })
}

/// Run ComBat batch correction.
///
/// Returns a PyCorrectedResult with the corrected expression matrix.
#[pyfunction]
fn combat(data: Vec<f64>, n_features: usize, batch: Vec<String>) -> PyResult<PyCorrectedResult> {
    let mut adata = build_adata(data, n_features)?;
    adata.add_obs("batch", batch).into_pyresult()?;
    let config = CombatConfig {
        batch_key: "batch".to_string(),
        parametric: true,
    };
    cyanea_omics::sc_integrate::combat(&mut adata, &config).into_pyresult()?;
    let (n_obs, n_vars) = adata.shape();
    Ok(PyCorrectedResult {
        corrected: adata.x().to_flat_row_major(),
        n_obs,
        n_vars,
    })
}

/// Run MNN batch correction.
///
/// Returns a PyCorrectedResult with the corrected expression matrix.
#[pyfunction]
#[pyo3(signature = (data, n_features, batch, k=20))]
fn mnn_correct(
    data: Vec<f64>,
    n_features: usize,
    batch: Vec<String>,
    k: usize,
) -> PyResult<PyCorrectedResult> {
    let mut adata = build_adata(data, n_features)?;
    adata.add_obs("batch", batch).into_pyresult()?;
    let config = MnnConfig {
        batch_key: "batch".to_string(),
        k,
        sigma: 1.0,
        cos_norm: true,
    };
    cyanea_omics::sc_integrate::mnn_correct(&mut adata, &config).into_pyresult()?;
    let (n_obs, n_vars) = adata.shape();
    Ok(PyCorrectedResult {
        corrected: adata.x().to_flat_row_major(),
        n_obs,
        n_vars,
    })
}

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "sc")?;
    // Result classes
    m.add_class::<PyHvgResult>()?;
    m.add_class::<PyNeighborsResult>()?;
    m.add_class::<PyClusterResult>()?;
    m.add_class::<PyDiffusionResult>()?;
    m.add_class::<PyDptResult>()?;
    m.add_class::<PyPagaResult>()?;
    m.add_class::<PyMarkerGene>()?;
    m.add_class::<PyMarkerResult>()?;
    m.add_class::<PyCorrectedResult>()?;
    // Preprocessing
    m.add_function(wrap_pyfunction!(normalize_total, &m)?)?;
    m.add_function(wrap_pyfunction!(hvg, &m)?)?;
    m.add_function(wrap_pyfunction!(regress_out, &m)?)?;
    m.add_function(wrap_pyfunction!(score_genes, &m)?)?;
    // Clustering
    m.add_function(wrap_pyfunction!(neighbors, &m)?)?;
    m.add_function(wrap_pyfunction!(leiden, &m)?)?;
    m.add_function(wrap_pyfunction!(louvain, &m)?)?;
    // Trajectory
    m.add_function(wrap_pyfunction!(diffusion_map, &m)?)?;
    m.add_function(wrap_pyfunction!(dpt, &m)?)?;
    m.add_function(wrap_pyfunction!(paga, &m)?)?;
    // Markers
    m.add_function(wrap_pyfunction!(rank_genes_groups, &m)?)?;
    m.add_function(wrap_pyfunction!(filter_markers, &m)?)?;
    // Integration
    m.add_function(wrap_pyfunction!(harmony, &m)?)?;
    m.add_function(wrap_pyfunction!(combat, &m)?)?;
    m.add_function(wrap_pyfunction!(mnn_correct, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
