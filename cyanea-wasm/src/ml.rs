//! ML primitive wrappers with JSON input/output.
//!
//! Provides k-mer counting and distance metric functions that accept strings
//! and return JSON, suitable for the WASM boundary.

use std::collections::HashMap;

use serde::Serialize;

use cyanea_ml::distance;
use cyanea_ml::kmer::KmerCounter;
use cyanea_ml::forest::{RandomForest, RandomForestConfig};
use cyanea_ml::gbdt::{GbdtConfig, GradientBoostedTrees};
use cyanea_ml::hmm::HmmModel;
use cyanea_ml::metrics::{self, ConfusionMatrix};
use cyanea_ml::cross_validation;
use cyanea_ml::feature_selection;

use crate::error::{wasm_err, wasm_ok};

#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

// ── Wrapper types ────────────────────────────────────────────────────────

/// Serializable k-mer counts with string keys (instead of `Vec<u8>`).
#[derive(Debug, Serialize)]
pub struct JsKmerCounts {
    pub k: usize,
    pub total: usize,
    pub distinct: usize,
    pub counts: HashMap<String, usize>,
}

// ── JSON boundary functions ──────────────────────────────────────────────

fn parse_f64_array(json: &str) -> Result<Vec<f64>, String> {
    serde_json::from_str::<Vec<f64>>(json).map_err(|e| format!("invalid JSON array: {e}"))
}

fn parse_usize_array(json: &str) -> Result<Vec<usize>, String> {
    serde_json::from_str::<Vec<usize>>(json).map_err(|e| format!("invalid JSON array: {e}"))
}

fn parse_bool_array(json: &str) -> Result<Vec<bool>, String> {
    serde_json::from_str::<Vec<bool>>(json).map_err(|e| format!("invalid JSON array: {e}"))
}

/// Count k-mers in a nucleotide/protein sequence string.
///
/// Returns JSON `JsKmerCounts` with string keys.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn kmer_count(seq: &str, k: usize) -> String {
    let counter = match KmerCounter::new(k) {
        Ok(c) => c,
        Err(e) => return wasm_err(e),
    };
    let counts = counter.count_sequence(seq.as_bytes());
    let js = JsKmerCounts {
        k,
        total: counts.total(),
        distinct: counts.distinct(),
        counts: counts
            .iter()
            .map(|(kmer, count)| (String::from_utf8_lossy(kmer).into_owned(), count))
            .collect(),
    };
    wasm_ok(&js)
}

/// Euclidean distance between two JSON arrays of numbers.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn euclidean_distance(a_json: &str, b_json: &str) -> String {
    let a = match parse_f64_array(a_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let b = match parse_f64_array(b_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match distance::euclidean(&a, &b) {
        Ok(d) => wasm_ok(&d),
        Err(e) => wasm_err(e),
    }
}

/// Manhattan distance between two JSON arrays of numbers.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn manhattan_distance(a_json: &str, b_json: &str) -> String {
    let a = match parse_f64_array(a_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let b = match parse_f64_array(b_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match distance::manhattan(&a, &b) {
        Ok(d) => wasm_ok(&d),
        Err(e) => wasm_err(e),
    }
}

/// Hamming distance between two raw strings (byte-level comparison).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn hamming_distance(a: &str, b: &str) -> String {
    match distance::hamming(a.as_bytes(), b.as_bytes()) {
        Ok(d) => wasm_ok(&d),
        Err(e) => wasm_err(e),
    }
}

/// Cosine similarity between two JSON arrays of numbers.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn cosine_similarity(a_json: &str, b_json: &str) -> String {
    let a = match parse_f64_array(a_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let b = match parse_f64_array(b_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match distance::cosine_similarity(&a, &b) {
        Ok(d) => wasm_ok(&d),
        Err(e) => wasm_err(e),
    }
}

// ── UMAP ────────────────────────────────────────────────────────────────

/// Serializable UMAP result.
#[derive(Debug, Serialize)]
pub struct JsUmapResult {
    pub embedding: Vec<f64>,
    pub n_samples: usize,
    pub n_components: usize,
    pub n_epochs: usize,
}

// ── PCA ──────────────────────────────────────────────────────────────────

/// Serializable PCA result.
#[derive(Debug, Serialize)]
pub struct JsPcaResult {
    pub components: Vec<f64>,
    pub explained_variance: Vec<f64>,
    pub explained_variance_ratio: Vec<f64>,
    pub transformed: Vec<f64>,
    pub mean: Vec<f64>,
    pub n_features: usize,
    pub n_components: usize,
}

// ── t-SNE ────────────────────────────────────────────────────────────────

/// Serializable t-SNE result.
#[derive(Debug, Serialize)]
pub struct JsTsneResult {
    pub embedding: Vec<f64>,
    pub n_samples: usize,
    pub n_components: usize,
    pub kl_divergence: f64,
}

// ── K-means ──────────────────────────────────────────────────────────────

/// Serializable K-means result.
#[derive(Debug, Serialize)]
pub struct JsKmeansResult {
    pub centroids: Vec<f64>,
    pub labels: Vec<usize>,
    pub inertia: f64,
    pub n_iter: usize,
    pub n_features: usize,
}

// ── Random Forest ────────────────────────────────────────────────────────

/// Serializable random forest classification result.
#[derive(Debug, Serialize)]
pub struct JsRandomForestResult {
    pub predictions: Vec<usize>,
    pub feature_importance: Vec<f64>,
    pub n_trees: usize,
    pub n_classes: usize,
}

// ── GBDT Regression ─────────────────────────────────────────────────────

/// Serializable GBDT regression result.
#[derive(Debug, Serialize)]
pub struct JsGbdtRegressionResult {
    pub predictions: Vec<f64>,
    pub feature_importance: Vec<f64>,
    pub n_estimators: usize,
}

// ── GBDT Classification ─────────────────────────────────────────────────

/// Serializable GBDT classification result.
#[derive(Debug, Serialize)]
pub struct JsGbdtClassifyResult {
    pub predictions: Vec<usize>,
    pub feature_importance: Vec<f64>,
    pub n_estimators: usize,
    pub n_classes: usize,
}

// ── HMM Viterbi ─────────────────────────────────────────────────────────

/// Serializable HMM Viterbi result.
#[derive(Debug, Serialize)]
pub struct JsHmmViterbiResult {
    pub path: Vec<usize>,
    pub log_probability: f64,
}

// ── Confusion Matrix ────────────────────────────────────────────────────

/// Serializable confusion matrix result.
#[derive(Debug, Serialize)]
pub struct JsConfusionMatrix {
    pub matrix: Vec<usize>,
    pub n_classes: usize,
    pub accuracy: f64,
}

// ── ROC Curve ───────────────────────────────────────────────────────────

/// Serializable ROC curve result.
#[derive(Debug, Serialize)]
pub struct JsRocCurve {
    pub points: Vec<JsRocPoint>,
    pub auc: f64,
}

/// A single point on the ROC curve.
#[derive(Debug, Serialize)]
pub struct JsRocPoint {
    pub threshold: f64,
    pub fpr: f64,
    pub tpr: f64,
}

// ── PR Curve ────────────────────────────────────────────────────────────

/// Serializable precision-recall curve result.
#[derive(Debug, Serialize)]
pub struct JsPrCurve {
    pub points: Vec<JsPrPoint>,
    pub auc: f64,
}

/// A single point on the PR curve.
#[derive(Debug, Serialize)]
pub struct JsPrPoint {
    pub threshold: f64,
    pub precision: f64,
    pub recall: f64,
}

// ── Cross-validation ────────────────────────────────────────────────────

/// Serializable cross-validation result.
#[derive(Debug, Serialize)]
pub struct JsCvResult {
    pub mean_score: f64,
    pub std_score: f64,
    pub fold_scores: Vec<f64>,
}

// ── Feature Selection ───────────────────────────────────────────────────

/// Serializable feature selection result.
#[derive(Debug, Serialize)]
pub struct JsFeatureSelection {
    pub selected: Vec<usize>,
    pub scores: Vec<f64>,
    pub n_features: usize,
}

/// Parse a distance metric string.
fn parse_metric(s: &str) -> Result<cyanea_ml::DistanceMetric, String> {
    match s {
        "euclidean" => Ok(cyanea_ml::DistanceMetric::Euclidean),
        "manhattan" => Ok(cyanea_ml::DistanceMetric::Manhattan),
        "cosine" => Ok(cyanea_ml::DistanceMetric::Cosine),
        _ => Err(format!(
            "unknown metric: {s} (expected euclidean, manhattan, or cosine)"
        )),
    }
}

/// UMAP dimensionality reduction.
///
/// `data_json`: JSON array of numbers (flat row-major matrix).
/// `n_features`: number of features per sample.
/// `n_components`: output dimensionality (default 2).
/// `n_neighbors`: number of nearest neighbors (default 15).
/// `min_dist`: minimum distance in embedding (default 0.1).
/// `n_epochs`: optimization epochs (default 200).
/// `metric`: distance metric ("euclidean", "manhattan", "cosine").
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn umap(
    data_json: &str,
    n_features: usize,
    n_components: usize,
    n_neighbors: usize,
    min_dist: f64,
    n_epochs: usize,
    metric: &str,
) -> String {
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let metric = match parse_metric(metric) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let config = cyanea_ml::UmapConfig {
        n_components,
        n_neighbors,
        min_dist,
        n_epochs,
        metric,
        ..Default::default()
    };
    match cyanea_ml::umap(&data, n_features, &config) {
        Ok(r) => {
            let js = JsUmapResult {
                embedding: r.embedding,
                n_samples: r.n_samples,
                n_components: r.n_components,
                n_epochs: r.n_epochs,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// PCA dimensionality reduction.
///
/// `data_json`: JSON array of numbers (flat row-major matrix).
/// `n_features`: number of features per sample.
/// `n_components`: output dimensionality.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn pca(data_json: &str, n_features: usize, n_components: usize) -> String {
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let config = cyanea_ml::PcaConfig {
        n_components,
        ..Default::default()
    };
    match cyanea_ml::pca(&data, n_features, &config) {
        Ok(r) => {
            let js = JsPcaResult {
                components: r.components,
                explained_variance: r.explained_variance,
                explained_variance_ratio: r.explained_variance_ratio,
                transformed: r.transformed,
                mean: r.mean,
                n_features: r.n_features,
                n_components: r.n_components,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// t-SNE dimensionality reduction.
///
/// `data_json`: JSON array of numbers (flat row-major matrix).
/// `n_features`: number of features per sample.
/// `n_components`: output dimensionality (typically 2 or 3).
/// `perplexity`: perplexity parameter (5-50 typical).
/// `learning_rate`: learning rate.
/// `n_iter`: number of iterations.
/// `seed`: random seed.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn tsne(
    data_json: &str,
    n_features: usize,
    n_components: usize,
    perplexity: f64,
    learning_rate: f64,
    n_iter: usize,
    seed: u64,
) -> String {
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let config = cyanea_ml::TsneConfig {
        n_components,
        perplexity,
        learning_rate,
        n_iter,
        seed,
    };
    match cyanea_ml::tsne(&data, n_features, &config) {
        Ok(r) => {
            let js = JsTsneResult {
                embedding: r.embedding,
                n_samples: r.n_samples,
                n_components: r.n_components,
                kl_divergence: r.kl_divergence,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// K-means clustering.
///
/// `data_json`: JSON array of numbers (flat row-major matrix).
/// `n_features`: number of features per sample.
/// `n_clusters`: number of clusters.
/// `max_iter`: maximum iterations.
/// `seed`: random seed.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn kmeans(
    data_json: &str,
    n_features: usize,
    n_clusters: usize,
    max_iter: usize,
    seed: u64,
) -> String {
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let config = cyanea_ml::KMeansConfig {
        n_clusters,
        max_iter,
        seed,
        ..Default::default()
    };
    let rows: Vec<&[f64]> = data.chunks_exact(n_features).collect();
    match cyanea_ml::kmeans(&rows, &config) {
        Ok(r) => {
            let js = JsKmeansResult {
                centroids: r.centroids,
                labels: r.labels,
                inertia: r.inertia,
                n_iter: r.n_iter,
                n_features: r.n_features,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

// ── Random Forest Classification ─────────────────────────────────────────

/// Random forest classification.
///
/// `data_json`: JSON array of numbers (flat row-major matrix).
/// `n_features`: number of features per sample.
/// `labels_json`: JSON array of class labels (usize).
/// `n_trees`: number of trees in the ensemble.
/// `max_depth`: maximum depth per tree.
/// `seed`: random seed.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn random_forest_classify(
    data_json: &str,
    n_features: usize,
    labels_json: &str,
    n_trees: usize,
    max_depth: usize,
    seed: u64,
) -> String {
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let labels = match parse_usize_array(labels_json) {
        Ok(l) => l,
        Err(e) => return wasm_err(e),
    };
    let config = RandomForestConfig {
        n_trees,
        max_depth,
        max_features: None,
        seed,
    };
    match RandomForest::fit(&data, n_features, &labels, &config) {
        Ok(rf) => {
            let predictions = rf.predict_batch(&data, n_features);
            let feature_importance = rf.feature_importance(n_features);
            let js = JsRandomForestResult {
                predictions,
                feature_importance,
                n_trees: rf.n_trees(),
                n_classes: rf.n_classes(),
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

// ── GBDT Regression ──────────────────────────────────────────────────────

/// Gradient boosted decision tree regression.
///
/// `data_json`: JSON array of numbers (flat row-major matrix).
/// `n_features`: number of features per sample.
/// `targets_json`: JSON array of target values (f64).
/// `n_estimators`: number of boosting rounds.
/// `learning_rate`: shrinkage factor.
/// `max_depth`: maximum depth per tree.
/// `seed`: random seed.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn gbdt_regression(
    data_json: &str,
    n_features: usize,
    targets_json: &str,
    n_estimators: usize,
    learning_rate: f64,
    max_depth: usize,
    seed: u64,
) -> String {
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let targets = match parse_f64_array(targets_json) {
        Ok(t) => t,
        Err(e) => return wasm_err(e),
    };
    let config = GbdtConfig {
        n_estimators,
        learning_rate,
        max_depth,
        seed,
        ..Default::default()
    };
    match GradientBoostedTrees::fit_regression(&data, n_features, &targets, &config) {
        Ok(model) => {
            let predictions = model.predict_batch(&data, n_features);
            let feature_importance = model.feature_importance();
            let js = JsGbdtRegressionResult {
                predictions,
                feature_importance,
                n_estimators: model.n_estimators(),
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

// ── GBDT Classification ─────────────────────────────────────────────────

/// Gradient boosted decision tree classification.
///
/// `data_json`: JSON array of numbers (flat row-major matrix).
/// `n_features`: number of features per sample.
/// `labels_json`: JSON array of class labels (usize).
/// `n_estimators`: number of boosting rounds.
/// `learning_rate`: shrinkage factor.
/// `max_depth`: maximum depth per tree.
/// `seed`: random seed.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn gbdt_classify(
    data_json: &str,
    n_features: usize,
    labels_json: &str,
    n_estimators: usize,
    learning_rate: f64,
    max_depth: usize,
    seed: u64,
) -> String {
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let labels = match parse_usize_array(labels_json) {
        Ok(l) => l,
        Err(e) => return wasm_err(e),
    };
    let config = GbdtConfig {
        n_estimators,
        learning_rate,
        max_depth,
        seed,
        ..Default::default()
    };
    match GradientBoostedTrees::fit_classification(&data, n_features, &labels, &config) {
        Ok(model) => {
            let predictions = model.predict_class_batch(&data, n_features);
            let feature_importance = model.feature_importance();
            let js = JsGbdtClassifyResult {
                predictions,
                feature_importance,
                n_estimators: model.n_estimators(),
                n_classes: model.n_classes(),
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

// ── HMM Viterbi ─────────────────────────────────────────────────────────

/// HMM Viterbi decoding.
///
/// `n_states`: number of hidden states.
/// `n_symbols`: number of observable symbols.
/// `initial_json`: JSON array of initial probabilities.
/// `transition_json`: JSON array of transition probabilities (row-major).
/// `emission_json`: JSON array of emission probabilities (row-major).
/// `observations_json`: JSON array of observation indices.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn hmm_viterbi(
    n_states: usize,
    n_symbols: usize,
    initial_json: &str,
    transition_json: &str,
    emission_json: &str,
    observations_json: &str,
) -> String {
    let initial = match parse_f64_array(initial_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let transition = match parse_f64_array(transition_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let emission = match parse_f64_array(emission_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let observations = match parse_usize_array(observations_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let model = match HmmModel::new(n_states, n_symbols, initial, transition, emission) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    match model.viterbi(&observations) {
        Ok((path, log_probability)) => {
            let js = JsHmmViterbiResult {
                path,
                log_probability,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

// ── HMM Likelihood ──────────────────────────────────────────────────────

/// HMM log-likelihood of an observation sequence.
///
/// `n_states`: number of hidden states.
/// `n_symbols`: number of observable symbols.
/// `initial_json`: JSON array of initial probabilities.
/// `transition_json`: JSON array of transition probabilities (row-major).
/// `emission_json`: JSON array of emission probabilities (row-major).
/// `observations_json`: JSON array of observation indices.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn hmm_likelihood(
    n_states: usize,
    n_symbols: usize,
    initial_json: &str,
    transition_json: &str,
    emission_json: &str,
    observations_json: &str,
) -> String {
    let initial = match parse_f64_array(initial_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let transition = match parse_f64_array(transition_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let emission = match parse_f64_array(emission_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let observations = match parse_usize_array(observations_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let model = match HmmModel::new(n_states, n_symbols, initial, transition, emission) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    match model.log_likelihood(&observations) {
        Ok(ll) => wasm_ok(&ll),
        Err(e) => wasm_err(e),
    }
}

// ── Confusion Matrix ────────────────────────────────────────────────────

/// Compute a confusion matrix from actual and predicted labels.
///
/// `actual_json`: JSON array of actual class labels (usize).
/// `predicted_json`: JSON array of predicted class labels (usize).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn confusion_matrix(actual_json: &str, predicted_json: &str) -> String {
    let actual = match parse_usize_array(actual_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let predicted = match parse_usize_array(predicted_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match ConfusionMatrix::from_labels(&actual, &predicted, None) {
        Ok(cm) => {
            let js = JsConfusionMatrix {
                matrix: cm.matrix.clone(),
                n_classes: cm.n_classes,
                accuracy: cm.accuracy(),
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

// ── ROC Curve ───────────────────────────────────────────────────────────

/// Compute the ROC curve from predicted scores and binary labels.
///
/// `scores_json`: JSON array of predicted scores (f64).
/// `labels_json`: JSON array of binary labels (bool).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn roc_curve(scores_json: &str, labels_json: &str) -> String {
    let scores = match parse_f64_array(scores_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let labels = match parse_bool_array(labels_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match metrics::roc_curve(&scores, &labels) {
        Ok(roc) => {
            let points: Vec<JsRocPoint> = roc
                .points
                .iter()
                .map(|p| JsRocPoint {
                    threshold: p.threshold,
                    fpr: p.fpr,
                    tpr: p.tpr,
                })
                .collect();
            let js = JsRocCurve {
                points,
                auc: roc.auc,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

// ── PR Curve ────────────────────────────────────────────────────────────

/// Compute the precision-recall curve from predicted scores and binary labels.
///
/// `scores_json`: JSON array of predicted scores (f64).
/// `labels_json`: JSON array of binary labels (bool).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn pr_curve(scores_json: &str, labels_json: &str) -> String {
    let scores = match parse_f64_array(scores_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let labels = match parse_bool_array(labels_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match metrics::pr_curve(&scores, &labels) {
        Ok(pr) => {
            let points: Vec<JsPrPoint> = pr
                .points
                .iter()
                .map(|p| JsPrPoint {
                    threshold: p.threshold,
                    precision: p.precision,
                    recall: p.recall,
                })
                .collect();
            let js = JsPrCurve {
                points,
                auc: pr.auc,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

// ── Cross-validation (RF) ───────────────────────────────────────────────

/// K-fold cross-validation with a random forest classifier.
///
/// `data_json`: JSON array of numbers (flat row-major matrix).
/// `n_features`: number of features per sample.
/// `labels_json`: JSON array of class labels (usize).
/// `k`: number of folds.
/// `seed`: random seed.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn cross_validate_rf(
    data_json: &str,
    n_features: usize,
    labels_json: &str,
    k: usize,
    seed: u64,
) -> String {
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let labels = match parse_usize_array(labels_json) {
        Ok(l) => l,
        Err(e) => return wasm_err(e),
    };
    let n_samples = data.len() / n_features;
    match cross_validation::cross_validate_kfold(n_samples, k, seed, |train_idx, test_idx| {
        // Build training data
        let train_data: Vec<f64> = train_idx
            .iter()
            .flat_map(|&i| data[i * n_features..(i + 1) * n_features].iter().copied())
            .collect();
        let train_labels: Vec<usize> = train_idx.iter().map(|&i| labels[i]).collect();

        let config = RandomForestConfig {
            n_trees: 20,
            max_depth: 5,
            max_features: None,
            seed,
        };
        let rf = RandomForest::fit(&train_data, n_features, &train_labels, &config)?;

        // Predict on test data and compute accuracy
        let correct: usize = test_idx
            .iter()
            .filter(|&&i| {
                let sample = &data[i * n_features..(i + 1) * n_features];
                rf.predict(sample) == labels[i]
            })
            .count();
        Ok(correct as f64 / test_idx.len() as f64)
    }) {
        Ok(cv) => {
            let fold_scores: Vec<f64> = cv.folds.iter().map(|f| f.score).collect();
            let js = JsCvResult {
                mean_score: cv.mean_score,
                std_score: cv.std_score,
                fold_scores,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

// ── Feature Importance (Variance Threshold) ─────────────────────────────

/// Feature selection via variance threshold.
///
/// `data_json`: JSON array of numbers (flat row-major matrix).
/// `n_features`: number of features per sample.
/// `threshold`: minimum variance (features with variance > threshold are kept).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn feature_importance_variance(
    data_json: &str,
    n_features: usize,
    threshold: f64,
) -> String {
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match feature_selection::variance_threshold(&data, n_features, threshold) {
        Ok(fs) => {
            let js = JsFeatureSelection {
                selected: fs.selected,
                scores: fs.scores,
                n_features: fs.n_features,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmer_count_k2() {
        let json = kmer_count("ACGT", 2);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["k"], 2);
        assert_eq!(obj["total"], 3);
        assert_eq!(obj["distinct"], 3);
        assert_eq!(obj["counts"]["AC"], 1);
        assert_eq!(obj["counts"]["CG"], 1);
        assert_eq!(obj["counts"]["GT"], 1);
    }

    #[test]
    fn kmer_count_invalid_k() {
        let json = kmer_count("ACGT", 0);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn euclidean_known() {
        let json = euclidean_distance("[0,0]", "[3,4]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!((v["ok"].as_f64().unwrap() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn manhattan_known() {
        let json = manhattan_distance("[0,0]", "[3,4]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!((v["ok"].as_f64().unwrap() - 7.0).abs() < 1e-10);
    }

    #[test]
    fn hamming_known() {
        let json = hamming_distance("ACGT", "ACGA");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"], 1);
    }

    #[test]
    fn cosine_orthogonal() {
        let json = cosine_similarity("[1,0]", "[0,1]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!((v["ok"].as_f64().unwrap()).abs() < 1e-10);
    }

    #[test]
    fn cosine_identical() {
        let json = cosine_similarity("[1,2,3]", "[1,2,3]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!((v["ok"].as_f64().unwrap() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn pca_basic() {
        // 4 samples, 3 features each -> reduce to 2 components
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0];
        let json = pca(&serde_json::to_string(&data).unwrap(), 3, 2);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["n_components"], 2);
        assert_eq!(obj["n_features"], 3);
        assert!(obj["components"].as_array().unwrap().len() > 0);
        assert!(obj["explained_variance"].as_array().unwrap().len() > 0);
        assert!(obj["transformed"].as_array().unwrap().len() > 0);
    }

    #[test]
    fn pca_invalid_data() {
        let json = pca("not json", 3, 2);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn tsne_basic() {
        // 4 samples, 2 features each
        let data = vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let json = tsne(&serde_json::to_string(&data).unwrap(), 2, 2, 2.0, 100.0, 100, 42);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["n_components"], 2);
        assert_eq!(obj["n_samples"], 4);
        assert!(obj["embedding"].as_array().unwrap().len() == 8); // 4 samples * 2 components
    }

    #[test]
    fn tsne_invalid_data() {
        let json = tsne("[]", 2, 2, 30.0, 200.0, 1000, 42);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn kmeans_basic() {
        // 6 points in 2D, should cluster into 2 groups
        let data = vec![0.0, 0.0, 0.1, 0.1, 0.2, 0.0, 10.0, 10.0, 10.1, 10.1, 10.2, 10.0];
        let json = kmeans(&serde_json::to_string(&data).unwrap(), 2, 2, 100, 42);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["labels"].as_array().unwrap().len(), 6);
        assert_eq!(obj["n_features"], 2);
        // First 3 should be in same cluster, last 3 in another
        let labels: Vec<usize> = obj["labels"].as_array().unwrap()
            .iter().map(|l| l.as_u64().unwrap() as usize).collect();
        assert_eq!(labels[0], labels[1]);
        assert_eq!(labels[0], labels[2]);
        assert_eq!(labels[3], labels[4]);
        assert_eq!(labels[3], labels[5]);
        assert_ne!(labels[0], labels[3]);
    }

    #[test]
    fn kmeans_invalid_data() {
        let json = kmeans("bad", 2, 2, 100, 42);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    // ── Random Forest tests ─────────────────────────────────────────────

    #[test]
    fn random_forest_basic() {
        // 6 samples, 2 features, 2 classes
        let data = vec![0.0, 0.0, 0.1, 0.1, 0.2, 0.0, 10.0, 10.0, 10.1, 10.1, 10.2, 10.0];
        let labels = vec![0, 0, 0, 1, 1, 1];
        let json = random_forest_classify(
            &serde_json::to_string(&data).unwrap(),
            2,
            &serde_json::to_string(&labels).unwrap(),
            10,
            5,
            42,
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["n_trees"], 10);
        assert_eq!(obj["n_classes"], 2);
        assert_eq!(obj["predictions"].as_array().unwrap().len(), 6);
    }

    #[test]
    fn random_forest_importance() {
        let data = vec![0.0, 0.0, 0.1, 0.1, 0.2, 0.0, 10.0, 10.0, 10.1, 10.1, 10.2, 10.0];
        let labels = vec![0, 0, 0, 1, 1, 1];
        let json = random_forest_classify(
            &serde_json::to_string(&data).unwrap(),
            2,
            &serde_json::to_string(&labels).unwrap(),
            10,
            5,
            42,
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["feature_importance"].as_array().unwrap().len(), 2);
    }

    // ── GBDT Regression tests ───────────────────────────────────────────

    #[test]
    fn gbdt_regression_basic() {
        // Simple regression: y = 2*x
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let targets = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        let json = gbdt_regression(
            &serde_json::to_string(&data).unwrap(),
            1,
            &serde_json::to_string(&targets).unwrap(),
            50,
            0.1,
            5,
            42,
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert!(obj["n_estimators"].as_u64().unwrap() > 0);
        assert_eq!(obj["predictions"].as_array().unwrap().len(), 5);
    }

    #[test]
    fn gbdt_regression_predictions() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let targets = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        let json = gbdt_regression(
            &serde_json::to_string(&data).unwrap(),
            1,
            &serde_json::to_string(&targets).unwrap(),
            50,
            0.1,
            5,
            42,
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        let preds = obj["predictions"].as_array().unwrap();
        assert!(!preds.is_empty());
        // Check feature importance exists
        let fi = obj["feature_importance"].as_array().unwrap();
        assert_eq!(fi.len(), 1);
    }

    // ── GBDT Classification tests ───────────────────────────────────────

    #[test]
    fn gbdt_classify_basic() {
        // 2-class classification
        let data = vec![0.0, 0.0, 0.1, 0.1, 0.2, 0.0, 10.0, 10.0, 10.1, 10.1, 10.2, 10.0];
        let labels = vec![0, 0, 0, 1, 1, 1];
        let json = gbdt_classify(
            &serde_json::to_string(&data).unwrap(),
            2,
            &serde_json::to_string(&labels).unwrap(),
            30,
            0.1,
            3,
            42,
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["n_classes"], 2);
        assert_eq!(obj["predictions"].as_array().unwrap().len(), 6);
    }

    #[test]
    fn gbdt_classify_predictions() {
        let data = vec![0.0, 0.0, 0.1, 0.1, 0.2, 0.0, 10.0, 10.0, 10.1, 10.1, 10.2, 10.0];
        let labels = vec![0, 0, 0, 1, 1, 1];
        let json = gbdt_classify(
            &serde_json::to_string(&data).unwrap(),
            2,
            &serde_json::to_string(&labels).unwrap(),
            30,
            0.1,
            3,
            42,
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        let preds = obj["predictions"].as_array().unwrap();
        assert_eq!(preds.len(), 6);
        let fi = obj["feature_importance"].as_array().unwrap();
        assert_eq!(fi.len(), 2);
    }

    // ── HMM tests ───────────────────────────────────────────────────────

    #[test]
    fn hmm_viterbi_basic() {
        // 2 states, 2 symbols
        let initial = vec![0.5, 0.5];
        let transition = vec![0.9, 0.1, 0.2, 0.8];
        let emission = vec![0.5, 0.5, 0.8, 0.2];
        let observations = vec![0, 0, 1, 0, 0];
        let json = hmm_viterbi(
            2,
            2,
            &serde_json::to_string(&initial).unwrap(),
            &serde_json::to_string(&transition).unwrap(),
            &serde_json::to_string(&emission).unwrap(),
            &serde_json::to_string(&observations).unwrap(),
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["path"].as_array().unwrap().len(), 5);
        assert!(obj["log_probability"].as_f64().unwrap().is_finite());
    }

    #[test]
    fn hmm_viterbi_path_valid() {
        let initial = vec![0.5, 0.5];
        let transition = vec![0.9, 0.1, 0.2, 0.8];
        let emission = vec![0.5, 0.5, 0.8, 0.2];
        let observations = vec![0, 1, 0, 1, 0];
        let json = hmm_viterbi(
            2,
            2,
            &serde_json::to_string(&initial).unwrap(),
            &serde_json::to_string(&transition).unwrap(),
            &serde_json::to_string(&emission).unwrap(),
            &serde_json::to_string(&observations).unwrap(),
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let path = v["ok"]["path"].as_array().unwrap();
        for p in path {
            let state = p.as_u64().unwrap() as usize;
            assert!(state < 2, "state {} should be < n_states (2)", state);
        }
    }

    #[test]
    fn hmm_likelihood_basic() {
        let initial = vec![0.5, 0.5];
        let transition = vec![0.9, 0.1, 0.2, 0.8];
        let emission = vec![0.5, 0.5, 0.8, 0.2];
        let observations = vec![0, 0, 1, 0, 0];
        let json = hmm_likelihood(
            2,
            2,
            &serde_json::to_string(&initial).unwrap(),
            &serde_json::to_string(&transition).unwrap(),
            &serde_json::to_string(&emission).unwrap(),
            &serde_json::to_string(&observations).unwrap(),
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let ll = v["ok"].as_f64().unwrap();
        assert!(ll.is_finite());
        assert!(ll < 0.0, "log-likelihood should be negative, got {}", ll);
    }

    // ── Confusion Matrix tests ──────────────────────────────────────────

    #[test]
    fn confusion_matrix_basic() {
        let actual = vec![0, 0, 1, 1];
        let predicted = vec![0, 1, 0, 1];
        let json = confusion_matrix(
            &serde_json::to_string(&actual).unwrap(),
            &serde_json::to_string(&predicted).unwrap(),
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert_eq!(obj["n_classes"], 2);
        assert!((obj["accuracy"].as_f64().unwrap() - 0.5).abs() < 1e-10);
        assert_eq!(obj["matrix"].as_array().unwrap().len(), 4);
    }

    // ── ROC Curve tests ─────────────────────────────────────────────────

    #[test]
    fn roc_curve_basic() {
        let scores = vec![0.9, 0.8, 0.3, 0.1];
        let labels = vec![true, true, false, false];
        let json = roc_curve(
            &serde_json::to_string(&scores).unwrap(),
            &serde_json::to_string(&labels).unwrap(),
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        let auc = obj["auc"].as_f64().unwrap();
        assert!((auc - 1.0).abs() < 1e-10, "perfect separation should give AUC=1.0, got {}", auc);
        assert!(!obj["points"].as_array().unwrap().is_empty());
    }

    // ── PR Curve tests ──────────────────────────────────────────────────

    #[test]
    fn pr_curve_basic() {
        let scores = vec![0.9, 0.8, 0.3, 0.1];
        let labels = vec![true, true, false, false];
        let json = pr_curve(
            &serde_json::to_string(&scores).unwrap(),
            &serde_json::to_string(&labels).unwrap(),
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        let auc = obj["auc"].as_f64().unwrap();
        assert!(auc > 0.0, "PR AUC should be > 0, got {}", auc);
        assert!(!obj["points"].as_array().unwrap().is_empty());
    }

    // ── Cross-validation tests ──────────────────────────────────────────

    #[test]
    fn cross_validate_rf_basic() {
        // Well-separated 2-class data, k=2 fold CV
        let data = vec![
            0.0, 0.0, 0.1, 0.1, 0.2, 0.0, 0.0, 0.2,
            10.0, 10.0, 10.1, 10.1, 10.2, 10.0, 10.0, 10.2,
        ];
        let labels = vec![0, 0, 0, 0, 1, 1, 1, 1];
        let json = cross_validate_rf(
            &serde_json::to_string(&data).unwrap(),
            2,
            &serde_json::to_string(&labels).unwrap(),
            2,
            42,
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        assert!(obj["mean_score"].as_f64().unwrap() >= 0.0);
        assert_eq!(obj["fold_scores"].as_array().unwrap().len(), 2);
    }

    // ── Feature Selection tests ─────────────────────────────────────────

    #[test]
    fn feature_importance_variance_basic() {
        // Feature 0: varies, Feature 1: constant, Feature 2: varies
        let data = vec![
            1.0, 5.0, 0.0,
            2.0, 5.0, 1.0,
            3.0, 5.0, 2.0,
            4.0, 5.0, 3.0,
        ];
        let json = feature_importance_variance(
            &serde_json::to_string(&data).unwrap(),
            3,
            0.0,
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let obj = &v["ok"];
        let selected: Vec<usize> = obj["selected"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_u64().unwrap() as usize)
            .collect();
        // Feature 1 (constant) should be removed
        assert!(!selected.contains(&1));
        // Features 0 and 2 should be selected
        assert!(selected.contains(&0));
        assert!(selected.contains(&2));
        assert_eq!(obj["n_features"], 3);
    }
}
