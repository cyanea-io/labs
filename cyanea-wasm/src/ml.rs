//! ML primitive wrappers with JSON input/output.
//!
//! Provides k-mer counting and distance metric functions that accept strings
//! and return JSON, suitable for the WASM boundary.

use std::collections::HashMap;

use serde::Serialize;

use cyanea_ml::distance;
use cyanea_ml::kmer::KmerCounter;

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
}
