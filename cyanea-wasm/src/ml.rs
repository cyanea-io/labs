//! ML primitive wrappers with JSON input/output.
//!
//! Provides k-mer counting and distance metric functions that accept strings
//! and return JSON, suitable for the WASM boundary.

use std::collections::HashMap;

use serde::Serialize;

use cyanea_ml::distance;
use cyanea_ml::kmer::KmerCounter;

use crate::error::{wasm_err, wasm_ok};

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
pub fn hamming_distance(a: &str, b: &str) -> String {
    match distance::hamming(a.as_bytes(), b.as_bytes()) {
        Ok(d) => wasm_ok(&d),
        Err(e) => wasm_err(e),
    }
}

/// Cosine similarity between two JSON arrays of numbers.
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
}
