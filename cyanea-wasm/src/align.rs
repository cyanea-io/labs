//! Sequence alignment wrappers with string-based API.
//!
//! Wraps [`cyanea_align::align`] so that all inputs are plain strings and the
//! output is JSON-serialized [`AlignmentResult`].

use serde::Deserialize;

use cyanea_align::{align, AlignmentMode, AlignmentResult, ScoringMatrix, ScoringScheme, SubstitutionMatrix};
use cyanea_core::CyaneaError;

use crate::error::{wasm_err, wasm_ok, wasm_result};

#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

/// Parse an alignment mode string ("local", "global") into [`AlignmentMode`].
fn parse_mode(mode: &str) -> cyanea_core::Result<AlignmentMode> {
    match mode.to_ascii_lowercase().as_str() {
        "local" => Ok(AlignmentMode::Local),
        "global" => Ok(AlignmentMode::Global),
        "semiglobal" | "semi-global" => Ok(AlignmentMode::SemiGlobal),
        _ => Err(CyaneaError::InvalidInput(format!(
            "unknown alignment mode: {mode:?} (expected \"local\" or \"global\")"
        ))),
    }
}

/// Parse a substitution matrix name into [`SubstitutionMatrix`].
fn parse_matrix(name: &str) -> cyanea_core::Result<SubstitutionMatrix> {
    match name.to_ascii_lowercase().as_str() {
        "blosum62" => Ok(SubstitutionMatrix::blosum62()),
        "blosum45" => Ok(SubstitutionMatrix::blosum45()),
        "blosum80" => Ok(SubstitutionMatrix::blosum80()),
        "pam250" => Ok(SubstitutionMatrix::pam250()),
        _ => Err(CyaneaError::InvalidInput(format!(
            "unknown substitution matrix: {name:?} \
             (expected \"blosum62\", \"blosum45\", \"blosum80\", or \"pam250\")"
        ))),
    }
}

fn do_align(
    query: &[u8],
    target: &[u8],
    mode: AlignmentMode,
    scoring: &ScoringScheme,
) -> cyanea_core::Result<AlignmentResult> {
    if query.is_empty() || target.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "query and target must be non-empty".into(),
        ));
    }
    align(query, target, mode, scoring)
}

// ── JSON boundary functions ──────────────────────────────────────────────

/// Align two DNA sequences with default scoring (+2/-1/-5/-2).
///
/// `mode` is `"local"`, `"global"`, or `"semiglobal"`. Returns JSON `AlignmentResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn align_dna(query: &str, target: &str, mode: &str) -> String {
    let m = match parse_mode(mode) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
    wasm_result(do_align(query.as_bytes(), target.as_bytes(), m, &scoring))
}

/// Align two DNA sequences with custom scoring parameters.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn align_dna_custom(
    query: &str,
    target: &str,
    mode: &str,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> String {
    let m = match parse_mode(mode) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let scoring = ScoringScheme::Simple(ScoringMatrix {
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
    });
    wasm_result(do_align(query.as_bytes(), target.as_bytes(), m, &scoring))
}

/// Align two protein sequences using a named substitution matrix.
///
/// `matrix` is one of `"blosum62"`, `"blosum45"`, `"blosum80"`, `"pam250"`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn align_protein(query: &str, target: &str, mode: &str, matrix: &str) -> String {
    let m = match parse_mode(mode) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let sub = match parse_matrix(matrix) {
        Ok(s) => s,
        Err(e) => return wasm_err(e),
    };
    let scoring = ScoringScheme::Substitution(sub);
    wasm_result(do_align(query.as_bytes(), target.as_bytes(), m, &scoring))
}

/// A query/target pair for batch alignment, deserialized from JSON.
#[derive(Deserialize)]
struct SeqPair {
    query: String,
    target: String,
}

/// Batch-align multiple sequence pairs with custom scoring.
///
/// `pairs_json` is a JSON array of `{"query": "...", "target": "..."}` objects.
/// `mode` is `"local"`, `"global"`, or `"semiglobal"`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn align_batch(
    pairs_json: &str,
    mode: &str,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> String {
    let m = match parse_mode(mode) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let pairs: Vec<SeqPair> = match serde_json::from_str(pairs_json) {
        Ok(p) => p,
        Err(e) => return wasm_err(format!("invalid JSON pairs: {e}")),
    };
    let scoring = ScoringScheme::Simple(ScoringMatrix {
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
    });
    let refs: Vec<(&[u8], &[u8])> = pairs
        .iter()
        .map(|p| (p.query.as_bytes(), p.target.as_bytes()))
        .collect();
    match cyanea_align::align_batch(&refs, m, &scoring) {
        Ok(results) => wasm_ok(&results),
        Err(e) => wasm_err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_global_perfect() {
        let json = align_dna("ACGT", "ACGT", "global");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["ok"]["score"].as_i64().unwrap() > 0);
    }

    #[test]
    fn dna_local_mismatch() {
        let json = align_dna("AAACGTAAA", "TTTCGTTTT", "local");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["ok"]["score"].as_i64().unwrap() > 0);
    }

    #[test]
    fn dna_custom_scoring() {
        let json = align_dna_custom("ACGT", "ACGT", "global", 3, -2, -7, -3);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["ok"]["score"].as_i64().unwrap() > 0);
    }

    #[test]
    fn protein_blosum62() {
        let json = align_protein("HEAGAWGHEE", "PAWHEAE", "global", "blosum62");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["ok"]["score"].is_number());
    }

    #[test]
    fn invalid_mode() {
        let json = align_dna("ACGT", "ACGT", "foobar");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].as_str().unwrap().contains("unknown alignment mode"));
    }

    #[test]
    fn empty_sequence() {
        let json = align_dna("", "ACGT", "global");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn invalid_matrix() {
        let json = align_protein("HEAG", "PAWH", "global", "unknown");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].as_str().unwrap().contains("unknown substitution matrix"));
    }

    #[test]
    fn batch_global() {
        let pairs = r#"[{"query":"ACGT","target":"ACGT"},{"query":"AAA","target":"TTT"}]"#;
        let json = align_batch(pairs, "global", 2, -1, -5, -2);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let results = v["ok"].as_array().unwrap();
        assert_eq!(results.len(), 2);
        assert_eq!(results[0]["score"], 8);
    }

    #[test]
    fn batch_invalid_json() {
        let json = align_batch("not json", "global", 2, -1, -5, -2);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn dna_semiglobal() {
        let json = align_dna("CGT", "AACGTAA", "semiglobal");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"]["score"].as_i64().unwrap(), 6);
    }
}
