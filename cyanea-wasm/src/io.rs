//! In-memory SAM pileup bindings for WASM environments.
//!
//! Parses SAM text, generates pileup, and computes depth statistics.

use std::collections::HashMap;

use serde::Serialize;

use crate::error::{wasm_err, wasm_ok};

#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

// ── Wrapper types ────────────────────────────────────────────────────────

/// A single pileup column for JSON serialization.
#[derive(Debug, Serialize)]
pub struct JsPileupColumn {
    pub pos: u64,
    pub ref_base: String,
    pub depth: u32,
    pub base_counts: HashMap<String, u32>,
}

/// A pileup for a single reference sequence.
#[derive(Debug, Serialize)]
pub struct JsPileup {
    pub rname: String,
    pub columns: Vec<JsPileupColumn>,
}

/// Depth statistics for a single reference sequence.
#[derive(Debug, Serialize)]
pub struct JsDepthStats {
    pub rname: String,
    pub length: u64,
    pub covered: u64,
    pub breadth: f64,
    pub min_depth: u32,
    pub max_depth: u32,
    pub mean_depth: f64,
    pub median_depth: f64,
}

// ── Helper ───────────────────────────────────────────────────────────────

/// Convert base_counts array [A, C, G, T, N, del] into a string-keyed map.
fn base_counts_map(counts: &[u32; 6]) -> HashMap<String, u32> {
    let mut map = HashMap::new();
    let labels = ["A", "C", "G", "T", "N", "del"];
    for (i, &count) in counts.iter().enumerate() {
        if count > 0 {
            map.insert(labels[i].to_string(), count);
        }
    }
    map
}

// ── JSON boundary functions ──────────────────────────────────────────────

/// Generate pileup from SAM text.
///
/// Parses SAM-formatted text and generates per-position pileup data.
/// Returns JSON array of pileups (one per reference sequence).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn pileup_from_sam(sam_text: &str) -> String {
    let records = match cyanea_io::parse_sam_str(sam_text) {
        Ok(r) => r,
        Err(e) => return wasm_err(e),
    };
    match cyanea_io::pileup(&records, None) {
        Ok(pileups) => {
            let js: Vec<JsPileup> = pileups
                .iter()
                .map(|p| JsPileup {
                    rname: p.rname.clone(),
                    columns: p
                        .columns
                        .iter()
                        .map(|c| JsPileupColumn {
                            pos: c.pos,
                            ref_base: String::from(c.ref_base as char),
                            depth: c.depth,
                            base_counts: base_counts_map(&c.base_counts),
                        })
                        .collect(),
                })
                .collect();
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Compute depth statistics from SAM text.
///
/// Parses SAM-formatted text and returns per-reference depth statistics.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn depth_stats_from_sam(sam_text: &str) -> String {
    let records = match cyanea_io::parse_sam_str(sam_text) {
        Ok(r) => r,
        Err(e) => return wasm_err(e),
    };
    match cyanea_io::pileup(&records, None) {
        Ok(pileups) => {
            let js: Vec<JsDepthStats> = pileups
                .iter()
                .map(|p| {
                    let ds = cyanea_io::depth_stats(p);
                    JsDepthStats {
                        rname: ds.rname,
                        length: ds.length,
                        covered: ds.covered,
                        breadth: ds.breadth,
                        min_depth: ds.min_depth,
                        max_depth: ds.max_depth,
                        mean_depth: ds.mean_depth,
                        median_depth: ds.median_depth,
                    }
                })
                .collect();
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Convert SAM text to mpileup format.
///
/// Parses SAM-formatted text and produces mpileup text output.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn pileup_to_mpileup_text(sam_text: &str) -> String {
    let records = match cyanea_io::parse_sam_str(sam_text) {
        Ok(r) => r,
        Err(e) => return wasm_err(e),
    };
    match cyanea_io::pileup(&records, None) {
        Ok(pileups) => {
            let mut result = String::new();
            for p in &pileups {
                let mpileup = cyanea_io::pileup_to_mpileup(p);
                if !result.is_empty() {
                    result.push('\n');
                }
                result.push_str(&mpileup);
            }
            wasm_ok(&result)
        }
        Err(e) => wasm_err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAM_TEXT: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*
r2\t0\tchr1\t3\t60\t8M\t*\t0\t0\tGTACGTAC\t*
";

    #[test]
    fn pileup_from_sam_basic() {
        let json = pileup_from_sam(SAM_TEXT);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let pileups = v["ok"].as_array().unwrap();
        assert_eq!(pileups.len(), 1);
        assert_eq!(pileups[0]["rname"], "chr1");
        assert!(pileups[0]["columns"].as_array().unwrap().len() > 0);
    }

    #[test]
    fn depth_stats_from_sam_basic() {
        let json = depth_stats_from_sam(SAM_TEXT);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let stats = v["ok"].as_array().unwrap();
        assert_eq!(stats.len(), 1);
        assert_eq!(stats[0]["rname"], "chr1");
        assert!(stats[0]["max_depth"].as_u64().unwrap() >= 1);
    }

    #[test]
    fn mpileup_text_from_sam() {
        let json = pileup_to_mpileup_text(SAM_TEXT);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let text = v["ok"].as_str().unwrap();
        assert!(text.contains("chr1"));
    }

    #[test]
    fn pileup_from_invalid_sam() {
        let json = pileup_from_sam("not\tvalid\tsam");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }
}
