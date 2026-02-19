//! Sequence alignment wrappers with string-based API.
//!
//! Wraps [`cyanea_align::align`] so that all inputs are plain strings and the
//! output is JSON-serialized [`AlignmentResult`].

use serde::Deserialize;

use cyanea_align::{align, AlignmentMode, AlignmentResult, ScoringMatrix, ScoringScheme, SubstitutionMatrix};
use cyanea_align::msa;
use cyanea_align::poa::{PoaGraph, PoaScoring};
use cyanea_align::simd::{banded_nw, banded_sw, banded_semi_global};
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

// ── CIGAR utility functions ──────────────────────────────────────────────

/// CIGAR statistics returned as a JSON object.
#[derive(serde::Serialize)]
struct CigarStats {
    cigar_string: String,
    reference_consumed: usize,
    query_consumed: usize,
    alignment_columns: usize,
    identity: f64,
    gap_count: usize,
    gap_bases: usize,
    soft_clipped: usize,
    hard_clipped: usize,
}

/// Parse a SAM CIGAR string and return the operations as JSON.
///
/// Returns a JSON array of CIGAR operations (e.g. `[{"AlnMatch":10},{"Insertion":3}]`).
/// Accepts the full SAM alphabet (M, I, D, N, S, H, P, =, X) and `*` for unavailable.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn parse_cigar(cigar: &str) -> String {
    wasm_result(cyanea_align::cigar::parse_cigar(cigar))
}

/// Validate a CIGAR string against SAM spec rules.
///
/// Returns `{"ok": true}` if valid, or `{"error": "..."}` describing the violation.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn validate_cigar(cigar: &str) -> String {
    let ops = match cyanea_align::cigar::parse_cigar(cigar) {
        Ok(ops) => ops,
        Err(e) => return wasm_err(e),
    };
    match cyanea_align::cigar::validate_cigar(&ops) {
        Ok(()) => wasm_ok(&true),
        Err(e) => wasm_err(e),
    }
}

/// Compute statistics from a CIGAR string.
///
/// Returns a JSON object with reference/query consumed, alignment columns,
/// identity, gap counts, and clipping totals.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn cigar_stats(cigar: &str) -> String {
    let ops = match cyanea_align::cigar::parse_cigar(cigar) {
        Ok(ops) => ops,
        Err(e) => return wasm_err(e),
    };
    let (soft, hard) = cyanea_align::cigar::clipped_bases(&ops);
    wasm_ok(&CigarStats {
        cigar_string: cyanea_align::cigar::cigar_string(&ops),
        reference_consumed: cyanea_align::cigar::reference_consumed(&ops),
        query_consumed: cyanea_align::cigar::query_consumed(&ops),
        alignment_columns: cyanea_align::cigar::alignment_columns(&ops),
        identity: cyanea_align::cigar::identity(&ops),
        gap_count: cyanea_align::cigar::gap_count(&ops),
        gap_bases: cyanea_align::cigar::gap_bases(&ops),
        soft_clipped: soft,
        hard_clipped: hard,
    })
}

/// Reconstruct gapped alignment from CIGAR + ungapped sequences.
///
/// Returns `{"ok": {"aligned_query": [...], "aligned_target": [...]}}`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn cigar_to_alignment(cigar: &str, query: &str, target: &str) -> String {
    let ops = match cyanea_align::cigar::parse_cigar(cigar) {
        Ok(ops) => ops,
        Err(e) => return wasm_err(e),
    };
    match cyanea_align::cigar::cigar_to_alignment(&ops, query.as_bytes(), target.as_bytes()) {
        Ok((aq, at)) => {
            #[derive(serde::Serialize)]
            struct GappedAlignment {
                aligned_query: Vec<u8>,
                aligned_target: Vec<u8>,
            }
            wasm_ok(&GappedAlignment {
                aligned_query: aq,
                aligned_target: at,
            })
        }
        Err(e) => wasm_err(e),
    }
}

/// Extract CIGAR from a gapped alignment (using =/X distinction).
///
/// Both `query` and `target` must be gapped strings (same length, `-` for gaps).
/// Returns a CIGAR string.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn alignment_to_cigar(query: &str, target: &str) -> String {
    match cyanea_align::cigar::alignment_to_cigar(query.as_bytes(), target.as_bytes()) {
        Ok(ops) => wasm_ok(&cyanea_align::cigar::cigar_string(&ops)),
        Err(e) => wasm_err(e),
    }
}

/// Generate a SAM MD:Z tag from CIGAR + ungapped sequences.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn generate_md_tag(cigar: &str, query: &str, reference: &str) -> String {
    let ops = match cyanea_align::cigar::parse_cigar(cigar) {
        Ok(ops) => ops,
        Err(e) => return wasm_err(e),
    };
    wasm_result(cyanea_align::cigar::generate_md_tag(
        &ops,
        query.as_bytes(),
        reference.as_bytes(),
    ))
}

/// Merge adjacent same-type CIGAR operations.
///
/// Returns the merged CIGAR string.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn merge_cigar(cigar: &str) -> String {
    let ops = match cyanea_align::cigar::parse_cigar(cigar) {
        Ok(ops) => ops,
        Err(e) => return wasm_err(e),
    };
    wasm_ok(&cyanea_align::cigar::cigar_string(
        &cyanea_align::cigar::merge_adjacent(&ops),
    ))
}

/// Reverse CIGAR operation order.
///
/// Returns the reversed CIGAR string.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn reverse_cigar(cigar: &str) -> String {
    let ops = match cyanea_align::cigar::parse_cigar(cigar) {
        Ok(ops) => ops,
        Err(e) => return wasm_err(e),
    };
    wasm_ok(&cyanea_align::cigar::cigar_string(
        &cyanea_align::cigar::reverse_cigar(&ops),
    ))
}

/// Collapse =/X operations into M (alignment match).
///
/// Returns the collapsed CIGAR string.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn collapse_cigar(cigar: &str) -> String {
    let ops = match cyanea_align::cigar::parse_cigar(cigar) {
        Ok(ops) => ops,
        Err(e) => return wasm_err(e),
    };
    wasm_ok(&cyanea_align::cigar::cigar_string(
        &cyanea_align::cigar::collapse_matches(&ops),
    ))
}

/// Convert hard clips (H) to soft clips (S).
///
/// Returns the converted CIGAR string.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn hard_clip_to_soft(cigar: &str) -> String {
    let ops = match cyanea_align::cigar::parse_cigar(cigar) {
        Ok(ops) => ops,
        Err(e) => return wasm_err(e),
    };
    wasm_ok(&cyanea_align::cigar::cigar_string(
        &cyanea_align::cigar::hard_clip_to_soft(&ops),
    ))
}

/// Split CIGAR at a reference coordinate, returning two CIGAR strings.
///
/// Returns `{"ok": {"left": "...", "right": "..."}}`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn split_cigar(cigar: &str, ref_pos: usize) -> String {
    let ops = match cyanea_align::cigar::parse_cigar(cigar) {
        Ok(ops) => ops,
        Err(e) => return wasm_err(e),
    };
    let (left, right) = cyanea_align::cigar::split_at_reference(&ops, ref_pos);
    #[derive(serde::Serialize)]
    struct SplitResult {
        left: String,
        right: String,
    }
    wasm_ok(&SplitResult {
        left: cyanea_align::cigar::cigar_string(&left),
        right: cyanea_align::cigar::cigar_string(&right),
    })
}

// ── MSA, POA, and banded alignment functions ─────────────────────────────

/// Result of a progressive multiple sequence alignment.
#[derive(Debug, serde::Serialize)]
pub struct JsMsaResult {
    pub aligned: Vec<String>,
    pub n_columns: usize,
    pub n_sequences: usize,
}

/// Result of a POA consensus computation.
#[derive(Debug, serde::Serialize)]
pub struct JsPoaConsensus {
    pub consensus: String,
    pub n_sequences: usize,
    pub n_nodes: usize,
}

/// Perform progressive multiple sequence alignment.
///
/// `seqs_json` is a JSON array of sequence strings (e.g. `["ACGT","ACTT"]`).
/// Returns a JSON object with `aligned`, `n_columns`, and `n_sequences`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn progressive_msa(
    seqs_json: &str,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> String {
    let seqs: Vec<String> = match serde_json::from_str(seqs_json) {
        Ok(s) => s,
        Err(e) => return wasm_err(format!("invalid JSON sequences: {e}")),
    };
    let scoring = ScoringScheme::Simple(ScoringMatrix {
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
    });
    let byte_seqs: Vec<&[u8]> = seqs.iter().map(|s| s.as_bytes()).collect();
    match msa::progressive_msa(&byte_seqs, &scoring) {
        Ok(result) => {
            let aligned: Vec<String> = result
                .aligned
                .iter()
                .map(|s| String::from_utf8_lossy(s).into_owned())
                .collect();
            let n_columns = result.n_columns;
            let n_sequences = result.n_sequences();
            wasm_ok(&JsMsaResult {
                aligned,
                n_columns,
                n_sequences,
            })
        }
        Err(e) => wasm_err(e),
    }
}

/// Compute a POA consensus from multiple sequences.
///
/// `seqs_json` is a JSON array of sequence strings. The first sequence
/// initializes the graph; subsequent sequences are aligned and integrated.
/// Returns a JSON object with `consensus`, `n_sequences`, and `n_nodes`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn poa_consensus(seqs_json: &str) -> String {
    let seqs: Vec<String> = match serde_json::from_str(seqs_json) {
        Ok(s) => s,
        Err(e) => return wasm_err(format!("invalid JSON sequences: {e}")),
    };
    if seqs.is_empty() {
        return wasm_err("need at least one sequence");
    }
    let scoring = PoaScoring {
        match_score: 2,
        mismatch_score: -1,
        gap_score: -2,
    };
    let mut graph = PoaGraph::from_sequence(seqs[0].as_bytes());
    for seq in &seqs[1..] {
        if let Err(e) = graph.add_sequence(seq.as_bytes(), &scoring) {
            return wasm_err(e);
        }
    }
    let consensus = graph.consensus();
    wasm_ok(&JsPoaConsensus {
        consensus: String::from_utf8_lossy(&consensus).into_owned(),
        n_sequences: graph.num_sequences(),
        n_nodes: graph.num_nodes(),
    })
}

/// Perform banded alignment between two sequences.
///
/// `mode` is `"global"`, `"local"`, or `"semiglobal"`.
/// `bandwidth` controls the diagonal band width (actual band is `2 * bandwidth + 1`).
/// Returns a JSON `AlignmentResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn align_banded(
    query: &str,
    target: &str,
    mode: &str,
    bandwidth: usize,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> String {
    let scoring = ScoringScheme::Simple(ScoringMatrix {
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
    });
    let result = match mode.to_ascii_lowercase().as_str() {
        "global" => banded_nw(query.as_bytes(), target.as_bytes(), &scoring, bandwidth),
        "local" => banded_sw(query.as_bytes(), target.as_bytes(), &scoring, bandwidth),
        "semiglobal" | "semi-global" => {
            banded_semi_global(query.as_bytes(), target.as_bytes(), &scoring, bandwidth)
        }
        _ => {
            return wasm_err(CyaneaError::InvalidInput(format!(
                "unknown alignment mode: {mode:?} (expected \"global\", \"local\", or \"semiglobal\")"
            )));
        }
    };
    wasm_result(result)
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

    // -- CIGAR utilities --

    #[test]
    fn parse_cigar_valid() {
        let json = parse_cigar("10M3I4D2S");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let ops = v["ok"].as_array().unwrap();
        assert_eq!(ops.len(), 4);
        assert_eq!(ops[0]["AlnMatch"], 10);
        assert_eq!(ops[1]["Insertion"], 3);
        assert_eq!(ops[2]["Deletion"], 4);
        assert_eq!(ops[3]["SoftClip"], 2);
    }

    #[test]
    fn parse_cigar_invalid() {
        let json = parse_cigar("10Z");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn validate_cigar_valid() {
        let json = validate_cigar("5H3S10M2I3M1D5M3S5H");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"], true);
    }

    #[test]
    fn validate_cigar_invalid() {
        let json = validate_cigar("5M3M"); // adjacent same-type
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn cigar_stats_basic() {
        let json = cigar_stats("10M3I4D2S5H");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let s = &v["ok"];
        assert_eq!(s["reference_consumed"], 14); // M(10) + D(4)
        assert_eq!(s["query_consumed"], 15);     // M(10) + I(3) + S(2)
        assert_eq!(s["alignment_columns"], 17);  // M(10) + I(3) + D(4)
        assert_eq!(s["gap_count"], 2);           // I, D
        assert_eq!(s["gap_bases"], 7);           // 3 + 4
        assert_eq!(s["soft_clipped"], 2);
        assert_eq!(s["hard_clipped"], 5);
    }

    #[test]
    fn cigar_to_alignment_roundtrip() {
        let json = cigar_to_alignment("3=1I2=1D1=", "ACGTACG", "ACGACGA");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let aq: Vec<u8> = v["ok"]["aligned_query"]
            .as_array().unwrap().iter()
            .map(|x| x.as_u64().unwrap() as u8).collect();
        let at: Vec<u8> = v["ok"]["aligned_target"]
            .as_array().unwrap().iter()
            .map(|x| x.as_u64().unwrap() as u8).collect();
        assert_eq!(aq, b"ACGTAC-G");
        assert_eq!(at, b"ACG-ACGA");
    }

    #[test]
    fn alignment_to_cigar_basic() {
        let json = alignment_to_cigar("ACGTAC-G", "ACG-ACGA");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"], "3=1I2=1D1X");
    }

    #[test]
    fn md_tag_generation() {
        let json = generate_md_tag("3=1X2=", "ACGAAC", "ACGTAC");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"], "3T2");
    }

    #[test]
    fn merge_cigar_combines() {
        let json = merge_cigar("3M2M1I1I");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"], "5M2I");
    }

    #[test]
    fn reverse_cigar_order() {
        let json = reverse_cigar("3=2I1D");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"], "1D2I3=");
    }

    #[test]
    fn collapse_cigar_eq_x_to_m() {
        let json = collapse_cigar("3=2X1=");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"], "6M");
    }

    #[test]
    fn hard_clip_to_soft_conversion() {
        let json = hard_clip_to_soft("5H10M3H");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"], "5S10M3S");
    }

    #[test]
    fn split_cigar_basic() {
        let json = split_cigar("10M", 4);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["ok"]["left"], "4M");
        assert_eq!(v["ok"]["right"], "6M");
    }

    // -- MSA, POA, and banded alignment --

    #[test]
    fn progressive_msa_basic() {
        let seqs = r#"["ACGT","ACGT","ACGT"]"#;
        let json = progressive_msa(seqs, 2, -1, -5, -2);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert_eq!(result["n_sequences"], 3);
        assert_eq!(result["n_columns"], 4);
        let aligned = result["aligned"].as_array().unwrap();
        assert_eq!(aligned.len(), 3);
        for a in aligned {
            assert_eq!(a.as_str().unwrap(), "ACGT");
        }
    }

    #[test]
    fn progressive_msa_different() {
        let seqs = r#"["ACGT","ACTT"]"#;
        let json = progressive_msa(seqs, 2, -1, -5, -2);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert_eq!(result["n_sequences"], 2);
        let aligned = result["aligned"].as_array().unwrap();
        assert_eq!(aligned.len(), 2);
        // Both should have the same length (n_columns)
        let len0 = aligned[0].as_str().unwrap().len();
        let len1 = aligned[1].as_str().unwrap().len();
        assert_eq!(len0, len1);
        assert!(result["n_columns"].as_u64().unwrap() >= 4);
    }

    #[test]
    fn poa_consensus_identical() {
        let seqs = r#"["ACGT","ACGT","ACGT"]"#;
        let json = poa_consensus(seqs);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert_eq!(result["consensus"].as_str().unwrap(), "ACGT");
        assert_eq!(result["n_sequences"], 3);
    }

    #[test]
    fn poa_consensus_majority() {
        let seqs = r#"["ACGT","ACGT","ACTT"]"#;
        let json = poa_consensus(seqs);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        // Majority of sequences (2 out of 3) have ACGT
        assert_eq!(result["consensus"].as_str().unwrap(), "ACGT");
    }

    #[test]
    fn align_banded_global() {
        let json = align_banded("ACGT", "ACGT", "global", 5, 2, -1, -5, -2);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["ok"]["score"].as_i64().unwrap() > 0);
        assert_eq!(v["ok"]["score"].as_i64().unwrap(), 8);
    }

    #[test]
    fn align_banded_local() {
        let json = align_banded("XXXACGTXXX", "ACGT", "local", 10, 2, -1, -5, -2);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["ok"]["score"].as_i64().unwrap() > 0);
    }
}
