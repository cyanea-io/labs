//! CIGAR string parsing, validation, arithmetic, and conversion utilities.
//!
//! Operates on the full SAM CIGAR alphabet (M, I, D, N, S, H, P, =, X).
//!
//! # Examples
//!
//! ```
//! use cyanea_align::cigar::{parse_cigar, cigar_string, reference_consumed, query_consumed};
//!
//! let ops = parse_cigar("10M3I4D2S").unwrap();
//! assert_eq!(cigar_string(&ops), "10M3I4D2S");
//! assert_eq!(reference_consumed(&ops), 14); // M + D
//! assert_eq!(query_consumed(&ops), 15);     // M + I + S
//! ```

use crate::types::CigarOp;
use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Parsing & formatting
// ---------------------------------------------------------------------------

/// Parse a SAM CIGAR string into a vector of operations.
///
/// Accepts the full SAM alphabet: `M`, `I`, `D`, `N`, `S`, `H`, `P`, `=`, `X`.
/// The special CIGAR `*` (unavailable) returns an empty vector.
///
/// # Errors
///
/// Returns [`CyaneaError::Parse`] on invalid characters or missing counts.
pub fn parse_cigar(s: &str) -> Result<Vec<CigarOp>> {
    if s == "*" {
        return Ok(Vec::new());
    }
    if s.is_empty() {
        return Err(CyaneaError::Parse("empty CIGAR string".into()));
    }

    let mut ops = Vec::new();
    let mut num_start: Option<usize> = None;

    for (i, c) in s.char_indices() {
        if c.is_ascii_digit() {
            if num_start.is_none() {
                num_start = Some(i);
            }
        } else {
            let start = num_start.ok_or_else(|| {
                CyaneaError::Parse(format!("CIGAR op '{}' at position {} has no count", c, i))
            })?;
            let n: usize = s[start..i].parse().map_err(|e| {
                CyaneaError::Parse(format!("invalid CIGAR count: {}", e))
            })?;
            let op = char_to_cigar_op(c, n).ok_or_else(|| {
                CyaneaError::Parse(format!("invalid CIGAR op character '{}'", c))
            })?;
            ops.push(op);
            num_start = None;
        }
    }

    if num_start.is_some() {
        return Err(CyaneaError::Parse(
            "CIGAR string ends with digits but no op character".into(),
        ));
    }

    Ok(ops)
}

/// Format CIGAR operations as a compact string (e.g. `"10M3I4D"`).
///
/// Returns `"*"` for an empty slice.
pub fn cigar_string(ops: &[CigarOp]) -> String {
    if ops.is_empty() {
        return "*".into();
    }
    let mut s = String::new();
    for op in ops {
        s.push_str(&format!("{}{}", op.len(), op.code()));
    }
    s
}

fn char_to_cigar_op(c: char, n: usize) -> Option<CigarOp> {
    match c {
        '=' => Some(CigarOp::Match(n)),
        'X' => Some(CigarOp::Mismatch(n)),
        'I' => Some(CigarOp::Insertion(n)),
        'D' => Some(CigarOp::Deletion(n)),
        'M' => Some(CigarOp::AlnMatch(n)),
        'N' => Some(CigarOp::Skip(n)),
        'S' => Some(CigarOp::SoftClip(n)),
        'H' => Some(CigarOp::HardClip(n)),
        'P' => Some(CigarOp::Padding(n)),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

/// Validate a CIGAR operation sequence according to SAM spec rules.
///
/// Checks:
/// - No zero-length operations
/// - No adjacent operations of the same type
/// - `H` (hard clip) only at the ends
/// - `S` (soft clip) adjacent to `H` or at the ends
///
/// # Errors
///
/// Returns [`CyaneaError::InvalidInput`] describing the first violation found.
pub fn validate_cigar(ops: &[CigarOp]) -> Result<()> {
    if ops.is_empty() {
        return Ok(());
    }

    for (i, op) in ops.iter().enumerate() {
        if op.is_empty() {
            return Err(CyaneaError::InvalidInput(format!(
                "zero-length CIGAR op at position {}", i
            )));
        }
    }

    // Adjacent same-type check
    for i in 1..ops.len() {
        if ops[i].code() == ops[i - 1].code() {
            return Err(CyaneaError::InvalidInput(format!(
                "adjacent CIGAR ops of same type '{}' at positions {} and {}",
                ops[i].code(), i - 1, i
            )));
        }
    }

    // H only at ends
    for (i, op) in ops.iter().enumerate() {
        if matches!(op, CigarOp::HardClip(_)) && i != 0 && i != ops.len() - 1 {
            return Err(CyaneaError::InvalidInput(format!(
                "hard clip (H) at interior position {}", i
            )));
        }
    }

    // S must be adjacent to H or at ends
    for (i, op) in ops.iter().enumerate() {
        if matches!(op, CigarOp::SoftClip(_)) {
            let at_start = i == 0 || (i == 1 && matches!(ops[0], CigarOp::HardClip(_)));
            let at_end = i == ops.len() - 1
                || (i == ops.len() - 2
                    && matches!(ops[ops.len() - 1], CigarOp::HardClip(_)));
            if !at_start && !at_end {
                return Err(CyaneaError::InvalidInput(format!(
                    "soft clip (S) at invalid position {}", i
                )));
            }
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Coordinate queries
// ---------------------------------------------------------------------------

/// Number of reference bases consumed (M, =, X, D, N).
pub fn reference_consumed(ops: &[CigarOp]) -> usize {
    ops.iter()
        .filter(|op| op.consumes_reference())
        .map(|op| op.len())
        .sum()
}

/// Number of query bases consumed (M, =, X, I, S).
pub fn query_consumed(ops: &[CigarOp]) -> usize {
    ops.iter()
        .filter(|op| op.consumes_query())
        .map(|op| op.len())
        .sum()
}

/// Number of aligned columns (M, =, X, I, D) — excludes clips, skips, padding.
pub fn alignment_columns(ops: &[CigarOp]) -> usize {
    ops.iter()
        .filter(|op| {
            matches!(
                op,
                CigarOp::Match(_)
                    | CigarOp::Mismatch(_)
                    | CigarOp::Insertion(_)
                    | CigarOp::Deletion(_)
                    | CigarOp::AlnMatch(_)
            )
        })
        .map(|op| op.len())
        .sum()
}

/// Total soft-clipped and hard-clipped bases as `(soft, hard)`.
pub fn clipped_bases(ops: &[CigarOp]) -> (usize, usize) {
    let mut soft = 0;
    let mut hard = 0;
    for op in ops {
        match op {
            CigarOp::SoftClip(n) => soft += n,
            CigarOp::HardClip(n) => hard += n,
            _ => {}
        }
    }
    (soft, hard)
}

// ---------------------------------------------------------------------------
// Statistics
// ---------------------------------------------------------------------------

/// Sequence identity: `=` bases / (M + = + X + I + D) columns.
///
/// Returns 0.0 for empty CIGAR. When only `M` ops are present (no `=`/`X`
/// distinction), returns 0.0 since exact matches cannot be determined.
pub fn identity(ops: &[CigarOp]) -> f64 {
    let cols = alignment_columns(ops);
    if cols == 0 {
        return 0.0;
    }
    let eq: usize = ops
        .iter()
        .filter_map(|op| match op {
            CigarOp::Match(n) => Some(n),
            _ => None,
        })
        .sum();
    eq as f64 / cols as f64
}

/// Fraction of reference covered: `reference_consumed(ops) / ref_len`.
///
/// Returns 0.0 if `ref_len` is 0.
pub fn coverage(ops: &[CigarOp], ref_len: usize) -> f64 {
    if ref_len == 0 {
        return 0.0;
    }
    reference_consumed(ops) as f64 / ref_len as f64
}

/// Number of gap openings (runs of I or D operations).
pub fn gap_count(ops: &[CigarOp]) -> usize {
    ops.iter()
        .filter(|op| matches!(op, CigarOp::Insertion(_) | CigarOp::Deletion(_)))
        .count()
}

/// Total bases in gaps (sum of I and D lengths).
pub fn gap_bases(ops: &[CigarOp]) -> usize {
    ops.iter()
        .filter_map(|op| match op {
            CigarOp::Insertion(n) | CigarOp::Deletion(n) => Some(n),
            _ => None,
        })
        .sum()
}

// ---------------------------------------------------------------------------
// Arithmetic
// ---------------------------------------------------------------------------

/// Merge consecutive operations of the same type.
pub fn merge_adjacent(ops: &[CigarOp]) -> Vec<CigarOp> {
    if ops.is_empty() {
        return Vec::new();
    }
    let mut merged = Vec::with_capacity(ops.len());
    merged.push(ops[0]);
    for &op in &ops[1..] {
        let last = merged.last_mut().unwrap();
        if last.code() == op.code() {
            *last = with_len(*last, last.len() + op.len());
        } else {
            merged.push(op);
        }
    }
    merged
}

/// Reverse the order of CIGAR operations.
pub fn reverse_cigar(ops: &[CigarOp]) -> Vec<CigarOp> {
    ops.iter().copied().rev().collect()
}

/// Split CIGAR at a reference coordinate, returning `(left, right)`.
///
/// `ref_pos` is 0-based; `left` covers `[0, ref_pos)`, `right` covers
/// `[ref_pos, end)`. Operations that don't consume reference (I, S, H, P)
/// are attached to whichever side they fall on during the left-to-right walk.
pub fn split_at_reference(ops: &[CigarOp], ref_pos: usize) -> (Vec<CigarOp>, Vec<CigarOp>) {
    let mut left = Vec::new();
    let mut right = Vec::new();
    let mut consumed = 0usize;
    let mut split = false;

    for &op in ops {
        if split {
            right.push(op);
            continue;
        }
        if !op.consumes_reference() {
            left.push(op);
            continue;
        }
        let remaining = ref_pos - consumed;
        if op.len() <= remaining {
            left.push(op);
            consumed += op.len();
            if consumed == ref_pos {
                split = true;
            }
        } else {
            // Split this op
            if remaining > 0 {
                left.push(with_len(op, remaining));
            }
            let rest = op.len() - remaining;
            right.push(with_len(op, rest));
            split = true;
        }
    }

    (left, right)
}

/// Convert hard clips (`H`) to soft clips (`S`).
pub fn hard_clip_to_soft(ops: &[CigarOp]) -> Vec<CigarOp> {
    ops.iter()
        .map(|op| match op {
            CigarOp::HardClip(n) => CigarOp::SoftClip(*n),
            other => *other,
        })
        .collect()
}

/// Collapse `=` and `X` operations into `M` (alignment match).
pub fn collapse_matches(ops: &[CigarOp]) -> Vec<CigarOp> {
    let collapsed: Vec<CigarOp> = ops
        .iter()
        .map(|op| match op {
            CigarOp::Match(n) | CigarOp::Mismatch(n) => CigarOp::AlnMatch(*n),
            other => *other,
        })
        .collect();
    merge_adjacent(&collapsed)
}

/// Create a CigarOp of the same variant but with a different length.
fn with_len(op: CigarOp, n: usize) -> CigarOp {
    match op {
        CigarOp::Match(_) => CigarOp::Match(n),
        CigarOp::Mismatch(_) => CigarOp::Mismatch(n),
        CigarOp::Insertion(_) => CigarOp::Insertion(n),
        CigarOp::Deletion(_) => CigarOp::Deletion(n),
        CigarOp::AlnMatch(_) => CigarOp::AlnMatch(n),
        CigarOp::Skip(_) => CigarOp::Skip(n),
        CigarOp::SoftClip(_) => CigarOp::SoftClip(n),
        CigarOp::HardClip(_) => CigarOp::HardClip(n),
        CigarOp::Padding(_) => CigarOp::Padding(n),
    }
}

// ---------------------------------------------------------------------------
// Conversion
// ---------------------------------------------------------------------------

/// Reconstruct gapped alignment sequences from CIGAR and ungapped input.
///
/// Walks the CIGAR operations against `query` and `target` (both ungapped),
/// producing two gapped sequences with `-` for gaps. Soft/hard clips, skips,
/// and padding are not included in the output alignment.
///
/// # Errors
///
/// Returns an error if the CIGAR consumes more bases than available in query
/// or target.
pub fn cigar_to_alignment(
    ops: &[CigarOp],
    query: &[u8],
    target: &[u8],
) -> Result<(Vec<u8>, Vec<u8>)> {
    let mut aligned_q = Vec::new();
    let mut aligned_t = Vec::new();
    let mut qi = 0usize;
    let mut ti = 0usize;

    for op in ops {
        match op {
            CigarOp::Match(n) | CigarOp::Mismatch(n) | CigarOp::AlnMatch(n) => {
                if qi + n > query.len() || ti + n > target.len() {
                    return Err(CyaneaError::InvalidInput(
                        "CIGAR consumes more bases than available in sequences".into(),
                    ));
                }
                aligned_q.extend_from_slice(&query[qi..qi + n]);
                aligned_t.extend_from_slice(&target[ti..ti + n]);
                qi += n;
                ti += n;
            }
            CigarOp::Insertion(n) => {
                if qi + n > query.len() {
                    return Err(CyaneaError::InvalidInput(
                        "CIGAR insertion consumes more query bases than available".into(),
                    ));
                }
                aligned_q.extend_from_slice(&query[qi..qi + n]);
                aligned_t.extend(core::iter::repeat(b'-').take(*n));
                qi += n;
            }
            CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                if ti + n > target.len() {
                    return Err(CyaneaError::InvalidInput(
                        "CIGAR deletion/skip consumes more target bases than available".into(),
                    ));
                }
                aligned_q.extend(core::iter::repeat(b'-').take(*n));
                aligned_t.extend_from_slice(&target[ti..ti + n]);
                ti += n;
            }
            CigarOp::SoftClip(n) => {
                qi += n;
            }
            CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
        }
    }

    Ok((aligned_q, aligned_t))
}

/// Extract CIGAR operations from a gapped alignment using `=`/`X` distinction.
///
/// Both `query` and `target` must be gapped (same length, with `-` for gaps).
///
/// # Errors
///
/// Returns an error if the sequences have different lengths.
pub fn alignment_to_cigar(query: &[u8], target: &[u8]) -> Result<Vec<CigarOp>> {
    if query.len() != target.len() {
        return Err(CyaneaError::InvalidInput(
            "gapped alignment sequences must have equal length".into(),
        ));
    }

    let mut ops: Vec<CigarOp> = Vec::new();

    for (q, t) in query.iter().zip(target.iter()) {
        let op = match (*q, *t) {
            (b'-', b'-') => continue,
            (b'-', _) => CigarOp::Deletion(1),
            (_, b'-') => CigarOp::Insertion(1),
            (a, b) if a == b => CigarOp::Match(1),
            _ => CigarOp::Mismatch(1),
        };
        // Merge with previous if same type
        if let Some(last) = ops.last_mut() {
            if last.code() == op.code() {
                *last = with_len(*last, last.len() + 1);
                continue;
            }
        }
        ops.push(op);
    }

    Ok(ops)
}

// ---------------------------------------------------------------------------
// MD tag
// ---------------------------------------------------------------------------

/// Generate a SAM MD:Z tag from CIGAR operations and ungapped sequences.
///
/// The MD tag records mismatched and deleted reference bases. It walks
/// reference-consuming CIGAR ops:
/// - `=` / match: increments match counter
/// - `X` / mismatch: emits match count + reference base
/// - `M`: compares query/reference and emits accordingly
/// - `D` / `N`: emits `^` + deleted reference bases
/// - `I`, `S`, `H`, `P`: skipped (don't consume reference)
///
/// # Errors
///
/// Returns an error if sequences are too short for the CIGAR.
pub fn generate_md_tag(
    ops: &[CigarOp],
    query: &[u8],
    reference: &[u8],
) -> Result<String> {
    let mut md = String::new();
    let mut match_count = 0usize;
    let mut qi = 0usize;
    let mut ri = 0usize;

    for op in ops {
        match op {
            CigarOp::Match(n) => {
                match_count += n;
                qi += n;
                ri += n;
            }
            CigarOp::Mismatch(n) => {
                if ri + n > reference.len() {
                    return Err(CyaneaError::InvalidInput(
                        "MD tag: reference too short for mismatch".into(),
                    ));
                }
                for _ in 0..*n {
                    md.push_str(&match_count.to_string());
                    match_count = 0;
                    md.push(reference[ri] as char);
                    qi += 1;
                    ri += 1;
                }
            }
            CigarOp::AlnMatch(n) => {
                if qi + n > query.len() || ri + n > reference.len() {
                    return Err(CyaneaError::InvalidInput(
                        "MD tag: sequences too short for M op".into(),
                    ));
                }
                for _ in 0..*n {
                    if query[qi] == reference[ri] {
                        match_count += 1;
                    } else {
                        md.push_str(&match_count.to_string());
                        match_count = 0;
                        md.push(reference[ri] as char);
                    }
                    qi += 1;
                    ri += 1;
                }
            }
            CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                if ri + n > reference.len() {
                    return Err(CyaneaError::InvalidInput(
                        "MD tag: reference too short for deletion".into(),
                    ));
                }
                md.push_str(&match_count.to_string());
                match_count = 0;
                md.push('^');
                for _ in 0..*n {
                    md.push(reference[ri] as char);
                    ri += 1;
                }
            }
            CigarOp::Insertion(n) => {
                qi += n;
            }
            CigarOp::SoftClip(n) => {
                qi += n;
            }
            CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
        }
    }

    // Emit trailing match count
    md.push_str(&match_count.to_string());

    Ok(md)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- Parsing --

    #[test]
    fn parse_simple() {
        let ops = parse_cigar("10M").unwrap();
        assert_eq!(ops, vec![CigarOp::AlnMatch(10)]);
    }

    #[test]
    fn parse_complex() {
        let ops = parse_cigar("5=2X3I4D1N2S3H1P").unwrap();
        assert_eq!(
            ops,
            vec![
                CigarOp::Match(5),
                CigarOp::Mismatch(2),
                CigarOp::Insertion(3),
                CigarOp::Deletion(4),
                CigarOp::Skip(1),
                CigarOp::SoftClip(2),
                CigarOp::HardClip(3),
                CigarOp::Padding(1),
            ]
        );
    }

    #[test]
    fn parse_all_nine_ops() {
        let ops = parse_cigar("1M1I1D1N1S1H1P1=1X").unwrap();
        assert_eq!(ops.len(), 9);
        let codes: String = ops.iter().map(|op| op.code()).collect();
        assert_eq!(codes, "MIDNSHP=X");
    }

    #[test]
    fn parse_star() {
        assert_eq!(parse_cigar("*").unwrap(), vec![]);
    }

    #[test]
    fn parse_empty_is_error() {
        assert!(parse_cigar("").is_err());
    }

    #[test]
    fn parse_bad_char() {
        assert!(parse_cigar("10Z").is_err());
    }

    #[test]
    fn parse_missing_count() {
        assert!(parse_cigar("M").is_err());
    }

    #[test]
    fn parse_trailing_digits() {
        assert!(parse_cigar("10M5").is_err());
    }

    #[test]
    fn format_round_trip() {
        let original = "3=1X2I1D4M2S5H100N1P";
        let ops = parse_cigar(original).unwrap();
        assert_eq!(cigar_string(&ops), original);
    }

    #[test]
    fn format_empty() {
        assert_eq!(cigar_string(&[]), "*");
    }

    // -- Validation --

    #[test]
    fn validate_valid() {
        let ops = parse_cigar("5H3S10M2I3M1D5M3S5H").unwrap();
        assert!(validate_cigar(&ops).is_ok());
    }

    #[test]
    fn validate_empty() {
        assert!(validate_cigar(&[]).is_ok());
    }

    #[test]
    fn validate_zero_length() {
        let ops = vec![CigarOp::AlnMatch(0)];
        assert!(validate_cigar(&ops).is_err());
    }

    #[test]
    fn validate_adjacent_duplicates() {
        let ops = vec![CigarOp::AlnMatch(3), CigarOp::AlnMatch(2)];
        assert!(validate_cigar(&ops).is_err());
    }

    #[test]
    fn validate_h_not_at_ends() {
        let ops = vec![
            CigarOp::AlnMatch(5),
            CigarOp::HardClip(3),
            CigarOp::AlnMatch(5),
        ];
        assert!(validate_cigar(&ops).is_err());
    }

    #[test]
    fn validate_s_at_invalid_position() {
        let ops = vec![
            CigarOp::AlnMatch(5),
            CigarOp::SoftClip(3),
            CigarOp::AlnMatch(5),
        ];
        assert!(validate_cigar(&ops).is_err());
    }

    #[test]
    fn validate_s_adjacent_to_h() {
        let ops = vec![
            CigarOp::HardClip(2),
            CigarOp::SoftClip(3),
            CigarOp::AlnMatch(10),
            CigarOp::SoftClip(2),
            CigarOp::HardClip(1),
        ];
        assert!(validate_cigar(&ops).is_ok());
    }

    // -- Coordinates --

    #[test]
    fn ref_query_consumed() {
        let ops = parse_cigar("3=2I4D1X2S3H").unwrap();
        // ref: = (3) + D (4) + X (1) = 8
        assert_eq!(reference_consumed(&ops), 8);
        // query: = (3) + I (2) + X (1) + S (2) = 8
        assert_eq!(query_consumed(&ops), 8);
    }

    #[test]
    fn alignment_columns_excludes_clips() {
        let ops = parse_cigar("3S10M2I3D100N5S").unwrap();
        // M (10) + I (2) + D (3) = 15 (N excluded)
        assert_eq!(alignment_columns(&ops), 15);
    }

    #[test]
    fn clipped_bases_totals() {
        let ops = parse_cigar("5H3S10M2S4H").unwrap();
        assert_eq!(clipped_bases(&ops), (5, 9));
    }

    // -- Statistics --

    #[test]
    fn identity_with_eq_x() {
        let ops = parse_cigar("8=2X").unwrap();
        assert!((identity(&ops) - 0.8).abs() < 1e-10);
    }

    #[test]
    fn identity_empty() {
        assert!((identity(&[]) - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn identity_m_only_is_zero() {
        let ops = parse_cigar("10M").unwrap();
        // M ops cannot distinguish = from X, so identity returns 0
        assert!((identity(&ops) - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn coverage_basic() {
        let ops = parse_cigar("50M").unwrap();
        assert!((coverage(&ops, 100) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn coverage_zero_ref() {
        assert!((coverage(&[], 0) - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn gap_count_and_bases() {
        let ops = parse_cigar("5=3I2=4D1=1I").unwrap();
        assert_eq!(gap_count(&ops), 3); // 3I, 4D, 1I
        assert_eq!(gap_bases(&ops), 8); // 3 + 4 + 1
    }

    // -- Arithmetic --

    #[test]
    fn merge_adjacent_combines() {
        let ops = vec![
            CigarOp::AlnMatch(3),
            CigarOp::AlnMatch(2),
            CigarOp::Insertion(1),
            CigarOp::Insertion(4),
        ];
        let merged = merge_adjacent(&ops);
        assert_eq!(merged, vec![CigarOp::AlnMatch(5), CigarOp::Insertion(5)]);
    }

    #[test]
    fn merge_adjacent_idempotent() {
        let ops = parse_cigar("10M3I4D").unwrap();
        assert_eq!(merge_adjacent(&ops), ops);
    }

    #[test]
    fn merge_adjacent_empty() {
        assert_eq!(merge_adjacent(&[]), vec![]);
    }

    #[test]
    fn reverse_cigar_order() {
        let ops = parse_cigar("3=2I1D").unwrap();
        let rev = reverse_cigar(&ops);
        assert_eq!(
            rev,
            vec![
                CigarOp::Deletion(1),
                CigarOp::Insertion(2),
                CigarOp::Match(3),
            ]
        );
    }

    #[test]
    fn split_at_reference_basic() {
        let ops = parse_cigar("10M").unwrap();
        let (left, right) = split_at_reference(&ops, 4);
        assert_eq!(left, vec![CigarOp::AlnMatch(4)]);
        assert_eq!(right, vec![CigarOp::AlnMatch(6)]);
    }

    #[test]
    fn split_at_reference_boundary() {
        let ops = parse_cigar("5=3I5=").unwrap();
        // Split at ref pos 5: first = consumes 5 ref bases → exact boundary
        // I doesn't consume ref, so it goes to right side
        let (left, right) = split_at_reference(&ops, 5);
        assert_eq!(left, vec![CigarOp::Match(5)]);
        assert_eq!(right, vec![CigarOp::Insertion(3), CigarOp::Match(5)]);
    }

    #[test]
    fn split_at_reference_zero() {
        let ops = parse_cigar("5M").unwrap();
        let (left, right) = split_at_reference(&ops, 0);
        assert!(left.is_empty());
        assert_eq!(right, vec![CigarOp::AlnMatch(5)]);
    }

    #[test]
    fn hard_clip_to_soft_conversion() {
        let ops = parse_cigar("5H10M3H").unwrap();
        let converted = hard_clip_to_soft(&ops);
        assert_eq!(
            converted,
            vec![
                CigarOp::SoftClip(5),
                CigarOp::AlnMatch(10),
                CigarOp::SoftClip(3),
            ]
        );
    }

    #[test]
    fn collapse_matches_eq_x_to_m() {
        let ops = parse_cigar("3=2X1=").unwrap();
        let collapsed = collapse_matches(&ops);
        assert_eq!(collapsed, vec![CigarOp::AlnMatch(6)]);
    }

    // -- Conversion --

    #[test]
    fn cigar_to_alignment_basic() {
        let ops = parse_cigar("3=1I2=1D1=").unwrap();
        let query = b"ACGTACG"; // 3 + 1 + 2 + 1 = 7
        let target = b"ACGACGA"; // 3 + 2 + 1 + 1 = 7
        let (aq, at) = cigar_to_alignment(&ops, query, target).unwrap();
        assert_eq!(aq, b"ACGTAC-G");
        assert_eq!(at, b"ACG-ACGA");
    }

    #[test]
    fn cigar_to_alignment_with_clips() {
        // Soft clip at start, query has extra bases
        let ops = vec![CigarOp::SoftClip(2), CigarOp::AlnMatch(3)];
        let query = b"NNACG";
        let target = b"ACG";
        let (aq, at) = cigar_to_alignment(&ops, query, target).unwrap();
        assert_eq!(aq, b"ACG");
        assert_eq!(at, b"ACG");
    }

    #[test]
    fn alignment_to_cigar_basic() {
        let query  = b"ACGTAC-G";
        let target = b"ACG-ACGA";
        let ops = alignment_to_cigar(query, target).unwrap();
        assert_eq!(
            ops,
            vec![
                CigarOp::Match(3),
                CigarOp::Insertion(1),
                CigarOp::Match(2),
                CigarOp::Deletion(1),
                CigarOp::Mismatch(1), // G vs A
            ]
        );
    }

    #[test]
    fn alignment_to_cigar_length_mismatch() {
        assert!(alignment_to_cigar(b"ACG", b"AC").is_err());
    }

    #[test]
    fn cigar_alignment_round_trip() {
        // Start with gapped alignment, extract CIGAR, reconstruct
        let query  = b"AC-GT";
        let target = b"ACAGT";
        let ops = alignment_to_cigar(query, target).unwrap();
        // ungapped sequences
        let uq = b"ACGT";
        let ut = b"ACAGT";
        let (aq, at) = cigar_to_alignment(&ops, uq, ut).unwrap();
        assert_eq!(&aq, query);
        assert_eq!(&at, target);
    }

    // -- MD tag --

    #[test]
    fn md_tag_perfect_match() {
        let ops = parse_cigar("10=").unwrap();
        let seq = b"ACGTACGTAC";
        let md = generate_md_tag(&ops, seq, seq).unwrap();
        assert_eq!(md, "10");
    }

    #[test]
    fn md_tag_with_mismatches() {
        let ops = parse_cigar("3=1X2=").unwrap();
        let query = b"ACGAAC";
        let refer = b"ACGTAC";
        let md = generate_md_tag(&ops, query, refer).unwrap();
        assert_eq!(md, "3T2");
    }

    #[test]
    fn md_tag_with_deletion() {
        let ops = parse_cigar("3=2D3=").unwrap();
        let query = b"ACGAAC";
        let refer = b"ACGTTAAC";
        let md = generate_md_tag(&ops, query, refer).unwrap();
        assert_eq!(md, "3^TT3");
    }

    #[test]
    fn md_tag_with_insertion() {
        // Insertions are skipped in MD tag
        let ops = parse_cigar("3=2I3=").unwrap();
        let query = b"ACGTTACG";
        let refer = b"ACGACG";
        let md = generate_md_tag(&ops, query, refer).unwrap();
        assert_eq!(md, "6");
    }

    #[test]
    fn md_tag_with_m_ops() {
        // M ops compare query vs reference
        let ops = parse_cigar("3M1I2M").unwrap();
        let query = b"ACGTAC";
        let refer = b"ACGAC";
        let md = generate_md_tag(&ops, query, refer).unwrap();
        assert_eq!(md, "5");
    }

    // -- Edge cases --

    #[test]
    fn single_op() {
        let ops = parse_cigar("1M").unwrap();
        assert_eq!(reference_consumed(&ops), 1);
        assert_eq!(query_consumed(&ops), 1);
        assert_eq!(alignment_columns(&ops), 1);
    }

    #[test]
    fn all_clips() {
        let ops = parse_cigar("5H3S").unwrap();
        assert_eq!(reference_consumed(&ops), 0);
        assert_eq!(query_consumed(&ops), 3);
        assert_eq!(alignment_columns(&ops), 0);
        assert_eq!(clipped_bases(&ops), (3, 5));
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    fn arb_cigar_op() -> impl Strategy<Value = CigarOp> {
        (1..100usize, 0..9u8).prop_map(|(n, kind)| match kind {
            0 => CigarOp::Match(n),
            1 => CigarOp::Mismatch(n),
            2 => CigarOp::Insertion(n),
            3 => CigarOp::Deletion(n),
            4 => CigarOp::AlnMatch(n),
            5 => CigarOp::Skip(n),
            6 => CigarOp::SoftClip(n),
            7 => CigarOp::HardClip(n),
            _ => CigarOp::Padding(n),
        })
    }

    proptest! {
        #[test]
        fn parse_format_round_trip(ops in proptest::collection::vec(arb_cigar_op(), 1..20)) {
            let merged = merge_adjacent(&ops);
            let s = cigar_string(&merged);
            let parsed = parse_cigar(&s).unwrap();
            prop_assert_eq!(parsed, merged);
        }

        #[test]
        fn ref_consumed_consistent(ops in proptest::collection::vec(arb_cigar_op(), 0..20)) {
            let total: usize = ops.iter()
                .filter(|op| op.consumes_reference())
                .map(|op| op.len())
                .sum();
            prop_assert_eq!(reference_consumed(&ops), total);
        }

        #[test]
        fn query_consumed_consistent(ops in proptest::collection::vec(arb_cigar_op(), 0..20)) {
            let total: usize = ops.iter()
                .filter(|op| op.consumes_query())
                .map(|op| op.len())
                .sum();
            prop_assert_eq!(query_consumed(&ops), total);
        }

        #[test]
        fn merge_preserves_lengths(ops in proptest::collection::vec(arb_cigar_op(), 0..20)) {
            let merged = merge_adjacent(&ops);
            prop_assert_eq!(reference_consumed(&ops), reference_consumed(&merged));
            prop_assert_eq!(query_consumed(&ops), query_consumed(&merged));
        }
    }
}
