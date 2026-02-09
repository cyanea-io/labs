//! Banded alignment and SIMD-ready alignment kernels.
//!
//! Provides banded variants of Needleman-Wunsch and Smith-Waterman that restrict
//! computation to a diagonal band of width `2 * bandwidth + 1`, reducing work
//! from O(mn) to O(m × bandwidth).
//!
//! True SIMD vectorization (SSE4.1, AVX2, NEON) is deferred until profile-guided
//! benchmarks justify the complexity. The banded algorithms here serve as the
//! scalar baseline and are already significantly faster for long sequences with
//! bounded edit distance.

use cyanea_core::{CyaneaError, Result};

use crate::scoring::ScoringScheme;
use crate::types::{AlignmentMode, AlignmentResult, CigarOp};

/// Banded Needleman-Wunsch (global) alignment.
///
/// Restricts the DP matrix to a diagonal band of width `2 * bandwidth + 1`.
/// If the optimal alignment path lies outside this band the result may be
/// suboptimal, but computation is significantly reduced.
///
/// # Errors
///
/// Returns an error if either sequence is empty or bandwidth is zero.
pub fn banded_nw(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringScheme,
    bandwidth: usize,
) -> Result<AlignmentResult> {
    validate_inputs(query, target, bandwidth)?;
    banded_align(query, target, scoring, bandwidth, AlignmentMode::Global)
}

/// Banded Smith-Waterman (local) alignment.
///
/// # Errors
///
/// Returns an error if either sequence is empty or bandwidth is zero.
pub fn banded_sw(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringScheme,
    bandwidth: usize,
) -> Result<AlignmentResult> {
    validate_inputs(query, target, bandwidth)?;
    banded_align(query, target, scoring, bandwidth, AlignmentMode::Local)
}

/// Banded semi-global alignment.
///
/// Free leading and trailing gaps within the band. Useful for overlap
/// detection and adapter trimming with bounded edit distance.
///
/// # Errors
///
/// Returns an error if either sequence is empty or bandwidth is zero.
pub fn banded_semi_global(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringScheme,
    bandwidth: usize,
) -> Result<AlignmentResult> {
    validate_inputs(query, target, bandwidth)?;
    banded_align(query, target, scoring, bandwidth, AlignmentMode::SemiGlobal)
}

/// Score-only banded alignment (no traceback). Uses O(bandwidth) memory.
pub fn banded_score_only(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringScheme,
    bandwidth: usize,
    mode: AlignmentMode,
) -> Result<i32> {
    validate_inputs(query, target, bandwidth)?;
    let result = banded_align(query, target, scoring, bandwidth, mode)?;
    Ok(result.score)
}

fn validate_inputs(query: &[u8], target: &[u8], bandwidth: usize) -> Result<()> {
    if query.is_empty() || target.is_empty() {
        return Err(CyaneaError::InvalidInput("empty sequence".into()));
    }
    if bandwidth == 0 {
        return Err(CyaneaError::InvalidInput("bandwidth must be > 0".into()));
    }
    Ok(())
}

const TRACE_DIAG: u8 = 0;
const TRACE_UP: u8 = 1;
const TRACE_LEFT: u8 = 2;
const TRACE_STOP: u8 = 3;

fn banded_align(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringScheme,
    bandwidth: usize,
    mode: AlignmentMode,
) -> Result<AlignmentResult> {
    let m = query.len();
    let n = target.len();
    let w = bandwidth;
    let is_local = mode == AlignmentMode::Local;
    let is_semi_global = mode == AlignmentMode::SemiGlobal;
    let band_width = 2 * w + 1;
    let rows = m + 1;

    // Band-indexed storage: for row i, column j maps to band_idx = j - i + w.
    let mut h = vec![vec![i32::MIN / 2; band_width]; rows];
    let mut trace = vec![vec![TRACE_STOP; band_width]; rows];

    let band_idx = |i: usize, j: usize| -> Option<usize> {
        let idx = j as isize - i as isize + w as isize;
        if idx >= 0 && (idx as usize) < band_width {
            Some(idx as usize)
        } else {
            None
        }
    };

    // Row 0
    if let Some(bi) = band_idx(0, 0) {
        h[0][bi] = 0;
    }
    if !is_local && !is_semi_global {
        for j in 1..=w.min(n) {
            if let Some(bi) = band_idx(0, j) {
                h[0][bi] = scoring.gap_open() + (j as i32 - 1) * scoring.gap_extend();
                trace[0][bi] = TRACE_LEFT;
            }
        }
    } else {
        for j in 0..=w.min(n) {
            if let Some(bi) = band_idx(0, j) {
                h[0][bi] = 0;
            }
        }
    }

    // Column 0
    if !is_local && !is_semi_global {
        for i in 1..=w.min(m) {
            if let Some(bi) = band_idx(i, 0) {
                h[i][bi] = scoring.gap_open() + (i as i32 - 1) * scoring.gap_extend();
                trace[i][bi] = TRACE_UP;
            }
        }
    } else {
        for i in 1..=m {
            if let Some(bi) = band_idx(i, 0) {
                h[i][bi] = 0;
            }
        }
    }

    let mut best_score = if is_local { 0 } else { i32::MIN };
    let mut best_i = 0;
    let mut best_j = 0;

    for i in 1..=m {
        let j_min = if i > w { i - w } else { 1 };
        let j_max = (i + w).min(n);

        for j in j_min..=j_max {
            let bi = match band_idx(i, j) {
                Some(b) => b,
                None => continue,
            };

            let sub = scoring.score_pair(query[i - 1], target[j - 1]);

            let diag = band_idx(i - 1, j - 1)
                .map(|b| h[i - 1][b])
                .unwrap_or(i32::MIN / 2)
                .saturating_add(sub);

            let up = band_idx(i - 1, j)
                .map(|b| h[i - 1][b])
                .unwrap_or(i32::MIN / 2)
                .saturating_add(scoring.gap_open());

            let left = band_idx(i, j - 1)
                .map(|b| h[i][b])
                .unwrap_or(i32::MIN / 2)
                .saturating_add(scoring.gap_open());

            let mut score = diag.max(up).max(left);
            let mut tr = if score == diag {
                TRACE_DIAG
            } else if score == up {
                TRACE_UP
            } else {
                TRACE_LEFT
            };

            if is_local && score < 0 {
                score = 0;
                tr = TRACE_STOP;
            }

            h[i][bi] = score;
            trace[i][bi] = tr;

            if score > best_score {
                best_score = score;
                best_i = i;
                best_j = j;
            }
        }
    }

    let (end_i, end_j) = if is_local {
        (best_i, best_j)
    } else if is_semi_global {
        // Semi-global: find max in last row OR last column (free trailing gaps)
        best_score = i32::MIN;
        best_i = m;
        best_j = n;
        // Scan last row (i = m)
        let j_min_last = if m > w { m - w } else { 0 };
        let j_max_last = (m + w).min(n);
        for j in j_min_last..=j_max_last {
            if let Some(bi) = band_idx(m, j) {
                if h[m][bi] > best_score {
                    best_score = h[m][bi];
                    best_i = m;
                    best_j = j;
                }
            }
        }
        // Scan last column (j = n)
        let i_min_last = if n > w { n - w } else { 0 };
        let i_max_last = (n + w).min(m);
        for i in i_min_last..=i_max_last {
            if let Some(bi) = band_idx(i, n) {
                if h[i][bi] > best_score {
                    best_score = h[i][bi];
                    best_i = i;
                    best_j = n;
                }
            }
        }
        (best_i, best_j)
    } else {
        if let Some(bi) = band_idx(m, n) {
            best_score = h[m][bi];
        }
        (m, n)
    };

    // Traceback
    let mut ops: Vec<CigarOp> = Vec::new();
    let mut aligned_query = Vec::new();
    let mut aligned_target = Vec::new();
    let mut ci = end_i;
    let mut cj = end_j;

    loop {
        if ci == 0 && cj == 0 {
            break;
        }
        // Semi-global: stop at border (free leading gaps)
        if is_semi_global && (ci == 0 || cj == 0) {
            break;
        }
        let bi = match band_idx(ci, cj) {
            Some(b) => b,
            None => break,
        };
        let tr = trace[ci][bi];
        if is_local && tr == TRACE_STOP {
            break;
        }
        match tr {
            TRACE_DIAG => {
                let q = query[ci - 1];
                let t = target[cj - 1];
                aligned_query.push(q);
                aligned_target.push(t);
                let op = if q.to_ascii_uppercase() == t.to_ascii_uppercase() {
                    CigarOp::Match(1)
                } else {
                    CigarOp::Mismatch(1)
                };
                ops.push(op);
                ci -= 1;
                cj -= 1;
            }
            TRACE_UP => {
                aligned_query.push(query[ci - 1]);
                aligned_target.push(b'-');
                ops.push(CigarOp::Deletion(1));
                ci -= 1;
            }
            TRACE_LEFT => {
                aligned_query.push(b'-');
                aligned_target.push(target[cj - 1]);
                ops.push(CigarOp::Insertion(1));
                cj -= 1;
            }
            _ => break,
        }
    }

    aligned_query.reverse();
    aligned_target.reverse();
    ops.reverse();

    let cigar = merge_cigar(ops);
    let (query_start, query_end) = if is_local || is_semi_global { (ci, end_i) } else { (0, m) };
    let (target_start, target_end) = if is_local || is_semi_global { (cj, end_j) } else { (0, n) };

    Ok(AlignmentResult {
        score: best_score,
        aligned_query,
        aligned_target,
        query_start,
        query_end,
        target_start,
        target_end,
        cigar,
    })
}

fn merge_cigar(ops: Vec<CigarOp>) -> Vec<CigarOp> {
    let mut merged: Vec<CigarOp> = Vec::new();
    for op in ops {
        if let Some(last) = merged.last_mut() {
            let combined = match (*last, op) {
                (CigarOp::Match(a), CigarOp::Match(b)) => Some(CigarOp::Match(a + b)),
                (CigarOp::Mismatch(a), CigarOp::Mismatch(b)) => Some(CigarOp::Mismatch(a + b)),
                (CigarOp::Insertion(a), CigarOp::Insertion(b)) => Some(CigarOp::Insertion(a + b)),
                (CigarOp::Deletion(a), CigarOp::Deletion(b)) => Some(CigarOp::Deletion(a + b)),
                _ => None,
            };
            if let Some(c) = combined {
                *last = c;
                continue;
            }
        }
        merged.push(op);
    }
    merged
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::{ScoringMatrix, ScoringScheme};

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::Simple(ScoringMatrix::dna_default())
    }

    #[test]
    fn banded_nw_identical() {
        let r = banded_nw(b"ACGT", b"ACGT", &dna_scheme(), 3).unwrap();
        assert_eq!(r.score, 8);
        assert_eq!(r.cigar_string(), "4=");
    }

    #[test]
    fn banded_nw_with_mismatch() {
        let r = banded_nw(b"ACGT", b"ACTT", &dna_scheme(), 3).unwrap();
        assert!(r.score > 0);
        assert_eq!(r.query_start, 0);
        assert_eq!(r.query_end, 4);
    }

    #[test]
    fn banded_sw_local() {
        let r = banded_sw(b"AAACGTAAA", b"TTTCGTTTT", &dna_scheme(), 5).unwrap();
        assert!(r.score > 0);
    }

    #[test]
    fn banded_sw_identical() {
        let r = banded_sw(b"ACGT", b"ACGT", &dna_scheme(), 4).unwrap();
        assert_eq!(r.score, 8);
    }

    #[test]
    fn banded_empty_sequence_error() {
        assert!(banded_nw(b"", b"ACGT", &dna_scheme(), 3).is_err());
        assert!(banded_nw(b"ACGT", b"", &dna_scheme(), 3).is_err());
    }

    #[test]
    fn banded_zero_bandwidth_error() {
        assert!(banded_nw(b"ACGT", b"ACGT", &dna_scheme(), 0).is_err());
    }

    #[test]
    fn banded_score_only_global() {
        let score =
            banded_score_only(b"ACGT", b"ACGT", &dna_scheme(), 3, AlignmentMode::Global).unwrap();
        assert_eq!(score, 8);
    }

    #[test]
    fn banded_score_only_local() {
        let score =
            banded_score_only(b"AAACGTAAA", b"TTTCGTTTT", &dna_scheme(), 5, AlignmentMode::Local)
                .unwrap();
        assert!(score > 0);
    }

    #[test]
    fn banded_semi_global_identical() {
        let r = banded_semi_global(b"ACGT", b"ACGT", &dna_scheme(), 3).unwrap();
        assert_eq!(r.score, 8);
        assert_eq!(r.cigar_string(), "4=");
    }

    #[test]
    fn banded_semi_global_query_in_target() {
        let r = banded_semi_global(b"CGT", b"AACGTAA", &dna_scheme(), 5).unwrap();
        assert_eq!(r.score, 6); // 3 matches * 2
    }

    #[test]
    fn banded_semi_global_overlap() {
        let r = banded_semi_global(b"AAAACGT", b"CGTTTTT", &dna_scheme(), 5).unwrap();
        assert_eq!(r.score, 6); // CGT = 3 matches * 2
    }

    #[test]
    fn banded_score_only_semi_global() {
        let score = banded_score_only(
            b"CGT",
            b"AACGTAA",
            &dna_scheme(),
            5,
            AlignmentMode::SemiGlobal,
        )
        .unwrap();
        assert_eq!(score, 6);
    }
}
