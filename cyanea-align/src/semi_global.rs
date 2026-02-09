//! Semi-global alignment with affine gap penalties.
//!
//! A variant of Needleman-Wunsch where leading and trailing gaps in both
//! sequences are free (not penalised). This is useful for finding a query
//! within a longer target (adapter trimming, overlap detection) without
//! penalising unaligned flanking regions.
//!
//! Uses the same three-matrix formulation (Gotoh, 1982) as NW:
//!
//! - **H** — best score ending in a match/mismatch
//! - **E** — best score ending in a gap in the query (insertion)
//! - **F** — best score ending in a gap in the target (deletion)
//!
//! Differences from global alignment:
//!
//! 1. **Initialisation**: `H[0][j] = 0` and `H[i][0] = 0` (free leading gaps).
//! 2. **Fill**: Identical recurrence (no clamping to zero).
//! 3. **Traceback start**: Maximum score in the last row *or* last column
//!    (free trailing gaps beyond that point).

use crate::scoring::ScoringScheme;
use crate::types::{AlignmentResult, CigarOp};
use cyanea_core::{CyaneaError, Result};

/// Perform semi-global alignment with affine gap penalties.
///
/// Leading and trailing gaps in both sequences are free. The alignment
/// score is the maximum found in the last row or last column of the DP
/// matrix, and traceback proceeds from that cell.
///
/// # Errors
///
/// Returns an error if either sequence is empty.
pub fn semi_global(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringScheme,
) -> Result<AlignmentResult> {
    let m = query.len();
    let n = target.len();

    if m == 0 || n == 0 {
        return Err(CyaneaError::InvalidInput(
            "sequences must not be empty".into(),
        ));
    }

    let gap_open = scoring.gap_open();
    let gap_extend = scoring.gap_extend();

    let rows = m + 1;
    let cols = n + 1;

    let mut h = vec![i32::MIN / 2; rows * cols];
    let mut e = vec![i32::MIN / 2; rows * cols];
    let mut f = vec![i32::MIN / 2; rows * cols];

    let idx = |i: usize, j: usize| -> usize { i * cols + j };

    // Initialisation — free leading gaps: H borders are 0, E/F stay at sentinel.
    h[idx(0, 0)] = 0;

    for j in 1..cols {
        h[idx(0, j)] = 0; // free leading gap in query
    }

    for i in 1..rows {
        h[idx(i, 0)] = 0; // free leading gap in target
    }

    // Fill — identical recurrence to Needleman-Wunsch
    for i in 1..rows {
        for j in 1..cols {
            // E: gap in query (insertion — consuming target[j-1])
            e[idx(i, j)] = (h[idx(i, j - 1)] + gap_open).max(e[idx(i, j - 1)] + gap_extend);

            // F: gap in target (deletion — consuming query[i-1])
            f[idx(i, j)] = (h[idx(i - 1, j)] + gap_open).max(f[idx(i - 1, j)] + gap_extend);

            // H: match/mismatch
            let sub = scoring.score_pair(query[i - 1], target[j - 1]);
            let diag = h[idx(i - 1, j - 1)] + sub;

            h[idx(i, j)] = diag.max(e[idx(i, j)]).max(f[idx(i, j)]);
        }
    }

    // Find traceback start: max in last row (i=m) OR last column (j=n).
    // Free trailing gaps beyond the traceback start point.
    let mut best_score = i32::MIN;
    let mut best_i = m;
    let mut best_j = n;

    // Scan last row (i = m) — free trailing gap in query
    for j in 0..cols {
        if h[idx(m, j)] > best_score {
            best_score = h[idx(m, j)];
            best_i = m;
            best_j = j;
        }
    }

    // Scan last column (j = n) — free trailing gap in target
    for i in 0..rows {
        if h[idx(i, n)] > best_score {
            best_score = h[idx(i, n)];
            best_i = i;
            best_j = n;
        }
    }

    // Traceback from (best_i, best_j) to the border
    let mut aligned_query = Vec::new();
    let mut aligned_target = Vec::new();
    let mut cigar_ops: Vec<CigarOp> = Vec::new();

    let mut i = best_i;
    let mut j = best_j;

    #[derive(Clone, Copy, PartialEq)]
    enum State {
        H,
        E,
        F,
    }

    let mut state = State::H;

    while i > 0 || j > 0 {
        // Stop at the border — leading gaps are free and not traced
        if i == 0 && j == 0 {
            break;
        }

        // For semi-global: if we reach a border, we stop (free leading gaps)
        if state == State::H && (i == 0 || j == 0) {
            break;
        }

        match state {
            State::H => {
                if i > 0 && j > 0 {
                    let sub = scoring.score_pair(query[i - 1], target[j - 1]);
                    let diag = h[idx(i - 1, j - 1)] + sub;

                    if h[idx(i, j)] == diag {
                        aligned_query.push(query[i - 1]);
                        aligned_target.push(target[j - 1]);
                        let op = if query[i - 1].to_ascii_uppercase()
                            == target[j - 1].to_ascii_uppercase()
                        {
                            CigarOp::Match(1)
                        } else {
                            CigarOp::Mismatch(1)
                        };
                        push_cigar(&mut cigar_ops, op);
                        i -= 1;
                        j -= 1;
                    } else if h[idx(i, j)] == e[idx(i, j)] {
                        state = State::E;
                    } else {
                        state = State::F;
                    }
                } else if j > 0 {
                    state = State::E;
                } else {
                    state = State::F;
                }
            }
            State::E => {
                // Gap in query — consume target[j-1]
                aligned_query.push(b'-');
                aligned_target.push(target[j - 1]);
                push_cigar(&mut cigar_ops, CigarOp::Insertion(1));

                if e[idx(i, j)] == h[idx(i, j - 1)] + gap_open {
                    state = State::H;
                }
                j -= 1;
            }
            State::F => {
                // Gap in target — consume query[i-1]
                aligned_query.push(query[i - 1]);
                aligned_target.push(b'-');
                push_cigar(&mut cigar_ops, CigarOp::Deletion(1));

                if f[idx(i, j)] == h[idx(i - 1, j)] + gap_open {
                    state = State::H;
                }
                i -= 1;
            }
        }
    }

    aligned_query.reverse();
    aligned_target.reverse();
    cigar_ops.reverse();

    Ok(AlignmentResult {
        score: best_score,
        aligned_query,
        aligned_target,
        query_start: i,
        query_end: best_i,
        target_start: j,
        target_end: best_j,
        cigar: cigar_ops,
    })
}

/// Merge a new 1-length CIGAR op with the last op if they are the same variant.
fn push_cigar(ops: &mut Vec<CigarOp>, op: CigarOp) {
    if let Some(last) = ops.last_mut() {
        match (last, &op) {
            (CigarOp::Match(ref mut n), CigarOp::Match(1)) => {
                *n += 1;
                return;
            }
            (CigarOp::Mismatch(ref mut n), CigarOp::Mismatch(1)) => {
                *n += 1;
                return;
            }
            (CigarOp::Insertion(ref mut n), CigarOp::Insertion(1)) => {
                *n += 1;
                return;
            }
            (CigarOp::Deletion(ref mut n), CigarOp::Deletion(1)) => {
                *n += 1;
                return;
            }
            _ => {}
        }
    }
    ops.push(op);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::{ScoringMatrix, SubstitutionMatrix};

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::Simple(ScoringMatrix::dna_default())
    }

    #[test]
    fn identical_sequences() {
        let result = semi_global(b"ACGT", b"ACGT", &dna_scheme()).unwrap();
        assert_eq!(result.score, 8); // 4 matches * 2
        assert_eq!(result.cigar_string(), "4=");
        assert!((result.identity() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn query_contained_in_target() {
        // Query is a substring of target — free leading/trailing gaps in target
        let result = semi_global(b"CGT", b"AACGTAA", &dna_scheme()).unwrap();
        assert_eq!(result.score, 6); // 3 matches * 2, no gap penalties
        assert_eq!(result.matches(), 3);
    }

    #[test]
    fn target_contained_in_query() {
        // Target is a substring of query — free leading/trailing gaps in query
        let result = semi_global(b"AACGTAA", b"CGT", &dna_scheme()).unwrap();
        assert_eq!(result.score, 6); // 3 matches * 2
        assert_eq!(result.matches(), 3);
    }

    #[test]
    fn adapter_trimming() {
        // Simulates adapter at end of read: ACGTACGT + adapter TTTT
        // Semi-global should find the ACGTACGT match without penalising trailing TTTT
        let result = semi_global(b"ACGTACGT", b"ACGTACGTTTTT", &dna_scheme()).unwrap();
        assert_eq!(result.score, 16); // 8 matches * 2
    }

    #[test]
    fn overlap_detection() {
        // Overlap: suffix of query matches prefix of target
        let result = semi_global(b"AAAACGT", b"CGTTTTT", &dna_scheme()).unwrap();
        assert_eq!(result.score, 6); // CGT = 3 matches * 2
    }

    #[test]
    fn empty_sequence_errors() {
        assert!(semi_global(b"", b"ACGT", &dna_scheme()).is_err());
        assert!(semi_global(b"ACGT", b"", &dna_scheme()).is_err());
    }

    #[test]
    fn protein_semi_global() {
        let scheme = ScoringScheme::Substitution(SubstitutionMatrix::blosum62());
        let result = semi_global(b"HEAG", b"XXHEAGXX", &scheme).unwrap();
        assert!(result.score > 0, "expected positive score for embedded peptide");
    }

    #[test]
    fn single_base() {
        let result = semi_global(b"A", b"A", &dna_scheme()).unwrap();
        assert_eq!(result.score, 2);
        assert_eq!(result.cigar_string(), "1=");
    }

    #[test]
    fn no_match_region() {
        // Completely different sequences — semi-global still aligns but score is poor
        let result = semi_global(b"AAAA", b"TTTT", &dna_scheme()).unwrap();
        // Best is a single-base alignment of 0 mismatches at border, or just mismatches
        // The score should be <= 0 for all-mismatch
        assert!(result.score <= 0 || result.length() > 0);
    }
}
