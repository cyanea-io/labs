//! Smith-Waterman local alignment with affine gap penalties.
//!
//! Uses the same three-matrix formulation as Needleman-Wunsch but clamps all
//! scores to zero (no negative scores), and traceback begins at the
//! highest-scoring cell rather than at `(m, n)`.

use crate::scoring::ScoringScheme;
use crate::types::{AlignmentResult, CigarOp};
use cyanea_core::{CyaneaError, Result};

/// Perform local (Smith-Waterman) alignment with affine gap penalties.
///
/// # Errors
///
/// Returns an error if either sequence is empty.
pub fn smith_waterman(
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

    let mut h = vec![0i32; rows * cols];
    let mut e = vec![0i32; rows * cols];
    let mut f = vec![0i32; rows * cols];

    let idx = |i: usize, j: usize| -> usize { i * cols + j };

    // First row and column are already 0 (local alignment)

    let mut max_score = 0i32;
    let mut max_i = 0usize;
    let mut max_j = 0usize;

    // Fill
    for i in 1..rows {
        for j in 1..cols {
            // E: gap in query (insertion — consuming target[j-1])
            e[idx(i, j)] = (h[idx(i, j - 1)] + gap_open)
                .max(e[idx(i, j - 1)] + gap_extend)
                .max(0);

            // F: gap in target (deletion — consuming query[i-1])
            f[idx(i, j)] = (h[idx(i - 1, j)] + gap_open)
                .max(f[idx(i - 1, j)] + gap_extend)
                .max(0);

            // H: match/mismatch or reset to 0
            let sub = scoring.score_pair(query[i - 1], target[j - 1]);
            let diag = h[idx(i - 1, j - 1)] + sub;

            h[idx(i, j)] = diag.max(e[idx(i, j)]).max(f[idx(i, j)]).max(0);

            if h[idx(i, j)] > max_score {
                max_score = h[idx(i, j)];
                max_i = i;
                max_j = j;
            }
        }
    }

    // No positive-scoring region found
    if max_score == 0 {
        return Ok(AlignmentResult {
            score: 0,
            aligned_query: Vec::new(),
            aligned_target: Vec::new(),
            query_start: 0,
            query_end: 0,
            target_start: 0,
            target_end: 0,
            cigar: Vec::new(),
        });
    }

    // Traceback from (max_i, max_j) until we hit a cell with score 0
    let mut aligned_query = Vec::new();
    let mut aligned_target = Vec::new();
    let mut cigar_ops: Vec<CigarOp> = Vec::new();

    let mut i = max_i;
    let mut j = max_j;

    #[derive(Clone, Copy, PartialEq)]
    enum State {
        H,
        E,
        F,
    }

    let mut state = State::H;

    while i > 0 || j > 0 {
        // In local alignment, stop when H reaches 0 (only in H state)
        if state == State::H && h[idx(i, j)] == 0 {
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
                aligned_query.push(b'-');
                aligned_target.push(target[j - 1]);
                push_cigar(&mut cigar_ops, CigarOp::Insertion(1));

                if e[idx(i, j)] == h[idx(i, j - 1)] + gap_open {
                    state = State::H;
                }
                j -= 1;
            }
            State::F => {
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

    // Reverse since we traced back
    aligned_query.reverse();
    aligned_target.reverse();
    cigar_ops.reverse();

    Ok(AlignmentResult {
        score: max_score,
        aligned_query,
        aligned_target,
        query_start: i,
        query_end: max_i,
        target_start: j,
        target_end: max_j,
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
    use crate::scoring::ScoringMatrix;

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::Simple(ScoringMatrix::dna_default())
    }

    #[test]
    fn local_match_in_flanks() {
        // CGT is the conserved region in unrelated flanking sequence
        let result = smith_waterman(b"AAACGTAAA", b"TTTCGTTTT", &dna_scheme()).unwrap();
        assert!(result.score > 0, "expected positive score");
        // The local alignment should find CGT
        let aligned_q = String::from_utf8_lossy(&result.aligned_query);
        let aligned_t = String::from_utf8_lossy(&result.aligned_target);
        assert!(
            aligned_q.contains("CGT"),
            "expected CGT in query alignment, got: {aligned_q}"
        );
        assert!(
            aligned_t.contains("CGT"),
            "expected CGT in target alignment, got: {aligned_t}"
        );
    }

    #[test]
    fn no_match_returns_zero() {
        // Use a scoring scheme where even single matches can't overcome the mismatches
        let scheme = ScoringScheme::Simple(ScoringMatrix::new(1, -4, -10, -5).unwrap());
        let result = smith_waterman(b"AAAA", b"CCCC", &scheme).unwrap();
        assert_eq!(result.score, 0);
        assert!(result.aligned_query.is_empty());
        assert!(result.cigar.is_empty());
    }

    #[test]
    fn full_match() {
        let result = smith_waterman(b"ACGT", b"ACGT", &dna_scheme()).unwrap();
        assert_eq!(result.score, 8); // 4 * 2
        assert_eq!(result.aligned_query, b"ACGT");
        assert_eq!(result.aligned_target, b"ACGT");
    }

    #[test]
    fn query_start_end_reflect_local_region() {
        let result = smith_waterman(b"AAACGTAAA", b"TTTCGTTTT", &dna_scheme()).unwrap();
        // CGT starts at index 3 in both
        assert!(result.query_start >= 3, "query_start should be >= 3");
        assert!(result.target_start >= 3, "target_start should be >= 3");
    }

    #[test]
    fn empty_sequence_errors() {
        assert!(smith_waterman(b"", b"ACGT", &dna_scheme()).is_err());
        assert!(smith_waterman(b"ACGT", b"", &dna_scheme()).is_err());
    }

    #[test]
    fn single_base_match() {
        let result = smith_waterman(b"A", b"A", &dna_scheme()).unwrap();
        assert_eq!(result.score, 2);
        assert_eq!(result.cigar_string(), "1=");
    }

    #[test]
    fn local_alignment_ignores_poor_flanks() {
        // GATTACA with a perfect match in the middle but bad flanks
        let result = smith_waterman(b"CCGATTACACC", b"TTGATTACATT", &dna_scheme()).unwrap();
        assert!(result.score >= 14, "GATTACA = 7 bases * 2 = 14");
    }
}
