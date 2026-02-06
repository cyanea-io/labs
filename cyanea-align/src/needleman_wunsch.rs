//! Needleman-Wunsch global alignment with affine gap penalties.
//!
//! Uses a three-matrix dynamic programming formulation (Gotoh, 1982):
//!
//! - **H** — best score ending in a match/mismatch
//! - **E** — best score ending in a gap in the query (insertion)
//! - **F** — best score ending in a gap in the target (deletion)

use crate::scoring::ScoringScheme;
use crate::types::{AlignmentResult, CigarOp};
use cyanea_core::{CyaneaError, Result};

/// Perform global (Needleman-Wunsch) alignment with affine gap penalties.
///
/// # Errors
///
/// Returns an error if either sequence is empty.
pub fn needleman_wunsch(
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

    // H[i][j]: best score for aligning query[..i] and target[..j]
    // E[i][j]: best score ending with a gap in the query (consuming target)
    // F[i][j]: best score ending with a gap in the target (consuming query)
    let mut h = vec![i32::MIN / 2; rows * cols];
    let mut e = vec![i32::MIN / 2; rows * cols];
    let mut f = vec![i32::MIN / 2; rows * cols];

    let idx = |i: usize, j: usize| -> usize { i * cols + j };

    // Initialization
    h[idx(0, 0)] = 0;

    for i in 1..rows {
        // Opening a gap in the target of length i
        h[idx(i, 0)] = gap_open + (i as i32 - 1) * gap_extend;
        f[idx(i, 0)] = h[idx(i, 0)];
    }

    for j in 1..cols {
        // Opening a gap in the query of length j
        h[idx(0, j)] = gap_open + (j as i32 - 1) * gap_extend;
        e[idx(0, j)] = h[idx(0, j)];
    }

    // Fill
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

    // Traceback from (m, n) to (0, 0)
    let mut aligned_query = Vec::new();
    let mut aligned_target = Vec::new();
    let mut cigar_ops: Vec<CigarOp> = Vec::new();

    let mut i = m;
    let mut j = n;

    // Track which matrix we're currently in for traceback
    #[derive(Clone, Copy, PartialEq)]
    enum State {
        H,
        E,
        F,
    }

    let mut state = State::H;

    while i > 0 || j > 0 {
        match state {
            State::H => {
                if i > 0 && j > 0 {
                    let sub = scoring.score_pair(query[i - 1], target[j - 1]);
                    let diag = h[idx(i - 1, j - 1)] + sub;

                    if h[idx(i, j)] == diag {
                        // Match or mismatch
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

                // Decide whether to stay in E or return to H
                if e[idx(i, j)] == h[idx(i, j - 1)] + gap_open {
                    state = State::H;
                }
                // else: stay in E (extending the gap)
                j -= 1;
            }
            State::F => {
                // Gap in target — consume query[i-1]
                aligned_query.push(query[i - 1]);
                aligned_target.push(b'-');
                push_cigar(&mut cigar_ops, CigarOp::Deletion(1));

                // Decide whether to stay in F or return to H
                if f[idx(i, j)] == h[idx(i - 1, j)] + gap_open {
                    state = State::H;
                }
                // else: stay in F (extending the gap)
                i -= 1;
            }
        }
    }

    // Reverse since we traced back from (m, n)
    aligned_query.reverse();
    aligned_target.reverse();
    cigar_ops.reverse();

    Ok(AlignmentResult {
        score: h[idx(m, n)],
        aligned_query,
        aligned_target,
        query_start: 0,
        query_end: m,
        target_start: 0,
        target_end: n,
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
        let result = needleman_wunsch(b"ACGT", b"ACGT", &dna_scheme()).unwrap();
        // 4 matches * 2 = 8
        assert_eq!(result.score, 8);
        assert_eq!(result.aligned_query, b"ACGT");
        assert_eq!(result.aligned_target, b"ACGT");
        assert_eq!(result.cigar_string(), "4=");
        assert!((result.identity() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn single_mismatch() {
        let result = needleman_wunsch(b"ACGT", b"ACAT", &dna_scheme()).unwrap();
        // Should align A=A, C=C, G!=A, T=T → 3 matches * 2 + 1 mismatch * -1 = 5
        assert_eq!(result.score, 5);
        assert_eq!(result.matches(), 3);
        assert_eq!(result.mismatches(), 1);
    }

    #[test]
    fn gap_insertion() {
        // ACGT vs ACT — should introduce one deletion in query (or insertion in target)
        let result = needleman_wunsch(b"ACGT", b"ACT", &dna_scheme()).unwrap();
        assert!(result.gaps() > 0, "expected at least one gap");
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 4);
        assert_eq!(result.target_start, 0);
        assert_eq!(result.target_end, 3);
    }

    #[test]
    fn protein_with_blosum62() {
        let scheme = ScoringScheme::Substitution(SubstitutionMatrix::blosum62());
        let result = needleman_wunsch(b"HEAGAWGHEE", b"PAWHEAE", &scheme).unwrap();
        // Just check it runs and produces a reasonable alignment
        assert!(result.score > 0, "expected positive score for related peptides");
        assert!(result.length() > 0);
    }

    #[test]
    fn empty_sequence_errors() {
        assert!(needleman_wunsch(b"", b"ACGT", &dna_scheme()).is_err());
        assert!(needleman_wunsch(b"ACGT", b"", &dna_scheme()).is_err());
    }

    #[test]
    fn single_base() {
        let result = needleman_wunsch(b"A", b"A", &dna_scheme()).unwrap();
        assert_eq!(result.score, 2);
        assert_eq!(result.cigar_string(), "1=");
    }

    #[test]
    fn completely_different() {
        let result = needleman_wunsch(b"AAAA", b"TTTT", &dna_scheme()).unwrap();
        // 4 mismatches * -1 = -4
        assert_eq!(result.score, -4);
        assert_eq!(result.mismatches(), 4);
    }
}
