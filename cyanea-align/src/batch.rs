//! Batch pairwise alignment over multiple sequence pairs.
//!
//! Provides a simple serial CPU implementation that dispatches each pair to
//! the appropriate algorithm based on [`AlignmentMode`].

use crate::needleman_wunsch::needleman_wunsch;
use crate::scoring::ScoringScheme;
use crate::semi_global::semi_global;
use crate::smith_waterman::smith_waterman;
use crate::types::{AlignmentMode, AlignmentResult};
use cyanea_core::Result;

/// Align a batch of sequence pairs using the specified mode and scoring scheme.
///
/// Each pair is aligned independently and results are returned in the same order.
///
/// # Errors
///
/// Returns an error if any individual alignment fails.
pub fn align_batch(
    pairs: &[(&[u8], &[u8])],
    mode: AlignmentMode,
    scoring: &ScoringScheme,
) -> Result<Vec<AlignmentResult>> {
    pairs
        .iter()
        .map(|(query, target)| match mode {
            AlignmentMode::Local => smith_waterman(query, target, scoring),
            AlignmentMode::Global => needleman_wunsch(query, target, scoring),
            AlignmentMode::SemiGlobal => semi_global(query, target, scoring),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::ScoringMatrix;

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::Simple(ScoringMatrix::dna_default())
    }

    #[test]
    fn batch_multiple_pairs() {
        let pairs: Vec<(&[u8], &[u8])> = vec![
            (b"ACGT", b"ACGT"),
            (b"AAAA", b"TTTT"),
            (b"ACGT", b"ACT"),
        ];
        let results = align_batch(&pairs, AlignmentMode::Global, &dna_scheme()).unwrap();
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].score, 8); // perfect match
        assert_eq!(results[1].score, -4); // all mismatches
    }

    #[test]
    fn batch_local_mode() {
        let pairs: Vec<(&[u8], &[u8])> = vec![(b"AAACGTAAA", b"TTTCGTTTT")];
        let results = align_batch(&pairs, AlignmentMode::Local, &dna_scheme()).unwrap();
        assert_eq!(results.len(), 1);
        assert!(results[0].score > 0);
    }

    #[test]
    fn empty_batch() {
        let pairs: Vec<(&[u8], &[u8])> = vec![];
        let results = align_batch(&pairs, AlignmentMode::Global, &dna_scheme()).unwrap();
        assert!(results.is_empty());
    }

    #[test]
    fn batch_semi_global() {
        let pairs: Vec<(&[u8], &[u8])> = vec![(b"ACGT", b"ACGT"), (b"CGT", b"AACGTAA")];
        let results = align_batch(&pairs, AlignmentMode::SemiGlobal, &dna_scheme()).unwrap();
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].score, 8);
        assert_eq!(results[1].score, 6);
    }
}
