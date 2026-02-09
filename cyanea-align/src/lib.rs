//! Sequence alignment algorithms for the Cyanea bioinformatics ecosystem.
//!
//! Provides pairwise alignment via Smith-Waterman (local) and Needleman-Wunsch
//! (global) with affine gap penalties, supporting both nucleotide and protein
//! scoring schemes (BLOSUM, PAM, custom).
//!
//! # Quick start
//!
//! ```
//! use cyanea_align::{align, AlignmentMode, ScoringMatrix, ScoringScheme};
//!
//! let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
//! let result = align(b"ACGT", b"ACGT", AlignmentMode::Global, &scoring).unwrap();
//! assert_eq!(result.score, 8);
//! ```

pub mod types;
pub mod scoring;
pub mod needleman_wunsch;
pub mod smith_waterman;
pub mod semi_global;
pub mod batch;
pub mod simd;
pub mod msa;

pub mod gpu;

pub use types::{AlignmentMode, AlignmentResult, CigarOp};
pub use scoring::{ScoringMatrix, ScoringScheme, SubstitutionMatrix};
pub use needleman_wunsch::needleman_wunsch;
pub use smith_waterman::smith_waterman;
pub use semi_global::semi_global;
pub use batch::align_batch;

/// Convenience function: align two sequences using the specified mode and scoring.
///
/// Dispatches to [`smith_waterman`] for [`AlignmentMode::Local`],
/// [`needleman_wunsch`] for [`AlignmentMode::Global`], or
/// [`semi_global`] for [`AlignmentMode::SemiGlobal`].
///
/// # Errors
///
/// Returns an error if either sequence is empty.
pub fn align(
    query: &[u8],
    target: &[u8],
    mode: AlignmentMode,
    scoring: &ScoringScheme,
) -> cyanea_core::Result<AlignmentResult> {
    match mode {
        AlignmentMode::Local => crate::smith_waterman::smith_waterman(query, target, scoring),
        AlignmentMode::Global => crate::needleman_wunsch::needleman_wunsch(query, target, scoring),
        AlignmentMode::SemiGlobal => crate::semi_global::semi_global(query, target, scoring),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::Simple(ScoringMatrix::dna_default())
    }

    #[test]
    fn align_global_end_to_end() {
        let result = align(b"ACGT", b"ACGT", AlignmentMode::Global, &dna_scheme()).unwrap();
        assert_eq!(result.score, 8);
        assert_eq!(result.cigar_string(), "4=");
    }

    #[test]
    fn align_local_end_to_end() {
        let result = align(
            b"AAACGTAAA",
            b"TTTCGTTTT",
            AlignmentMode::Local,
            &dna_scheme(),
        )
        .unwrap();
        assert!(result.score > 0);
    }

    #[test]
    fn align_semi_global() {
        let result = align(b"ACGT", b"ACGT", AlignmentMode::SemiGlobal, &dna_scheme()).unwrap();
        assert_eq!(result.score, 8);
    }

    #[test]
    fn align_protein_global() {
        let scoring = ScoringScheme::Substitution(SubstitutionMatrix::blosum62());
        let result = align(b"HEAGAWGHEE", b"PAWHEAE", AlignmentMode::Global, &scoring).unwrap();
        assert!(result.score > 0);
        assert!(result.length() > 0);
    }
}
