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
pub mod simd_sw;
pub mod msa;
pub mod poa;
pub mod minimizers;
pub mod seed_extend;
pub mod wfa;
pub mod lcsk;
pub mod pair_hmm;

pub mod gpu;

pub use types::{AlignmentMode, AlignmentResult, CigarOp};
pub use scoring::{ScoringMatrix, ScoringScheme, SubstitutionMatrix};
pub use needleman_wunsch::needleman_wunsch;
pub use smith_waterman::smith_waterman;
pub use semi_global::semi_global;
pub use batch::align_batch;
pub use simd_sw::{sw_simd_score, sw_scalar_score};
pub use minimizers::{minimizers, find_seed_matches, Minimizer};
pub use seed_extend::{seed_extend_align, chain_seeds, Seed, SeedChain};
pub use wfa::wfa_align;
pub use lcsk::{sparse_align, find_kmer_matches, lcsk_plusplus, SparseAlignResult};
pub use poa::{PoaGraph, PoaScoring};
pub use pair_hmm::{
    pair_hmm_forward, pair_hmm_viterbi, PairHmmAlignment, PairHmmParams, PairHmmState,
};

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

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    fn dna_seq(max_len: usize) -> impl Strategy<Value = Vec<u8>> {
        proptest::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..=max_len)
    }

    proptest! {
        #[test]
        fn alignment_score_is_deterministic(
            q in dna_seq(50),
            t in dna_seq(50),
        ) {
            let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
            let r1 = align(&q, &t, AlignmentMode::Global, &scoring).unwrap();
            let r2 = align(&q, &t, AlignmentMode::Global, &scoring).unwrap();
            prop_assert_eq!(r1.score, r2.score);
        }

        #[test]
        fn identity_in_unit_interval(
            q in dna_seq(50),
            t in dna_seq(50),
        ) {
            let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
            let result = align(&q, &t, AlignmentMode::Global, &scoring).unwrap();
            let id = result.identity();
            prop_assert!(id >= 0.0 && id <= 1.0, "identity={} out of [0,1]", id);
        }

        #[test]
        fn identical_sequences_perfect_identity(seq in dna_seq(50)) {
            let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
            let result = align(&seq, &seq, AlignmentMode::Global, &scoring).unwrap();
            prop_assert!((result.identity() - 1.0).abs() < 1e-10,
                "identical seqs should have identity=1.0, got {}", result.identity());
        }

        #[test]
        fn local_score_nonnegative(
            q in dna_seq(50),
            t in dna_seq(50),
        ) {
            let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
            let result = align(&q, &t, AlignmentMode::Local, &scoring).unwrap();
            prop_assert!(result.score >= 0, "SW score should be >= 0, got {}", result.score);
        }
    }
}
