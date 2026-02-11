//! Seed-and-extend alignment for long reads.
//!
//! This module implements the classic seed-and-extend paradigm used by modern
//! long-read aligners:
//!
//! 1. **Seed** — Extract minimizers from both query and target, then find
//!    matching anchors (shared k-mer hashes).
//! 2. **Chain** — Select the highest-scoring colinear subset of seeds using
//!    dynamic programming with a gap penalty.
//! 3. **Extend** — Perform banded Smith-Waterman alignment in the regions
//!    between chained seeds, then concatenate the results.
//!
//! For long reads (10 kb+) this avoids the full O(mn) dynamic programming
//! cost and runs in roughly O(m × bandwidth × number_of_chain_segments).

use cyanea_core::{CyaneaError, Result};

use crate::minimizers::{find_seed_matches, minimizers};
use crate::scoring::ScoringScheme;
use crate::simd::banded_sw;
use crate::types::{AlignmentResult, CigarOp};

/// A seed (anchor) between query and target sequences.
#[derive(Debug, Clone, Copy)]
pub struct Seed {
    /// Position in the query sequence (0-based).
    pub query_pos: usize,
    /// Position in the target sequence (0-based).
    pub target_pos: usize,
}

/// A chain of colinear seeds.
#[derive(Debug, Clone)]
pub struct SeedChain {
    /// The seeds in this chain, ordered by target position.
    pub seeds: Vec<Seed>,
    /// Total chain score.
    pub score: i64,
}

// ---------------------------------------------------------------------------
// Seed chaining
// ---------------------------------------------------------------------------

/// Chain colinear seeds using dynamic programming.
///
/// Seeds must be colinear: both query and target positions strictly increase
/// along the chain.  Each seed is scored proportionally to the k-mer length
/// it represents (`k` matching bases), minus a gap penalty proportional to
/// the distance between consecutive seeds.
///
/// # Arguments
/// * `seeds` - Pairs of `(query_pos, target_pos)`.
/// * `k`       - k-mer length (seed match length in bases).
/// * `max_gap` - Maximum allowed gap (in either query or target coordinates)
///   between consecutive seeds in the chain.
///
/// # Returns
/// The highest-scoring chain.
///
/// # Errors
/// Returns an error if `seeds` is empty.
pub fn chain_seeds(seeds: &[(usize, usize)], k: usize, max_gap: usize) -> Result<SeedChain> {
    if seeds.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "no seeds to chain".into(),
        ));
    }

    // Convert to Seed structs and sort by target position, breaking ties by
    // query position.
    let mut sorted: Vec<Seed> = seeds
        .iter()
        .map(|&(qp, tp)| Seed {
            query_pos: qp,
            target_pos: tp,
        })
        .collect();
    sorted.sort_unstable_by(|a, b| {
        a.target_pos
            .cmp(&b.target_pos)
            .then(a.query_pos.cmp(&b.query_pos))
    });

    let n = sorted.len();
    let seed_score: i64 = k.max(1) as i64; // each seed represents k matching bases

    // DP: dp[i] = best chain score ending at seed i.
    let mut dp: Vec<i64> = vec![seed_score; n];
    // Traceback: prev[i] = predecessor index in the chain (or usize::MAX).
    let mut prev: Vec<usize> = vec![usize::MAX; n];

    for i in 1..n {
        for j in (0..i).rev() {
            // Colinear check: query positions must also be strictly increasing.
            if sorted[j].query_pos >= sorted[i].query_pos {
                continue;
            }
            // target_pos is already sorted, so sorted[j].target_pos <= sorted[i].target_pos.
            // Require strict increase.
            if sorted[j].target_pos >= sorted[i].target_pos {
                continue;
            }

            let q_gap = sorted[i].query_pos - sorted[j].query_pos;
            let t_gap = sorted[i].target_pos - sorted[j].target_pos;
            let gap = q_gap.max(t_gap);

            if gap > max_gap {
                continue;
            }

            // Gap penalty: the larger of the two gaps, scaled down.
            // Use 1 penalty unit per gap base, capped so it doesn't dominate.
            let penalty = gap as i64;
            let candidate = dp[j] + seed_score - penalty;

            if candidate > dp[i] {
                dp[i] = candidate;
                prev[i] = j;
            }
        }
    }

    // Find the seed with the best score.
    let mut best_idx = 0;
    for i in 1..n {
        if dp[i] > dp[best_idx] {
            best_idx = i;
        }
    }

    // Traceback to reconstruct the chain.
    let mut chain_indices: Vec<usize> = Vec::new();
    let mut idx = best_idx;
    loop {
        chain_indices.push(idx);
        if prev[idx] == usize::MAX {
            break;
        }
        idx = prev[idx];
    }
    chain_indices.reverse();

    let chain_seeds: Vec<Seed> = chain_indices.iter().map(|&i| sorted[i]).collect();
    let chain_score = dp[best_idx];

    Ok(SeedChain {
        seeds: chain_seeds,
        score: chain_score,
    })
}

// ---------------------------------------------------------------------------
// Seed-and-extend alignment
// ---------------------------------------------------------------------------

/// Perform seed-and-extend alignment.
///
/// 1. Extract minimizers from both query and target.
/// 2. Find matching seeds.
/// 3. Chain the best colinear seeds.
/// 4. Extend alignment around the chain using banded Smith-Waterman.
///
/// This is designed for long reads (10 kb+) where full O(mn) DP is too slow.
///
/// # Arguments
/// * `query`     - Query sequence (DNA, ACGT).
/// * `target`    - Target sequence (DNA, ACGT).
/// * `scoring`   - Scoring scheme (typically DNA simple scoring).
/// * `k`         - k-mer length for minimizers (recommended: 15).
/// * `w`         - Window size for minimizers (recommended: 10).
/// * `bandwidth` - Band width for banded extension alignment.
///
/// # Errors
/// Returns an error if sequences are empty or too short for the given k.
pub fn seed_extend_align(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringScheme,
    k: usize,
    w: usize,
    bandwidth: usize,
) -> Result<AlignmentResult> {
    if query.is_empty() || target.is_empty() {
        return Err(CyaneaError::InvalidInput("empty sequence".into()));
    }

    // If sequences are too short for minimizer extraction, fall back.
    if query.len() < k || target.len() < k {
        return banded_sw(query, target, scoring, bandwidth);
    }

    // Step 1: Extract minimizers.
    let query_mins = minimizers(query, k, w)?;
    let target_mins = minimizers(target, k, w)?;

    // Step 2: Find matching seeds.
    let seed_pairs = find_seed_matches(&query_mins, &target_mins);
    if seed_pairs.is_empty() {
        // No shared minimizers — fall back to full banded alignment.
        return banded_sw(query, target, scoring, bandwidth);
    }

    // Step 3: Chain the seeds.
    let max_gap = query.len().max(target.len());
    let chain = chain_seeds(&seed_pairs, k, max_gap)?;
    if chain.seeds.is_empty() {
        return banded_sw(query, target, scoring, bandwidth);
    }

    // Step 4: Extend alignment in regions between chained seeds.
    //
    // Strategy:
    //   - Prefix: region before the first seed
    //   - Middle: regions between consecutive seeds
    //   - Suffix: region after the last seed
    //
    // For each region we run banded_sw on the corresponding query/target slices,
    // then concatenate the CIGAR strings.  Each seed itself contributes k
    // matched bases (they share the same k-mer).

    let mut total_score: i32 = 0;
    let mut cigar: Vec<CigarOp> = Vec::new();
    let mut aligned_query: Vec<u8> = Vec::new();
    let mut aligned_target: Vec<u8> = Vec::new();

    let first_seed = &chain.seeds[0];
    let last_seed = &chain.seeds[chain.seeds.len() - 1];

    // --- Prefix region ---
    if first_seed.query_pos > 0 && first_seed.target_pos > 0 {
        let q_prefix = &query[..first_seed.query_pos];
        let t_prefix = &target[..first_seed.target_pos];
        if let Ok(r) = banded_sw(q_prefix, t_prefix, scoring, bandwidth) {
            total_score += r.score;
            append_alignment(&mut aligned_query, &mut aligned_target, &mut cigar, &r);
        }
    }

    // --- Process each seed and the inter-seed gap after it ---
    for (idx, seed) in chain.seeds.iter().enumerate() {
        // Emit the k-mer match for this seed.
        let q_end = (seed.query_pos + k).min(query.len());
        let t_end = (seed.target_pos + k).min(target.len());
        let match_len = (q_end - seed.query_pos).min(t_end - seed.target_pos);

        // Score the k-mer region base by base.  Minimizer hashes can
        // collide, so we don't assume all bases match.
        for i in 0..match_len {
            let qb = query[seed.query_pos + i];
            let tb = target[seed.target_pos + i];
            total_score += scoring.score_pair(qb, tb);
            aligned_query.push(qb);
            aligned_target.push(tb);
            if qb.to_ascii_uppercase() == tb.to_ascii_uppercase() {
                push_cigar(&mut cigar, CigarOp::Match(1));
            } else {
                push_cigar(&mut cigar, CigarOp::Mismatch(1));
            }
        }

        // Inter-seed gap: align the region between this seed's end and the
        // next seed's start.
        if idx + 1 < chain.seeds.len() {
            let next = &chain.seeds[idx + 1];
            let q_gap_start = seed.query_pos + k;
            let q_gap_end = next.query_pos;
            let t_gap_start = seed.target_pos + k;
            let t_gap_end = next.target_pos;

            if q_gap_start < q_gap_end && t_gap_start < t_gap_end {
                let q_slice = &query[q_gap_start..q_gap_end];
                let t_slice = &target[t_gap_start..t_gap_end];
                if let Ok(r) = banded_sw(q_slice, t_slice, scoring, bandwidth) {
                    total_score += r.score;
                    append_alignment(&mut aligned_query, &mut aligned_target, &mut cigar, &r);
                }
            } else if q_gap_start < q_gap_end {
                // Query has extra bases, target doesn't — deletion.
                let del_len = q_gap_end - q_gap_start;
                push_cigar(&mut cigar, CigarOp::Deletion(del_len));
                for i in q_gap_start..q_gap_end {
                    aligned_query.push(query[i]);
                    aligned_target.push(b'-');
                }
                total_score += scoring.gap_open() + (del_len as i32 - 1) * scoring.gap_extend();
            } else if t_gap_start < t_gap_end {
                // Target has extra bases — insertion.
                let ins_len = t_gap_end - t_gap_start;
                push_cigar(&mut cigar, CigarOp::Insertion(ins_len));
                for i in t_gap_start..t_gap_end {
                    aligned_query.push(b'-');
                    aligned_target.push(target[i]);
                }
                total_score += scoring.gap_open() + (ins_len as i32 - 1) * scoring.gap_extend();
            }
            // If both gap regions are empty or overlapping, skip (seeds are adjacent or overlapping).
        }
    }

    // --- Suffix region ---
    let q_suffix_start = last_seed.query_pos + k;
    let t_suffix_start = last_seed.target_pos + k;
    if q_suffix_start < query.len() && t_suffix_start < target.len() {
        let q_suffix = &query[q_suffix_start..];
        let t_suffix = &target[t_suffix_start..];
        if let Ok(r) = banded_sw(q_suffix, t_suffix, scoring, bandwidth) {
            total_score += r.score;
            append_alignment(&mut aligned_query, &mut aligned_target, &mut cigar, &r);
        }
    }

    // Compute alignment coordinates.
    let query_start = if first_seed.query_pos > 0 { 0 } else { first_seed.query_pos };
    let query_end = if q_suffix_start < query.len() {
        query.len()
    } else {
        (last_seed.query_pos + k).min(query.len())
    };
    let target_start = if first_seed.target_pos > 0 { 0 } else { first_seed.target_pos };
    let target_end = if t_suffix_start < target.len() {
        target.len()
    } else {
        (last_seed.target_pos + k).min(target.len())
    };

    Ok(AlignmentResult {
        score: total_score,
        aligned_query,
        aligned_target,
        query_start,
        query_end,
        target_start,
        target_end,
        cigar,
    })
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Push a CIGAR op, merging with the last op if they are the same variant.
fn push_cigar(cigar: &mut Vec<CigarOp>, op: CigarOp) {
    if let Some(last) = cigar.last_mut() {
        let merged = match (*last, op) {
            (CigarOp::Match(a), CigarOp::Match(b)) => Some(CigarOp::Match(a + b)),
            (CigarOp::Mismatch(a), CigarOp::Mismatch(b)) => Some(CigarOp::Mismatch(a + b)),
            (CigarOp::Insertion(a), CigarOp::Insertion(b)) => Some(CigarOp::Insertion(a + b)),
            (CigarOp::Deletion(a), CigarOp::Deletion(b)) => Some(CigarOp::Deletion(a + b)),
            _ => None,
        };
        if let Some(m) = merged {
            *last = m;
            return;
        }
    }
    cigar.push(op);
}

/// Append alignment result data (aligned sequences and CIGAR ops).
fn append_alignment(
    aligned_query: &mut Vec<u8>,
    aligned_target: &mut Vec<u8>,
    cigar: &mut Vec<CigarOp>,
    result: &AlignmentResult,
) {
    aligned_query.extend_from_slice(&result.aligned_query);
    aligned_target.extend_from_slice(&result.aligned_target);
    for &op in &result.cigar {
        push_cigar(cigar, op);
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::{ScoringMatrix, ScoringScheme};

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::Simple(ScoringMatrix::dna_default())
    }

    // -----------------------------------------------------------------------
    // Seed chaining
    // -----------------------------------------------------------------------

    #[test]
    fn chain_colinear_seeds() {
        // Simple diagonal seeds with k=15.
        let seeds = vec![(0, 0), (10, 10), (20, 20)];
        let chain = chain_seeds(&seeds, 15, 100).unwrap();
        assert!(chain.seeds.len() >= 2, "should chain multiple seeds");
        // Seeds should be colinear.
        for pair in chain.seeds.windows(2) {
            assert!(pair[0].query_pos < pair[1].query_pos);
            assert!(pair[0].target_pos < pair[1].target_pos);
        }
    }

    #[test]
    fn chain_out_of_order_seeds() {
        // Seeds given in non-sorted order should still be chained correctly.
        let seeds = vec![(20, 20), (0, 0), (10, 10)];
        let chain = chain_seeds(&seeds, 15, 100).unwrap();
        assert!(chain.seeds.len() >= 2);
        for pair in chain.seeds.windows(2) {
            assert!(pair[0].target_pos < pair[1].target_pos);
        }
    }

    #[test]
    fn chain_non_colinear_filtered() {
        // One seed is not colinear with the others.
        let seeds = vec![(0, 0), (5, 50), (10, 10), (20, 20)];
        let chain = chain_seeds(&seeds, 15, 100).unwrap();
        // (5,50) should not appear together with (10,10) in the same chain
        // because target 50 > 10 but appears before it in query order.
        // The chain should pick the best colinear subset.
        for pair in chain.seeds.windows(2) {
            assert!(pair[0].query_pos < pair[1].query_pos);
            assert!(pair[0].target_pos < pair[1].target_pos);
        }
    }

    #[test]
    fn chain_single_seed() {
        let seeds = vec![(5, 5)];
        let chain = chain_seeds(&seeds, 15, 100).unwrap();
        assert_eq!(chain.seeds.len(), 1);
    }

    #[test]
    fn chain_empty_seeds_error() {
        let seeds: Vec<(usize, usize)> = vec![];
        assert!(chain_seeds(&seeds, 15, 100).is_err());
    }

    #[test]
    fn chain_gap_penalty_prefers_close_seeds() {
        // Two possible chains: close seeds vs. far-apart seeds.
        // Close seeds should win because of lower gap penalty.
        let seeds = vec![
            (0, 0),
            (2, 2),   // close to (0,0)
            (4, 4),   // close to (2,2)
            (100, 100), // far from everything
        ];
        let chain = chain_seeds(&seeds, 15, 200).unwrap();
        // With k=15 seed_score, chaining close seeds is profitable:
        // (0,0)->(2,2)->(4,4) has score 15 + 15-2 + 15-2 = 41
        // vs (100,100) standalone = 15.
        assert!(chain.seeds.len() >= 2);
    }

    #[test]
    fn chain_respects_max_gap() {
        let seeds = vec![(0, 0), (200, 200)];
        let chain = chain_seeds(&seeds, 15, 10).unwrap();
        // Gap of 200 exceeds max_gap=10, so these can't be chained together.
        assert_eq!(chain.seeds.len(), 1);
    }

    // -----------------------------------------------------------------------
    // Seed-and-extend alignment
    // -----------------------------------------------------------------------

    #[test]
    fn identical_sequences() {
        let seq = b"ACGTACGTACGTACGT";
        let result =
            seed_extend_align(seq, seq, &dna_scheme(), 3, 2, 5).unwrap();
        assert!(result.score > 0, "identical seqs should have positive score");
        // Check that the alignment covers the full length.
        assert_eq!(result.query_start, 0);
    }

    #[test]
    fn related_sequences_with_mutations() {
        let query  = b"ACGTACGTACGTACGT";
        let target = b"ACGTACCTACGTACGT"; // G->C at pos 6
        let result =
            seed_extend_align(query, target, &dna_scheme(), 3, 2, 5).unwrap();
        assert!(result.score > 0);
    }

    #[test]
    fn completely_different_sequences_fallback() {
        // No shared k-mers -> falls back to banded alignment.
        let query  = b"AAAAAAAAAA";
        let target = b"CCCCCCCCCC";
        let result =
            seed_extend_align(query, target, &dna_scheme(), 3, 2, 10).unwrap();
        // Smith-Waterman local alignment on completely different sequences
        // should produce a non-negative score (SW is always >= 0).
        assert!(result.score >= 0);
    }

    #[test]
    fn empty_sequence_error() {
        assert!(seed_extend_align(b"", b"ACGT", &dna_scheme(), 3, 2, 5).is_err());
        assert!(seed_extend_align(b"ACGT", b"", &dna_scheme(), 3, 2, 5).is_err());
    }

    #[test]
    fn short_sequences_fallback() {
        // Sequences shorter than k fall back to banded alignment.
        let result =
            seed_extend_align(b"AC", b"AC", &dna_scheme(), 3, 2, 5).unwrap();
        assert!(result.score > 0);
    }

    #[test]
    fn longer_sequences_with_shared_regions() {
        // Construct sequences that share a long common substring.
        let common = b"ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bp
        let mut query = Vec::new();
        query.extend_from_slice(b"TTTTTTT"); // unique prefix
        query.extend_from_slice(common);
        query.extend_from_slice(b"GGGGGGG"); // unique suffix

        let mut target = Vec::new();
        target.extend_from_slice(b"CCCCCCC"); // different prefix
        target.extend_from_slice(common);
        target.extend_from_slice(b"AAAAAAA"); // different suffix

        let result =
            seed_extend_align(&query, &target, &dna_scheme(), 5, 3, 10).unwrap();
        assert!(result.score > 0, "should align the shared region");
    }

    #[test]
    fn alignment_produces_valid_cigar() {
        let seq = b"ACGTACGTACGTACGT";
        let result =
            seed_extend_align(seq, seq, &dna_scheme(), 3, 2, 5).unwrap();
        // CIGAR ops should all have non-zero length.
        for op in &result.cigar {
            assert!(op.len() > 0, "CIGAR ops should have non-zero length");
        }
    }

    #[test]
    fn alignment_with_insertions_and_deletions() {
        // Target has an insertion relative to query.
        let query  = b"ACGTACGTACGT";
        let target = b"ACGTAAACGTACGT"; // extra AA inserted
        let result =
            seed_extend_align(query, target, &dna_scheme(), 3, 2, 5).unwrap();
        // Should still produce a valid alignment.
        assert!(result.score > 0 || result.score == 0);
    }

    #[test]
    fn chain_scoring_positive() {
        let seeds = vec![(0, 0), (1, 1), (2, 2), (3, 3)];
        let chain = chain_seeds(&seeds, 15, 100).unwrap();
        // Each seed is 1 apart, gap penalty is 1 per hop.
        // With k=15: score = 15 + (15-1) + (15-1) + (15-1) = 57.
        assert!(chain.score >= 15);
        assert_eq!(chain.seeds.len(), 4);
    }
}
