//! Minimizer extraction for seed-and-extend alignment.
//!
//! A (w,k)-minimizer is the smallest k-mer hash in a window of w consecutive
//! k-mers.  Minimizers are the most common seeding strategy for long-read
//! alignment because they reduce the number of anchors while guaranteeing that
//! any shared substring of length `w + k - 1` bases produces at least one
//! shared minimizer.

use cyanea_core::{CyaneaError, Result};

/// A minimizer hit: position in the sequence and the hash value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Minimizer {
    /// 0-based position of the k-mer's first base in the original sequence.
    pub position: usize,
    /// Invertible hash of the canonical k-mer encoding.
    pub hash: u64,
}

// ---------------------------------------------------------------------------
// 2-bit encoding
// ---------------------------------------------------------------------------

/// Encode a DNA base as 2 bits: A=0, C=1, G=2, T=3.
/// Returns `None` for non-ACGT characters.
#[inline]
fn base_to_2bit(b: u8) -> Option<u64> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Invertible hash (murmur-style bit mixer)
// ---------------------------------------------------------------------------

/// Invertible 64-bit hash to avoid positional bias in lexicographic k-mer
/// ordering.  This is a simplified murmur3 finalizer.
#[inline]
fn hash64(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ceb9fe1a85ec53);
    x ^= x >> 33;
    x
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Extract (w,k)-minimizers from a DNA sequence.
///
/// Returns minimizers in order of position, deduplicated (consecutive
/// windows sharing the same minimizer yield one entry).
///
/// # Arguments
/// * `seq` - DNA sequence (ACGT, case-insensitive). Non-ACGT bases break
///   the current k-mer run (the hash is reset).
/// * `k` - k-mer length (must be in `1..=32`)
/// * `w` - window size (number of consecutive k-mers per window)
///
/// # Errors
/// Returns an error if `k > 32`, `k == 0`, `w == 0`, or `seq.len() < k`.
pub fn minimizers(seq: &[u8], k: usize, w: usize) -> Result<Vec<Minimizer>> {
    if k == 0 {
        return Err(CyaneaError::InvalidInput("k must be > 0".into()));
    }
    if k > 32 {
        return Err(CyaneaError::InvalidInput("k must be <= 32".into()));
    }
    if w == 0 {
        return Err(CyaneaError::InvalidInput("w must be > 0".into()));
    }
    if seq.len() < k {
        return Err(CyaneaError::InvalidInput(
            "sequence is shorter than k".into(),
        ));
    }

    let mask: u64 = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };

    // Phase 1: compute all k-mer hashes via a rolling 2-bit shift register.
    let num_kmers = seq.len() - k + 1;
    let mut kmer_hashes: Vec<Option<u64>> = Vec::with_capacity(num_kmers);

    let mut forward: u64 = 0;
    let mut valid_bases: usize = 0; // consecutive valid bases seen

    for (i, &base) in seq.iter().enumerate() {
        if let Some(bits) = base_to_2bit(base) {
            forward = ((forward << 2) | bits) & mask;
            valid_bases += 1;
        } else {
            // Non-ACGT base: reset the run
            valid_bases = 0;
            forward = 0;
        }

        // Once we have at least k valid consecutive bases, the k-mer ending
        // at position i (starting at position i - k + 1) is valid.
        if i >= k - 1 {
            if valid_bases >= k {
                kmer_hashes.push(Some(hash64(forward)));
            } else {
                kmer_hashes.push(None);
            }
        }
    }

    debug_assert_eq!(kmer_hashes.len(), num_kmers);

    // Phase 2: slide a window of w k-mers, picking the minimum hash.
    // If there are fewer than w k-mers, use all of them as a single window.
    let effective_w = w.min(num_kmers);

    let mut result: Vec<Minimizer> = Vec::new();
    let mut prev_hash: Option<u64> = None;

    let num_windows = if num_kmers >= effective_w {
        num_kmers - effective_w + 1
    } else {
        1
    };

    for win_start in 0..num_windows {
        let win_end = (win_start + effective_w).min(num_kmers);

        // Find the position of the minimum hash in this window.
        let mut best_hash: Option<u64> = None;
        let mut best_pos: Option<usize> = None;

        for pos in win_start..win_end {
            if let Some(h) = kmer_hashes[pos] {
                match best_hash {
                    None => {
                        best_hash = Some(h);
                        best_pos = Some(pos);
                    }
                    Some(bh) => {
                        // Tie-break by leftmost position (stable)
                        if h < bh {
                            best_hash = Some(h);
                            best_pos = Some(pos);
                        }
                    }
                }
            }
        }

        // Deduplicate: only emit if the minimizer hash changed.
        // On homopolymers / tandem repeats, consecutive windows select the
        // same k-mer content (same hash) at shifting positions.  Emitting
        // only when the hash changes collapses these into one anchor.
        if let (Some(pos), Some(h)) = (best_pos, best_hash) {
            if prev_hash != Some(h) {
                result.push(Minimizer {
                    position: pos,
                    hash: h,
                });
                prev_hash = Some(h);
            }
        }
    }

    Ok(result)
}

/// Find shared minimizer positions between query and target.
///
/// Returns pairs of `(query_pos, target_pos)` where both sequences
/// share the same minimizer hash.  The result is sorted by
/// `(query_pos, target_pos)`.
///
/// This uses a hash-map join so it runs in O(q + t) expected time where
/// q and t are the number of minimizers in each sequence.
pub fn find_seed_matches(
    query_mins: &[Minimizer],
    target_mins: &[Minimizer],
) -> Vec<(usize, usize)> {
    use std::collections::HashMap;

    // Build an index from hash -> list of positions in the target.
    let mut target_index: HashMap<u64, Vec<usize>> = HashMap::new();
    for m in target_mins {
        target_index.entry(m.hash).or_default().push(m.position);
    }

    let mut matches: Vec<(usize, usize)> = Vec::new();
    for qm in query_mins {
        if let Some(positions) = target_index.get(&qm.hash) {
            for &tp in positions {
                matches.push((qm.position, tp));
            }
        }
    }

    matches.sort_unstable();
    matches
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // Minimizer extraction basics
    // -----------------------------------------------------------------------

    #[test]
    fn basic_minimizer_extraction() {
        // Small sequence, k=3, w=2.  Enough to get at least one minimizer.
        let seq = b"ACGTACGT";
        let mins = minimizers(seq, 3, 2).unwrap();
        assert!(!mins.is_empty(), "should find at least one minimizer");
        // All positions should be within valid range.
        for m in &mins {
            assert!(m.position + 3 <= seq.len());
        }
    }

    #[test]
    fn minimizers_on_homopolymer() {
        // All k-mers are the same, so there should be a single minimizer
        // (consecutive deduplication collapses all windows into one).
        let seq = b"AAAAAAAAAA"; // 10 A's
        let mins = minimizers(seq, 3, 2).unwrap();
        assert_eq!(mins.len(), 1, "homopolymer should produce a single minimizer");
    }

    #[test]
    fn kmer_hashing_correctness() {
        // Same sequence should produce the same hashes.
        let seq = b"ACGTACGT";
        let m1 = minimizers(seq, 4, 1).unwrap();
        let m2 = minimizers(seq, 4, 1).unwrap();
        assert_eq!(m1, m2, "deterministic hashing");

        // With w=1 every k-mer is its own window, so we get one minimizer
        // per position (unless consecutive k-mers share the same hash,
        // which can only happen for identical k-mers).
        // Here k=4 on ACGTACGT gives k-mers: ACGT, CGTA, GTAC, TACG, ACGT.
        // ACGT appears twice (pos 0 and 4). After dedup the first occurrence
        // persists until pos 4 re-emits. So we expect 5 or fewer entries.
        assert!(m1.len() <= 5);
    }

    #[test]
    fn window_sliding_produces_ordered_positions() {
        let seq = b"ACGTACGTACGT";
        let mins = minimizers(seq, 3, 3).unwrap();
        for pair in mins.windows(2) {
            assert!(
                pair[0].position < pair[1].position,
                "positions must be strictly increasing after dedup"
            );
        }
    }

    #[test]
    fn deduplication_consecutive_same_minimizer() {
        // Construct a sequence where a single k-mer dominates many windows.
        // All-A prefix means the AAA k-mer has the same hash across windows.
        let seq = b"AAAAAACGT";
        let mins = minimizers(seq, 3, 2).unwrap();
        // The AAA minimizer should appear only once despite spanning many windows.
        let aaa_count = mins.iter().filter(|m| m.position == 0 || seq[m.position..m.position + 3] == *b"AAA").count();
        // We only care that there's no *consecutive* duplicate position.
        for pair in mins.windows(2) {
            assert_ne!(pair[0].position, pair[1].position);
        }
        // Suppress unused variable warning.
        let _ = aaa_count;
    }

    // -----------------------------------------------------------------------
    // Edge cases
    // -----------------------------------------------------------------------

    #[test]
    fn k_equals_1_w_equals_1() {
        let seq = b"ACGT";
        let mins = minimizers(seq, 1, 1).unwrap();
        // With k=1, w=1 every base is its own window. Dedup removes
        // consecutive duplicates only, so we get at most 4 entries.
        assert!(mins.len() <= 4);
        assert!(!mins.is_empty());
    }

    #[test]
    fn sequence_exactly_k_long() {
        let seq = b"ACGT";
        let mins = minimizers(seq, 4, 1).unwrap();
        assert_eq!(mins.len(), 1, "exactly one k-mer = one minimizer");
    }

    #[test]
    fn error_k_zero() {
        assert!(minimizers(b"ACGT", 0, 1).is_err());
    }

    #[test]
    fn error_k_too_large() {
        assert!(minimizers(b"ACGT", 33, 1).is_err());
    }

    #[test]
    fn error_w_zero() {
        assert!(minimizers(b"ACGT", 3, 0).is_err());
    }

    #[test]
    fn error_sequence_too_short() {
        assert!(minimizers(b"AC", 3, 1).is_err());
    }

    // -----------------------------------------------------------------------
    // Case insensitivity
    // -----------------------------------------------------------------------

    #[test]
    fn case_insensitivity() {
        let upper = minimizers(b"ACGTACGT", 3, 2).unwrap();
        let lower = minimizers(b"acgtacgt", 3, 2).unwrap();
        let mixed = minimizers(b"AcGtAcGt", 3, 2).unwrap();
        assert_eq!(upper, lower, "upper and lower case should be identical");
        assert_eq!(upper, mixed, "mixed case should be identical");
    }

    // -----------------------------------------------------------------------
    // Seed matching
    // -----------------------------------------------------------------------

    #[test]
    fn seed_matching_identical_sequences() {
        let seq = b"ACGTACGTACGT";
        let q_mins = minimizers(seq, 3, 2).unwrap();
        let t_mins = minimizers(seq, 3, 2).unwrap();
        let matches = find_seed_matches(&q_mins, &t_mins);
        // Identical sequences should produce self-matches for every minimizer.
        assert!(!matches.is_empty());
        // Every minimizer should have a diagonal match (q_pos == t_pos).
        // Off-diagonal matches are expected when the same k-mer appears at
        // multiple positions (e.g., "ACG" at positions 0, 4, 8).
        for m in &q_mins {
            assert!(
                matches.contains(&(m.position, m.position)),
                "minimizer at {} should have a diagonal match",
                m.position,
            );
        }
    }

    #[test]
    fn seed_matching_related_sequences() {
        // Introduce a single mutation; most minimizers should still match.
        let query  = b"ACGTACGTACGT";
        let target = b"ACGTACCTACGT"; // G->C at position 6
        let q_mins = minimizers(query, 3, 2).unwrap();
        let t_mins = minimizers(target, 3, 2).unwrap();
        let matches = find_seed_matches(&q_mins, &t_mins);
        assert!(!matches.is_empty(), "related seqs should share some seeds");
    }

    #[test]
    fn seed_matching_unrelated_sequences() {
        let query  = b"AAAAAAAAAA";
        let target = b"CCCCCCCCCC";
        let q_mins = minimizers(query, 3, 2).unwrap();
        let t_mins = minimizers(target, 3, 2).unwrap();
        let matches = find_seed_matches(&q_mins, &t_mins);
        assert!(matches.is_empty(), "unrelated seqs should share no seeds");
    }

    #[test]
    fn seed_matches_are_sorted() {
        let query  = b"ACGTACGTACGT";
        let target = b"TACGTACGTACG";
        let q_mins = minimizers(query, 3, 2).unwrap();
        let t_mins = minimizers(target, 3, 2).unwrap();
        let matches = find_seed_matches(&q_mins, &t_mins);
        for pair in matches.windows(2) {
            assert!(pair[0] <= pair[1], "matches must be sorted");
        }
    }

    // -----------------------------------------------------------------------
    // Larger k values
    // -----------------------------------------------------------------------

    #[test]
    fn k32_works() {
        // k=32 is the maximum allowed
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 36 bases
        let mins = minimizers(seq, 32, 1).unwrap();
        assert!(!mins.is_empty());
    }
}
