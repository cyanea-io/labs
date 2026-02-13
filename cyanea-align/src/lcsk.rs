//! LCSk++ sparse alignment algorithm.
//!
//! Computes the longest common subsequence in k-mer space (LCSk++) between two
//! sequences.  Rather than aligning every base pair via O(mn) dynamic
//! programming, this approach:
//!
//! 1. **Finds k-mer matches** — positions where an identical k-mer occurs in
//!    both sequences.
//! 2. **Chains matches** — selects the highest-scoring colinear subset using a
//!    Fenwick-tree-accelerated DP in O(n log n) time.
//!
//! The LCSk++ score counts the number of matching bases covered by the optimal
//! chain, allowing adjacent anchors to merge (the "++" extension over plain
//! LCSk).
//!
//! Reference: Filip Pavetić et al., "Fast and simple algorithms for computing
//! both LCSk and LCSk+", *Information Processing Letters*, 2017.
//!
//! # Example
//!
//! ```
//! use cyanea_align::lcsk::{sparse_align, SparseAlignResult};
//!
//! let result = sparse_align(b"ACGTACGT", b"ACGTACGT", 4).unwrap();
//! assert_eq!(result.score, 8);
//! assert_eq!(result.k, 4);
//! ```

use cyanea_core::{CyaneaError, Result};
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Result of a sparse alignment.
#[derive(Debug, Clone)]
pub struct SparseAlignResult {
    /// LCSk++ score (number of matching bases in the chain).
    pub score: usize,
    /// Anchored match positions: (pos_in_a, pos_in_b) for each k-mer anchor in the chain.
    pub anchors: Vec<(usize, usize)>,
    /// k-mer length used.
    pub k: usize,
}

// ---------------------------------------------------------------------------
// Fenwick tree for maximum prefix queries
// ---------------------------------------------------------------------------

/// Fenwick tree (Binary Indexed Tree) for maximum prefix queries.
///
/// Each node stores a `(value, backpointer)` pair.  `update` sets position `i`
/// to value `v` (if it improves the current maximum), and `query` returns the
/// maximum `(value, backpointer)` over the prefix `[0, i]`.
struct FenwickMax {
    tree: Vec<(usize, i32)>,
    n: usize,
}

impl FenwickMax {
    fn new(n: usize) -> Self {
        Self {
            tree: vec![(0, -1); n + 1],
            n,
        }
    }

    /// Update position `i` with value `v` and backpointer `bp`.
    ///
    /// Only updates if `v` strictly exceeds the current stored value at each
    /// affected tree node.
    fn update(&mut self, mut i: usize, v: usize, bp: i32) {
        i += 1; // 1-indexed
        while i <= self.n {
            if v > self.tree[i].0 {
                self.tree[i] = (v, bp);
            }
            i += i & i.wrapping_neg();
        }
    }

    /// Query maximum value in `[0, i]`, returning `(value, backpointer)`.
    fn query(&self, mut i: usize) -> (usize, i32) {
        i += 1; // 1-indexed
        let mut best = (0usize, -1i32);
        while i > 0 {
            if self.tree[i].0 > best.0 {
                best = self.tree[i];
            }
            i -= i & i.wrapping_neg();
        }
        best
    }
}

// ---------------------------------------------------------------------------
// K-mer matching
// ---------------------------------------------------------------------------

/// Find all positions where k-mers match between two sequences.
///
/// Returns a vector of `(pos_in_a, pos_in_b)` pairs for every occurrence of an
/// identical k-mer in both `seq_a` and `seq_b`.  The pairs are not sorted in
/// any particular order.
///
/// # Panics
///
/// Callers must ensure `k > 0` and that both sequences are at least `k` bases
/// long.  The public [`sparse_align`] function validates these preconditions.
pub fn find_kmer_matches(seq_a: &[u8], seq_b: &[u8], k: usize) -> Vec<(usize, usize)> {
    if seq_a.len() < k || seq_b.len() < k {
        return Vec::new();
    }

    // Index all k-mers in seq_a by their byte content.
    let mut index: HashMap<&[u8], Vec<usize>> = HashMap::new();
    for i in 0..=seq_a.len() - k {
        index.entry(&seq_a[i..i + k]).or_default().push(i);
    }

    // Scan seq_b and collect matches.
    let mut matches = Vec::new();
    for j in 0..=seq_b.len() - k {
        if let Some(positions) = index.get(&seq_b[j..j + k]) {
            for &i in positions {
                matches.push((i, j));
            }
        }
    }

    matches
}

// ---------------------------------------------------------------------------
// Core LCSk++ algorithm
// ---------------------------------------------------------------------------

/// Core LCSk++ algorithm using Fenwick-tree-accelerated DP.
///
/// Takes a list of k-mer match positions `(pos_in_a, pos_in_b)` and returns the
/// optimal chain score together with the chain anchors.
///
/// # Algorithm
///
/// Matches are sorted by `(pos_b, pos_a)`.  For each match `(i, j)`, the DP
/// recurrence is:
///
/// ```text
/// dp[i][j] = k + max { dp[i'][j'] : i' + k <= i and j' + k <= j }
/// ```
///
/// The "++" extension additionally checks whether the previous anchor is
/// immediately adjacent (i.e. `i' + k == i` and `j' + k == j`), in which case a
/// single base is added instead of a full k-mer to avoid double-counting
/// overlapping coverage.
///
/// Coordinate compression and a Fenwick tree reduce the per-match work to
/// O(log n).
///
/// Returns `(score, chain)` where `chain` contains the anchors in order.
pub fn lcsk_plusplus(matches: &[(usize, usize)], k: usize) -> (usize, Vec<(usize, usize)>) {
    if matches.is_empty() || k == 0 {
        return (0, Vec::new());
    }

    let n = matches.len();

    // Sort matches by (pos_b, pos_a).
    let mut sorted: Vec<(usize, usize)> = matches.to_vec();
    sorted.sort_unstable_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));

    // Coordinate compression on the "a" axis.
    // We need to query: for all previous matches with a' + k <= a_current,
    // what is the maximum DP value?  We compress the "end-of-anchor" positions
    // (a + k - 1) into a contiguous range for the Fenwick tree.
    let mut a_coords: Vec<usize> = sorted.iter().map(|&(a, _)| a).collect();
    a_coords.sort_unstable();
    a_coords.dedup();
    let compress = |val: usize| -> usize {
        a_coords.binary_search(&val).unwrap_or_else(|x| x)
    };

    // DP arrays.
    let mut dp: Vec<usize> = vec![0; n];
    let mut prev: Vec<i32> = vec![-1; n];

    // Fenwick tree over compressed a-coordinates.
    // We process matches in order of increasing pos_b.  When we process match
    // at (a, b), we need predecessors with b' + k <= b (non-overlapping on the
    // b axis) and a' + k <= a (non-overlapping on the a axis).
    //
    // Strategy: process in sweepline order over b.  Group matches by b value.
    // Before processing group at b_cur, insert all matches with b' + k <= b_cur
    // into the Fenwick tree.
    let mut fenwick = FenwickMax::new(a_coords.len());

    // Track which matches have been inserted into the Fenwick tree.
    // We use a pointer into the sorted array, advancing as b increases.
    // Matches are sorted by b, so we can insert in order.
    let mut insert_ptr: usize = 0;

    // We also need to handle the "++" adjacency check.  For that we maintain a
    // map from (a + k, b + k) -> (dp_value, match_index) for the immediately-
    // preceding adjacent extension.
    let mut adjacent: HashMap<(usize, usize), (usize, i32)> = HashMap::new();

    // Process each match.
    for idx in 0..n {
        let (a_i, b_i) = sorted[idx];

        // Insert all matches whose b-end (b' + k) <= b_i into the Fenwick tree.
        // Since sorted is by b, we advance insert_ptr while sorted[insert_ptr].1 + k <= b_i.
        while insert_ptr < n && sorted[insert_ptr].1 + k <= b_i {
            let (a_p, _) = sorted[insert_ptr];
            let ca = compress(a_p);
            fenwick.update(ca, dp[insert_ptr], insert_ptr as i32);
            insert_ptr += 1;
        }

        // Base case: this anchor starts a new chain.
        dp[idx] = k;
        prev[idx] = -1;

        // Case 1: non-overlapping predecessor via Fenwick tree.
        // We need a' + k <= a_i, which means a' <= a_i - k.
        if a_i >= k {
            // Find the largest compressed index with a_coords[c] <= a_i - k.
            let upper = a_i - k;
            let ca_upper = match a_coords.binary_search(&upper) {
                Ok(pos) => pos,
                Err(pos) => {
                    if pos == 0 {
                        usize::MAX // no valid predecessor
                    } else {
                        pos - 1
                    }
                }
            };
            if ca_upper != usize::MAX {
                let (best_val, best_bp) = fenwick.query(ca_upper);
                if best_val > 0 {
                    let candidate = best_val + k;
                    if candidate > dp[idx] {
                        dp[idx] = candidate;
                        prev[idx] = best_bp;
                    }
                }
            }
        }

        // Case 2: "++" adjacent extension.
        // Check if there is a predecessor at (a_i - 1, b_i - 1), meaning the
        // previous anchor ends exactly where this one starts (shifted by 1 base).
        if a_i > 0 && b_i > 0 {
            if let Some(&(adj_val, adj_idx)) = adjacent.get(&(a_i, b_i)) {
                // Adjacent extension adds 1 base (not k, to avoid double-counting).
                let candidate = adj_val + 1;
                if candidate > dp[idx] {
                    dp[idx] = candidate;
                    prev[idx] = adj_idx;
                }
            }
        }

        // Record this match for future "++" adjacency checks.
        // A match at (a_i, b_i) with k-mer length k ends at (a_i + k, b_i + k).
        // The next adjacent anchor would start at (a_i + 1, b_i + 1), so the
        // adjacency key is (a_i + 1, b_i + 1).
        let adj_key = (a_i + 1, b_i + 1);
        match adjacent.get(&adj_key) {
            Some(&(existing_val, _)) if existing_val >= dp[idx] => {}
            _ => {
                adjacent.insert(adj_key, (dp[idx], idx as i32));
            }
        }
    }

    // Find the global best.
    let mut best_idx = 0;
    for i in 1..n {
        if dp[i] > dp[best_idx] {
            best_idx = i;
        }
    }

    let score = dp[best_idx];

    // Traceback to recover the chain.
    let mut chain_indices: Vec<usize> = Vec::new();
    let mut idx = best_idx as i32;
    while idx >= 0 {
        chain_indices.push(idx as usize);
        idx = prev[idx as usize];
    }
    chain_indices.reverse();

    let chain: Vec<(usize, usize)> = chain_indices.iter().map(|&i| sorted[i]).collect();

    (score, chain)
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Perform sparse alignment between two sequences using LCSk++.
///
/// This is the full pipeline: find k-mer matches, run the LCSk++ DP, and
/// return a [`SparseAlignResult`].
///
/// # Arguments
/// * `seq_a` - First sequence (arbitrary alphabet).
/// * `seq_b` - Second sequence.
/// * `k`     - k-mer length.  Must be > 0 and both sequences must be at least
///   `k` bases long.
///
/// # Errors
///
/// Returns [`CyaneaError::InvalidInput`] if `k == 0` or either sequence is
/// shorter than `k`.
///
/// # Example
///
/// ```
/// use cyanea_align::lcsk::sparse_align;
///
/// let result = sparse_align(b"ACGTACGT", b"ACGTACGT", 4).unwrap();
/// assert_eq!(result.score, 8);
/// ```
pub fn sparse_align(seq_a: &[u8], seq_b: &[u8], k: usize) -> Result<SparseAlignResult> {
    if k == 0 {
        return Err(CyaneaError::InvalidInput(
            "k must be greater than 0".into(),
        ));
    }
    if seq_a.len() < k {
        return Err(CyaneaError::InvalidInput(format!(
            "sequence A length ({}) is shorter than k ({})",
            seq_a.len(),
            k,
        )));
    }
    if seq_b.len() < k {
        return Err(CyaneaError::InvalidInput(format!(
            "sequence B length ({}) is shorter than k ({})",
            seq_b.len(),
            k,
        )));
    }

    let matches = find_kmer_matches(seq_a, seq_b, k);
    let (score, anchors) = lcsk_plusplus(&matches, k);

    Ok(SparseAlignResult { score, anchors, k })
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // Fenwick tree
    // -----------------------------------------------------------------------

    #[test]
    fn fenwick_basic_update_query() {
        let mut fw = FenwickMax::new(10);
        fw.update(3, 5, 0);
        fw.update(7, 3, 1);
        assert_eq!(fw.query(3).0, 5);
        assert_eq!(fw.query(7).0, 5); // prefix max includes index 3
        fw.update(5, 10, 2);
        assert_eq!(fw.query(5).0, 10);
        assert_eq!(fw.query(9).0, 10);
    }

    #[test]
    fn fenwick_empty() {
        let fw = FenwickMax::new(5);
        assert_eq!(fw.query(4), (0, -1));
    }

    // -----------------------------------------------------------------------
    // find_kmer_matches
    // -----------------------------------------------------------------------

    #[test]
    fn kmer_matches_identical() {
        let seq = b"ACGT";
        let matches = find_kmer_matches(seq, seq, 2);
        // k-mers: AC, CG, GT — each matches at the same position.
        assert_eq!(matches.len(), 3);
        for &(a, b) in &matches {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn kmer_matches_no_overlap() {
        let a = b"AAAA";
        let b = b"CCCC";
        let matches = find_kmer_matches(a, b, 2);
        assert!(matches.is_empty());
    }

    #[test]
    fn kmer_matches_short_sequence() {
        let matches = find_kmer_matches(b"AC", b"ACGT", 3);
        assert!(matches.is_empty());
    }

    // -----------------------------------------------------------------------
    // lcsk_plusplus
    // -----------------------------------------------------------------------

    #[test]
    fn lcsk_empty_matches() {
        let (score, chain) = lcsk_plusplus(&[], 4);
        assert_eq!(score, 0);
        assert!(chain.is_empty());
    }

    #[test]
    fn lcsk_single_match() {
        let matches = vec![(0, 0)];
        let (score, chain) = lcsk_plusplus(&matches, 4);
        assert_eq!(score, 4);
        assert_eq!(chain.len(), 1);
        assert_eq!(chain[0], (0, 0));
    }

    #[test]
    fn lcsk_non_overlapping_chain() {
        // Two non-overlapping anchors with k=3: (0,0) and (5,5).
        // They don't overlap (0+3=3 <= 5), so score = 3 + 3 = 6.
        let matches = vec![(0, 0), (5, 5)];
        let (score, chain) = lcsk_plusplus(&matches, 3);
        assert_eq!(score, 6);
        assert_eq!(chain.len(), 2);
    }

    #[test]
    fn lcsk_adjacent_plusplus() {
        // Adjacent anchors: (0,0), (1,1), (2,2) with k=3.
        // With "++" extension: first anchor contributes 3, each subsequent
        // adjacent one contributes 1 → score = 3 + 1 + 1 = 5.
        // This corresponds to covering positions 0..5 in both sequences.
        let matches = vec![(0, 0), (1, 1), (2, 2)];
        let (score, chain) = lcsk_plusplus(&matches, 3);
        assert_eq!(score, 5);
        assert_eq!(chain.len(), 3);
    }

    // -----------------------------------------------------------------------
    // sparse_align — full pipeline
    // -----------------------------------------------------------------------

    #[test]
    fn identical_sequences() {
        let seq = b"ACGTACGTACGT";
        let result = sparse_align(seq, seq, 4).unwrap();
        assert_eq!(result.score, seq.len());
        assert_eq!(result.k, 4);
        assert!(!result.anchors.is_empty());
    }

    #[test]
    fn completely_different_sequences() {
        let a = b"AAAAAAAAAA";
        let b = b"CCCCCCCCCC";
        let result = sparse_align(a, b, 4).unwrap();
        assert_eq!(result.score, 0);
        assert!(result.anchors.is_empty());
    }

    #[test]
    fn known_lcs_case() {
        // "ACGTACGT" vs "ACGTXXACGT" with k=4.
        // Shared 4-mers from seq_a: ACGT(0), CGTA(1), GTAC(2), TACG(3), ACGT(4).
        // In seq_b: ACGT at pos 0 and 6, ACGT at pos 6.
        // Best chain should cover at least 8 matching bases.
        let a = b"ACGTACGT";
        let b = b"ACGTXXACGT";
        let result = sparse_align(a, b, 4).unwrap();
        assert!(
            result.score >= 8,
            "expected score >= 8, got {}",
            result.score
        );
    }

    #[test]
    fn single_kmer_match() {
        // Only one shared k-mer.
        let a = b"AAAACGTAAAA";
        let b = b"TTTTCGTTTTT";
        let result = sparse_align(a, b, 3).unwrap();
        assert!(result.score >= 3, "expected at least one k-mer match");
        assert!(!result.anchors.is_empty());
    }

    #[test]
    fn repetitive_sequences() {
        // Repetitive k-mers produce many match pairs; the algorithm should
        // still find a valid colinear chain.
        let a = b"ACACACACACAC";
        let b = b"ACACACACACAC";
        let result = sparse_align(a, b, 2).unwrap();
        assert_eq!(result.score, a.len());
        // All anchors should be colinear.
        for pair in result.anchors.windows(2) {
            assert!(pair[0].0 < pair[1].0, "anchors must be colinear on a-axis");
            assert!(pair[0].1 < pair[1].1, "anchors must be colinear on b-axis");
        }
    }

    #[test]
    fn error_on_k_zero() {
        let result = sparse_align(b"ACGT", b"ACGT", 0);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("k must be greater than 0"));
    }

    #[test]
    fn error_on_sequence_shorter_than_k() {
        // seq_a too short.
        let result = sparse_align(b"AC", b"ACGTACGT", 4);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("shorter than k"));

        // seq_b too short.
        let result = sparse_align(b"ACGTACGT", b"AC", 4);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("shorter than k"));
    }

    #[test]
    fn empty_matches_gives_zero() {
        // Sequences with no shared k-mers.
        let a = b"AAAA";
        let b = b"CCCC";
        let result = sparse_align(a, b, 2).unwrap();
        assert_eq!(result.score, 0);
        assert!(result.anchors.is_empty());
    }

    #[test]
    fn chain_is_colinear() {
        // Verify colinearity invariant on a non-trivial case.
        let a = b"ACGTACGTACGTACGT";
        let b = b"XXACGTXXACGTXX";
        let result = sparse_align(a, b, 4).unwrap();
        for pair in result.anchors.windows(2) {
            assert!(
                pair[0].0 < pair[1].0,
                "a-positions must strictly increase: {} < {}",
                pair[0].0,
                pair[1].0
            );
            assert!(
                pair[0].1 < pair[1].1,
                "b-positions must strictly increase: {} < {}",
                pair[0].1,
                pair[1].1
            );
        }
    }

    #[test]
    fn score_does_not_exceed_shorter_sequence() {
        let a = b"ACGTACGT";
        let b = b"ACGTACGTACGTACGT";
        let result = sparse_align(a, b, 3).unwrap();
        assert!(
            result.score <= a.len(),
            "score {} should not exceed shorter sequence length {}",
            result.score,
            a.len()
        );
    }

    #[test]
    fn k_equals_sequence_length() {
        let seq = b"ACGT";
        let result = sparse_align(seq, seq, 4).unwrap();
        assert_eq!(result.score, 4);
        assert_eq!(result.anchors.len(), 1);
        assert_eq!(result.anchors[0], (0, 0));
    }
}
