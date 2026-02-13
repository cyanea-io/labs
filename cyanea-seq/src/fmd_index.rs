//! Bidirectional FM-Index (FMD-Index) for strand-aware DNA search.
//!
//! The FMD-Index indexes the concatenation of a DNA sequence with its reverse
//! complement, enabling efficient bidirectional extension and super-maximal
//! exact match (SMEM) enumeration. Both strands are searched simultaneously.
//!
//! This is a self-contained implementation that does not depend on the
//! [`crate::fm_index`] module.

/// Map a DNA base to its index in the occurrence/C tables.
///
/// A=0, C=1, G=2, T=3. Returns `None` for non-ACGT characters.
#[inline]
fn base_to_idx(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Return the complement of a DNA base.
#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

/// Compute the reverse complement of a DNA sequence.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Build a suffix array by direct sorting of suffix indices.
fn build_sa(text: &[u8]) -> Vec<usize> {
    let mut sa: Vec<usize> = (0..text.len()).collect();
    sa.sort_by(|&a, &b| text[a..].cmp(&text[b..]));
    sa
}

/// Build all FM-index components from an augmented text (must already contain sentinel).
///
/// Returns (bwt, sa, occ, c_table).
fn build_fm(text: &[u8]) -> (Vec<u8>, Vec<usize>, Vec<[usize; 4]>, [usize; 4]) {
    let n = text.len();
    let sa = build_sa(text);

    // BWT: bwt[i] = text[sa[i] - 1], wrapping around for position 0.
    let mut bwt = Vec::with_capacity(n);
    for &pos in &sa {
        if pos == 0 {
            bwt.push(text[n - 1]);
        } else {
            bwt.push(text[pos - 1]);
        }
    }

    // Occurrence table: occ[i][c] = count of character c in bwt[0..=i].
    let mut occ = Vec::with_capacity(n);
    let mut counts = [0usize; 4];
    for &b in &bwt {
        if let Some(idx) = base_to_idx(b) {
            counts[idx] += 1;
        }
        occ.push(counts);
    }

    // C table: c_table[c] = number of characters lexicographically smaller than c.
    // Lexicographic order with sentinels: '#' < '$' < 'A' < 'C' < 'G' < 'T'
    // We count how many sentinel characters appear (they sort before A).
    let num_sentinels = n - counts[0] - counts[1] - counts[2] - counts[3];
    let c_table = [
        num_sentinels,
        num_sentinels + counts[0],
        num_sentinels + counts[0] + counts[1],
        num_sentinels + counts[0] + counts[1] + counts[2],
    ];

    (bwt, sa, occ, c_table)
}

/// Count occurrences of character `c` in `bwt[0..pos]`.
#[inline]
fn occ_count(occ: &[[usize; 4]], pos: usize, c: usize) -> usize {
    if pos == 0 { 0 } else { occ[pos - 1][c] }
}

/// Perform a single backward-search step on an FM-index.
///
/// Given a current SA interval [lo, hi) and a character, compute the new interval.
/// Returns an empty interval (lo >= hi) if there are no matches.
#[inline]
fn lf_step(
    c_table: &[usize; 4],
    occ: &[[usize; 4]],
    lo: usize,
    hi: usize,
    c: usize,
) -> (usize, usize) {
    if lo >= hi {
        return (0, 0);
    }
    let occ_lo = occ_count(occ, lo, c);
    let occ_hi = occ_count(occ, hi, c);
    (c_table[c] + occ_lo, c_table[c] + occ_hi)
}

/// Bidirectional interval in the FM-index.
///
/// Represents matching intervals in both the forward and reverse BWT simultaneously.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BiInterval {
    /// Lower bound in the forward index.
    pub lower: usize,
    /// Size of the interval (same in both directions).
    pub size: usize,
    /// Lower bound in the reverse index.
    pub lower_rev: usize,
}

impl BiInterval {
    /// Returns `true` if the interval is empty (matches nothing).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }
}

/// Bidirectional FM-Index (FMD-Index) for strand-aware DNA search.
///
/// Indexes the concatenation `text#reverse_complement(text)$` to enable
/// efficient bidirectional extension for finding super-maximal exact matches (SMEMs).
#[derive(Debug, Clone)]
pub struct FmdIndex {
    /// Forward FM-index components (built on the concatenated text).
    bwt: Vec<u8>,
    sa: Vec<usize>,
    occ: Vec<[usize; 4]>,
    c_table: [usize; 4],
    /// Reverse FM-index components (built on the reversed concatenated text).
    _bwt_rev: Vec<u8>,
    _sa_rev: Vec<usize>,
    occ_rev: Vec<[usize; 4]>,
    c_table_rev: [usize; 4],
    /// Length of original text (without sentinels or reverse complement).
    text_len: usize,
}

impl FmdIndex {
    /// Build an FMD-Index from a DNA sequence.
    ///
    /// The input `seq` should contain only A, C, G, T characters.
    /// Internally constructs the concatenation `seq#rev_comp(seq)$` and builds
    /// both the forward and reverse FM-indexes.
    ///
    /// # Example
    ///
    /// ```
    /// use cyanea_seq::fmd_index::FmdIndex;
    ///
    /// let fmd = FmdIndex::new(b"ACGT");
    /// let interval = fmd.backward_search(b"ACG");
    /// assert!(!interval.is_empty());
    /// ```
    pub fn new(seq: &[u8]) -> Self {
        let text_len = seq.len();
        let rc = reverse_complement(seq);

        // Forward text: seq # rev_comp $
        let mut forward = Vec::with_capacity(seq.len() + rc.len() + 2);
        forward.extend_from_slice(seq);
        forward.push(b'#');
        forward.extend_from_slice(&rc);
        forward.push(b'$');

        // Reverse text: reverse of the forward concatenation
        let reversed: Vec<u8> = forward.iter().rev().copied().collect();

        let (bwt, sa, occ, c_table) = build_fm(&forward);
        let (bwt_rev, _sa_rev, occ_rev, c_table_rev) = build_fm(&reversed);

        Self {
            bwt,
            sa,
            occ,
            c_table,
            _bwt_rev: bwt_rev,
            _sa_rev,
            occ_rev,
            c_table_rev,
            text_len,
        }
    }

    /// Initialize a bi-interval for a single character.
    ///
    /// Returns the interval representing all suffixes starting with `c` in both
    /// the forward and reverse indexes.
    pub fn init_interval(&self, c: u8) -> BiInterval {
        let ci = match base_to_idx(c) {
            Some(idx) => idx,
            None => return BiInterval { lower: 0, size: 0, lower_rev: 0 },
        };

        let n = self.bwt.len();

        // Forward interval for character c.
        let lo = self.c_table[ci];
        let hi = if ci < 3 { self.c_table[ci + 1] } else { n };
        let size = hi - lo;

        if size == 0 {
            return BiInterval { lower: 0, size: 0, lower_rev: 0 };
        }

        // Reverse interval for the same character.
        let lo_rev = self.c_table_rev[ci];

        BiInterval {
            lower: lo,
            size,
            lower_rev: lo_rev,
        }
    }

    /// Extend the bi-interval by prepending character `c` (backward extension).
    ///
    /// This performs a standard backward search step on the forward index and
    /// updates the reverse interval accordingly.
    pub fn extend_backward(&self, interval: &BiInterval, c: u8) -> BiInterval {
        if interval.is_empty() {
            return *interval;
        }
        let ci = match base_to_idx(c) {
            Some(idx) => idx,
            None => return BiInterval { lower: 0, size: 0, lower_rev: 0 },
        };

        let lo = interval.lower;
        let hi = lo + interval.size;

        // Backward step on forward index.
        let (new_lo, new_hi) = lf_step(
            &self.c_table, &self.occ, lo, hi, ci,
        );
        let new_size = if new_hi > new_lo { new_hi - new_lo } else { 0 };

        if new_size == 0 {
            return BiInterval { lower: 0, size: 0, lower_rev: 0 };
        }

        // Update reverse interval: count occurrences of characters ranked below ci
        // in the forward BWT within [lo, hi).
        let lo_rev = interval.lower_rev;
        let mut offset = 0;
        for a in 0..ci {
            offset += occ_count(&self.occ, hi, a) - occ_count(&self.occ, lo, a);
        }
        let new_lo_rev = lo_rev + offset;

        BiInterval {
            lower: new_lo,
            size: new_size,
            lower_rev: new_lo_rev,
        }
    }

    /// Extend the bi-interval by appending character `c` (forward extension).
    ///
    /// This performs a backward search step on the reverse index and updates
    /// the forward interval accordingly.
    pub fn extend_forward(&self, interval: &BiInterval, c: u8) -> BiInterval {
        if interval.is_empty() {
            return *interval;
        }
        let ci = match base_to_idx(c) {
            Some(idx) => idx,
            None => return BiInterval { lower: 0, size: 0, lower_rev: 0 },
        };

        let lo_rev = interval.lower_rev;
        let hi_rev = lo_rev + interval.size;

        // Backward step on reverse index.
        let (new_lo_rev, new_hi_rev) = lf_step(
            &self.c_table_rev, &self.occ_rev, lo_rev, hi_rev, ci,
        );
        let new_size = if new_hi_rev > new_lo_rev {
            new_hi_rev - new_lo_rev
        } else {
            0
        };

        if new_size == 0 {
            return BiInterval { lower: 0, size: 0, lower_rev: 0 };
        }

        // Update forward interval: count occurrences of characters ranked below ci
        // in the reverse BWT within [lo_rev, hi_rev).
        let lo = interval.lower;
        let mut offset = 0;
        for a in 0..ci {
            offset += occ_count(&self.occ_rev, hi_rev, a) - occ_count(&self.occ_rev, lo_rev, a);
        }
        let new_lo = lo + offset;

        BiInterval {
            lower: new_lo,
            size: new_size,
            lower_rev: new_lo_rev,
        }
    }

    /// Retrieve text positions from a bi-interval.
    ///
    /// Returns all suffix array positions within `[0, text_len)` — i.e., positions
    /// in the original sequence (not its reverse complement).
    pub fn locate(&self, interval: &BiInterval) -> Vec<usize> {
        if interval.is_empty() {
            return vec![];
        }
        let lo = interval.lower;
        let hi = lo + interval.size;
        let mut positions: Vec<usize> = (lo..hi)
            .map(|i| self.sa[i])
            .filter(|&pos| pos < self.text_len)
            .collect();
        positions.sort_unstable();
        positions
    }

    /// Full backward search for an exact pattern.
    ///
    /// Processes the pattern from right to left using the forward index, returning
    /// the final bi-interval. The interval may be empty if the pattern is not found.
    pub fn backward_search(&self, pattern: &[u8]) -> BiInterval {
        if pattern.is_empty() {
            return BiInterval { lower: 0, size: 0, lower_rev: 0 };
        }

        let mut interval = self.init_interval(pattern[pattern.len() - 1]);
        if interval.is_empty() {
            return interval;
        }

        for i in (0..pattern.len() - 1).rev() {
            interval = self.extend_backward(&interval, pattern[i]);
            if interval.is_empty() {
                return interval;
            }
        }

        interval
    }

    /// Find all Super-Maximal Exact Matches (SMEMs) of `query` against the index.
    ///
    /// Returns a vector of `(query_start, query_end, interval)` triples, where
    /// `query_start..query_end` is the matching substring of the query. An SMEM is
    /// a maximal exact match that is not contained within any longer match at the
    /// same query position.
    ///
    /// Only matches of length `>= min_len` are returned.
    ///
    /// # Algorithm
    ///
    /// For each starting position in the query, extend backward as far as possible,
    /// recording the maximal intervals. Then filter to retain only super-maximal
    /// matches — those not properly contained in another match.
    pub fn smems(
        &self,
        query: &[u8],
        min_len: usize,
    ) -> Vec<(usize, usize, BiInterval)> {
        if query.is_empty() || min_len == 0 {
            return vec![];
        }

        let qlen = query.len();
        let mut mems: Vec<(usize, usize, BiInterval)> = Vec::new();

        // For each position, find the longest match ending at or beyond that position
        // by extending forward from a backward-maximal seed.
        let mut pos = 0;
        while pos < qlen {
            // Start with the character at `pos` and extend forward.
            let mut interval = self.init_interval(query[pos]);
            if interval.is_empty() {
                pos += 1;
                continue;
            }

            // Extend forward as far as possible.
            let mut end = pos + 1;
            while end < qlen {
                let next = self.extend_forward(&interval, query[end]);
                if next.is_empty() {
                    break;
                }
                interval = next;
                end += 1;
            }

            // Now try extending backward from [pos, end) to find longer matches.
            let start = pos;
            let mut best_interval = interval;
            let mut best_start = start;

            if start > 0 {
                let mut back_interval = interval;
                let mut s = start;
                while s > 0 {
                    let prev = self.extend_backward(&back_interval, query[s - 1]);
                    if prev.is_empty() {
                        break;
                    }
                    back_interval = prev;
                    s -= 1;
                }
                if s < best_start {
                    best_start = s;
                    best_interval = back_interval;
                }
            }

            let match_len = end - best_start;
            if match_len >= min_len {
                mems.push((best_start, end, best_interval));
            }

            // Advance past the current forward-maximal position.
            pos = end;
        }

        // Filter to super-maximal: remove any match contained in another.
        Self::filter_supermaximal(&mut mems);
        mems
    }

    /// Remove matches that are properly contained within another match.
    fn filter_supermaximal(mems: &mut Vec<(usize, usize, BiInterval)>) {
        if mems.len() <= 1 {
            return;
        }
        // Sort by start ascending, then by end descending (longer first).
        mems.sort_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)));

        let mut keep = vec![true; mems.len()];
        let mut max_end = 0;
        for i in 0..mems.len() {
            if mems[i].1 <= max_end {
                // This match is contained within a previously seen one.
                keep[i] = false;
            }
            if mems[i].1 > max_end {
                max_end = mems[i].1;
            }
        }

        let mut write = 0;
        for read in 0..mems.len() {
            if keep[read] {
                mems[write] = mems[read];
                write += 1;
            }
        }
        mems.truncate(write);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_bases() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'G'), b'C');
    }

    #[test]
    fn revcomp() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAC"), b"GTTT");
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn build_simple() {
        let fmd = FmdIndex::new(b"ACGT");
        assert_eq!(fmd.text_len, 4);
        // Concatenated text: ACGT#ACGT$ (len 10)
        assert_eq!(fmd.bwt.len(), 10);
        assert_eq!(fmd._bwt_rev.len(), 10);
    }

    #[test]
    fn init_interval_valid() {
        let fmd = FmdIndex::new(b"ACGTACGT");
        for &base in &[b'A', b'C', b'G', b'T'] {
            let iv = fmd.init_interval(base);
            assert!(!iv.is_empty(), "init_interval for {} should be non-empty", base as char);
            // Each base appears at least twice in the original + twice in rev comp.
            assert!(iv.size >= 4);
        }
    }

    #[test]
    fn init_interval_invalid() {
        let fmd = FmdIndex::new(b"ACGT");
        let iv = fmd.init_interval(b'N');
        assert!(iv.is_empty());
    }

    #[test]
    fn backward_search_exact_match() {
        let fmd = FmdIndex::new(b"ACGTACGT");
        let iv = fmd.backward_search(b"ACGT");
        assert!(!iv.is_empty());
        let positions = fmd.locate(&iv);
        // "ACGT" appears at 0 and 4 in the original text.
        assert_eq!(positions, vec![0, 4]);
    }

    #[test]
    fn backward_search_no_match() {
        let fmd = FmdIndex::new(b"ACGTACGT");
        let iv = fmd.backward_search(b"AAA");
        assert!(iv.is_empty());
        assert!(fmd.locate(&iv).is_empty());
    }

    #[test]
    fn backward_search_single_base() {
        let fmd = FmdIndex::new(b"ACGTACGT");
        let iv = fmd.backward_search(b"A");
        assert!(!iv.is_empty());
        let positions = fmd.locate(&iv);
        assert_eq!(positions, vec![0, 4]);
    }

    #[test]
    fn backward_search_full_text() {
        let fmd = FmdIndex::new(b"ACGTACGT");
        let iv = fmd.backward_search(b"ACGTACGT");
        assert!(!iv.is_empty());
        let positions = fmd.locate(&iv);
        assert_eq!(positions, vec![0]);
    }

    #[test]
    fn backward_search_empty_pattern() {
        let fmd = FmdIndex::new(b"ACGT");
        let iv = fmd.backward_search(b"");
        assert!(iv.is_empty());
    }

    #[test]
    fn forward_extension() {
        let fmd = FmdIndex::new(b"ACGTACGT");
        // Start with 'A', then extend forward with 'C', 'G', 'T'.
        let iv = fmd.init_interval(b'A');
        assert!(!iv.is_empty());

        let iv2 = fmd.extend_forward(&iv, b'C');
        assert!(!iv2.is_empty());

        let iv3 = fmd.extend_forward(&iv2, b'G');
        assert!(!iv3.is_empty());

        let iv4 = fmd.extend_forward(&iv3, b'T');
        assert!(!iv4.is_empty());

        // "ACGT" should match at positions 0 and 4.
        let positions = fmd.locate(&iv4);
        assert_eq!(positions, vec![0, 4]);
    }

    #[test]
    fn forward_and_backward_agree() {
        let fmd = FmdIndex::new(b"ACGTACGT");

        // Build "ACG" via backward search.
        let iv_back = fmd.backward_search(b"ACG");

        // Build "ACG" via forward extension: A -> AC -> ACG.
        let iv_fwd = fmd.init_interval(b'A');
        let iv_fwd = fmd.extend_forward(&iv_fwd, b'C');
        let iv_fwd = fmd.extend_forward(&iv_fwd, b'G');

        // Both should locate the same positions.
        let pos_back = fmd.locate(&iv_back);
        let pos_fwd = fmd.locate(&iv_fwd);
        assert_eq!(pos_back, pos_fwd);
    }

    #[test]
    fn locate_filters_to_original_text() {
        // Positions in the reverse complement half should be filtered out.
        let fmd = FmdIndex::new(b"AACC");
        let iv = fmd.backward_search(b"AA");
        let positions = fmd.locate(&iv);
        // "AA" appears at position 0 in the original text.
        // rev_comp("AACC") = "GGTT", which does not contain "AA".
        assert_eq!(positions, vec![0]);
    }

    #[test]
    fn reverse_complement_symmetry() {
        // Searching for a pattern should also find its reverse complement
        // in the reverse complement half (but locate filters to original text).
        let seq = b"ACGTAAAA";
        let fmd = FmdIndex::new(seq);

        // "ACGT" appears at position 0 in the original text.
        let iv = fmd.backward_search(b"ACGT");
        let positions = fmd.locate(&iv);
        assert!(positions.contains(&0));

        // The reverse complement of the sequence contains "ACGT" at the end:
        // rev_comp("ACGTAAAA") = "TTTTACGT"
        // So "ACGT" appears in the concatenated text at position in the rev_comp half too,
        // but locate only returns positions < text_len.
        // The full interval should have size >= 2 (one from each strand).
        assert!(iv.size >= 2);
    }

    #[test]
    fn smems_simple() {
        // Index: "ACGTACGT", query a substring.
        let fmd = FmdIndex::new(b"ACGTACGT");
        let smems = fmd.smems(b"ACGT", 2);
        assert!(!smems.is_empty());
        // The entire query "ACGT" should be one SMEM.
        let (start, end, _) = smems[0];
        assert_eq!(start, 0);
        assert_eq!(end, 4);
    }

    #[test]
    fn smems_min_len_filter() {
        let fmd = FmdIndex::new(b"ACGTACGT");
        // With min_len = 10, no SMEM of length >= 10 can exist in a query of length 4.
        let smems = fmd.smems(b"ACGT", 10);
        assert!(smems.is_empty());
    }

    #[test]
    fn smems_multiple_matches() {
        // Construct text with two distinct regions.
        let fmd = FmdIndex::new(b"ACGTTTTTGGGG");
        // Query that spans two regions.
        let smems = fmd.smems(b"ACGTTTTT", 3);
        assert!(!smems.is_empty());
        // Should find at least one SMEM covering most of the query.
        let max_len = smems.iter().map(|(s, e, _)| e - s).max().unwrap();
        assert!(max_len >= 4);
    }

    #[test]
    fn smems_empty_query() {
        let fmd = FmdIndex::new(b"ACGT");
        let smems = fmd.smems(b"", 1);
        assert!(smems.is_empty());
    }

    #[test]
    fn smems_no_match_in_query() {
        // If the query has no match at all (impossible for single bases in DNA,
        // but a non-DNA character won't match).
        let fmd = FmdIndex::new(b"ACGT");
        let smems = fmd.smems(b"N", 1);
        assert!(smems.is_empty());
    }

    #[test]
    fn smems_supermaximal_property() {
        // Verify that returned SMEMs are not contained in each other.
        let fmd = FmdIndex::new(b"ACGTACGTACGT");
        let smems = fmd.smems(b"ACGTACGT", 2);
        for i in 0..smems.len() {
            for j in 0..smems.len() {
                if i == j {
                    continue;
                }
                let (si, ei, _) = smems[i];
                let (sj, ej, _) = smems[j];
                // No SMEM should be properly contained in another.
                assert!(!(si >= sj && ei <= ej && (si, ei) != (sj, ej)),
                    "SMEM ({}, {}) is contained in ({}, {})", si, ei, sj, ej);
            }
        }
    }

    #[test]
    fn locate_repeated_sequence() {
        let fmd = FmdIndex::new(b"AAAA");
        let iv = fmd.backward_search(b"AA");
        let positions = fmd.locate(&iv);
        // "AA" appears at positions 0, 1, 2 in "AAAA".
        assert_eq!(positions, vec![0, 1, 2]);
    }

    #[test]
    fn backward_extension() {
        let fmd = FmdIndex::new(b"ACGTACGT");
        // Start with 'T', extend backward with G, C, A to spell "ACGT".
        let iv = fmd.init_interval(b'T');
        let iv = fmd.extend_backward(&iv, b'G');
        let iv = fmd.extend_backward(&iv, b'C');
        let iv = fmd.extend_backward(&iv, b'A');
        let positions = fmd.locate(&iv);
        assert_eq!(positions, vec![0, 4]);
    }

    #[test]
    fn bi_interval_empty() {
        let iv = BiInterval { lower: 0, size: 0, lower_rev: 0 };
        assert!(iv.is_empty());
        let iv2 = BiInterval { lower: 5, size: 3, lower_rev: 2 };
        assert!(!iv2.is_empty());
    }
}
