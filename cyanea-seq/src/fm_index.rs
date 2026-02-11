//! FM-Index for fast substring search on DNA sequences.
//!
//! Built on top of the suffix array, the FM-Index provides O(m) pattern
//! matching via backward search using the Burrows-Wheeler Transform (BWT),
//! occurrence table, and C table.
//!
//! Only supports the DNA alphabet (A, C, G, T) plus the sentinel character '$'.

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

/// FM-Index for DNA sequences.
///
/// Supports efficient exact pattern matching using backward search.
///
/// # Fields
///
/// - `bwt` — the Burrows-Wheeler Transform of the text (with sentinel '$')
/// - `sa` — the full suffix array (for recovering positions)
/// - `occ` — occurrence table: `occ[i][c]` = number of occurrences of character c
///   in `bwt[0..=i]`
/// - `c_table` — `c_table[c]` = number of characters in the text that are
///   lexicographically smaller than character c (among A, C, G, T)
#[derive(Debug, Clone)]
pub struct FmIndex {
    bwt: Vec<u8>,
    sa: Vec<usize>,
    occ: Vec<[usize; 4]>,
    c_table: [usize; 4],
}

impl FmIndex {
    /// Build an FM-Index from a DNA sequence.
    ///
    /// The input `text` should contain only A, C, G, T characters.
    /// A sentinel '$' is appended internally.
    ///
    /// # Example
    ///
    /// ```
    /// use cyanea_seq::fm_index::FmIndex;
    ///
    /// let fm = FmIndex::build(b"ACGTACGT");
    /// let positions = fm.search(b"ACG");
    /// assert_eq!(positions.len(), 2);
    /// ```
    pub fn build(text: &[u8]) -> Self {
        // Append sentinel
        let mut augmented = Vec::with_capacity(text.len() + 1);
        augmented.extend_from_slice(text);
        augmented.push(b'$');

        let n = augmented.len();

        // Build suffix array for the augmented text (with '$' sentinel).
        // We use a direct sort since we need the raw SA positions for the
        // augmented text with the '$' character as the actual sentinel.
        let sa_raw = build_sa_for_fm(&augmented);

        // Compute BWT from SA: bwt[i] = augmented[sa[i] - 1] or sentinel if sa[i] == 0
        let mut bwt = Vec::with_capacity(n);
        for &pos in &sa_raw {
            if pos == 0 {
                bwt.push(augmented[n - 1]); // wrap around: character before position 0 is the last char
            } else {
                bwt.push(augmented[pos - 1]);
            }
        }

        // Build occurrence table
        let mut occ = Vec::with_capacity(n);
        let mut counts = [0usize; 4];
        for &b in &bwt {
            if let Some(idx) = base_to_idx(b) {
                counts[idx] += 1;
            }
            occ.push(counts);
        }

        // Build C table: c_table[c] = number of chars lexicographically < c
        // In our alphabet ordering: $ < A < C < G < T
        // So c_table[A] = 1 (just $), c_table[C] = 1 + count(A), etc.
        let total_a = counts[0]; // total A's in BWT = total A's in text
        let total_c = counts[1];
        let total_g = counts[2];
        // total_t = counts[3]

        let c_table = [
            1,                                  // A: just '$' is smaller
            1 + total_a,                        // C: '$' + all A's
            1 + total_a + total_c,              // G: '$' + A's + C's
            1 + total_a + total_c + total_g,    // T: '$' + A's + C's + G's
        ];

        Self {
            bwt,
            sa: sa_raw,
            occ,
            c_table,
        }
    }

    /// Search for all occurrences of `pattern` in the indexed text.
    ///
    /// Returns the starting positions (0-indexed, in the original text without
    /// sentinel) of each occurrence. Positions are returned in sorted order.
    ///
    /// Uses backward search: O(m) for the search phase, where m is the pattern
    /// length, plus O(k) to look up the k matching positions from the SA.
    pub fn search(&self, pattern: &[u8]) -> Vec<usize> {
        if pattern.is_empty() {
            return vec![];
        }

        let (top, bottom) = self.backward_search(pattern);

        if top > bottom {
            return vec![];
        }

        let text_len = self.bwt.len() - 1; // exclude sentinel
        let mut positions: Vec<usize> = (top..=bottom)
            .map(|i| self.sa[i])
            .filter(|&pos| pos < text_len)
            .collect();
        positions.sort_unstable();
        positions
    }

    /// Count the number of occurrences of `pattern` without returning positions.
    ///
    /// More efficient than `search` when only the count is needed, as it
    /// avoids the SA lookups.
    pub fn count(&self, pattern: &[u8]) -> usize {
        if pattern.is_empty() {
            return 0;
        }

        let (top, bottom) = self.backward_search(pattern);

        if top > bottom {
            0
        } else {
            bottom - top + 1
        }
    }

    /// Perform backward search, returning the SA interval [top, bottom].
    ///
    /// If no match is found, returns a state where top > bottom.
    fn backward_search(&self, pattern: &[u8]) -> (usize, usize) {
        let n = self.bwt.len();

        // Start with the last character of the pattern
        let last = pattern[pattern.len() - 1];
        let c_idx = match base_to_idx(last) {
            Some(idx) => idx,
            None => return (1, 0), // non-DNA char, no match
        };

        let mut top = self.c_table[c_idx];
        let mut bottom = if c_idx < 3 {
            self.c_table[c_idx + 1] - 1
        } else {
            n - 1
        };

        if top > bottom {
            return (1, 0);
        }

        // Extend backward through the pattern
        for i in (0..pattern.len() - 1).rev() {
            let c = pattern[i];
            let c_idx = match base_to_idx(c) {
                Some(idx) => idx,
                None => return (1, 0),
            };

            let occ_before_top = if top == 0 { 0 } else { self.occ[top - 1][c_idx] };
            let occ_at_bottom = self.occ[bottom][c_idx];

            if occ_before_top >= occ_at_bottom {
                return (1, 0); // no match
            }

            top = self.c_table[c_idx] + occ_before_top;
            bottom = self.c_table[c_idx] + occ_at_bottom - 1;
        }

        (top, bottom)
    }
}

/// Build a suffix array for the FM-Index (direct construction on augmented text).
///
/// This builds a simple SA by sorting suffix indices. For the FM-Index we need
/// the SA of the exact augmented text with '$' as the actual sentinel character.
fn build_sa_for_fm(text: &[u8]) -> Vec<usize> {
    let mut sa: Vec<usize> = (0..text.len()).collect();
    sa.sort_by(|&a, &b| text[a..].cmp(&text[b..]));
    sa
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_simple() {
        let fm = FmIndex::build(b"ACGT");
        assert_eq!(fm.bwt.len(), 5); // 4 + sentinel
    }

    #[test]
    fn search_single_occurrence() {
        let text = b"ACGTACGT";
        let fm = FmIndex::build(text);
        let positions = fm.search(b"TACG");
        assert_eq!(positions, vec![3]);
    }

    #[test]
    fn search_multiple_occurrences() {
        let text = b"ACGTACGT";
        let fm = FmIndex::build(text);
        let mut positions = fm.search(b"ACG");
        positions.sort();
        assert_eq!(positions, vec![0, 4]);
    }

    #[test]
    fn search_full_text() {
        let text = b"ACGTACGT";
        let fm = FmIndex::build(text);
        let positions = fm.search(b"ACGTACGT");
        assert_eq!(positions, vec![0]);
    }

    #[test]
    fn search_single_base() {
        let text = b"ACGTACGT";
        let fm = FmIndex::build(text);
        let mut positions = fm.search(b"A");
        positions.sort();
        assert_eq!(positions, vec![0, 4]);
    }

    #[test]
    fn search_no_match() {
        let text = b"ACGTACGT";
        let fm = FmIndex::build(text);
        let positions = fm.search(b"AAA");
        assert!(positions.is_empty());
    }

    #[test]
    fn search_empty_pattern() {
        let text = b"ACGTACGT";
        let fm = FmIndex::build(text);
        let positions = fm.search(b"");
        assert!(positions.is_empty());
    }

    #[test]
    fn count_matches_search() {
        let text = b"ACGTACGT";
        let fm = FmIndex::build(text);

        assert_eq!(fm.count(b"ACG"), fm.search(b"ACG").len());
        assert_eq!(fm.count(b"ACG"), 2);
        assert_eq!(fm.count(b"TACG"), 1);
        assert_eq!(fm.count(b"AAA"), 0);
    }

    #[test]
    fn count_empty_pattern() {
        let text = b"ACGTACGT";
        let fm = FmIndex::build(text);
        assert_eq!(fm.count(b""), 0);
    }

    #[test]
    fn search_non_existent() {
        let text = b"AAACCCGGGTTT";
        let fm = FmIndex::build(text);
        let positions = fm.search(b"ACGT");
        assert!(positions.is_empty());
    }

    #[test]
    fn search_at_boundaries() {
        let text = b"ACGTTTTT";
        let fm = FmIndex::build(text);

        // Pattern at the start
        let positions = fm.search(b"ACG");
        assert_eq!(positions, vec![0]);

        // Pattern at the end (overlapping matches at positions 3 and 4)
        let positions = fm.search(b"TTTT");
        assert_eq!(positions, vec![3, 4]);
    }

    #[test]
    fn build_single_base() {
        let text = b"A";
        let fm = FmIndex::build(text);
        let positions = fm.search(b"A");
        assert_eq!(positions, vec![0]);
        assert_eq!(fm.count(b"A"), 1);
    }

    #[test]
    fn search_repeated_sequence() {
        let text = b"AAAA";
        let fm = FmIndex::build(text);
        let mut positions = fm.search(b"AA");
        positions.sort();
        assert_eq!(positions, vec![0, 1, 2]);
        assert_eq!(fm.count(b"AA"), 3);
    }

    #[test]
    fn search_longer_dna() {
        let text = b"ACGTACGTACGTACGT";
        let fm = FmIndex::build(text);

        let mut positions = fm.search(b"ACGT");
        positions.sort();
        assert_eq!(positions, vec![0, 4, 8, 12]);
        assert_eq!(fm.count(b"ACGT"), 4);
    }

    #[test]
    fn count_non_dna_pattern() {
        let text = b"ACGTACGT";
        let fm = FmIndex::build(text);
        // 'N' is not in the DNA alphabet for FM-Index
        assert_eq!(fm.count(b"N"), 0);
        assert_eq!(fm.count(b"ACN"), 0);
    }
}
