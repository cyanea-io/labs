//! Suffix array construction via the SA-IS algorithm.
//!
//! Implements the induced-sorting algorithm of Nong, Zhang & Chan (2009)
//! for O(n) suffix array construction. Provides binary search for
//! pattern matching in O(m log n).

/// A suffix array built from a text.
///
/// The suffix array stores the starting positions of all suffixes of the
/// input text, sorted in lexicographic order. This enables efficient
/// substring search via binary search.
#[derive(Debug, Clone)]
pub struct SuffixArray {
    sa: Vec<usize>,
}

// ---------------------------------------------------------------------------
// SA-IS implementation
// ---------------------------------------------------------------------------

/// Suffix type: S-type or L-type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SuffixType {
    S,
    L,
}

/// Classify each position as S-type or L-type.
///
/// A suffix at position i is S-type if text[i] < text[i+1], or
/// text[i] == text[i+1] and suffix[i+1] is S-type. The last
/// position (sentinel) is always S-type.
fn classify_suffixes(text: &[usize]) -> Vec<SuffixType> {
    let n = text.len();
    if n == 0 {
        return vec![];
    }
    let mut types = vec![SuffixType::S; n];
    // Last character (sentinel) is S-type
    types[n - 1] = SuffixType::S;

    if n > 1 {
        // Scan right to left
        for i in (0..n - 1).rev() {
            if text[i] > text[i + 1] {
                types[i] = SuffixType::L;
            } else if text[i] < text[i + 1] {
                types[i] = SuffixType::S;
            } else {
                types[i] = types[i + 1];
            }
        }
    }

    types
}

/// Check if position i is an LMS (leftmost S-type) suffix.
#[inline]
fn is_lms(types: &[SuffixType], i: usize) -> bool {
    i > 0 && types[i] == SuffixType::S && types[i - 1] == SuffixType::L
}

/// Get the bucket sizes for each character in the alphabet.
fn get_bucket_sizes(text: &[usize], alphabet_size: usize) -> Vec<usize> {
    let mut buckets = vec![0usize; alphabet_size];
    for &c in text {
        buckets[c] += 1;
    }
    buckets
}

/// Get the start positions of each bucket.
fn get_bucket_starts(bucket_sizes: &[usize]) -> Vec<usize> {
    let mut starts = vec![0usize; bucket_sizes.len()];
    let mut sum = 0;
    for (i, &size) in bucket_sizes.iter().enumerate() {
        starts[i] = sum;
        sum += size;
    }
    starts
}

/// Get the end positions (exclusive) of each bucket.
fn get_bucket_ends(bucket_sizes: &[usize]) -> Vec<usize> {
    let mut ends = vec![0usize; bucket_sizes.len()];
    let mut sum = 0;
    for (i, &size) in bucket_sizes.iter().enumerate() {
        sum += size;
        ends[i] = sum;
    }
    ends
}

/// Place LMS suffixes into their bucket tails.
fn place_lms(
    sa: &mut [usize],
    text: &[usize],
    _types: &[SuffixType],
    bucket_sizes: &[usize],
    lms_positions: &[usize],
) {
    let n = sa.len();
    let sentinel = n; // use n as "empty" marker
    for slot in sa.iter_mut() {
        *slot = sentinel;
    }

    let mut tails = get_bucket_ends(bucket_sizes);

    // Place LMS suffixes at the end of their buckets, in reverse order
    for &pos in lms_positions.iter().rev() {
        let c = text[pos];
        tails[c] -= 1;
        sa[tails[c]] = pos;
    }
}

/// Induce L-type suffixes from left to right.
fn induce_l(
    sa: &mut [usize],
    text: &[usize],
    types: &[SuffixType],
    bucket_sizes: &[usize],
) {
    let n = sa.len();
    let sentinel = n;
    let mut heads = get_bucket_starts(bucket_sizes);

    for i in 0..n {
        if sa[i] == sentinel || sa[i] == 0 {
            continue;
        }
        let j = sa[i] - 1;
        if types[j] == SuffixType::L {
            let c = text[j];
            sa[heads[c]] = j;
            heads[c] += 1;
        }
    }
}

/// Induce S-type suffixes from right to left.
fn induce_s(
    sa: &mut [usize],
    text: &[usize],
    types: &[SuffixType],
    bucket_sizes: &[usize],
) {
    let n = sa.len();
    let sentinel = n;
    let mut tails = get_bucket_ends(bucket_sizes);

    for i in (0..n).rev() {
        if sa[i] == sentinel || sa[i] == 0 {
            continue;
        }
        let j = sa[i] - 1;
        if types[j] == SuffixType::S {
            let c = text[j];
            tails[c] -= 1;
            sa[tails[c]] = j;
        }
    }
}

/// Check if two LMS substrings are equal.
fn lms_substrings_equal(
    text: &[usize],
    types: &[SuffixType],
    pos1: usize,
    pos2: usize,
) -> bool {
    let n = text.len();

    // Both must be valid positions
    if pos1 >= n || pos2 >= n {
        return false;
    }

    let mut i = 0;
    loop {
        let end1 = i > 0 && is_lms(types, pos1 + i);
        let end2 = i > 0 && is_lms(types, pos2 + i);

        if end1 && end2 {
            return true; // Both ended at same relative position
        }
        if end1 != end2 {
            return false; // One ended before the other
        }
        if text[pos1 + i] != text[pos2 + i] {
            return false;
        }
        if types[pos1 + i] != types[pos2 + i] {
            return false;
        }

        i += 1;

        if pos1 + i >= n || pos2 + i >= n {
            return false;
        }
    }
}

/// Core SA-IS algorithm on integer alphabet.
fn sais(text: &[usize], alphabet_size: usize) -> Vec<usize> {
    let n = text.len();

    if n == 0 {
        return vec![];
    }
    if n == 1 {
        return vec![0];
    }
    if n == 2 {
        if text[0] <= text[1] {
            return vec![0, 1];
        } else {
            return vec![1, 0];
        }
    }

    // Step 1: Classify suffixes
    let types = classify_suffixes(text);

    // Step 2: Find LMS positions
    let lms_positions: Vec<usize> = (0..n).filter(|&i| is_lms(&types, i)).collect();

    // Step 3: Get bucket sizes
    let bucket_sizes = get_bucket_sizes(text, alphabet_size);

    // Step 4: Initial placement and induced sorting to sort LMS substrings
    let mut sa = vec![0usize; n];
    place_lms(&mut sa, text, &types, &bucket_sizes, &lms_positions);
    induce_l(&mut sa, text, &types, &bucket_sizes);
    induce_s(&mut sa, text, &types, &bucket_sizes);

    // Step 5: Compact sorted LMS suffixes and assign names
    let sentinel = n;
    // Collect the sorted LMS suffixes from sa
    let sorted_lms: Vec<usize> = sa.iter().copied().filter(|&x| x != sentinel && is_lms(&types, x)).collect();

    // Assign names (ranks) to LMS substrings
    let mut names = vec![sentinel; n];
    let mut current_name = 0;

    if !sorted_lms.is_empty() {
        names[sorted_lms[0]] = current_name;
        for i in 1..sorted_lms.len() {
            if !lms_substrings_equal(text, &types, sorted_lms[i - 1], sorted_lms[i]) {
                current_name += 1;
            }
            names[sorted_lms[i]] = current_name;
        }
    }

    let num_names = current_name + 1;

    // Step 6: If names are not unique, recurse
    let lms_count = lms_positions.len();
    let reduced: Vec<usize> = lms_positions.iter().map(|&pos| names[pos]).collect();

    let sorted_lms_indices = if num_names < lms_count {
        // Recurse
        let sub_sa = sais(&reduced, num_names);
        sub_sa
    } else {
        // Names are unique, directly compute order
        let mut order = vec![0usize; lms_count];
        for (i, &name) in reduced.iter().enumerate() {
            order[name] = i;
        }
        order
    };

    // Step 7: Final induced sort using correctly ordered LMS suffixes
    let sorted_lms_positions: Vec<usize> = sorted_lms_indices
        .iter()
        .map(|&i| lms_positions[i])
        .collect();

    place_lms(&mut sa, text, &types, &bucket_sizes, &sorted_lms_positions);
    induce_l(&mut sa, text, &types, &bucket_sizes);
    induce_s(&mut sa, text, &types, &bucket_sizes);

    sa
}

impl SuffixArray {
    /// Build a suffix array from the given text using the SA-IS algorithm.
    ///
    /// Appends a sentinel character (value 0, lexicographically smallest)
    /// to ensure correct suffix ordering.
    ///
    /// Time complexity: O(n) where n is the length of the text.
    ///
    /// # Example
    ///
    /// ```
    /// use cyanea_seq::suffix::SuffixArray;
    ///
    /// let sa = SuffixArray::build(b"banana");
    /// let positions = sa.search(b"banana", b"ana");
    /// assert_eq!(positions.len(), 2);
    /// ```
    pub fn build(text: &[u8]) -> Self {
        if text.is_empty() {
            return Self { sa: vec![0] }; // Just the sentinel
        }

        // Map bytes to integer alphabet with sentinel = 0
        // We need to map the byte values to a compact integer alphabet.
        // Sentinel gets value 0, all other characters get values 1..=alphabet_size
        let mut chars_present = [false; 256];
        for &b in text {
            chars_present[b as usize] = true;
        }

        // Create mapping: byte value -> compact integer (1-indexed, 0 is sentinel)
        let mut mapping = [0usize; 256];
        let mut next_id = 1usize;
        for i in 0..256 {
            if chars_present[i] {
                mapping[i] = next_id;
                next_id += 1;
            }
        }
        let alphabet_size = next_id; // includes sentinel

        // Build integer text with sentinel appended
        let mut int_text: Vec<usize> = Vec::with_capacity(text.len() + 1);
        for &b in text {
            int_text.push(mapping[b as usize]);
        }
        int_text.push(0); // sentinel

        let sa = sais(&int_text, alphabet_size);

        Self { sa }
    }

    /// Search for all occurrences of `pattern` in `text`.
    ///
    /// The `text` must be the same text used to build the suffix array
    /// (without the sentinel). Returns positions in arbitrary order.
    ///
    /// Time complexity: O(m log n) where m is the pattern length and n is the text length.
    pub fn search(&self, text: &[u8], pattern: &[u8]) -> Vec<usize> {
        if pattern.is_empty() || text.is_empty() {
            return vec![];
        }

        let n = text.len();

        // The SA includes the sentinel position, so sa.len() == text.len() + 1
        // We need to compare suffixes against the pattern.

        // Binary search for the leftmost occurrence
        let left = self.lower_bound(text, pattern);
        let right = self.upper_bound(text, pattern);

        if left >= right {
            return vec![];
        }

        let mut positions: Vec<usize> = (left..right)
            .map(|i| self.sa[i])
            .filter(|&pos| pos < n) // exclude sentinel position
            .collect();
        positions.sort_unstable();
        positions
    }

    /// Find the leftmost SA index where the suffix >= pattern.
    fn lower_bound(&self, text: &[u8], pattern: &[u8]) -> usize {
        let mut lo = 0;
        let mut hi = self.sa.len();

        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let pos = self.sa[mid];
            let suffix = if pos < text.len() {
                &text[pos..]
            } else {
                &[] // sentinel
            };

            if compare_prefix(suffix, pattern) == std::cmp::Ordering::Less {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }

        lo
    }

    /// Find the first SA index where the suffix > pattern (as a prefix).
    fn upper_bound(&self, text: &[u8], pattern: &[u8]) -> usize {
        let mut lo = 0;
        let mut hi = self.sa.len();

        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let pos = self.sa[mid];
            let suffix = if pos < text.len() {
                &text[pos..]
            } else {
                &[] // sentinel
            };

            if compare_prefix(suffix, pattern) == std::cmp::Ordering::Greater {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }

        lo
    }

    /// Number of entries in the suffix array (text length + 1 for sentinel).
    pub fn len(&self) -> usize {
        self.sa.len()
    }
}

/// Compare a suffix against a pattern, treating the pattern as a prefix.
///
/// Returns Ordering::Equal if the suffix starts with the pattern.
fn compare_prefix(suffix: &[u8], pattern: &[u8]) -> std::cmp::Ordering {
    let len = suffix.len().min(pattern.len());
    match suffix[..len].cmp(&pattern[..len]) {
        std::cmp::Ordering::Equal => {
            if suffix.len() >= pattern.len() {
                std::cmp::Ordering::Equal
            } else {
                std::cmp::Ordering::Less
            }
        }
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_banana() {
        let text = b"banana";
        let sa = SuffixArray::build(text);
        // SA should have 7 entries (6 chars + sentinel)
        assert_eq!(sa.len(), 7);

        // Verify suffixes are sorted
        let suffixes: Vec<&[u8]> = sa
            .sa
            .iter()
            .map(|&pos| {
                if pos < text.len() {
                    &text[pos..]
                } else {
                    b"" as &[u8] // sentinel
                }
            })
            .collect();

        for i in 1..suffixes.len() {
            assert!(
                suffixes[i - 1] <= suffixes[i],
                "suffixes not sorted at {}: {:?} > {:?}",
                i,
                std::str::from_utf8(suffixes[i - 1]),
                std::str::from_utf8(suffixes[i])
            );
        }
    }

    #[test]
    fn search_banana_ana() {
        let text = b"banana";
        let sa = SuffixArray::build(text);
        let mut positions = sa.search(text, b"ana");
        positions.sort();
        assert_eq!(positions, vec![1, 3]);
    }

    #[test]
    fn search_banana_ban() {
        let text = b"banana";
        let sa = SuffixArray::build(text);
        let positions = sa.search(text, b"ban");
        assert_eq!(positions, vec![0]);
    }

    #[test]
    fn search_banana_a() {
        let text = b"banana";
        let sa = SuffixArray::build(text);
        let mut positions = sa.search(text, b"a");
        positions.sort();
        assert_eq!(positions, vec![1, 3, 5]);
    }

    #[test]
    fn search_no_match() {
        let text = b"banana";
        let sa = SuffixArray::build(text);
        let positions = sa.search(text, b"xyz");
        assert!(positions.is_empty());
    }

    #[test]
    fn search_full_text() {
        let text = b"banana";
        let sa = SuffixArray::build(text);
        let positions = sa.search(text, b"banana");
        assert_eq!(positions, vec![0]);
    }

    #[test]
    fn search_empty_pattern() {
        let text = b"banana";
        let sa = SuffixArray::build(text);
        let positions = sa.search(text, b"");
        assert!(positions.is_empty());
    }

    #[test]
    fn build_single_char() {
        let text = b"a";
        let sa = SuffixArray::build(text);
        assert_eq!(sa.len(), 2); // 'a' + sentinel
        let positions = sa.search(text, b"a");
        assert_eq!(positions, vec![0]);
    }

    #[test]
    fn build_repeated_chars() {
        let text = b"aaaa";
        let sa = SuffixArray::build(text);
        assert_eq!(sa.len(), 5);

        // Verify sorted
        let suffixes: Vec<&[u8]> = sa
            .sa
            .iter()
            .map(|&pos| {
                if pos < text.len() {
                    &text[pos..]
                } else {
                    b"" as &[u8]
                }
            })
            .collect();

        for i in 1..suffixes.len() {
            assert!(suffixes[i - 1] <= suffixes[i]);
        }

        let mut positions = sa.search(text, b"aa");
        positions.sort();
        assert_eq!(positions, vec![0, 1, 2]);
    }

    #[test]
    fn build_dna_sequence() {
        let text = b"ACGTACGT";
        let sa = SuffixArray::build(text);

        let mut positions = sa.search(text, b"ACG");
        positions.sort();
        assert_eq!(positions, vec![0, 4]);

        let positions = sa.search(text, b"TACG");
        assert_eq!(positions, vec![3]);
    }

    #[test]
    fn build_empty_text() {
        let text = b"";
        let sa = SuffixArray::build(text);
        assert_eq!(sa.len(), 1); // just sentinel
        let positions = sa.search(text, b"a");
        assert!(positions.is_empty());
    }

    #[test]
    fn search_longer_than_text() {
        let text = b"abc";
        let sa = SuffixArray::build(text);
        let positions = sa.search(text, b"abcdefg");
        assert!(positions.is_empty());
    }

    #[test]
    fn build_mississippi() {
        let text = b"mississippi";
        let sa = SuffixArray::build(text);

        // Verify sorted
        let suffixes: Vec<&[u8]> = sa
            .sa
            .iter()
            .map(|&pos| {
                if pos < text.len() {
                    &text[pos..]
                } else {
                    b"" as &[u8]
                }
            })
            .collect();
        for i in 1..suffixes.len() {
            assert!(
                suffixes[i - 1] <= suffixes[i],
                "not sorted at {}: {:?} vs {:?}",
                i,
                std::str::from_utf8(suffixes[i - 1]),
                std::str::from_utf8(suffixes[i])
            );
        }

        let mut positions = sa.search(text, b"issi");
        positions.sort();
        assert_eq!(positions, vec![1, 4]);

        let mut positions = sa.search(text, b"ss");
        positions.sort();
        assert_eq!(positions, vec![2, 5]);
    }
}
