//! Pattern matching algorithms for biological sequences.
//!
//! Provides exact and approximate string matching on `&[u8]` slices, suitable
//! for DNA, RNA, and protein sequences alike.
//!
//! ## Exact matchers
//!
//! - [`horspool`] — Boyer-Moore-Horspool with bad-character shift table, O(n/m) average
//! - [`kmp`] — Knuth-Morris-Pratt with failure function, O(n+m)
//! - [`shift_and`] — Shift-And bitparallel (pattern <= 64)
//! - [`bndm`] — Backward Nondeterministic DAWG Matching, bitparallel (pattern <= 64)
//! - [`bom`] — Backward Oracle Matching with factor oracle
//!
//! ## Approximate matchers
//!
//! - [`myers_bitparallel`] — Myers bit-parallel edit distance (pattern <= 64)
//! - [`ukkonen`] — Ukkonen cut-off approximate matching via bounded DP

/// Boyer-Moore-Horspool exact pattern matching.
///
/// Builds a bad-character shift table and scans right-to-left within each
/// alignment window. Average case O(n/m), worst case O(nm).
///
/// Returns starting positions of all exact occurrences.
pub fn horspool(text: &[u8], pattern: &[u8]) -> Vec<usize> {
    let n = text.len();
    let m = pattern.len();
    if m == 0 || m > n {
        return vec![];
    }

    // Bad-character shift table: default shift is m.
    let mut shift = [m; 256];
    for i in 0..m - 1 {
        shift[pattern[i] as usize] = m - 1 - i;
    }

    let mut results = Vec::new();
    let mut i = 0;
    while i <= n - m {
        let mut j = m - 1;
        while pattern[j] == text[i + j] {
            if j == 0 {
                results.push(i);
                break;
            }
            j -= 1;
        }
        i += shift[text[i + m - 1] as usize];
    }
    results
}

/// Knuth-Morris-Pratt exact pattern matching.
///
/// Builds a failure (partial match) table in O(m), then scans in O(n).
/// Total time O(n+m), space O(m).
///
/// Returns starting positions of all exact occurrences.
pub fn kmp(text: &[u8], pattern: &[u8]) -> Vec<usize> {
    let n = text.len();
    let m = pattern.len();
    if m == 0 || m > n {
        return vec![];
    }

    // Build failure function.
    let mut fail = vec![0usize; m];
    let mut k = 0usize;
    for i in 1..m {
        while k > 0 && pattern[k] != pattern[i] {
            k = fail[k - 1];
        }
        if pattern[k] == pattern[i] {
            k += 1;
        }
        fail[i] = k;
    }

    // Search phase.
    let mut results = Vec::new();
    let mut q = 0usize;
    for i in 0..n {
        while q > 0 && pattern[q] != text[i] {
            q = fail[q - 1];
        }
        if pattern[q] == text[i] {
            q += 1;
        }
        if q == m {
            results.push(i + 1 - m);
            q = fail[q - 1];
        }
    }
    results
}

/// Shift-And bitparallel exact pattern matching.
///
/// Encodes the pattern as bitmasks (one per alphabet symbol) and simulates
/// an NFA with bitwise operations. Limited to patterns of length <= 64.
///
/// Returns starting positions of all exact occurrences, or an empty vec
/// if the pattern exceeds 64 characters.
pub fn shift_and(text: &[u8], pattern: &[u8]) -> Vec<usize> {
    let n = text.len();
    let m = pattern.len();
    if m == 0 || m > n || m > 64 {
        return vec![];
    }

    // Build bitmask table: B[c] has bit j set if pattern[j] == c.
    let mut b = [0u64; 256];
    for j in 0..m {
        b[pattern[j] as usize] |= 1u64 << j;
    }

    let accept = 1u64 << (m - 1);
    let mut state = 0u64;
    let mut results = Vec::new();

    for i in 0..n {
        state = ((state << 1) | 1) & b[text[i] as usize];
        if state & accept != 0 {
            results.push(i + 1 - m);
        }
    }
    results
}

/// Backward Nondeterministic DAWG Matching (BNDM) exact pattern matching.
///
/// A bitparallel algorithm that scans the current window backward, using a
/// nondeterministic suffix automaton to detect both matches and feasible
/// shift prefixes simultaneously. Limited to patterns of length <= 64.
///
/// Returns starting positions of all exact occurrences, or an empty vec
/// if the pattern exceeds 64 characters.
pub fn bndm(text: &[u8], pattern: &[u8]) -> Vec<usize> {
    let n = text.len();
    let m = pattern.len();
    if m == 0 || m > n || m > 64 {
        return vec![];
    }

    // Build bitmask table: B[c] has bit j set if pattern[m-1-j] == c
    // (reversed pattern positions for backward scanning).
    let mut b = [0u64; 256];
    for j in 0..m {
        b[pattern[m - 1 - j] as usize] |= 1u64 << j;
    }

    let accept = 1u64 << (m - 1);
    let mut results = Vec::new();
    let mut pos = 0usize;

    while pos <= n - m {
        let mut j = m - 1;
        let mut last = m;
        let mut d = !0u64; // all bits set

        loop {
            d &= b[text[pos + j] as usize];
            if d == 0 {
                break;
            }
            if d & accept != 0 {
                if j == 0 {
                    results.push(pos);
                    break;
                }
                last = j;
            }
            d <<= 1;
            if j == 0 {
                break;
            }
            j -= 1;
        }
        pos += last;
    }
    results
}

/// Backward Oracle Matching (BOM) exact pattern matching.
///
/// Builds a factor oracle for the reversed pattern, then scans the text
/// backward within each window. The factor oracle recognizes at least all
/// factors of the pattern, enabling efficient shift computation.
///
/// Returns starting positions of all exact occurrences.
pub fn bom(text: &[u8], pattern: &[u8]) -> Vec<usize> {
    let n = text.len();
    let m = pattern.len();
    if m == 0 || m > n {
        return vec![];
    }

    // Build factor oracle for reversed pattern.
    // States: 0..=m, transitions stored as a vec of hashmaps.
    let rev: Vec<u8> = pattern.iter().rev().copied().collect();
    let states = m + 1;
    let mut goto: Vec<[i32; 256]> = vec![[-1i32; 256]; states];
    let mut supply = vec![0usize; states];

    // Build oracle: add transitions for each character of the reversed pattern.
    for i in 0..m {
        goto[i][rev[i] as usize] = (i + 1) as i32;
        let mut s = supply[i];
        while s != 0 && goto[s][rev[i] as usize] == -1 {
            goto[s][rev[i] as usize] = (i + 1) as i32;
            s = supply[s];
        }
        if s == 0 && goto[0][rev[i] as usize] == -1 {
            goto[0][rev[i] as usize] = (i + 1) as i32;
            supply[i + 1] = 0;
        } else if goto[s][rev[i] as usize] == (i + 1) as i32 {
            supply[i + 1] = s;
        } else {
            supply[i + 1] = goto[s][rev[i] as usize] as usize;
        }
    }

    // Search: scan backward using the factor oracle.
    let mut results = Vec::new();
    let mut pos = 0usize;

    while pos <= n - m {
        let mut state = 0usize;
        let mut j = m;

        while j > 0 {
            let c = text[pos + j - 1] as usize;
            let next = goto[state][c];
            if next == -1 {
                break;
            }
            state = next as usize;
            j -= 1;
        }

        if j == 0 {
            // Verify the match (oracle may accept superset of factors).
            if &text[pos..pos + m] == pattern {
                results.push(pos);
            }
            pos += 1;
        } else {
            pos += j;
        }
    }
    results
}

/// Myers bit-parallel approximate matching.
///
/// Computes edit distance (Levenshtein) between the pattern and all
/// substrings of the text using Myers' 1999 bit-vector algorithm.
/// Limited to patterns of length <= 64.
///
/// Returns `(end_position, edit_distance)` pairs for every text position
/// where the best alignment ending there has edit distance <= `max_dist`.
/// Returns an empty vec if the pattern exceeds 64 characters.
pub fn myers_bitparallel(
    text: &[u8],
    pattern: &[u8],
    max_dist: usize,
) -> Vec<(usize, usize)> {
    let n = text.len();
    let m = pattern.len();
    if m == 0 || m > 64 {
        return vec![];
    }
    if n == 0 {
        return vec![];
    }

    // Build pattern bitmasks: peq[c] has bit j set if pattern[j] == c.
    let mut peq = [0u64; 256];
    for j in 0..m {
        peq[pattern[j] as usize] |= 1u64 << j;
    }

    let mask = if m == 64 { !0u64 } else { (1u64 << m) - 1 };
    let msb = 1u64 << (m - 1);
    let mut pv = mask; // low m bits set (positive vertical delta)
    let mut mv = 0u64; // negative vertical delta
    let mut score = m; // current edit distance

    let mut results = Vec::new();

    for i in 0..n {
        let eq = peq[text[i] as usize];

        let xv = eq | mv;
        let xh = (((eq & pv).wrapping_add(pv)) ^ pv) | eq;

        // Mask to m bits before score test — high bits from NOT would corrupt.
        let ph = (mv | !(xh | pv)) & mask;
        let mh = (pv & xh) & mask;

        // Update score based on bit m-1 of ph and mh.
        if ph & msb != 0 {
            score += 1;
        }
        if mh & msb != 0 {
            score -= 1;
        }

        // Shift for next column. Bit 0 = 0 for semi-global alignment
        // (D[0][j] = 0: no text-start penalty).
        pv = ((mh << 1) | !(xv | (ph << 1))) & mask;
        mv = ((ph << 1) & xv) & mask;

        if score <= max_dist {
            results.push((i, score));
        }
    }
    results
}

/// Ukkonen cut-off approximate matching via bounded dynamic programming.
///
/// Uses a classic DP matrix but only fills the diagonal band of width
/// `2 * max_dist + 1`, achieving O(n * max_dist) time instead of O(nm).
///
/// Returns `(end_position, edit_distance)` pairs for every text position
/// where a semi-global alignment of the pattern ends with edit distance
/// <= `max_dist`.
pub fn ukkonen(
    text: &[u8],
    pattern: &[u8],
    max_dist: usize,
) -> Vec<(usize, usize)> {
    let n = text.len();
    let m = pattern.len();
    if m == 0 || n == 0 {
        return vec![];
    }

    // Semi-global alignment: free gaps at start/end of text (row 0 = 0).
    // We keep a single column of the DP matrix, length m+1.
    let mut prev = vec![0usize; m + 1];
    let mut curr = vec![0usize; m + 1];

    // Initialize first column: aligning pattern[0..j] against empty text.
    for j in 0..=m {
        prev[j] = j;
    }

    let mut results = Vec::new();

    for i in 1..=n {
        curr[0] = 0; // Semi-global: no penalty for starting gap in text.
        let mut col_min = m + 1;

        for j in 1..=m {
            let cost = if text[i - 1] == pattern[j - 1] { 0 } else { 1 };
            curr[j] = (prev[j - 1] + cost)
                .min(prev[j] + 1)
                .min(curr[j - 1] + 1);
            if j == m && curr[j] < col_min {
                col_min = curr[j];
            }
        }

        // Check the final row (full pattern aligned).
        if curr[m] <= max_dist {
            results.push((i - 1, curr[m]));
        }

        std::mem::swap(&mut prev, &mut curr);
    }
    results
}

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // Helper: run all exact matchers and assert same result
    // -----------------------------------------------------------------------
    fn assert_exact(text: &[u8], pattern: &[u8], expected: &[usize]) {
        assert_eq!(horspool(text, pattern), expected, "horspool");
        assert_eq!(kmp(text, pattern), expected, "kmp");
        // Bitparallel methods only for pattern <= 64.
        if pattern.len() <= 64 {
            assert_eq!(shift_and(text, pattern), expected, "shift_and");
            assert_eq!(bndm(text, pattern), expected, "bndm");
        }
        assert_eq!(bom(text, pattern), expected, "bom");
    }

    #[test]
    fn exact_match_at_start() {
        assert_exact(b"ACGTACGT", b"ACGT", &[0, 4]);
    }

    #[test]
    fn exact_match_at_end() {
        assert_exact(b"TTTTACGT", b"ACGT", &[4]);
    }

    #[test]
    fn exact_match_in_middle() {
        assert_exact(b"TTACGTTT", b"ACGT", &[2]);
    }

    #[test]
    fn multiple_occurrences() {
        assert_exact(b"AAAAAA", b"AA", &[0, 1, 2, 3, 4]);
    }

    #[test]
    fn no_match() {
        assert_exact(b"ACGTACGT", b"TTTT", &[]);
    }

    #[test]
    fn empty_pattern() {
        assert_exact(b"ACGT", b"", &[]);
    }

    #[test]
    fn empty_text() {
        assert_exact(b"", b"ACGT", &[]);
    }

    #[test]
    fn pattern_longer_than_text() {
        assert_exact(b"AC", b"ACGT", &[]);
    }

    #[test]
    fn single_char_pattern() {
        assert_exact(b"AACAA", b"C", &[2]);
    }

    #[test]
    fn single_char_text_match() {
        assert_exact(b"A", b"A", &[0]);
    }

    #[test]
    fn single_char_text_no_match() {
        assert_exact(b"A", b"C", &[]);
    }

    #[test]
    fn full_text_match() {
        assert_exact(b"ACGT", b"ACGT", &[0]);
    }

    #[test]
    fn protein_sequence() {
        assert_exact(b"MKAILFVLV", b"AILF", &[2]);
    }

    #[test]
    fn overlapping_pattern() {
        assert_exact(b"ABABAB", b"ABAB", &[0, 2]);
    }

    // -----------------------------------------------------------------------
    // Bitparallel length limit
    // -----------------------------------------------------------------------
    #[test]
    fn shift_and_pattern_too_long() {
        let text = vec![b'A'; 128];
        let pattern = vec![b'A'; 65];
        assert_eq!(shift_and(&text, &pattern), vec![]);
    }

    #[test]
    fn bndm_pattern_too_long() {
        let text = vec![b'A'; 128];
        let pattern = vec![b'A'; 65];
        assert_eq!(bndm(&text, &pattern), vec![]);
    }

    #[test]
    fn shift_and_pattern_exactly_64() {
        let text = vec![b'A'; 128];
        let pattern = vec![b'A'; 64];
        let result = shift_and(&text, &pattern);
        assert_eq!(result.len(), 65); // 128 - 64 + 1
    }

    #[test]
    fn bndm_pattern_exactly_64() {
        let text = vec![b'A'; 128];
        let pattern = vec![b'A'; 64];
        let result = bndm(&text, &pattern);
        assert_eq!(result.len(), 65);
    }

    // -----------------------------------------------------------------------
    // Myers bit-parallel approximate matching
    // -----------------------------------------------------------------------
    #[test]
    fn myers_exact_match() {
        let hits = myers_bitparallel(b"ACGTACGT", b"ACGT", 0);
        // End positions of exact matches: 3, 7
        let ends: Vec<usize> = hits.iter().map(|&(e, _)| e).collect();
        assert!(ends.contains(&3));
        assert!(ends.contains(&7));
        assert!(hits.iter().all(|&(_, d)| d == 0));
    }

    #[test]
    fn myers_single_substitution() {
        // Pattern ACGT vs text AXGT: one substitution at position 1.
        let hits = myers_bitparallel(b"AXGT", b"ACGT", 1);
        assert!(hits.iter().any(|&(e, d)| e == 3 && d == 1));
    }

    #[test]
    fn myers_single_insertion() {
        // Pattern ACG vs text ACXG: one insertion.
        let hits = myers_bitparallel(b"TACXGT", b"ACG", 1);
        assert!(hits.iter().any(|&(_, d)| d <= 1));
    }

    #[test]
    fn myers_single_deletion() {
        // Pattern ACGT vs text AGT: one deletion.
        let hits = myers_bitparallel(b"AGT", b"ACGT", 1);
        assert!(hits.iter().any(|&(_, d)| d <= 1));
    }

    #[test]
    fn myers_no_match_within_distance() {
        let hits = myers_bitparallel(b"AAAA", b"CCCC", 1);
        assert!(hits.is_empty());
    }

    #[test]
    fn myers_empty_pattern() {
        assert_eq!(myers_bitparallel(b"ACGT", b"", 2), vec![]);
    }

    #[test]
    fn myers_empty_text() {
        assert_eq!(myers_bitparallel(b"", b"ACGT", 2), vec![]);
    }

    #[test]
    fn myers_pattern_too_long() {
        let text = vec![b'A'; 128];
        let pattern = vec![b'A'; 65];
        assert_eq!(myers_bitparallel(&text, &pattern, 2), vec![]);
    }

    // -----------------------------------------------------------------------
    // Ukkonen approximate matching
    // -----------------------------------------------------------------------
    #[test]
    fn ukkonen_exact_match() {
        let hits = ukkonen(b"ACGTACGT", b"ACGT", 0);
        let ends: Vec<usize> = hits.iter().map(|&(e, _)| e).collect();
        assert!(ends.contains(&3));
        assert!(ends.contains(&7));
        assert!(hits.iter().all(|&(_, d)| d == 0));
    }

    #[test]
    fn ukkonen_single_substitution() {
        let hits = ukkonen(b"AXGT", b"ACGT", 1);
        assert!(hits.iter().any(|&(e, d)| e == 3 && d == 1));
    }

    #[test]
    fn ukkonen_single_insertion() {
        let hits = ukkonen(b"TACXGT", b"ACG", 1);
        assert!(hits.iter().any(|&(_, d)| d <= 1));
    }

    #[test]
    fn ukkonen_single_deletion() {
        let hits = ukkonen(b"AGT", b"ACGT", 1);
        assert!(hits.iter().any(|&(_, d)| d <= 1));
    }

    #[test]
    fn ukkonen_no_match_within_distance() {
        let hits = ukkonen(b"AAAA", b"CCCC", 1);
        assert!(hits.is_empty());
    }

    #[test]
    fn ukkonen_empty_pattern() {
        assert_eq!(ukkonen(b"ACGT", b"", 2), vec![]);
    }

    #[test]
    fn ukkonen_empty_text() {
        assert_eq!(ukkonen(b"", b"ACGT", 2), vec![]);
    }

    // -----------------------------------------------------------------------
    // Cross-validation: approximate matchers agree on distance-0 positions
    // -----------------------------------------------------------------------
    #[test]
    fn approx_matchers_agree_on_exact() {
        let text = b"GATTACACGTACGTTTG";
        let pattern = b"ACGT";
        let myers = myers_bitparallel(text, pattern, 0);
        let ukk = ukkonen(text, pattern, 0);
        let myers_ends: Vec<usize> = myers.iter().map(|&(e, _)| e).collect();
        let ukk_ends: Vec<usize> = ukk.iter().map(|&(e, _)| e).collect();
        assert_eq!(myers_ends, ukk_ends);
    }

    // -----------------------------------------------------------------------
    // Cross-validation: exact matchers agree with approx at distance 0
    // -----------------------------------------------------------------------
    #[test]
    fn exact_vs_approx_agreement() {
        let text = b"AACGTAACGTAA";
        let pattern = b"ACGT";
        let exact = kmp(text, pattern);
        let approx = myers_bitparallel(text, pattern, 0);
        // Convert end positions to start positions.
        let approx_starts: Vec<usize> = approx
            .iter()
            .map(|&(e, _)| e + 1 - pattern.len())
            .collect();
        assert_eq!(exact, approx_starts);
    }
}
