//! Burrows-Wheeler Transform (BWT) for arbitrary byte sequences.
//!
//! Provides forward BWT construction via suffix array and inverse BWT
//! via LF-mapping reconstruction.

/// Burrows-Wheeler Transform of a byte string.
///
/// Constructed by appending a sentinel `$` to the input, building the suffix
/// array, and reading off the last column of the sorted rotation matrix.
#[derive(Debug, Clone)]
pub struct Bwt {
    /// The BWT string (same length as text + sentinel).
    bwt: Vec<u8>,
    /// Position of the sentinel character in the BWT.
    primary_index: usize,
}

impl Bwt {
    /// Build the BWT of `text`.
    ///
    /// A sentinel `$` (0x00) is appended internally and is lexicographically
    /// smaller than every input byte.
    ///
    /// # Example
    ///
    /// ```
    /// use cyanea_seq::bwt::Bwt;
    ///
    /// let bwt = Bwt::build(b"banana");
    /// assert_eq!(bwt.len(), 7); // 6 + sentinel
    /// ```
    pub fn build(text: &[u8]) -> Self {
        if text.is_empty() {
            return Self {
                bwt: vec![0u8], // just the sentinel
                primary_index: 0,
            };
        }

        // Append sentinel (0x00, sorts before everything)
        let mut augmented = Vec::with_capacity(text.len() + 1);
        augmented.extend_from_slice(text);
        augmented.push(0u8);

        let n = augmented.len();

        // Build suffix array via sort
        let mut sa: Vec<usize> = (0..n).collect();
        sa.sort_by(|&a, &b| augmented[a..].cmp(&augmented[b..]));

        // BWT: last column of sorted rotations = text[(sa[i] - 1) % n]
        let mut bwt = Vec::with_capacity(n);
        let mut primary_index = 0;
        for (i, &pos) in sa.iter().enumerate() {
            if pos == 0 {
                bwt.push(augmented[n - 1]);
                primary_index = i;
            } else {
                bwt.push(augmented[pos - 1]);
            }
        }

        Self { bwt, primary_index }
    }

    /// The BWT string as a byte slice.
    pub fn as_bytes(&self) -> &[u8] {
        &self.bwt
    }

    /// Position of the sentinel character (`$`) in the BWT.
    pub fn primary_index(&self) -> usize {
        self.primary_index
    }

    /// Length of the BWT (including sentinel).
    pub fn len(&self) -> usize {
        self.bwt.len()
    }

    /// Whether the BWT is empty (only possible for empty input + sentinel).
    pub fn is_empty(&self) -> bool {
        self.bwt.len() <= 1
    }

    /// Reconstruct the original text (without sentinel) from the BWT.
    ///
    /// Uses the LF-mapping: build C table (cumulative counts of characters
    /// smaller than c) and occurrence counts, then walk from the row
    /// containing the sentinel to reconstruct the text in reverse.
    ///
    /// # Example
    ///
    /// ```
    /// use cyanea_seq::bwt::Bwt;
    ///
    /// let text = b"banana";
    /// let bwt = Bwt::build(text);
    /// let recovered = bwt.invert();
    /// assert_eq!(recovered, text);
    /// ```
    pub fn invert(&self) -> Vec<u8> {
        let n = self.bwt.len();
        if n <= 1 {
            return vec![];
        }

        // Count occurrences of each byte value in the BWT
        let mut counts = [0usize; 256];
        for &b in &self.bwt {
            counts[b as usize] += 1;
        }

        // Build C table: c_table[c] = number of characters in BWT that are < c
        let mut c_table = [0usize; 256];
        let mut cumulative = 0;
        for i in 0..256 {
            c_table[i] = cumulative;
            cumulative += counts[i];
        }

        // Build occurrence array: occ[i] = number of occurrences of bwt[i]
        // in bwt[0..i] (exclusive of i, so occ[i] is the rank of bwt[i] among
        // its equal characters up to position i)
        let mut running = [0usize; 256];
        let mut occ = Vec::with_capacity(n);
        for &b in &self.bwt {
            occ.push(running[b as usize]);
            running[b as usize] += 1;
        }

        // LF-mapping: for row i, lf(i) = c_table[bwt[i]] + occ[i]
        // Walk from the sentinel row (primary_index) to reconstruct text
        let text_len = n - 1; // exclude sentinel
        let mut result = vec![0u8; text_len];
        let mut idx = self.primary_index;
        for i in (0..text_len).rev() {
            idx = c_table[self.bwt[idx] as usize] + occ[idx];
            result[i] = self.bwt[idx];
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_banana() {
        let text = b"banana";
        let bwt = Bwt::build(text);
        assert_eq!(bwt.invert(), text);
    }

    #[test]
    fn roundtrip_empty() {
        let bwt = Bwt::build(b"");
        assert!(bwt.is_empty());
        assert_eq!(bwt.len(), 1); // just sentinel
        assert_eq!(bwt.invert(), b"");
    }

    #[test]
    fn roundtrip_single_char() {
        let text = b"A";
        let bwt = Bwt::build(text);
        assert_eq!(bwt.len(), 2);
        assert_eq!(bwt.invert(), text);
    }

    #[test]
    fn roundtrip_repeated_chars() {
        let text = b"aaaa";
        let bwt = Bwt::build(text);
        assert_eq!(bwt.invert(), text);
    }

    #[test]
    fn roundtrip_dna() {
        let text = b"ACGTACGT";
        let bwt = Bwt::build(text);
        assert_eq!(bwt.invert(), text);
    }

    #[test]
    fn roundtrip_mississippi() {
        let text = b"mississippi";
        let bwt = Bwt::build(text);
        assert_eq!(bwt.invert(), text);
    }

    #[test]
    fn roundtrip_all_different() {
        let text = b"abcdefgh";
        let bwt = Bwt::build(text);
        assert_eq!(bwt.invert(), text);
    }

    #[test]
    fn primary_index_is_valid() {
        let text = b"banana";
        let bwt = Bwt::build(text);
        // The sentinel position should be within bounds
        assert!(bwt.primary_index() < bwt.len());
        // The character at primary_index in the BWT is the last char of the original text
        // (wrapping: the char before the sentinel in the rotation)
    }

    #[test]
    fn bwt_length_matches() {
        let text = b"ACGTACGT";
        let bwt = Bwt::build(text);
        assert_eq!(bwt.len(), text.len() + 1);
        assert_eq!(bwt.as_bytes().len(), text.len() + 1);
    }

    #[test]
    fn roundtrip_long_sequence() {
        let text = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let bwt = Bwt::build(text);
        assert_eq!(bwt.invert(), text);
    }
}
