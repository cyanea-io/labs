//! 2-bit DNA encoding for compact storage.
//!
//! Packs DNA sequences (A, C, G, T only) into 2 bits per base,
//! achieving 4x compression over ASCII representation.
//!
//! Encoding: A=00, C=01, G=10, T=11

use cyanea_core::{CyaneaError, Result};

/// A DNA sequence stored in 2-bit packed representation.
///
/// Each byte holds 4 bases. The encoding is:
/// - A = 0b00
/// - C = 0b01
/// - G = 0b10
/// - T = 0b11
///
/// Bases are packed from the most significant bits: the first base occupies
/// bits 7..6, the second bits 5..4, and so on.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TwoBitSequence {
    data: Vec<u8>,
    len: usize,
}

/// Encode a single ASCII base to its 2-bit representation.
#[inline]
fn encode_base(b: u8) -> Result<u8> {
    match b {
        b'A' | b'a' => Ok(0b00),
        b'C' | b'c' => Ok(0b01),
        b'G' | b'g' => Ok(0b10),
        b'T' | b't' => Ok(0b11),
        _ => Err(CyaneaError::InvalidInput(format!(
            "invalid DNA base for 2-bit encoding: '{}' (0x{:02X}). Only A, C, G, T are supported",
            b as char, b
        ))),
    }
}

/// Decode a 2-bit value back to ASCII.
#[inline]
fn decode_base(bits: u8) -> u8 {
    match bits & 0b11 {
        0b00 => b'A',
        0b01 => b'C',
        0b10 => b'G',
        0b11 => b'T',
        _ => unreachable!(),
    }
}

impl TwoBitSequence {
    /// Encode an ASCII DNA sequence into 2-bit packed representation.
    ///
    /// Only unambiguous bases (A, C, G, T) are supported. Case-insensitive.
    ///
    /// # Errors
    ///
    /// Returns `CyaneaError::InvalidInput` if any base is not A, C, G, or T.
    ///
    /// # Example
    ///
    /// ```
    /// use cyanea_seq::twobit::TwoBitSequence;
    ///
    /// let seq = TwoBitSequence::encode(b"ACGT").unwrap();
    /// assert_eq!(seq.len(), 4);
    /// assert_eq!(seq.decode(), b"ACGT");
    /// ```
    pub fn encode(seq: &[u8]) -> Result<Self> {
        let len = seq.len();
        let num_bytes = (len + 3) / 4;
        let mut data = vec![0u8; num_bytes];

        for (i, &base) in seq.iter().enumerate() {
            let bits = encode_base(base)?;
            let byte_idx = i / 4;
            let bit_offset = 6 - (i % 4) * 2; // 6, 4, 2, 0
            data[byte_idx] |= bits << bit_offset;
        }

        Ok(Self { data, len })
    }

    /// Decode back to ASCII DNA bytes.
    ///
    /// Always produces uppercase A, C, G, T.
    pub fn decode(&self) -> Vec<u8> {
        let mut result = Vec::with_capacity(self.len);
        for i in 0..self.len {
            let byte_idx = i / 4;
            let bit_offset = 6 - (i % 4) * 2;
            let bits = (self.data[byte_idx] >> bit_offset) & 0b11;
            result.push(decode_base(bits));
        }
        result
    }

    /// Get the base at a specific position.
    ///
    /// Returns `None` if `index >= len`.
    pub fn get(&self, index: usize) -> Option<u8> {
        if index >= self.len {
            return None;
        }
        let byte_idx = index / 4;
        let bit_offset = 6 - (index % 4) * 2;
        let bits = (self.data[byte_idx] >> bit_offset) & 0b11;
        Some(decode_base(bits))
    }

    /// Number of bases in the sequence.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Extract a k-mer starting at `pos` as an integer encoding.
    ///
    /// Each base occupies 2 bits in the returned `u64`, with the first base
    /// in the most significant position. Maximum k is 32 (64 bits / 2 bits per base).
    ///
    /// Returns `None` if `pos + k > len` or `k > 32` or `k == 0`.
    pub fn kmer(&self, pos: usize, k: usize) -> Option<u64> {
        if k == 0 || k > 32 || pos + k > self.len {
            return None;
        }

        let mut value: u64 = 0;
        for i in 0..k {
            let idx = pos + i;
            let byte_idx = idx / 4;
            let bit_offset = 6 - (idx % 4) * 2;
            let bits = ((self.data[byte_idx] >> bit_offset) & 0b11) as u64;
            value = (value << 2) | bits;
        }

        Some(value)
    }

    /// Compute the bitwise complement (A<->T, C<->G).
    ///
    /// In 2-bit encoding, complement is a simple XOR with 0b11:
    /// - A (00) -> T (11)
    /// - C (01) -> G (10)
    /// - G (10) -> C (01)
    /// - T (11) -> A (00)
    pub fn complement(&self) -> Self {
        let mut data = self.data.clone();

        // XOR all full bytes
        for byte in &mut data {
            *byte ^= 0xFF;
        }

        // If the last byte has unused bits (padding), clear them back to zero
        let remainder = self.len % 4;
        if remainder != 0 && !data.is_empty() {
            let last = data.len() - 1;
            // Clear the unused low bits
            let used_bits = remainder * 2;
            let mask = !0u8 << (8 - used_bits);
            data[last] &= mask;
        }

        Self {
            data,
            len: self.len,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_decode_roundtrip() {
        let original = b"ACGTACGT";
        let encoded = TwoBitSequence::encode(original).unwrap();
        assert_eq!(encoded.decode(), original);
    }

    #[test]
    fn encode_decode_all_bases() {
        let original = b"AAAA";
        assert_eq!(TwoBitSequence::encode(original).unwrap().decode(), original.to_vec());

        let original = b"CCCC";
        assert_eq!(TwoBitSequence::encode(original).unwrap().decode(), original.to_vec());

        let original = b"GGGG";
        assert_eq!(TwoBitSequence::encode(original).unwrap().decode(), original.to_vec());

        let original = b"TTTT";
        assert_eq!(TwoBitSequence::encode(original).unwrap().decode(), original.to_vec());
    }

    #[test]
    fn encode_decode_non_multiple_of_four() {
        // 1 base
        let seq = TwoBitSequence::encode(b"A").unwrap();
        assert_eq!(seq.len(), 1);
        assert_eq!(seq.decode(), b"A");

        // 2 bases
        let seq = TwoBitSequence::encode(b"CG").unwrap();
        assert_eq!(seq.len(), 2);
        assert_eq!(seq.decode(), b"CG");

        // 3 bases
        let seq = TwoBitSequence::encode(b"ACT").unwrap();
        assert_eq!(seq.len(), 3);
        assert_eq!(seq.decode(), b"ACT");

        // 5 bases
        let seq = TwoBitSequence::encode(b"ACGTA").unwrap();
        assert_eq!(seq.len(), 5);
        assert_eq!(seq.decode(), b"ACGTA");
    }

    #[test]
    fn encode_case_insensitive() {
        let upper = TwoBitSequence::encode(b"ACGT").unwrap();
        let lower = TwoBitSequence::encode(b"acgt").unwrap();
        assert_eq!(upper, lower);
        assert_eq!(lower.decode(), b"ACGT");
    }

    #[test]
    fn encode_non_dna_error() {
        assert!(TwoBitSequence::encode(b"ACGN").is_err());
        assert!(TwoBitSequence::encode(b"ACGR").is_err());
        assert!(TwoBitSequence::encode(b"ACGU").is_err());
        assert!(TwoBitSequence::encode(b"XYZ").is_err());
    }

    #[test]
    fn empty_sequence() {
        let seq = TwoBitSequence::encode(b"").unwrap();
        assert_eq!(seq.len(), 0);
        assert!(seq.is_empty());
        assert_eq!(seq.decode(), Vec::<u8>::new());
        assert_eq!(seq.get(0), None);
    }

    #[test]
    fn get_individual_bases() {
        let seq = TwoBitSequence::encode(b"ACGT").unwrap();
        assert_eq!(seq.get(0), Some(b'A'));
        assert_eq!(seq.get(1), Some(b'C'));
        assert_eq!(seq.get(2), Some(b'G'));
        assert_eq!(seq.get(3), Some(b'T'));
        assert_eq!(seq.get(4), None);
    }

    #[test]
    fn get_bases_non_aligned() {
        let seq = TwoBitSequence::encode(b"TAGCAA").unwrap();
        assert_eq!(seq.get(0), Some(b'T'));
        assert_eq!(seq.get(1), Some(b'A'));
        assert_eq!(seq.get(2), Some(b'G'));
        assert_eq!(seq.get(3), Some(b'C'));
        assert_eq!(seq.get(4), Some(b'A'));
        assert_eq!(seq.get(5), Some(b'A'));
    }

    #[test]
    fn kmer_extraction() {
        // ACGT: A=00, C=01, G=10, T=11
        let seq = TwoBitSequence::encode(b"ACGT").unwrap();

        // 2-mer at pos 0: AC = 0b0001 = 1
        assert_eq!(seq.kmer(0, 2), Some(0b0001));
        // 2-mer at pos 1: CG = 0b0110 = 6
        assert_eq!(seq.kmer(1, 2), Some(0b0110));
        // 2-mer at pos 2: GT = 0b1011 = 11
        assert_eq!(seq.kmer(2, 2), Some(0b1011));

        // 4-mer at pos 0: ACGT = 0b00011011 = 27
        assert_eq!(seq.kmer(0, 4), Some(0b00011011));
    }

    #[test]
    fn kmer_edge_cases() {
        let seq = TwoBitSequence::encode(b"ACGT").unwrap();

        // k = 0
        assert_eq!(seq.kmer(0, 0), None);
        // k > 32
        assert_eq!(seq.kmer(0, 33), None);
        // pos + k > len
        assert_eq!(seq.kmer(3, 2), None);
        // k = 1
        assert_eq!(seq.kmer(0, 1), Some(0b00)); // A
        assert_eq!(seq.kmer(3, 1), Some(0b11)); // T
    }

    #[test]
    fn kmer_max_k32() {
        // 32 A's should give 0 (all 00 bits)
        let seq = TwoBitSequence::encode(&vec![b'A'; 32]).unwrap();
        assert_eq!(seq.kmer(0, 32), Some(0u64));

        // 32 T's should give all 1s in 64 bits
        let seq = TwoBitSequence::encode(&vec![b'T'; 32]).unwrap();
        assert_eq!(seq.kmer(0, 32), Some(u64::MAX));
    }

    #[test]
    fn complement_basic() {
        let seq = TwoBitSequence::encode(b"ACGT").unwrap();
        let comp = seq.complement();
        assert_eq!(comp.decode(), b"TGCA");
    }

    #[test]
    fn complement_all_same() {
        let seq = TwoBitSequence::encode(b"AAAA").unwrap();
        assert_eq!(seq.complement().decode(), b"TTTT");

        let seq = TwoBitSequence::encode(b"CCCC").unwrap();
        assert_eq!(seq.complement().decode(), b"GGGG");
    }

    #[test]
    fn complement_non_aligned() {
        let seq = TwoBitSequence::encode(b"ACG").unwrap();
        let comp = seq.complement();
        assert_eq!(comp.decode(), b"TGC");
        assert_eq!(comp.len(), 3);
    }

    #[test]
    fn complement_involution() {
        let seq = TwoBitSequence::encode(b"ACGTACGTAA").unwrap();
        let double_comp = seq.complement().complement();
        assert_eq!(seq.decode(), double_comp.decode());
    }

    #[test]
    fn complement_empty() {
        let seq = TwoBitSequence::encode(b"").unwrap();
        let comp = seq.complement();
        assert!(comp.is_empty());
        assert_eq!(comp.decode(), Vec::<u8>::new());
    }

    #[test]
    fn compact_storage() {
        // 8 bases should use exactly 2 bytes
        let seq = TwoBitSequence::encode(b"ACGTACGT").unwrap();
        assert_eq!(seq.data.len(), 2);

        // 9 bases should use 3 bytes
        let seq = TwoBitSequence::encode(b"ACGTACGTA").unwrap();
        assert_eq!(seq.data.len(), 3);
    }

    #[test]
    fn long_sequence_roundtrip() {
        let bases = b"ACGT";
        let long: Vec<u8> = (0..1000).map(|i| bases[i % 4]).collect();
        let seq = TwoBitSequence::encode(&long).unwrap();
        assert_eq!(seq.len(), 1000);
        assert_eq!(seq.decode(), long);
    }
}
