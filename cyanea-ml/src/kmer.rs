//! K-mer feature extraction from biological sequences.

use cyanea_core::{CyaneaError, Result, Summarizable};

use crate::encoding::Alphabet;

use std::collections::HashMap;

/// K-mer counter that extracts fixed-length subsequence counts from sequences.
#[derive(Debug, Clone)]
pub struct KmerCounter {
    k: usize,
}

impl KmerCounter {
    /// Create a new counter for the given k-mer size (1..=12).
    pub fn new(k: usize) -> Result<Self> {
        if k == 0 || k > 12 {
            return Err(CyaneaError::InvalidInput(
                "k must be between 1 and 12".into(),
            ));
        }
        Ok(Self { k })
    }

    /// Count all k-mers in the sequence.
    pub fn count_sequence(&self, seq: &[u8]) -> KmerCounts {
        let mut counts = HashMap::new();
        if seq.len() >= self.k {
            for window in seq.windows(self.k) {
                let kmer: Vec<u8> = window.iter().map(|b| b.to_ascii_uppercase()).collect();
                *counts.entry(kmer).or_insert(0usize) += 1;
            }
        }
        KmerCounts {
            counts,
            k: self.k,
        }
    }

    /// Count k-mers and return a normalized frequency vector for the given alphabet.
    ///
    /// The vector has length `alphabet.size().pow(k)` in lexicographic order.
    /// Each entry is count / total.
    pub fn count_normalized(&self, seq: &[u8], alphabet: Alphabet) -> Vec<f64> {
        let counts = self.count_sequence(seq);
        let total = counts.total();
        let mut freq = counts.to_frequency_vector(alphabet);
        if total > 0 {
            let t = total as f64;
            freq.iter_mut().for_each(|v| *v /= t);
        }
        freq
    }
}

/// A collection of k-mer counts from a sequence.
#[derive(Debug, Clone)]
pub struct KmerCounts {
    counts: HashMap<Vec<u8>, usize>,
    k: usize,
}

impl KmerCounts {
    /// Get the count for a specific k-mer.
    pub fn get(&self, kmer: &[u8]) -> usize {
        let upper: Vec<u8> = kmer.iter().map(|b| b.to_ascii_uppercase()).collect();
        self.counts.get(&upper).copied().unwrap_or(0)
    }

    /// Total number of k-mer occurrences.
    pub fn total(&self) -> usize {
        self.counts.values().sum()
    }

    /// Number of distinct k-mers observed.
    pub fn distinct(&self) -> usize {
        self.counts.len()
    }

    /// The k-mer size.
    pub fn k(&self) -> usize {
        self.k
    }

    /// Convert to a fixed-length frequency vector in lexicographic order.
    ///
    /// The vector has `alphabet.size().pow(k)` entries containing raw counts.
    pub fn to_frequency_vector(&self, alphabet: Alphabet) -> Vec<f64> {
        let syms = alphabet.symbols();
        let alpha_size = syms.len();
        let vec_len = alpha_size.pow(self.k as u32);
        let mut freq = vec![0.0; vec_len];

        for (kmer, &count) in &self.counts {
            if let Some(idx) = kmer_to_index(kmer, syms) {
                if idx < vec_len {
                    freq[idx] = count as f64;
                }
            }
        }
        freq
    }

    /// Iterate over observed (k-mer, count) pairs.
    pub fn iter(&self) -> impl Iterator<Item = (&[u8], usize)> {
        self.counts.iter().map(|(k, &v)| (k.as_slice(), v))
    }
}

impl Summarizable for KmerCounts {
    fn summary(&self) -> String {
        format!(
            "KmerCounts: k={}, {} distinct k-mers, {} total",
            self.k,
            self.distinct(),
            self.total(),
        )
    }
}

/// Map a k-mer to its lexicographic index within the given alphabet.
fn kmer_to_index(kmer: &[u8], symbols: &[u8]) -> Option<usize> {
    let base = symbols.len();
    let mut index = 0;
    for &b in kmer {
        let pos = symbols.iter().position(|&s| s == b)?;
        index = index * base + pos;
    }
    Some(index)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counter_new_valid() {
        assert!(KmerCounter::new(1).is_ok());
        assert!(KmerCounter::new(12).is_ok());
    }

    #[test]
    fn counter_new_invalid() {
        assert!(KmerCounter::new(0).is_err());
        assert!(KmerCounter::new(13).is_err());
    }

    #[test]
    fn count_sequence_basic() {
        let counter = KmerCounter::new(2).unwrap();
        let counts = counter.count_sequence(b"ACGT");
        assert_eq!(counts.get(b"AC"), 1);
        assert_eq!(counts.get(b"CG"), 1);
        assert_eq!(counts.get(b"GT"), 1);
        assert_eq!(counts.get(b"AA"), 0);
    }

    #[test]
    fn count_sequence_repeated() {
        let counter = KmerCounter::new(2).unwrap();
        let counts = counter.count_sequence(b"AAAA");
        assert_eq!(counts.get(b"AA"), 3);
        assert_eq!(counts.total(), 3);
        assert_eq!(counts.distinct(), 1);
    }

    #[test]
    fn count_sequence_short() {
        let counter = KmerCounter::new(3).unwrap();
        let counts = counter.count_sequence(b"AC");
        assert_eq!(counts.total(), 0);
        assert_eq!(counts.distinct(), 0);
    }

    #[test]
    fn count_sequence_lowercase() {
        let counter = KmerCounter::new(2).unwrap();
        let counts = counter.count_sequence(b"acgt");
        assert_eq!(counts.get(b"AC"), 1);
        assert_eq!(counts.get(b"ac"), 1); // get also uppercases
    }

    #[test]
    fn frequency_vector_length() {
        let counter = KmerCounter::new(2).unwrap();
        let counts = counter.count_sequence(b"ACGT");
        let freq = counts.to_frequency_vector(Alphabet::Dna);
        assert_eq!(freq.len(), 16); // 4^2
    }

    #[test]
    fn frequency_vector_values() {
        let counter = KmerCounter::new(1).unwrap();
        let counts = counter.count_sequence(b"AACG");
        let freq = counts.to_frequency_vector(Alphabet::Dna);
        assert_eq!(freq.len(), 4);
        assert_eq!(freq[0], 2.0); // A
        assert_eq!(freq[1], 1.0); // C
        assert_eq!(freq[2], 1.0); // G
        assert_eq!(freq[3], 0.0); // T
    }

    #[test]
    fn count_normalized_sums_to_one() {
        let counter = KmerCounter::new(2).unwrap();
        let freq = counter.count_normalized(b"ACGTACGT", Alphabet::Dna);
        let sum: f64 = freq.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn summary_format() {
        let counter = KmerCounter::new(3).unwrap();
        let counts = counter.count_sequence(b"ACGTACGT");
        let s = counts.summary();
        assert!(s.starts_with("KmerCounts: k=3"));
    }

    #[test]
    fn iter_pairs() {
        let counter = KmerCounter::new(1).unwrap();
        let counts = counter.count_sequence(b"AA");
        let pairs: Vec<_> = counts.iter().collect();
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].1, 2);
    }
}
