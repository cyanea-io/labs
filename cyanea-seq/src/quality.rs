//! Phred quality scores for sequencing reads.
//!
//! [`QualityScores`] stores decoded Phred quality values (not raw ASCII).
//! Supports both Phred+33 (Illumina 1.8+) and Phred+64 (older Illumina)
//! encodings via [`PhredEncoding`].

use cyanea_core::{CyaneaError, Result, Scored};

/// Quality score encoding scheme.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum PhredEncoding {
    /// Phred+33 (Sanger / Illumina 1.8+). Most common modern encoding.
    Phred33,
    /// Phred+64 (Illumina 1.3–1.7).
    Phred64,
}

impl PhredEncoding {
    fn offset(self) -> u8 {
        match self {
            PhredEncoding::Phred33 => 33,
            PhredEncoding::Phred64 => 64,
        }
    }
}

/// Decoded Phred quality scores.
///
/// Stores raw Phred values (typically 0–41), not ASCII-encoded bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct QualityScores {
    scores: Vec<u8>,
}

impl QualityScores {
    /// Create from raw Phred quality values (already decoded).
    pub fn from_raw(scores: Vec<u8>) -> Self {
        Self { scores }
    }

    /// Decode from ASCII-encoded quality bytes.
    pub fn from_ascii(ascii: &[u8], encoding: PhredEncoding) -> Result<Self> {
        let offset = encoding.offset();
        let mut scores = Vec::with_capacity(ascii.len());
        for (i, &b) in ascii.iter().enumerate() {
            if b < offset {
                return Err(CyaneaError::InvalidInput(format!(
                    "quality byte {} (0x{:02X}) at position {} is below {:?} offset {}",
                    b as char, b, i, encoding, offset
                )));
            }
            scores.push(b - offset);
        }
        Ok(Self { scores })
    }

    /// Encode back to ASCII using the given encoding.
    pub fn to_ascii(&self, encoding: PhredEncoding) -> Vec<u8> {
        let offset = encoding.offset();
        self.scores.iter().map(|&q| q + offset).collect()
    }

    /// Number of quality scores.
    pub fn len(&self) -> usize {
        self.scores.len()
    }

    /// Whether the quality scores are empty.
    pub fn is_empty(&self) -> bool {
        self.scores.is_empty()
    }

    /// Raw quality scores as a slice.
    pub fn as_slice(&self) -> &[u8] {
        &self.scores
    }

    /// Mean quality score.
    pub fn mean(&self) -> f64 {
        if self.scores.is_empty() {
            return 0.0;
        }
        let sum: u64 = self.scores.iter().map(|&q| q as u64).sum();
        sum as f64 / self.scores.len() as f64
    }

    /// Minimum quality score.
    pub fn min(&self) -> Option<u8> {
        self.scores.iter().copied().min()
    }

    /// Maximum quality score.
    pub fn max(&self) -> Option<u8> {
        self.scores.iter().copied().max()
    }

    /// Fraction of scores at or above the given threshold.
    pub fn fraction_above(&self, threshold: u8) -> f64 {
        if self.scores.is_empty() {
            return 0.0;
        }
        let count = self.scores.iter().filter(|&&q| q >= threshold).count();
        count as f64 / self.scores.len() as f64
    }

    /// Error probability for a single Phred quality score: 10^(-Q/10).
    pub fn error_probability(phred: u8) -> f64 {
        10.0_f64.powf(-(phred as f64) / 10.0)
    }
}

impl Scored for QualityScores {
    fn score(&self) -> f64 {
        self.mean()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_ascii_phred33() {
        // '!' = 33 → Q0, 'I' = 73 → Q40
        let ascii = b"!I";
        let q = QualityScores::from_ascii(ascii, PhredEncoding::Phred33).unwrap();
        assert_eq!(q.as_slice(), &[0, 40]);
    }

    #[test]
    fn roundtrip_ascii() {
        let original = vec![0, 10, 20, 30, 40];
        let q = QualityScores::from_raw(original.clone());
        let ascii = q.to_ascii(PhredEncoding::Phred33);
        let q2 = QualityScores::from_ascii(&ascii, PhredEncoding::Phred33).unwrap();
        assert_eq!(q2.as_slice(), &original);
    }

    #[test]
    fn stats() {
        let q = QualityScores::from_raw(vec![10, 20, 30, 40]);
        assert!((q.mean() - 25.0).abs() < 1e-10);
        assert_eq!(q.min(), Some(10));
        assert_eq!(q.max(), Some(40));
        assert!((q.fraction_above(20) - 0.75).abs() < 1e-10);
        assert!((q.fraction_above(30) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn error_probability() {
        let p = QualityScores::error_probability(10);
        assert!((p - 0.1).abs() < 1e-10);
        let p = QualityScores::error_probability(20);
        assert!((p - 0.01).abs() < 1e-10);
    }

    #[test]
    fn scored_trait() {
        let q = QualityScores::from_raw(vec![20, 30]);
        assert!((q.score() - 25.0).abs() < 1e-10);
    }
}
