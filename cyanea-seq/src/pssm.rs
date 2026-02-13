//! Position-Specific Scoring Matrices (PSSMs) for motif representation and scanning.
//!
//! Build a PSSM from a count matrix, then score or scan sequences against it.
//! Supports any fixed-alphabet size via const generics, with [`PssmDna`] and
//! [`PssmProtein`] type aliases for common use cases.

use cyanea_core::{CyaneaError, Result};

/// A Position-Specific Scoring Matrix (PSSM) with a fixed alphabet size.
///
/// Scores are stored as log-odds ratios relative to background frequencies.
#[derive(Debug, Clone)]
pub struct Pssm<const N: usize> {
    /// Log-odds scores: `scores[pos][symbol]` for each position.
    scores: Vec<[f64; N]>,
    /// Background frequencies for each symbol.
    background: [f64; N],
}

/// PSSM for DNA sequences (A=0, C=1, G=2, T=3).
pub type PssmDna = Pssm<4>;

/// PSSM for protein sequences (20 standard amino acids).
pub type PssmProtein = Pssm<20>;

impl<const N: usize> Pssm<N> {
    /// Build a PSSM from a count matrix.
    ///
    /// Each row of `counts` represents one position; each column is a symbol.
    /// A `pseudocount` is added to every cell before converting to frequencies.
    /// Scores are computed as `ln(freq / background)`.
    ///
    /// # Errors
    ///
    /// Returns an error if `counts` is empty or any `background` entry is zero.
    pub fn from_counts(
        counts: &[[f64; N]],
        pseudocount: f64,
        background: [f64; N],
    ) -> Result<Self> {
        if counts.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "count matrix must have at least one position".into(),
            ));
        }
        for (i, &bg) in background.iter().enumerate() {
            if bg <= 0.0 {
                return Err(CyaneaError::InvalidInput(format!(
                    "background frequency at index {} must be positive, got {}",
                    i, bg
                )));
            }
        }

        let mut scores = Vec::with_capacity(counts.len());
        for row in counts {
            let total: f64 = row.iter().sum::<f64>() + pseudocount * N as f64;
            let mut log_odds = [0.0f64; N];
            for j in 0..N {
                let freq = (row[j] + pseudocount) / total;
                log_odds[j] = (freq / background[j]).ln();
            }
            scores.push(log_odds);
        }

        Ok(Self { scores, background })
    }

    /// Number of positions in the motif.
    pub fn len(&self) -> usize {
        self.scores.len()
    }

    /// Returns `true` if the PSSM has zero positions.
    pub fn is_empty(&self) -> bool {
        self.scores.is_empty()
    }

    /// Score a sequence window against this PSSM.
    ///
    /// `seq` must be exactly [`self.len()`] bytes. The `mapping` function
    /// converts each byte to an alphabet index in `0..N`.
    ///
    /// # Errors
    ///
    /// Returns an error if `seq.len() != self.len()` or if `mapping` returns
    /// `None` for any byte.
    pub fn score(&self, seq: &[u8], mapping: &dyn Fn(u8) -> Option<usize>) -> Result<f64> {
        if seq.len() != self.len() {
            return Err(CyaneaError::InvalidInput(format!(
                "sequence length {} does not match PSSM length {}",
                seq.len(),
                self.len()
            )));
        }
        let mut total = 0.0;
        for (i, &base) in seq.iter().enumerate() {
            let idx = mapping(base).ok_or_else(|| {
                CyaneaError::InvalidInput(format!(
                    "unmapped character '{}' at position {}",
                    base as char, i
                ))
            })?;
            total += self.scores[i][idx];
        }
        Ok(total)
    }

    /// Slide the PSSM across `seq` and return all hits at or above `threshold`.
    ///
    /// Returns `(position, score)` pairs. Positions with unmapped characters
    /// are silently skipped.
    pub fn scan(
        &self,
        seq: &[u8],
        threshold: f64,
        mapping: &dyn Fn(u8) -> Option<usize>,
    ) -> Vec<(usize, f64)> {
        let motif_len = self.len();
        if seq.len() < motif_len {
            return Vec::new();
        }
        let mut hits = Vec::new();
        for start in 0..=seq.len() - motif_len {
            if let Ok(s) = self.score(&seq[start..start + motif_len], mapping) {
                if s >= threshold {
                    hits.push((start, s));
                }
            }
        }
        hits
    }

    /// Information content (bits) at each position.
    ///
    /// `IC_j = sum_c freq_c * log2(freq_c / bg_c)` where frequencies are
    /// recovered from the stored log-odds scores.
    pub fn information_content(&self) -> Vec<f64> {
        self.scores
            .iter()
            .map(|row| {
                let mut ic = 0.0;
                for j in 0..N {
                    // score = ln(freq / bg), so freq = bg * exp(score)
                    let freq = self.background[j] * row[j].exp();
                    if freq > 0.0 {
                        ic += freq * (freq / self.background[j]).log2();
                    }
                }
                ic
            })
            .collect()
    }

    /// Maximum possible score (sum of the best symbol at each position).
    pub fn max_score(&self) -> f64 {
        self.scores
            .iter()
            .map(|row| row.iter().cloned().fold(f64::NEG_INFINITY, f64::max))
            .sum()
    }

    /// Minimum possible score (sum of the worst symbol at each position).
    pub fn min_score(&self) -> f64 {
        self.scores
            .iter()
            .map(|row| row.iter().cloned().fold(f64::INFINITY, f64::min))
            .sum()
    }
}

/// Map a DNA base to index (A=0, C=1, G=2, T=3).
pub fn dna_mapping(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Map an amino acid to index (0--19, alphabetical order:
/// A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y).
pub fn protein_mapping(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'D' | b'd' => Some(2),
        b'E' | b'e' => Some(3),
        b'F' | b'f' => Some(4),
        b'G' | b'g' => Some(5),
        b'H' | b'h' => Some(6),
        b'I' | b'i' => Some(7),
        b'K' | b'k' => Some(8),
        b'L' | b'l' => Some(9),
        b'M' | b'm' => Some(10),
        b'N' | b'n' => Some(11),
        b'P' | b'p' => Some(12),
        b'Q' | b'q' => Some(13),
        b'R' | b'r' => Some(14),
        b'S' | b's' => Some(15),
        b'T' | b't' => Some(16),
        b'V' | b'v' => Some(17),
        b'W' | b'w' => Some(18),
        b'Y' | b'y' => Some(19),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn uniform_bg() -> [f64; 4] {
        [0.25; 4]
    }

    #[test]
    fn uniform_counts_scores_near_zero() {
        let counts = vec![[10.0, 10.0, 10.0, 10.0]; 3];
        let pssm = PssmDna::from_counts(&counts, 0.0, uniform_bg()).unwrap();
        assert_eq!(pssm.len(), 3);
        let s = pssm.score(b"ACG", &dna_mapping).unwrap();
        assert!(s.abs() < 1e-10, "expected ~0, got {}", s);
    }

    #[test]
    fn biased_counts_high_score_for_consensus() {
        // Position 0 strongly favors A, position 1 favors C, position 2 favors G
        let counts = vec![
            [100.0, 1.0, 1.0, 1.0],
            [1.0, 100.0, 1.0, 1.0],
            [1.0, 1.0, 100.0, 1.0],
        ];
        let pssm = PssmDna::from_counts(&counts, 0.0, uniform_bg()).unwrap();
        let consensus = pssm.score(b"ACG", &dna_mapping).unwrap();
        let mismatch = pssm.score(b"TTA", &dna_mapping).unwrap();
        assert!(consensus > mismatch, "consensus {} should beat mismatch {}", consensus, mismatch);
        assert!(consensus > 0.0);
    }

    #[test]
    fn score_known_motif() {
        let counts = vec![
            [50.0, 0.0, 0.0, 0.0],
            [0.0, 50.0, 0.0, 0.0],
        ];
        let pssm = PssmDna::from_counts(&counts, 1.0, uniform_bg()).unwrap();
        let s = pssm.score(b"AC", &dna_mapping).unwrap();
        // freq = 51/54 ≈ 0.944, ln(0.944/0.25) ≈ 1.329 per position
        assert!(s > 2.0, "expected score > 2.0, got {}", s);
    }

    #[test]
    fn scan_finds_positions() {
        let counts = vec![
            [100.0, 0.0, 0.0, 0.0],
            [0.0, 100.0, 0.0, 0.0],
            [0.0, 0.0, 100.0, 0.0],
        ];
        let pssm = PssmDna::from_counts(&counts, 1.0, uniform_bg()).unwrap();
        let seq = b"TTACGTTACGTT";
        let hits = pssm.scan(seq, 0.0, &dna_mapping);
        // "ACG" appears at positions 2 and 7
        let positions: Vec<usize> = hits.iter().map(|&(p, _)| p).collect();
        assert!(positions.contains(&2), "expected hit at 2, got {:?}", positions);
        assert!(positions.contains(&7), "expected hit at 7, got {:?}", positions);
    }

    #[test]
    fn information_content_uniform_is_zero() {
        let counts = vec![[25.0, 25.0, 25.0, 25.0]];
        let pssm = PssmDna::from_counts(&counts, 0.0, uniform_bg()).unwrap();
        let ic = pssm.information_content();
        assert!(ic[0].abs() < 1e-10, "uniform IC should be ~0, got {}", ic[0]);
    }

    #[test]
    fn information_content_conserved_is_two_bits() {
        // Perfectly conserved A (with tiny pseudocount to avoid log(0))
        let counts = vec![[1000.0, 0.0, 0.0, 0.0]];
        let pssm = PssmDna::from_counts(&counts, 0.01, uniform_bg()).unwrap();
        let ic = pssm.information_content();
        assert!((ic[0] - 2.0).abs() < 0.05, "conserved IC should be ~2 bits, got {}", ic[0]);
    }

    #[test]
    fn error_empty_counts() {
        let counts: Vec<[f64; 4]> = vec![];
        let result = PssmDna::from_counts(&counts, 1.0, uniform_bg());
        assert!(result.is_err());
    }

    #[test]
    fn error_zero_background() {
        let counts = vec![[10.0; 4]];
        let result = PssmDna::from_counts(&counts, 1.0, [0.25, 0.0, 0.25, 0.25]);
        assert!(result.is_err());
    }

    #[test]
    fn error_wrong_seq_length() {
        let counts = vec![[10.0; 4]; 3];
        let pssm = PssmDna::from_counts(&counts, 1.0, uniform_bg()).unwrap();
        let result = pssm.score(b"AC", &dna_mapping);
        assert!(result.is_err());
    }

    #[test]
    fn error_unmapped_character() {
        let counts = vec![[10.0; 4]];
        let pssm = PssmDna::from_counts(&counts, 1.0, uniform_bg()).unwrap();
        let result = pssm.score(b"X", &dna_mapping);
        assert!(result.is_err());
    }

    #[test]
    fn min_max_score_bounds() {
        let counts = vec![
            [100.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 100.0],
        ];
        let pssm = PssmDna::from_counts(&counts, 0.0, uniform_bg()).unwrap();
        let best = pssm.score(b"AT", &dna_mapping).unwrap();
        let worst = pssm.score(b"TA", &dna_mapping).unwrap();
        assert!((best - pssm.max_score()).abs() < 1e-10);
        assert!((worst - pssm.min_score()).abs() < 1e-10);
        assert!(pssm.max_score() > pssm.min_score());
    }

    #[test]
    fn dna_mapping_cases() {
        assert_eq!(dna_mapping(b'A'), Some(0));
        assert_eq!(dna_mapping(b'a'), Some(0));
        assert_eq!(dna_mapping(b'C'), Some(1));
        assert_eq!(dna_mapping(b'G'), Some(2));
        assert_eq!(dna_mapping(b'T'), Some(3));
        assert_eq!(dna_mapping(b't'), Some(3));
        assert_eq!(dna_mapping(b'N'), None);
        assert_eq!(dna_mapping(b'X'), None);
    }

    #[test]
    fn protein_mapping_cases() {
        assert_eq!(protein_mapping(b'A'), Some(0));
        assert_eq!(protein_mapping(b'Y'), Some(19));
        assert_eq!(protein_mapping(b'w'), Some(18));
        assert_eq!(protein_mapping(b'K'), Some(8));
        assert_eq!(protein_mapping(b'X'), None);
        assert_eq!(protein_mapping(b'B'), None);
    }

    #[test]
    fn scan_short_seq_returns_empty() {
        let counts = vec![[10.0; 4]; 5];
        let pssm = PssmDna::from_counts(&counts, 1.0, uniform_bg()).unwrap();
        let hits = pssm.scan(b"ACG", 0.0, &dna_mapping);
        assert!(hits.is_empty());
    }

    #[test]
    fn case_insensitive_scoring() {
        let counts = vec![[100.0, 0.0, 0.0, 0.0]];
        let pssm = PssmDna::from_counts(&counts, 1.0, uniform_bg()).unwrap();
        let upper = pssm.score(b"A", &dna_mapping).unwrap();
        let lower = pssm.score(b"a", &dna_mapping).unwrap();
        assert!((upper - lower).abs() < 1e-10);
    }
}
