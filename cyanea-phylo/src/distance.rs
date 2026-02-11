//! Evolutionary distance models for sequence data.
//!
//! Provides p-distance, Jukes-Cantor, and Kimura 2-parameter models for
//! computing evolutionary distances between aligned sequences.

use cyanea_core::{CyaneaError, Result};

/// Proportion of differing sites between two aligned sequences.
///
/// Both sequences must be the same length and non-empty.
pub fn p_distance(a: &[u8], b: &[u8]) -> Result<f64> {
    validate_sequences(a, b)?;
    let diffs = a.iter().zip(b).filter(|(x, y)| x != y).count();
    Ok(diffs as f64 / a.len() as f64)
}

/// Jukes-Cantor distance from an observed proportion of differences `p`.
///
/// Formula: d = -3/4 * ln(1 - 4p/3)
///
/// Errors if `p < 0` or `p >= 0.75` (saturation).
pub fn jukes_cantor(p: f64) -> Result<f64> {
    if p < 0.0 {
        return Err(CyaneaError::InvalidInput(format!(
            "p-distance must be non-negative, got {}",
            p
        )));
    }
    if p >= 0.75 {
        return Err(CyaneaError::InvalidInput(format!(
            "p-distance {} >= 0.75 (saturation limit for Jukes-Cantor)",
            p
        )));
    }
    if p == 0.0 {
        return Ok(0.0);
    }
    Ok(-0.75 * (1.0 - 4.0 * p / 3.0).ln())
}

/// Kimura 2-parameter distance from observed transition and transversion proportions.
///
/// Formula: d = -1/2 * ln((1 - 2S - V) * sqrt(1 - 2V))
///
/// where `S` is the proportion of transitions and `V` the proportion of transversions.
pub fn kimura_2p(transitions: f64, transversions: f64) -> Result<f64> {
    if transitions < 0.0 || transversions < 0.0 {
        return Err(CyaneaError::InvalidInput(
            "transition/transversion proportions must be non-negative".into(),
        ));
    }
    let s = transitions;
    let v = transversions;
    let a = 1.0 - 2.0 * s - v;
    let b = 1.0 - 2.0 * v;
    if a <= 0.0 || b <= 0.0 {
        return Err(CyaneaError::InvalidInput(format!(
            "substitution rates too high for Kimura 2-parameter model (S={}, V={})",
            s, v
        )));
    }
    if s == 0.0 && v == 0.0 {
        return Ok(0.0);
    }
    Ok(-0.5 * (a * b.sqrt()).ln())
}

/// Classify a substitution between two nucleotide bases.
///
/// Returns `(is_transition, is_transversion)`. Both are false if bases are equal
/// or not standard DNA bases (A, C, G, T).
#[cfg_attr(not(feature = "ml"), allow(dead_code))]
fn classify_substitution(a: u8, b: u8) -> (bool, bool) {
    let a = a.to_ascii_uppercase();
    let b = b.to_ascii_uppercase();
    if a == b {
        return (false, false);
    }
    match (a, b) {
        // Transitions: purine<->purine or pyrimidine<->pyrimidine
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C') => (true, false),
        // Transversions: purine<->pyrimidine
        (b'A', b'C')
        | (b'C', b'A')
        | (b'A', b'T')
        | (b'T', b'A')
        | (b'G', b'C')
        | (b'C', b'G')
        | (b'G', b'T')
        | (b'T', b'G') => (false, true),
        _ => (false, false),
    }
}

/// Evolutionary distance model for sequence comparison.
#[cfg(feature = "ml")]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistanceModel {
    /// Raw proportion of differences.
    P,
    /// Jukes-Cantor correction.
    JukesCantor,
    /// Kimura 2-parameter correction.
    Kimura2P,
}

/// Build a pairwise distance matrix from aligned sequences using the given model.
///
/// All sequences must be the same length. Returns a `cyanea_ml::DistanceMatrix`.
#[cfg(feature = "ml")]
pub fn sequence_distance_matrix(
    seqs: &[&[u8]],
    model: DistanceModel,
) -> Result<cyanea_ml::DistanceMatrix> {
    let n = seqs.len();
    if n < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 sequences".into(),
        ));
    }
    let len = seqs[0].len();
    if len == 0 {
        return Err(CyaneaError::InvalidInput("empty sequences".into()));
    }
    for (i, seq) in seqs.iter().enumerate() {
        if seq.len() != len {
            return Err(CyaneaError::InvalidInput(format!(
                "sequence {} has length {}, expected {}",
                i,
                seq.len(),
                len
            )));
        }
    }

    #[cfg(feature = "parallel")]
    let condensed = {
        use rayon::prelude::*;
        (0..n)
            .into_par_iter()
            .map(|i| {
                ((i + 1)..n)
                    .map(|j| match model {
                        DistanceModel::P => p_distance(seqs[i], seqs[j]),
                        DistanceModel::JukesCantor => {
                            let p = p_distance(seqs[i], seqs[j])?;
                            jukes_cantor(p)
                        }
                        DistanceModel::Kimura2P => {
                            let (ts, tv) = count_substitutions(seqs[i], seqs[j]);
                            let total = seqs[i].len() as f64;
                            kimura_2p(ts as f64 / total, tv as f64 / total)
                        }
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect::<Vec<f64>>()
    };

    #[cfg(not(feature = "parallel"))]
    let condensed = {
        let size = n * (n - 1) / 2;
        let mut condensed = Vec::with_capacity(size);
        for i in 0..n {
            for j in (i + 1)..n {
                let d = match model {
                    DistanceModel::P => p_distance(seqs[i], seqs[j])?,
                    DistanceModel::JukesCantor => {
                        let p = p_distance(seqs[i], seqs[j])?;
                        jukes_cantor(p)?
                    }
                    DistanceModel::Kimura2P => {
                        let (ts, tv) = count_substitutions(seqs[i], seqs[j]);
                        let total = seqs[i].len() as f64;
                        kimura_2p(ts as f64 / total, tv as f64 / total)?
                    }
                };
                condensed.push(d);
            }
        }
        condensed
    };

    cyanea_ml::DistanceMatrix::from_condensed(condensed, n)
}

/// Count transitions and transversions between two equal-length sequences.
#[cfg(feature = "ml")]
fn count_substitutions(a: &[u8], b: &[u8]) -> (usize, usize) {
    let mut transitions = 0;
    let mut transversions = 0;
    for (&x, &y) in a.iter().zip(b) {
        let (ts, tv) = classify_substitution(x, y);
        if ts {
            transitions += 1;
        }
        if tv {
            transversions += 1;
        }
    }
    (transitions, transversions)
}

fn validate_sequences(a: &[u8], b: &[u8]) -> Result<()> {
    if a.is_empty() {
        return Err(CyaneaError::InvalidInput("empty sequences".into()));
    }
    if a.len() != b.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "sequence length mismatch: {} vs {}",
            a.len(),
            b.len()
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn p_distance_identical() {
        assert_eq!(p_distance(b"ACGT", b"ACGT").unwrap(), 0.0);
    }

    #[test]
    fn p_distance_all_different() {
        assert_eq!(p_distance(b"AAAA", b"TTTT").unwrap(), 1.0);
    }

    #[test]
    fn p_distance_half() {
        let p = p_distance(b"AABB", b"AACC").unwrap();
        assert!((p - 0.5).abs() < 1e-12);
    }

    #[test]
    fn p_distance_empty_error() {
        assert!(p_distance(b"", b"").is_err());
    }

    #[test]
    fn p_distance_length_mismatch() {
        assert!(p_distance(b"AC", b"ACG").is_err());
    }

    #[test]
    fn jc_zero() {
        assert_eq!(jukes_cantor(0.0).unwrap(), 0.0);
    }

    #[test]
    fn jc_known_value() {
        // p = 0.1 -> d = -0.75 * ln(1 - 4*0.1/3) = -0.75 * ln(0.8667)
        let d = jukes_cantor(0.1).unwrap();
        let expected = -0.75 * (1.0 - 4.0 * 0.1 / 3.0_f64).ln();
        assert!((d - expected).abs() < 1e-12);
    }

    #[test]
    fn jc_saturation_error() {
        assert!(jukes_cantor(0.75).is_err());
        assert!(jukes_cantor(0.9).is_err());
    }

    #[test]
    fn jc_negative_error() {
        assert!(jukes_cantor(-0.1).is_err());
    }

    #[test]
    fn k2p_zero() {
        assert_eq!(kimura_2p(0.0, 0.0).unwrap(), 0.0);
    }

    #[test]
    fn k2p_known_value() {
        let d = kimura_2p(0.05, 0.03).unwrap();
        let a: f64 = 1.0 - 2.0 * 0.05 - 0.03;
        let b: f64 = 1.0 - 2.0 * 0.03;
        let expected = -0.5 * (a * b.sqrt()).ln();
        assert!((d - expected).abs() < 1e-12);
    }

    #[test]
    fn k2p_invalid_rates() {
        assert!(kimura_2p(0.4, 0.2).is_err()); // 2*0.4 + 0.2 = 1.0, a = 0
    }

    #[test]
    fn k2p_negative_error() {
        assert!(kimura_2p(-0.1, 0.0).is_err());
        assert!(kimura_2p(0.0, -0.1).is_err());
    }

    #[test]
    fn classify_transition() {
        assert_eq!(classify_substitution(b'A', b'G'), (true, false));
        assert_eq!(classify_substitution(b'C', b'T'), (true, false));
    }

    #[test]
    fn classify_transversion() {
        assert_eq!(classify_substitution(b'A', b'C'), (false, true));
        assert_eq!(classify_substitution(b'G', b'T'), (false, true));
    }

    #[test]
    fn classify_same_base() {
        assert_eq!(classify_substitution(b'A', b'A'), (false, false));
    }

    #[test]
    fn classify_lowercase() {
        assert_eq!(classify_substitution(b'a', b'g'), (true, false));
    }
}
