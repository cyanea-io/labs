//! False discovery rate (FDR) estimation via target-decoy approach.
//!
//! Implements the target-decoy competition strategy and q-value calculation
//! for controlling FDR in peptide/protein identifications.

use crate::error::{ProteomicsError, Result};
use crate::search::Psm;

/// FDR estimation result for a single PSM.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FdrResult {
    /// Spectrum identifier.
    pub spectrum_id: String,
    /// Peptide sequence.
    pub peptide_sequence: String,
    /// Original score.
    pub score: f64,
    /// Estimated q-value (minimum FDR at which this PSM would be accepted).
    pub q_value: f64,
    /// Whether this is a decoy hit.
    pub is_decoy: bool,
    /// Passes the specified FDR threshold.
    pub passes: bool,
}

/// Configuration for FDR estimation.
#[derive(Debug, Clone)]
pub struct FdrConfig {
    /// FDR threshold (e.g., 0.01 for 1% FDR).
    pub threshold: f64,
    /// Score field to use for ranking.
    pub score_type: ScoreType,
}

/// Which score to use for FDR ranking.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoreType {
    Hyperscore,
    Xcorr,
}

impl Default for FdrConfig {
    fn default() -> Self {
        Self {
            threshold: 0.01,
            score_type: ScoreType::Hyperscore,
        }
    }
}

/// Estimate FDR using target-decoy competition.
///
/// # Algorithm
///
/// 1. Combine target and decoy PSMs
/// 2. Sort by score (descending)
/// 3. At each score threshold, FDR = (decoy count) / (target count)
/// 4. Compute q-values as the minimum FDR at or below each PSM
/// 5. Mark PSMs passing the FDR threshold
pub fn estimate_fdr(
    target_psms: &[Psm],
    decoy_psms: &[Psm],
    config: &FdrConfig,
) -> Result<Vec<FdrResult>> {
    if target_psms.is_empty() {
        return Err(ProteomicsError::Fdr("no target PSMs".into()));
    }

    // Combine and sort by score descending
    let mut all: Vec<(&Psm, bool)> = target_psms.iter().map(|p| (p, false))
        .chain(decoy_psms.iter().map(|p| (p, true)))
        .collect();

    all.sort_by(|a, b| {
        let score_a = match config.score_type {
            ScoreType::Hyperscore => a.0.hyperscore,
            ScoreType::Xcorr => a.0.xcorr,
        };
        let score_b = match config.score_type {
            ScoreType::Hyperscore => b.0.hyperscore,
            ScoreType::Xcorr => b.0.xcorr,
        };
        score_b.partial_cmp(&score_a).unwrap_or(std::cmp::Ordering::Equal)
    });

    // Compute FDR at each rank
    let mut results = Vec::with_capacity(all.len());
    let mut target_count = 0usize;
    let mut decoy_count = 0usize;

    for &(psm, is_decoy) in &all {
        if is_decoy {
            decoy_count += 1;
        } else {
            target_count += 1;
        }

        let fdr = if target_count > 0 {
            decoy_count as f64 / target_count as f64
        } else {
            1.0
        };

        let score = match config.score_type {
            ScoreType::Hyperscore => psm.hyperscore,
            ScoreType::Xcorr => psm.xcorr,
        };

        results.push(FdrResult {
            spectrum_id: psm.spectrum_id.clone(),
            peptide_sequence: psm.peptide_sequence.clone(),
            score,
            q_value: fdr, // Will be corrected to q-value below
            is_decoy,
            passes: false, // Will be set below
        });
    }

    // Convert FDR to q-values (monotonically increasing from bottom)
    let n = results.len();
    if n > 1 {
        let mut min_fdr = results[n - 1].q_value;
        for i in (0..n).rev() {
            if results[i].q_value < min_fdr {
                min_fdr = results[i].q_value;
            }
            results[i].q_value = min_fdr;
        }
    }

    // Mark passing PSMs
    for r in &mut results {
        r.passes = !r.is_decoy && r.q_value <= config.threshold;
    }

    Ok(results)
}

/// Filter PSMs to a specified FDR level.
///
/// Returns only target PSMs passing the FDR threshold.
pub fn filter_fdr(
    target_psms: &[Psm],
    decoy_psms: &[Psm],
    fdr_threshold: f64,
) -> Result<Vec<Psm>> {
    let config = FdrConfig {
        threshold: fdr_threshold,
        ..Default::default()
    };

    let results = estimate_fdr(target_psms, decoy_psms, &config)?;

    let passing: Vec<Psm> = results.iter()
        .filter(|r| r.passes)
        .filter_map(|r| {
            target_psms.iter()
                .find(|p| p.spectrum_id == r.spectrum_id && p.peptide_sequence == r.peptide_sequence)
                .cloned()
        })
        .collect();

    Ok(passing)
}

/// Summary statistics for FDR analysis.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FdrSummary {
    /// Total target PSMs.
    pub total_targets: usize,
    /// Total decoy PSMs.
    pub total_decoys: usize,
    /// PSMs passing at 1% FDR.
    pub passing_1pct: usize,
    /// PSMs passing at 5% FDR.
    pub passing_5pct: usize,
    /// Score at 1% FDR cutoff.
    pub score_at_1pct: f64,
    /// Score at 5% FDR cutoff.
    pub score_at_5pct: f64,
}

/// Compute FDR summary statistics.
pub fn fdr_summary(
    target_psms: &[Psm],
    decoy_psms: &[Psm],
) -> Result<FdrSummary> {
    let config_1 = FdrConfig { threshold: 0.01, ..Default::default() };
    let results_1 = estimate_fdr(target_psms, decoy_psms, &config_1)?;

    let config_5 = FdrConfig { threshold: 0.05, ..Default::default() };
    let results_5 = estimate_fdr(target_psms, decoy_psms, &config_5)?;

    let passing_1: Vec<&FdrResult> = results_1.iter().filter(|r| r.passes).collect();
    let passing_5: Vec<&FdrResult> = results_5.iter().filter(|r| r.passes).collect();

    Ok(FdrSummary {
        total_targets: target_psms.len(),
        total_decoys: decoy_psms.len(),
        passing_1pct: passing_1.len(),
        passing_5pct: passing_5.len(),
        score_at_1pct: passing_1.last().map(|r| r.score).unwrap_or(f64::MAX),
        score_at_5pct: passing_5.last().map(|r| r.score).unwrap_or(f64::MAX),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_psm(id: &str, seq: &str, score: f64) -> Psm {
        Psm {
            spectrum_id: id.to_string(),
            peptide_sequence: seq.to_string(),
            xcorr: score * 0.1,
            hyperscore: score,
            matched_b: 5,
            matched_y: 5,
            total_b: 6,
            total_y: 6,
            delta_mass: 0.0,
            charge: 2,
            is_decoy: false,
        }
    }

    #[test]
    fn test_basic_fdr() {
        let targets = vec![
            make_psm("s1", "AAAK", 50.0),
            make_psm("s2", "BBBK", 40.0),
            make_psm("s3", "CCCK", 30.0),
            make_psm("s4", "DDDK", 20.0),
            make_psm("s5", "EEEK", 10.0),
        ];
        let decoys = vec![
            make_psm("d1", "KAAA", 25.0),
        ];

        let config = FdrConfig { threshold: 0.05, ..Default::default() };
        let results = estimate_fdr(&targets, &decoys, &config).unwrap();

        // At score 50, 40, 30: 0 decoys / 3 targets = 0% FDR
        // At score 25 (decoy): 1 decoy / 3 targets = 33% FDR
        // At score 20: 1 decoy / 4 targets = 25% FDR
        // First 3 targets should pass at 5% FDR
        let passing: Vec<&FdrResult> = results.iter().filter(|r| r.passes).collect();
        assert_eq!(passing.len(), 3);
    }

    #[test]
    fn test_q_values_monotonic() {
        let targets = vec![
            make_psm("s1", "AAAK", 50.0),
            make_psm("s2", "BBBK", 40.0),
            make_psm("s3", "CCCK", 30.0),
        ];
        let decoys = vec![
            make_psm("d1", "KAAA", 35.0),
        ];

        let config = FdrConfig::default();
        let results = estimate_fdr(&targets, &decoys, &config).unwrap();

        // q-values should be monotonically non-decreasing (from top to bottom by score)
        // Actually our q-values go from bottom up as minimum, so they should be non-decreasing
        for i in 1..results.len() {
            assert!(results[i].q_value >= results[i - 1].q_value - 1e-10,
                "q-values should be non-decreasing: {} < {}", results[i].q_value, results[i-1].q_value);
        }
    }

    #[test]
    fn test_filter_fdr() {
        let targets = vec![
            make_psm("s1", "AAAK", 50.0),
            make_psm("s2", "BBBK", 40.0),
            make_psm("s3", "CCCK", 5.0),
        ];
        let decoys = vec![
            make_psm("d1", "KAAA", 8.0),
        ];

        let passing = filter_fdr(&targets, &decoys, 0.01).unwrap();
        // Only top 2 should pass (decoy at 8.0 means FDR jumps for score<=8)
        assert_eq!(passing.len(), 2);
    }

    #[test]
    fn test_fdr_empty_targets() {
        let result = estimate_fdr(&[], &[], &FdrConfig::default());
        assert!(result.is_err());
    }

    #[test]
    fn test_fdr_no_decoys() {
        let targets = vec![
            make_psm("s1", "AAAK", 50.0),
            make_psm("s2", "BBBK", 40.0),
        ];

        let config = FdrConfig { threshold: 0.01, ..Default::default() };
        let results = estimate_fdr(&targets, &[], &config).unwrap();
        // No decoys means FDR=0 everywhere, all should pass
        let passing: Vec<&FdrResult> = results.iter().filter(|r| r.passes).collect();
        assert_eq!(passing.len(), 2);
    }

    #[test]
    fn test_fdr_summary() {
        let targets = vec![
            make_psm("s1", "AAAK", 50.0),
            make_psm("s2", "BBBK", 40.0),
            make_psm("s3", "CCCK", 30.0),
            make_psm("s4", "DDDK", 20.0),
        ];
        let decoys = vec![
            make_psm("d1", "KAAA", 25.0),
        ];

        let summary = fdr_summary(&targets, &decoys).unwrap();
        assert_eq!(summary.total_targets, 4);
        assert_eq!(summary.total_decoys, 1);
        assert!(summary.passing_5pct <= summary.total_targets);
    }

    #[test]
    fn test_xcorr_scoring() {
        let targets = vec![
            make_psm("s1", "AAAK", 50.0),
            make_psm("s2", "BBBK", 40.0),
        ];
        let decoys = vec![];

        let config = FdrConfig {
            threshold: 0.05,
            score_type: ScoreType::Xcorr,
        };
        let results = estimate_fdr(&targets, &decoys, &config).unwrap();
        // Should use xcorr scores (which are 0.1 * hyperscore in our test)
        assert!((results[0].score - 5.0).abs() < 1e-10);
    }
}
