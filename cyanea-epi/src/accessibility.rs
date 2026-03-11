//! ATAC-seq specific analysis and quality control metrics.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::pileup::TagPileup;

/// Insert size metrics for ATAC-seq QC.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct InsertSizeMetrics {
    /// Nucleosome-free region ratio (< 150 bp / 150-300 bp).
    pub nfr_ratio: f64,
    /// Nucleosomal periodicity score (0-1).
    pub periodicity: f64,
    /// Mean insert size.
    pub mean_insert_size: f64,
    /// Median insert size.
    pub median_insert_size: f64,
}

/// Comprehensive ATAC-seq QC results.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AtacQcResult {
    /// TSS enrichment score (ENCODE metric).
    pub tss_enrichment: f64,
    /// Fraction of reads in peaks.
    pub frip: f64,
    /// Nucleosome-free region ratio.
    pub nfr_ratio: f64,
    /// Nucleosomal periodicity.
    pub periodicity: f64,
}

/// Compute TSS enrichment score (ENCODE QC metric).
///
/// Aggregate signal in windows around transcription start sites
/// and compare to flanking regions.
pub fn tss_enrichment(
    pileup: &TagPileup,
    tss_positions: &[(String, u64)],
    up_window: u64,
    down_window: u64,
) -> f64 {
    if tss_positions.is_empty() {
        return 0.0;
    }

    let mut core_signal = 0.0;
    let mut flank_signal = 0.0;
    let mut count = 0;

    for (chrom, tss) in tss_positions {
        if let Some(cov) = pileup.coverage.get(chrom) {
            let core_start = tss.saturating_sub(up_window) as usize;
            let core_end = (*tss + down_window) as usize;
            let flank_start = tss.saturating_sub(up_window * 2) as usize;
            let flank_end = (*tss + down_window * 2) as usize;

            // Sum core region
            for i in core_start..core_end.min(cov.len()) {
                core_signal += cov[i] as f64;
            }

            // Sum flanking region (excluding core)
            for i in flank_start..flank_end.min(cov.len()) {
                if i < core_start || i >= core_end {
                    flank_signal += cov[i] as f64;
                }
            }

            count += 1;
        }
    }

    if count == 0 || flank_signal == 0.0 {
        return 0.0;
    }

    let avg_core = core_signal / (count as f64);
    let avg_flank = flank_signal / (count as f64 * 4.0);

    avg_core / avg_flank.max(1.0)
}

/// Compute fragment size distribution histogram.
///
/// Takes fragment lengths (insert sizes) and returns counts per bin.
pub fn fragment_size_distribution(fragment_sizes: &[u64], bin_size: u64) -> Vec<(u64, u64)> {
    if fragment_sizes.is_empty() {
        return Vec::new();
    }

    let max_size = *fragment_sizes.iter().max().unwrap_or(&0);
    let n_bins = ((max_size + bin_size - 1) / bin_size) as usize;

    let mut bins = vec![0u64; n_bins];

    for &size in fragment_sizes {
        let bin = (size / bin_size) as usize;
        if bin < n_bins {
            bins[bin] += 1;
        }
    }

    bins.into_iter()
        .enumerate()
        .map(|(i, count)| ((i as u64) * bin_size, count))
        .filter(|(_, count)| *count > 0)
        .collect()
}

/// Compute insert size metrics from fragment sizes.
pub fn insert_size_metrics(fragment_sizes: &[u64]) -> InsertSizeMetrics {
    if fragment_sizes.is_empty() {
        return InsertSizeMetrics {
            nfr_ratio: 0.0,
            periodicity: 0.0,
            mean_insert_size: 0.0,
            median_insert_size: 0.0,
        };
    }

    // Basic statistics
    let mean = (fragment_sizes.iter().sum::<u64>() as f64) / (fragment_sizes.len() as f64);

    let mut sorted = fragment_sizes.to_vec();
    sorted.sort_unstable();
    let median = sorted[sorted.len() / 2] as f64;

    // NFR ratio: fragments < 150 bp vs 150-300 bp
    let nfr_count = fragment_sizes.iter().filter(|&&s| s < 150).count() as f64;
    let mono_count = fragment_sizes
        .iter()
        .filter(|&&s| s >= 150 && s <= 300)
        .count() as f64;

    let nfr_ratio = if mono_count > 0.0 {
        nfr_count / mono_count
    } else {
        0.0
    };

    // Periodicity: check for 147 bp repeats
    let mut periodicity = 0.0;
    if fragment_sizes.len() > 10 {
        let hist = fragment_size_distribution(fragment_sizes, 1);
        if hist.len() > 147 {
            let mut corr = 0.0;
            let mut count = 0;
            for i in 0..(hist.len().saturating_sub(147)) {
                if hist[i].1 > 0 {
                    corr += (hist[i].1 as f64) * (hist[i + 147].1 as f64);
                    count += 1;
                }
            }
            if count > 0 {
                periodicity = (corr / (count as f64)).min(1.0);
            }
        }
    }

    InsertSizeMetrics {
        nfr_ratio,
        periodicity,
        mean_insert_size: mean,
        median_insert_size: median,
    }
}

/// Compute fraction of reads in peaks (FRiP).
pub fn frip(reads_in_peaks: u64, total_reads: u64) -> f64 {
    if total_reads == 0 {
        0.0
    } else {
        (reads_in_peaks as f64) / (total_reads as f64)
    }
}

/// Compute comprehensive ATAC-seq QC metrics.
pub fn atacqc(
    pileup: &TagPileup,
    tss_positions: &[(String, u64)],
    fragment_sizes: &[u64],
    reads_in_peaks: u64,
    total_reads: u64,
) -> AtacQcResult {
    let tss_enrichment = tss_enrichment(pileup, tss_positions, 500, 500);
    let frip_val = frip(reads_in_peaks, total_reads);
    let metrics = insert_size_metrics(fragment_sizes);

    AtacQcResult {
        tss_enrichment,
        frip: frip_val,
        nfr_ratio: metrics.nfr_ratio,
        periodicity: metrics.periodicity,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    #[test]
    fn test_tss_enrichment() {
        let mut coverage = BTreeMap::new();
        let mut cov_vec = vec![1u32; 2000];

        // Create enrichment around position 1000
        for i in 900..1100 {
            cov_vec[i] = 10;
        }

        coverage.insert("chr1".to_string(), cov_vec);
        let pileup = TagPileup { coverage };

        let tss_positions = vec![("chr1".to_string(), 1000)];
        let enrichment = tss_enrichment(&pileup, &tss_positions, 500, 500);

        assert!(enrichment > 1.0);
    }

    #[test]
    fn test_fragment_size_distribution() {
        let fragments = vec![80, 90, 100, 150, 160, 170, 200, 300];
        let dist = fragment_size_distribution(&fragments, 50);

        assert!(!dist.is_empty());
    }

    #[test]
    fn test_insert_size_metrics() {
        let fragments = vec![80, 90, 100, 100, 100, 200, 250, 280];
        let metrics = insert_size_metrics(&fragments);

        assert!(metrics.mean_insert_size > 0.0);
        assert!(metrics.median_insert_size > 0.0);
        assert!(metrics.nfr_ratio > 0.0);
    }

    #[test]
    fn test_frip() {
        let frip_val = frip(500, 1000);
        assert!((frip_val - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_atacqc() {
        let mut coverage = BTreeMap::new();
        coverage.insert("chr1".to_string(), vec![1u32; 1000]);

        let pileup = TagPileup { coverage };
        let tss_positions = vec![("chr1".to_string(), 500)];
        let fragments = vec![80, 100, 150, 200];

        let qc = atacqc(&pileup, &tss_positions, &fragments, 100, 200);

        assert!(qc.tss_enrichment >= 0.0);
        assert!(qc.frip >= 0.0 && qc.frip <= 1.0);
        assert!(qc.nfr_ratio >= 0.0);
        assert!(qc.periodicity >= 0.0);
    }
}
