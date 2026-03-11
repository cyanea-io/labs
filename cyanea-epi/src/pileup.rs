//! Signal track operations and pileup analysis.
//!
//! This module provides tools for building coverage pileups from aligned reads,
//! normalizing signal, smoothing, and computing QC metrics like correlation and enrichment.

use std::collections::BTreeMap;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::error::{EpiError, Result};

/// Per-chromosome coverage vectors (0-based positions).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TagPileup {
    /// Coverage per chromosome.
    pub coverage: BTreeMap<String, Vec<u32>>,
}

impl TagPileup {
    /// Create a new empty pileup.
    pub fn new() -> Self {
        Self {
            coverage: BTreeMap::new(),
        }
    }

    /// Get coverage at a specific position.
    pub fn get(&self, chrom: &str, position: u64) -> u32 {
        self.coverage
            .get(chrom)
            .and_then(|cov| cov.get(position as usize).copied())
            .unwrap_or(0)
    }

    /// Get total bases with coverage > 0.
    pub fn coverage_bases(&self) -> u64 {
        self.coverage
            .values()
            .map(|c| c.iter().filter(|&&v| v > 0).count() as u64)
            .sum()
    }

    /// Get mean coverage.
    pub fn mean_coverage(&self) -> f64 {
        let total_bases = self.coverage.values().map(|c| c.len()).sum::<usize>();
        if total_bases == 0 {
            return 0.0;
        }

        let total_cov: u64 = self
            .coverage
            .values()
            .flat_map(|c| c.iter().map(|&v| v as u64))
            .sum();

        (total_cov as f64) / (total_bases as f64)
    }
}

impl Default for TagPileup {
    fn default() -> Self {
        Self::new()
    }
}

/// Build a tag pileup from aligned reads.
///
/// Each read is extended to `fragment_size` and coverage is incremented across
/// the extended interval.
///
/// # Arguments
///
/// * `reads` — (chrom, start, original_length) tuples
/// * `fragment_size` — length to extend reads to (e.g., 200 for ChIP-seq)
pub fn build_pileup(reads: &[(String, u64, u64)], fragment_size: u64) -> TagPileup {
    let mut pileup = BTreeMap::new();

    for (chrom, start, read_len) in reads {
        let extended_len = (*read_len).max(fragment_size);
        let end = start + extended_len;

        let cov = pileup.entry(chrom.clone()).or_insert_with(Vec::new);

        // Extend coverage array if needed
        if end as usize >= cov.len() {
            cov.resize((end as usize) + 1, 0);
        }

        for i in *start as usize..=((end - 1) as usize).min(cov.len() - 1) {
            cov[i] = (cov[i] as u32).saturating_add(1);
        }
    }

    TagPileup {
        coverage: pileup,
    }
}

/// Normalize pileup coverage.
///
/// Normalization methods:
/// - `"cpm"` — counts per million mapped reads
/// - `"rpkm"` — reads per kilobase per million mapped reads
#[allow(clippy::result_large_error_types)]
pub fn normalize_pileup(pileup: &TagPileup, method: &str) -> Result<TagPileup> {
    let total_reads: u64 = pileup
        .coverage
        .values()
        .flat_map(|c| c.iter().map(|&v| v as u64))
        .sum();

    if total_reads == 0 {
        return Err(EpiError::InsufficientData("No coverage in pileup".to_string()));
    }

    let mut normalized = BTreeMap::new();

    for (chrom, cov) in &pileup.coverage {
        let norm_cov = match method {
            "cpm" => {
                let scale = 1e6 / (total_reads as f64);
                cov.iter().map(|&c| (c as f64 * scale) as u32).collect()
            }
            "rpkm" => {
                let scale = 1e9 / ((total_reads as f64) * (cov.len() as f64).max(1.0));
                cov.iter().map(|&c| (c as f64 * scale) as u32).collect()
            }
            _ => {
                return Err(EpiError::InvalidInput(format!(
                    "Unknown normalization method: {}",
                    method
                )))
            }
        };

        normalized.insert(chrom.clone(), norm_cov);
    }

    Ok(TagPileup {
        coverage: normalized,
    })
}

/// Smooth pileup using Gaussian kernel.
///
/// Applies a Gaussian smoothing filter with the specified bandwidth (standard deviation).
pub fn smooth_pileup(pileup: &TagPileup, bandwidth: f64) -> TagPileup {
    let mut smoothed = BTreeMap::new();

    for (chrom, cov) in &pileup.coverage {
        let mut smooth_cov = vec![0.0; cov.len()];

        // Pre-compute Gaussian kernel
        let kernel_radius = (3.0 * bandwidth) as i32;
        let mut kernel = Vec::new();
        let mut kernel_sum = 0.0;

        for i in (-kernel_radius)..=(kernel_radius) {
            let v = (-((i as f64) * (i as f64)) / (2.0 * bandwidth * bandwidth)).exp();
            kernel.push(v);
            kernel_sum += v;
        }

        // Normalize kernel
        kernel.iter_mut().for_each(|k| *k /= kernel_sum);

        // Apply kernel to each position
        for pos in 0..cov.len() {
            let mut weighted_sum = 0.0;

            for (k_idx, &k_val) in kernel.iter().enumerate() {
                let offset = k_idx as i32 - kernel_radius;
                let neighbor = (pos as i32 + offset) as usize;

                if neighbor < cov.len() {
                    weighted_sum += (cov[neighbor] as f64) * k_val;
                }
            }

            smooth_cov[pos] = weighted_sum;
        }

        smoothed.insert(chrom.clone(), smooth_cov.iter().map(|&v| v as u32).collect());
    }

    TagPileup { coverage: smoothed }
}

/// Compute Pearson correlation between two pileups (for replicate QC).
///
/// Compares coverage at all positions in both pileups on the same chromosome.
/// Returns correlation coefficient or None if insufficient data.
pub fn pileup_correlation(pileup1: &TagPileup, pileup2: &TagPileup) -> Option<f64> {
    let mut x_vals = Vec::new();
    let mut y_vals = Vec::new();

    for chrom in pileup1.coverage.keys() {
        if let Some(cov1) = pileup1.coverage.get(chrom) {
            if let Some(cov2) = pileup2.coverage.get(chrom) {
                for i in 0..cov1.len().min(cov2.len()) {
                    x_vals.push(cov1[i] as f64);
                    y_vals.push(cov2[i] as f64);
                }
            }
        }
    }

    if x_vals.len() < 2 {
        return None;
    }

    let n = x_vals.len() as f64;
    let mean_x = x_vals.iter().sum::<f64>() / n;
    let mean_y = y_vals.iter().sum::<f64>() / n;

    let mut cov = 0.0;
    let mut var_x = 0.0;
    let mut var_y = 0.0;

    for (xi, yi) in x_vals.iter().zip(y_vals.iter()) {
        let dx = xi - mean_x;
        let dy = yi - mean_y;
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }

    let denom = (var_x * var_y).sqrt();
    if denom == 0.0 {
        return Some(0.0);
    }

    Some(cov / denom)
}

/// Compute cumulative enrichment plot (deepTools-style fingerprint).
///
/// Returns (cumulative_bases, cumulative_coverage) pairs sorted by cumulative rank.
/// Useful for assessing quality of ChIP-seq data.
pub fn fingerprint(pileup: &TagPileup) -> Vec<(f64, f64)> {
    let mut all_cov: Vec<u32> = pileup
        .coverage
        .values()
        .flat_map(|c| c.iter().copied())
        .collect();

    if all_cov.is_empty() {
        return Vec::new();
    }

    // Sort by coverage (ascending)
    all_cov.sort_unstable();

    let total_cov: u64 = all_cov.iter().map(|&c| c as u64).sum();
    let max_bases = all_cov.len() as f64;

    let mut result = Vec::new();
    let mut cumulative_cov = 0u64;

    for (i, &cov) in all_cov.iter().enumerate() {
        cumulative_cov += cov as u64;
        let cum_bases = ((i + 1) as f64 / max_bases) * 100.0;
        let cum_cov = (cumulative_cov as f64 / (total_cov as f64).max(1.0)) * 100.0;
        result.push((cum_bases, cum_cov));
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_pileup() {
        let reads = vec![
            ("chr1".to_string(), 100, 50),
            ("chr1".to_string(), 150, 50),
        ];

        let pileup = build_pileup(&reads, 200);

        assert!(pileup.coverage.contains_key("chr1"));
        let chr1_cov = &pileup.coverage["chr1"];
        assert!(chr1_cov[100] > 0);
        assert!(chr1_cov[150] > 0);
    }

    #[test]
    fn test_pileup_get() {
        let reads = vec![("chr1".to_string(), 100, 50)];
        let pileup = build_pileup(&reads, 200);

        assert!(pileup.get("chr1", 100) > 0);
        assert_eq!(pileup.get("chr2", 100), 0);
    }

    #[test]
    fn test_pileup_mean_coverage() {
        let reads = vec![
            ("chr1".to_string(), 0, 50),
            ("chr1".to_string(), 100, 50),
        ];

        let pileup = build_pileup(&reads, 200);
        let mean = pileup.mean_coverage();

        assert!(mean > 0.0);
    }

    #[test]
    fn test_normalize_pileup_cpm() {
        let reads = vec![
            ("chr1".to_string(), 0, 50),
            ("chr1".to_string(), 100, 50),
        ];

        let pileup = build_pileup(&reads, 200);
        let normalized = normalize_pileup(&pileup, "cpm").unwrap();

        assert!(normalized.coverage.contains_key("chr1"));
    }

    #[test]
    fn test_smooth_pileup() {
        let reads = vec![
            ("chr1".to_string(), 100, 50),
            ("chr1".to_string(), 101, 50),
            ("chr1".to_string(), 102, 50),
        ];

        let pileup = build_pileup(&reads, 200);
        let smoothed = smooth_pileup(&pileup, 5.0);

        assert!(smoothed.coverage.contains_key("chr1"));
    }

    #[test]
    fn test_pileup_correlation() {
        let reads1 = vec![
            ("chr1".to_string(), 100, 50),
            ("chr1".to_string(), 200, 50),
        ];

        let reads2 = vec![
            ("chr1".to_string(), 100, 50),
            ("chr1".to_string(), 200, 50),
        ];

        let pileup1 = build_pileup(&reads1, 200);
        let pileup2 = build_pileup(&reads2, 200);

        let corr = pileup_correlation(&pileup1, &pileup2);
        assert!(corr.is_some());
    }

    #[test]
    fn test_fingerprint() {
        let reads = vec![
            ("chr1".to_string(), 100, 50),
            ("chr1".to_string(), 150, 50),
        ];

        let pileup = build_pileup(&reads, 200);
        let fp = fingerprint(&pileup);

        assert!(!fp.is_empty());
        assert!(fp[0].0 <= 100.0);
    }
}
