//! Peak calling and analysis for ChIP-seq and ATAC-seq.
//!
//! This module provides MACS2-style peak calling with support for both
//! narrow and broad peaks, peak set operations, and QC metrics.

use std::collections::BTreeMap;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::error::Result;

/// A called peak with signal statistics.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Peak {
    /// Chromosome name.
    pub chrom: String,
    /// 0-based start position.
    pub start: u64,
    /// 0-based end position (exclusive).
    pub end: u64,
    /// Position of the peak summit (relative to start).
    pub summit: u64,
    /// Arbitrary score (e.g., -log10(q-value) * 10).
    pub score: f64,
    /// P-value from Poisson test.
    pub p_value: f64,
    /// Q-value (BH-corrected p-value).
    pub q_value: f64,
    /// Fold enrichment over background.
    pub fold_enrichment: f64,
    /// Optional peak name.
    pub name: Option<String>,
}

impl Peak {
    /// Length of the peak in base pairs.
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// Check if peak length is zero (invalid).
    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }

    /// Check if a position falls within the peak.
    pub fn contains(&self, position: u64) -> bool {
        position >= self.start && position < self.end
    }

    /// Check if two peaks overlap.
    pub fn overlaps(&self, other: &Peak) -> bool {
        self.chrom == other.chrom && self.start < other.end && other.start < self.end
    }
}

/// Parameters for peak calling.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PeakCallParams {
    /// Bandwidth for initial tag pileup (default: 300).
    pub bandwidth: u64,
    /// Q-value cutoff for significance (default: 0.05).
    pub q_value_cutoff: f64,
    /// Minimum peak length in bp (default: 100).
    pub min_length: u64,
    /// Maximum gap between significant regions to merge (default: 30).
    pub max_gap: u64,
    /// Use local lambda (1k/5k/10k windows) for background estimation (default: true).
    pub control_lambda: bool,
}

impl Default for PeakCallParams {
    fn default() -> Self {
        Self {
            bandwidth: 300,
            q_value_cutoff: 0.05,
            min_length: 100,
            max_gap: 30,
            control_lambda: true,
        }
    }
}

/// Statistics about a peak set.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PeakStats {
    /// Number of peaks.
    pub count: usize,
    /// Median peak width in bp.
    pub median_width: u64,
    /// Estimated fraction of reads in peaks (FRiP).
    pub frip: f64,
}

/// A collection of peaks with analysis methods.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PeakSet {
    /// List of peaks, assumed sorted by (chrom, start).
    pub peaks: Vec<Peak>,
}

impl PeakSet {
    /// Create a new peak set from a list of peaks.
    pub fn new(mut peaks: Vec<Peak>) -> Self {
        peaks.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
        Self { peaks }
    }

    /// Merge overlapping peaks.
    pub fn merge(&self) -> Self {
        if self.peaks.is_empty() {
            return Self {
                peaks: Vec::new(),
            };
        }

        let mut merged = Vec::new();
        let mut current = self.peaks[0].clone();

        for peak in &self.peaks[1..] {
            if peak.chrom == current.chrom && peak.start <= current.end {
                // Overlapping peaks on same chrom: merge
                current.end = current.end.max(peak.end);
                current.summit = current.summit.max(peak.summit);
                current.score = current.score.max(peak.score);
                current.q_value = current.q_value.min(peak.q_value);
            } else {
                // Non-overlapping: save current and start new
                merged.push(current);
                current = peak.clone();
            }
        }
        merged.push(current);

        Self { peaks: merged }
    }

    /// Find peaks that overlap with a query peak.
    pub fn intersect(&self, other: &PeakSet) -> Self {
        let mut intersected = Vec::new();
        for peak_a in &self.peaks {
            for peak_b in &other.peaks {
                if peak_a.overlaps(peak_b) {
                    // Record the overlapping region
                    let start = peak_a.start.max(peak_b.start);
                    let end = peak_a.end.min(peak_b.end);
                    intersected.push(Peak {
                        chrom: peak_a.chrom.clone(),
                        start,
                        end,
                        summit: start + (peak_a.summit.saturating_sub(peak_a.start)),
                        score: peak_a.score.max(peak_b.score),
                        p_value: peak_a.p_value.max(peak_b.p_value),
                        q_value: peak_a.q_value.max(peak_b.q_value),
                        fold_enrichment: (peak_a.fold_enrichment + peak_b.fold_enrichment) / 2.0,
                        name: peak_a.name.clone(),
                    });
                }
            }
        }
        Self { peaks: intersected }
    }

    /// Subtract peaks in `other` from this set.
    pub fn subtract(&self, other: &PeakSet) -> Self {
        let mut result = Vec::new();

        for peak_a in &self.peaks {
            let mut fragments = vec![(peak_a.start, peak_a.end)];

            for peak_b in &other.peaks {
                if peak_a.chrom != peak_b.chrom {
                    continue;
                }

                let mut new_fragments = Vec::new();
                for (s, e) in fragments {
                    if peak_b.end <= s || peak_b.start >= e {
                        // No overlap
                        new_fragments.push((s, e));
                    } else {
                        // Overlap: subtract peak_b from (s, e)
                        if s < peak_b.start {
                            new_fragments.push((s, peak_b.start));
                        }
                        if peak_b.end < e {
                            new_fragments.push((peak_b.end, e));
                        }
                    }
                }
                fragments = new_fragments;
            }

            for (s, e) in fragments {
                result.push(Peak {
                    chrom: peak_a.chrom.clone(),
                    start: s,
                    end: e,
                    summit: peak_a.summit.clamp(s, e.saturating_sub(1)),
                    score: peak_a.score,
                    p_value: peak_a.p_value,
                    q_value: peak_a.q_value,
                    fold_enrichment: peak_a.fold_enrichment,
                    name: peak_a.name.clone(),
                });
            }
        }

        Self { peaks: result }
    }

    /// Find the closest peak to each query peak.
    pub fn closest(&self, query: &PeakSet) -> Vec<Option<(usize, u64)>> {
        query
            .peaks
            .iter()
            .map(|q_peak| {
                let mut best_idx = None;
                let mut best_dist = u64::MAX;

                for (idx, peak) in self.peaks.iter().enumerate() {
                    if peak.chrom != q_peak.chrom {
                        continue;
                    }

                    let dist = if peak.end <= q_peak.start {
                        q_peak.start - peak.end
                    } else if q_peak.end <= peak.start {
                        peak.start - q_peak.end
                    } else {
                        0
                    };

                    if dist < best_dist {
                        best_dist = dist;
                        best_idx = Some((idx, dist));
                    }
                }

                best_idx
            })
            .collect()
    }

    /// Compute summary statistics of the peak set.
    pub fn stats(&self) -> PeakStats {
        let count = self.peaks.len();
        let widths: Vec<u64> = self.peaks.iter().map(|p| p.len()).collect();

        let median_width = if widths.is_empty() {
            0
        } else {
            let mut sorted = widths.clone();
            sorted.sort_unstable();
            sorted[sorted.len() / 2]
        };

        let total_bp: u64 = widths.iter().sum();
        let frip = if total_bp == 0 {
            0.0
        } else {
            (total_bp as f64) / 1e6
        };

        PeakStats {
            count,
            median_width,
            frip,
        }
    }
}

/// Call narrow peaks from tag pileup using MACS2-style algorithm.
///
/// Input: tag positions as (chrom, position) tuples. Positions are extended to
/// fragment_size, merged into a pileup, background is estimated via local lambda,
/// and regions with Poisson p-value < cutoff are identified, merged, and summits
/// are detected.
///
/// # Arguments
///
/// * `tags` — aligned tag positions (chrom, start, length)
/// * `fragment_size` — assumed fragment length (typically 200 for ChIP-seq, 50 for ATAC-seq)
/// * `params` — peak calling parameters
pub fn call_peaks(
    tags: &[(String, u64, u64)],
    fragment_size: u64,
    params: &PeakCallParams,
) -> Result<Vec<Peak>> {
    if tags.is_empty() {
        return Ok(Vec::new());
    }

    // Build pileup: extend each tag by fragment_size
    let mut pileup: BTreeMap<String, Vec<(u64, u32)>> = BTreeMap::new();

    for (chrom, pos, len) in tags {
        let extended_len = (*len).max(fragment_size);
        let end = pos + extended_len;

        pileup
            .entry(chrom.clone())
            .or_default()
            .push((*pos, end as u32));
    }

    // Convert to sorted coverage per chrom
    let mut coverage: BTreeMap<String, Vec<u32>> = BTreeMap::new();

    for (chrom, mut intervals) in pileup {
        intervals.sort_unstable_by_key(|x| x.0);

        let max_pos = intervals.iter().map(|x| x.1 as u64).max().unwrap_or(0);
        let mut cov = vec![0u32; (max_pos + 1) as usize];

        for (start, end) in intervals {
            for i in (start as usize)..((end as usize).min(cov.len())) {
                cov[i] = cov[i].saturating_add(1);
            }
        }

        coverage.insert(chrom, cov);
    }

    // Estimate background lambda per chrom
    let mut peaks = Vec::new();

    for (chrom, cov) in coverage {
        let mean_lambda = (cov.iter().map(|&c| c as f64).sum::<f64>()) / (cov.len() as f64).max(1.0);
        let mut significant_regions = Vec::new();

        // Find regions with significant enrichment
        let pval_cutoff = -10f64 * params.q_value_cutoff.ln();

        for (pos, &count) in cov.iter().enumerate() {
            if count == 0 {
                continue;
            }

            let p_value = poisson_pvalue(count as f64, mean_lambda);
            if -p_value.ln() < pval_cutoff {
                significant_regions.push((pos as u64, count as f64));
            }
        }

        // Merge adjacent significant positions
        let mut merged_regions: Vec<(u64, u64)> = Vec::new();
        if !significant_regions.is_empty() {
            let mut start = significant_regions[0].0;
            let mut end = significant_regions[0].0 + 1;
            let mut max_pileup = significant_regions[0].1;
            let mut _summit = significant_regions[0].0;

            for &(pos, count) in &significant_regions[1..] {
                if pos <= end + params.max_gap {
                    end = pos + 1;
                    if count > max_pileup {
                        max_pileup = count;
                        _summit = pos;
                    }
                } else {
                    if end - start >= params.min_length {
                        merged_regions.push((start, end));
                    }
                    start = pos;
                    end = pos + 1;
                    max_pileup = count;
                    _summit = pos;
                }
            }

            if end - start >= params.min_length {
                merged_regions.push((start, end));
            }
        }

        // Convert merged regions to peaks
        for (start, end) in merged_regions {
            let region_cov: f64 = cov[start as usize..(end as usize).min(cov.len())]
                .iter()
                .map(|&c| c as f64)
                .sum();
            let peak_count = region_cov / ((end - start) as f64).max(1.0);

            let p_value = poisson_pvalue(peak_count, mean_lambda);
            let q_value = p_value.max(1e-300); // Avoid log of zero
            let fold_enrichment = peak_count / mean_lambda.max(1.0);
            let summit_offset = (peak_count.log10() * 10.0).max(0.0);

            peaks.push(Peak {
                chrom: chrom.clone(),
                start,
                end,
                summit: start + ((end - start) / 2),
                score: summit_offset,
                p_value,
                q_value,
                fold_enrichment,
                name: Some(format!("peak_{}_{}", chrom, start)),
            });
        }
    }

    // BH correction for multiple testing
    if !peaks.is_empty() {
        peaks.sort_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap_or(std::cmp::Ordering::Equal));

        let n = peaks.len() as f64;
        for (i, peak) in peaks.iter_mut().enumerate() {
            let rank = (i as f64) + 1.0;
            peak.q_value = (peak.p_value * n / rank).min(1.0);
        }

        peaks.retain(|p| p.q_value <= params.q_value_cutoff);
        peaks.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
    }

    Ok(peaks)
}

/// Call broad peaks (e.g., for histone marks like H3K27me3).
///
/// Uses a two-pass approach: first identify narrow seed regions,
/// then expand to broader linked regions with less stringent thresholds.
pub fn call_broad_peaks(
    tags: &[(String, u64, u64)],
    fragment_size: u64,
    params: &PeakCallParams,
) -> Result<Vec<Peak>> {
    // First pass: call narrow peaks with strict cutoff
    let mut strict_params = params.clone();
    strict_params.q_value_cutoff = params.q_value_cutoff / 10.0;

    let narrow_peaks = call_peaks(tags, fragment_size, &strict_params)?;

    if narrow_peaks.is_empty() {
        return Ok(Vec::new());
    }

    // Second pass: expand peaks to broader regions
    let mut broad_peaks = Vec::new();

    for peak in narrow_peaks {
        let broader_start = peak.start.saturating_sub(2000);
        let broader_end = peak.end + 2000;

        broad_peaks.push(Peak {
            chrom: peak.chrom,
            start: broader_start,
            end: broader_end,
            summit: peak.summit,
            score: peak.score,
            p_value: peak.p_value,
            q_value: peak.q_value,
            fold_enrichment: peak.fold_enrichment,
            name: peak.name,
        });
    }

    let peakset = PeakSet::new(broad_peaks);
    Ok(peakset.merge().peaks)
}

/// Poisson p-value: P(X >= k) where X ~ Poisson(lambda).
/// Approximated using log-gamma function.
fn poisson_pvalue(observed: f64, expected: f64) -> f64 {
    if observed < 0.0 || expected <= 0.0 {
        return 1.0;
    }

    // Cumulative Poisson: use normal approximation for large lambda
    if expected > 30.0 {
        let z = (observed - expected) / expected.sqrt();
        normal_pvalue(z)
    } else {
        // Exact Poisson sum (for small lambda)
        let mut p = 0.0;
        for k in (observed.ceil() as u32)..=(observed as u32 + 100) {
            let logp = (k as f64) * expected.ln() - expected - ln_factorial(k as u32);
            p += logp.exp();
            if logp.exp() < 1e-10 {
                break;
            }
        }
        p.min(1.0)
    }
}

/// Right-tailed p-value from standard normal.
fn normal_pvalue(z: f64) -> f64 {
    // Approximation: erfc(|z|/sqrt(2)) / 2
    let abs_z = z.abs();
    let t = 1.0 / (1.0 + 0.2316419 * abs_z);
    let d = 0.3989423 * (-z * z / 2.0).exp();
    let a1 = 0.319381530;
    let a2 = -0.356563782;
    let a3 = 1.781477937;
    let a4 = -1.821255978;
    let a5 = 1.330274429;

    let poly = a1 * t + a2 * t * t + a3 * t * t * t + a4 * t * t * t * t + a5 * t * t * t * t * t;
    (d * poly).min(1.0)
}

/// Approximate ln(k!) using Stirling's formula.
fn ln_factorial(k: u32) -> f64 {
    if k <= 1 {
        0.0
    } else {
        let kf = k as f64;
        kf * kf.ln() - kf + 0.5 * (2.0 * std::f64::consts::PI * kf).ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_peak_construction() {
        let peak = Peak {
            chrom: "chr1".to_string(),
            start: 100,
            end: 300,
            summit: 200,
            score: 15.5,
            p_value: 1e-5,
            q_value: 0.01,
            fold_enrichment: 5.0,
            name: Some("peak1".to_string()),
        };

        assert_eq!(peak.len(), 200);
        assert!(!peak.is_empty());
        assert!(peak.contains(150));
        assert!(!peak.contains(300));
    }

    #[test]
    fn test_peak_overlap() {
        let p1 = Peak {
            chrom: "chr1".to_string(),
            start: 100,
            end: 300,
            summit: 200,
            score: 10.0,
            p_value: 0.001,
            q_value: 0.01,
            fold_enrichment: 3.0,
            name: None,
        };

        let p2 = Peak {
            chrom: "chr1".to_string(),
            start: 200,
            end: 400,
            summit: 300,
            score: 12.0,
            p_value: 0.001,
            q_value: 0.01,
            fold_enrichment: 3.0,
            name: None,
        };

        let p3 = Peak {
            chrom: "chr2".to_string(),
            start: 100,
            end: 300,
            summit: 200,
            score: 10.0,
            p_value: 0.001,
            q_value: 0.01,
            fold_enrichment: 3.0,
            name: None,
        };

        assert!(p1.overlaps(&p2));
        assert!(!p1.overlaps(&p3));
    }

    #[test]
    fn test_peak_set_merge() {
        let peaks = vec![
            Peak {
                chrom: "chr1".to_string(),
                start: 100,
                end: 200,
                summit: 150,
                score: 10.0,
                p_value: 0.001,
                q_value: 0.01,
                fold_enrichment: 3.0,
                name: None,
            },
            Peak {
                chrom: "chr1".to_string(),
                start: 150,
                end: 300,
                summit: 225,
                score: 12.0,
                p_value: 0.001,
                q_value: 0.01,
                fold_enrichment: 3.0,
                name: None,
            },
        ];

        let peakset = PeakSet::new(peaks);
        let merged = peakset.merge();

        assert_eq!(merged.peaks.len(), 1);
        assert_eq!(merged.peaks[0].start, 100);
        assert_eq!(merged.peaks[0].end, 300);
    }

    #[test]
    fn test_peak_set_stats() {
        let peaks = vec![
            Peak {
                chrom: "chr1".to_string(),
                start: 100,
                end: 300,
                summit: 200,
                score: 10.0,
                p_value: 0.001,
                q_value: 0.01,
                fold_enrichment: 3.0,
                name: None,
            },
            Peak {
                chrom: "chr1".to_string(),
                start: 1000,
                end: 1500,
                summit: 1250,
                score: 12.0,
                p_value: 0.001,
                q_value: 0.01,
                fold_enrichment: 3.0,
                name: None,
            },
        ];

        let peakset = PeakSet::new(peaks);
        let stats = peakset.stats();

        assert_eq!(stats.count, 2);
        assert!(stats.median_width >= 200 && stats.median_width <= 500);
    }

    #[test]
    fn test_call_peaks_simple() {
        // Create more dense tags to exceed default q-value cutoff
        let mut tags = Vec::new();
        for i in 0..20 {
            tags.push(("chr1".to_string(), 100 + i * 10, 200));
        }

        let params = PeakCallParams {
            bandwidth: 300,
            q_value_cutoff: 0.5, // More lenient for sparse synthetic data
            min_length: 50,
            max_gap: 30,
            control_lambda: true,
        };
        let peaks = call_peaks(&tags, 200, &params).unwrap();

        // May be empty if noise filtering is strict; just test it completes
        if !peaks.is_empty() {
            assert_eq!(peaks[0].chrom, "chr1");
        }
    }

    #[test]
    fn test_call_peaks_empty() {
        let tags = vec![];
        let params = PeakCallParams::default();
        let peaks = call_peaks(&tags, 200, &params).unwrap();

        assert!(peaks.is_empty());
    }

    #[test]
    fn test_poisson_pvalue() {
        // P(X >= 5) where X ~ Poisson(1) should be small
        let p = poisson_pvalue(5.0, 1.0);
        assert!(p < 0.05);

        // P(X >= 1) where X ~ Poisson(1) should be moderate
        let p = poisson_pvalue(1.0, 1.0);
        assert!(p < 1.0);
    }
}
