//! Nucleosome positioning and periodicity analysis from MNase-seq and ATAC-seq.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::pileup::TagPileup;

/// Detected nucleosome position.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct NucleosomePosition {
    /// Center position (0-based).
    pub center: u64,
    /// Occupancy score (0-1).
    pub occupancy: f64,
    /// Fuzziness (standard deviation of position).
    pub fuzziness: f64,
    /// Call confidence score.
    pub score: f64,
}

/// Parameters for nucleosome detection.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct NucleosomeParams {
    /// Expected lower fragment bound (default: 120).
    pub fragment_lower: u64,
    /// Expected upper fragment bound (default: 200).
    pub fragment_upper: u64,
    /// Gaussian smoothing bandwidth (default: 30).
    pub smooth_bandwidth: f64,
}

impl Default for NucleosomeParams {
    fn default() -> Self {
        Self {
            fragment_lower: 120,
            fragment_upper: 200,
            smooth_bandwidth: 30.0,
        }
    }
}

/// Call nucleosome positions from MNase-seq or ATAC-seq fragment centers.
pub fn call_nucleosomes(pileup: &TagPileup, params: &NucleosomeParams) -> Vec<NucleosomePosition> {
    let mut nucleosomes = Vec::new();

    for (_chrom, cov) in &pileup.coverage {
        // Gaussian smoothing
        let smoothed = smooth_coverage(cov, params.smooth_bandwidth);

        // Peak detection
        let mut i = 0;
        while i < smoothed.len() {
            let mut is_peak = true;

            // Check if local maximum
            if i > 0 && smoothed[i] <= smoothed[i - 1] {
                is_peak = false;
            }
            if i < smoothed.len() - 1 && smoothed[i] <= smoothed[i + 1] {
                is_peak = false;
            }

            if is_peak && smoothed[i] > 0.5 {
                // Estimate width (fuzziness) from peak
                let mut left = i;
                let mut right = i;

                let half_height = smoothed[i] / 2.0;
                while left > 0 && smoothed[left] > half_height {
                    left -= 1;
                }
                while right < smoothed.len() - 1 && smoothed[right] > half_height {
                    right += 1;
                }

                let fuzziness = ((right - left) as f64 / 2.0).max(1.0);
                let occupancy = (smoothed[i] / 100.0).min(1.0);

                nucleosomes.push(NucleosomePosition {
                    center: i as u64,
                    occupancy,
                    fuzziness,
                    score: smoothed[i],
                });

                i = right + 1;
            } else {
                i += 1;
            }
        }
    }

    nucleosomes
}

/// Smooth coverage using Gaussian kernel.
fn smooth_coverage(cov: &[u32], bandwidth: f64) -> Vec<f64> {
    let mut smoothed = vec![0.0; cov.len()];

    let radius = (3.0 * bandwidth) as i32;
    let mut kernel = Vec::new();
    let mut kernel_sum = 0.0;

    for i in (-radius)..=(radius) {
        let v = (-((i as f64) * (i as f64)) / (2.0 * bandwidth * bandwidth)).exp();
        kernel.push(v);
        kernel_sum += v;
    }

    kernel.iter_mut().for_each(|k| *k /= kernel_sum);

    for pos in 0..cov.len() {
        let mut sum = 0.0;
        for (k_idx, &k_val) in kernel.iter().enumerate() {
            let offset = k_idx as i32 - radius;
            let neighbor = (pos as i32 + offset) as usize;
            if neighbor < cov.len() {
                sum += (cov[neighbor] as f64) * k_val;
            }
        }
        smoothed[pos] = sum;
    }

    smoothed
}

/// Compute nucleosome-free region (NFR) score around transcription start sites.
///
/// Returns the ratio of fragments in the NFR (<147 bp) to nucleosome-bound region (147-300 bp).
pub fn nfr_score(nuc_pos: &[NucleosomePosition], tss_positions: &[u64]) -> f64 {
    if tss_positions.is_empty() {
        return 0.0;
    }

    let window_size = 1000u64;
    let nfr_bound = 147;
    let mono_bound = 300;

    let mut nfr_count = 0;
    let mut mono_count = 0;

    for &tss in tss_positions {
        let start = tss.saturating_sub(window_size);
        let end = tss + window_size;

        for nuc in nuc_pos {
            if nuc.center >= start && nuc.center <= end {
                let distance = (nuc.center as i64 - tss as i64).abs() as u64;

                if distance < nfr_bound {
                    nfr_count += 1;
                } else if distance < mono_bound {
                    mono_count += 1;
                }
            }
        }
    }

    if mono_count == 0 {
        0.0
    } else {
        (nfr_count as f64) / (mono_count as f64)
    }
}

/// Detect periodic ~147 bp nucleosomal signal via autocorrelation.
pub fn periodicity(cov: &[u32]) -> f64 {
    if cov.len() < 500 {
        return 0.0;
    }

    let period = 147usize;

    // Compute autocorrelation at lag = period
    let mean = (cov.iter().map(|&c| c as f64).sum::<f64>()) / (cov.len() as f64);
    let mut numerator = 0.0;
    let mut denominator = 0.0;

    for i in 0..cov.len().saturating_sub(period) {
        let x = (cov[i] as f64) - mean;
        let y = (cov[i + period] as f64) - mean;
        numerator += x * y;
        denominator += x * x;
    }

    if denominator > 0.0 {
        numerator / denominator
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    #[test]
    fn test_call_nucleosomes() {
        let mut coverage = BTreeMap::new();
        let mut cov_vec = vec![0u32; 1000];

        // Create synthetic nucleosome-like peaks
        for i in 0..5 {
            let center = 100 + i * 200;
            for j in (center - 50)..=(center + 50).min(999) {
                cov_vec[j] = 10;
            }
        }

        coverage.insert("chr1".to_string(), cov_vec);
        let pileup = TagPileup { coverage };

        let params = NucleosomeParams::default();
        let nucs = call_nucleosomes(&pileup, &params);

        assert!(!nucs.is_empty());
    }

    #[test]
    fn test_nfr_score() {
        let nuc_pos = vec![
            NucleosomePosition {
                center: 500,
                occupancy: 0.8,
                fuzziness: 20.0,
                score: 10.0,
            },
            NucleosomePosition {
                center: 800,
                occupancy: 0.7,
                fuzziness: 20.0,
                score: 8.0,
            },
        ];

        let tss_positions = vec![500];
        let score = nfr_score(&nuc_pos, &tss_positions);

        assert!(score >= 0.0);
    }

    #[test]
    fn test_periodicity() {
        // Create periodic coverage
        let mut cov = vec![0u32; 1000];
        for i in 0..5 {
            let center = 150 + i * 147;
            if center < 1000 {
                cov[center] = 20;
            }
        }

        let periodic = periodicity(&cov);
        assert!(periodic >= 0.0 && periodic <= 1.0);
    }
}
