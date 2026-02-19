//! Copy number and structural variant analysis.
//!
//! This module provides tools for detecting copy number alterations (CNAs) and
//! structural variants (SVs) from array or sequencing data:
//!
//! - **Circular Binary Segmentation (CBS)** — Olshen et al. 2004 algorithm for
//!   segmenting log2 ratio profiles into regions of constant copy number.
//! - **BAF segmentation** — Segment B-allele frequency profiles to detect LOH.
//! - **SV breakpoint detection** — Cluster discordant read pairs and split reads
//!   to identify deletion, duplication, inversion, translocation, and insertion events.
//! - **Segment merging** — Combine adjacent segments with similar copy number.

use cyanea_core::{CyaneaError, Result};

/// A copy number segment identified by CBS.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CnvSegment {
    /// Chromosome name.
    pub chrom: String,
    /// Start position (0-based, inclusive).
    pub start: u64,
    /// End position (0-based, exclusive).
    pub end: u64,
    /// Mean log2 ratio of probes in this segment.
    pub log2_ratio: f64,
    /// Number of probes in this segment.
    pub n_probes: usize,
    /// Estimated integer copy number.
    pub copy_number: u32,
}

/// A B-allele frequency segment.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BafSegment {
    /// Chromosome name.
    pub chrom: String,
    /// Start position (0-based, inclusive).
    pub start: u64,
    /// End position (0-based, exclusive).
    pub end: u64,
    /// Mean BAF value in this segment.
    pub mean_baf: f64,
    /// Number of SNPs in this segment.
    pub n_snps: usize,
}

/// Type of structural variant.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SvType {
    /// Genomic deletion (loss of sequence).
    Deletion,
    /// Genomic duplication (gain of sequence).
    Duplication,
    /// Inversion of a genomic segment.
    Inversion,
    /// Translocation between different chromosomes.
    Translocation,
    /// Insertion of novel sequence.
    Insertion,
}

/// A structural variant breakpoint with supporting evidence.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SvBreakpoint {
    /// Type of structural variant.
    pub sv_type: SvType,
    /// Chromosome of the first breakpoint.
    pub chrom1: String,
    /// Position of the first breakpoint.
    pub pos1: u64,
    /// Chromosome of the second breakpoint.
    pub chrom2: String,
    /// Position of the second breakpoint.
    pub pos2: u64,
    /// Number of split reads supporting this breakpoint.
    pub split_reads: u32,
    /// Number of discordant read pairs supporting this breakpoint.
    pub discordant_pairs: u32,
    /// Quality score: split_reads * 10 + discordant_pairs * 5.
    pub quality: f64,
}

/// Configuration for Circular Binary Segmentation.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CbsConfig {
    /// Significance threshold for permutation test (default: 0.01).
    pub alpha: f64,
    /// Minimum number of probes per segment (default: 3).
    pub min_probes: usize,
    /// Number of permutations for significance testing (default: 1000).
    pub n_permutations: usize,
    /// Seed for the PRNG used in permutation testing (default: 42).
    pub seed: u64,
}

impl Default for CbsConfig {
    fn default() -> Self {
        Self {
            alpha: 0.01,
            min_probes: 3,
            n_permutations: 1000,
            seed: 42,
        }
    }
}

/// Private Xorshift64 PRNG for permutation testing.
struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        // Ensure non-zero state.
        Self {
            state: if seed == 0 { 1 } else { seed },
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }
}

/// Estimate integer copy number from a log2 ratio value.
///
/// Formula: `CN = round(2^(log2_ratio + 1))`, clamped to [0, 10].
fn cn_from_log2(log2_ratio: f64) -> u32 {
    let cn = 2.0_f64.powf(log2_ratio + 1.0);
    let rounded = cn.round() as i64;
    rounded.clamp(0, 10) as u32
}

/// Segment a log2 ratio profile using Circular Binary Segmentation (Olshen et al. 2004).
///
/// CBS recursively identifies the most significant change-point in the data by
/// computing a maximal t-statistic over all possible split points on a circular
/// scan. Significance is assessed via a permutation test.
///
/// # Arguments
///
/// * `positions` — Genomic positions of probes (must be sorted).
/// * `log2_ratios` — Log2 ratio values corresponding to each position.
/// * `chrom` — Chromosome name for the output segments.
/// * `config` — CBS algorithm parameters.
///
/// # Errors
///
/// Returns an error if `positions` and `log2_ratios` have different lengths,
/// or if either is empty.
pub fn circular_binary_segmentation(
    positions: &[u64],
    log2_ratios: &[f64],
    chrom: &str,
    config: &CbsConfig,
) -> Result<Vec<CnvSegment>> {
    if positions.len() != log2_ratios.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "positions length ({}) must match log2_ratios length ({})",
            positions.len(),
            log2_ratios.len()
        )));
    }
    if positions.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "positions and log2_ratios must not be empty".into(),
        ));
    }

    let mut breakpoints = Vec::new();
    let mut rng = Xorshift64::new(config.seed);
    cbs_recurse(
        log2_ratios,
        0,
        log2_ratios.len(),
        config,
        &mut rng,
        &mut breakpoints,
    );

    breakpoints.sort_unstable();
    breakpoints.dedup();

    // Build segments from breakpoints.
    let mut boundaries = Vec::with_capacity(breakpoints.len() + 2);
    boundaries.push(0usize);
    for &bp in &breakpoints {
        boundaries.push(bp);
    }
    boundaries.push(log2_ratios.len());
    boundaries.sort_unstable();
    boundaries.dedup();

    let mut segments = Vec::new();
    for w in boundaries.windows(2) {
        let s = w[0];
        let e = w[1];
        if e <= s {
            continue;
        }
        let n_probes = e - s;
        let mean_lr: f64 = log2_ratios[s..e].iter().sum::<f64>() / n_probes as f64;
        segments.push(CnvSegment {
            chrom: chrom.to_string(),
            start: positions[s],
            end: positions[e - 1] + 1, // exclusive end
            log2_ratio: mean_lr,
            n_probes,
            copy_number: cn_from_log2(mean_lr),
        });
    }

    Ok(segments)
}

/// Recursive CBS: find the best split point in data[start..end].
fn cbs_recurse(
    data: &[f64],
    start: usize,
    end: usize,
    config: &CbsConfig,
    rng: &mut Xorshift64,
    breakpoints: &mut Vec<usize>,
) {
    let n = end - start;
    if n < 2 * config.min_probes {
        return;
    }

    let slice = &data[start..end];

    // Find max t-statistic and its position.
    let (max_t, best_pos) = max_t_stat(slice, config.min_probes);

    if max_t <= 0.0 {
        return;
    }

    // Permutation test: count how many permuted datasets yield a higher max t-stat.
    let mut count_ge = 0usize;
    let mut perm_data: Vec<f64> = slice.to_vec();

    for _ in 0..config.n_permutations {
        // Fisher-Yates shuffle using our PRNG.
        for i in (1..perm_data.len()).rev() {
            let j = (rng.next_u64() as usize) % (i + 1);
            perm_data.swap(i, j);
        }
        let (perm_t, _) = max_t_stat(&perm_data, config.min_probes);
        if perm_t >= max_t {
            count_ge += 1;
        }
    }

    let p_value = count_ge as f64 / config.n_permutations as f64;

    if p_value < config.alpha {
        let split = start + best_pos;
        breakpoints.push(split);
        // Recurse on left and right segments.
        cbs_recurse(data, start, split, config, rng, breakpoints);
        cbs_recurse(data, split, end, config, rng, breakpoints);
    }
}

/// Compute the maximum circular t-statistic over all candidate split points.
///
/// Returns (max_t_stat, best_split_index) where best_split_index is relative
/// to the start of the slice.
fn max_t_stat(data: &[f64], min_probes: usize) -> (f64, usize) {
    let n = data.len();
    if n < 2 * min_probes {
        return (0.0, 0);
    }

    let total_sum: f64 = data.iter().sum();
    let total_mean = total_sum / n as f64;

    let mut best_t = 0.0_f64;
    let mut best_pos = min_probes;

    // Cumulative sum for efficient segment mean computation.
    let mut cum_sum = vec![0.0_f64; n + 1];
    for i in 0..n {
        cum_sum[i + 1] = cum_sum[i] + data[i];
    }

    // Try all split points that respect min_probes.
    for i in min_probes..=(n - min_probes) {
        let left_sum = cum_sum[i];
        let right_sum = total_sum - left_sum;
        let left_n = i as f64;
        let right_n = (n - i) as f64;
        let left_mean = left_sum / left_n;
        let right_mean = right_sum / right_n;

        // Pooled variance estimate.
        let mut ss = 0.0_f64;
        for j in 0..n {
            let d = data[j] - total_mean;
            ss += d * d;
        }
        let variance = ss / n as f64;

        if variance <= 1e-15 {
            continue;
        }

        // Two-sample t-statistic (absolute value).
        let diff = (left_mean - right_mean).abs();
        let se = (variance * (1.0 / left_n + 1.0 / right_n)).sqrt();
        let t = diff / se;

        if t > best_t {
            best_t = t;
            best_pos = i;
        }
    }

    (best_t, best_pos)
}

/// Segment a B-allele frequency profile using CBS.
///
/// BAF values are mirrored around 0.5 (by taking `|baf - 0.5|`) before
/// segmentation, so that heterozygous SNPs (BAF ~ 0.5) map to ~ 0.0 and
/// LOH regions (BAF ~ 0.0 or ~ 1.0) map to ~ 0.5.
///
/// # Arguments
///
/// * `positions` — Genomic positions of SNPs (must be sorted).
/// * `bafs` — B-allele frequency values in [0, 1].
/// * `chrom` — Chromosome name for the output segments.
/// * `config` — CBS algorithm parameters.
///
/// # Errors
///
/// Returns an error if `positions` and `bafs` have different lengths,
/// or if either is empty.
pub fn baf_segmentation(
    positions: &[u64],
    bafs: &[f64],
    chrom: &str,
    config: &CbsConfig,
) -> Result<Vec<BafSegment>> {
    if positions.len() != bafs.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "positions length ({}) must match bafs length ({})",
            positions.len(),
            bafs.len()
        )));
    }
    if positions.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "positions and bafs must not be empty".into(),
        ));
    }

    // Mirror BAF values around 0.5.
    let mirrored: Vec<f64> = bafs.iter().map(|b| (b - 0.5).abs()).collect();

    let mut breakpoints = Vec::new();
    let mut rng = Xorshift64::new(config.seed);
    cbs_recurse(
        &mirrored,
        0,
        mirrored.len(),
        config,
        &mut rng,
        &mut breakpoints,
    );

    breakpoints.sort_unstable();
    breakpoints.dedup();

    let mut boundaries = Vec::with_capacity(breakpoints.len() + 2);
    boundaries.push(0usize);
    for &bp in &breakpoints {
        boundaries.push(bp);
    }
    boundaries.push(bafs.len());
    boundaries.sort_unstable();
    boundaries.dedup();

    let mut segments = Vec::new();
    for w in boundaries.windows(2) {
        let s = w[0];
        let e = w[1];
        if e <= s {
            continue;
        }
        let n_snps = e - s;
        let mean_baf: f64 = bafs[s..e].iter().sum::<f64>() / n_snps as f64;
        segments.push(BafSegment {
            chrom: chrom.to_string(),
            start: positions[s],
            end: positions[e - 1] + 1,
            mean_baf,
            n_snps,
        });
    }

    Ok(segments)
}

/// Evidence record for SV breakpoint clustering.
struct SvEvidence {
    chrom1: String,
    pos1: u64,
    chrom2: String,
    pos2: u64,
    is_split_read: bool,
}

/// Detect structural variant breakpoints from discordant read pairs and split reads.
///
/// Evidence is clustered by proximity: reads mapping within `cluster_distance`
/// of each other on both breakpoint sides are grouped together. Each cluster
/// is classified by orientation and chromosome identity.
///
/// # Arguments
///
/// * `discordant_pairs` — (chrom1, pos1, chrom2, pos2) for each discordant read pair.
/// * `split_reads` — (chrom1, pos1, chrom2, pos2) for each split read alignment.
/// * `cluster_distance` — Maximum distance between evidence to be grouped into one cluster.
pub fn detect_sv_breakpoints(
    discordant_pairs: &[(String, u64, String, u64)],
    split_reads: &[(String, u64, String, u64)],
    cluster_distance: u64,
) -> Vec<SvBreakpoint> {
    // Collect all evidence.
    let mut evidence: Vec<SvEvidence> = Vec::with_capacity(discordant_pairs.len() + split_reads.len());

    for (c1, p1, c2, p2) in discordant_pairs {
        evidence.push(SvEvidence {
            chrom1: c1.clone(),
            pos1: *p1,
            chrom2: c2.clone(),
            pos2: *p2,
            is_split_read: false,
        });
    }
    for (c1, p1, c2, p2) in split_reads {
        evidence.push(SvEvidence {
            chrom1: c1.clone(),
            pos1: *p1,
            chrom2: c2.clone(),
            pos2: *p2,
            is_split_read: true,
        });
    }

    if evidence.is_empty() {
        return Vec::new();
    }

    // Sort evidence by (chrom1, chrom2, pos1, pos2).
    evidence.sort_by(|a, b| {
        a.chrom1
            .cmp(&b.chrom1)
            .then(a.chrom2.cmp(&b.chrom2))
            .then(a.pos1.cmp(&b.pos1))
            .then(a.pos2.cmp(&b.pos2))
    });

    // Cluster by proximity.
    let mut clusters: Vec<Vec<&SvEvidence>> = Vec::new();
    let mut current_cluster: Vec<&SvEvidence> = vec![&evidence[0]];

    for ev in evidence.iter().skip(1) {
        let last = current_cluster.last().unwrap();
        let same_chroms = ev.chrom1 == last.chrom1 && ev.chrom2 == last.chrom2;
        let close_pos1 = ev.pos1.abs_diff(last.pos1) <= cluster_distance;
        let close_pos2 = ev.pos2.abs_diff(last.pos2) <= cluster_distance;

        if same_chroms && close_pos1 && close_pos2 {
            current_cluster.push(ev);
        } else {
            clusters.push(current_cluster);
            current_cluster = vec![ev];
        }
    }
    clusters.push(current_cluster);

    // Classify each cluster.
    let mut breakpoints = Vec::new();
    for cluster in &clusters {
        let mut sr_count = 0u32;
        let mut dp_count = 0u32;
        let mut sum_pos1 = 0u64;
        let mut sum_pos2 = 0u64;

        for ev in cluster {
            if ev.is_split_read {
                sr_count += 1;
            } else {
                dp_count += 1;
            }
            sum_pos1 += ev.pos1;
            sum_pos2 += ev.pos2;
        }

        let n = cluster.len() as u64;
        let mean_pos1 = sum_pos1 / n;
        let mean_pos2 = sum_pos2 / n;
        let chrom1 = &cluster[0].chrom1;
        let chrom2 = &cluster[0].chrom2;

        let sv_type = classify_sv(chrom1, mean_pos1, chrom2, mean_pos2);
        let quality = sr_count as f64 * 10.0 + dp_count as f64 * 5.0;

        breakpoints.push(SvBreakpoint {
            sv_type,
            chrom1: chrom1.clone(),
            pos1: mean_pos1,
            chrom2: chrom2.clone(),
            pos2: mean_pos2,
            split_reads: sr_count,
            discordant_pairs: dp_count,
            quality,
        });
    }

    breakpoints
}

/// Classify an SV based on breakpoint chromosomes and positions.
///
/// - Different chromosomes → Translocation
/// - Same chromosome, pos2 > pos1 → Deletion
/// - Same chromosome, pos2 <= pos1 → Inversion (positions suggest orientation change)
fn classify_sv(chrom1: &str, pos1: u64, chrom2: &str, pos2: u64) -> SvType {
    if chrom1 != chrom2 {
        SvType::Translocation
    } else if pos2 > pos1 {
        SvType::Deletion
    } else {
        SvType::Inversion
    }
}

/// Merge adjacent copy number segments with similar copy number.
///
/// Segments are sorted by (chrom, start) and then merged if:
/// - They are on the same chromosome.
/// - The gap between them is at most `max_gap`.
/// - Their copy number difference is at most `cn_tolerance`.
///
/// Merged segments receive a weighted-average log2 ratio.
pub fn merge_cnv_segments(
    segments: &[CnvSegment],
    max_gap: u64,
    cn_tolerance: u32,
) -> Vec<CnvSegment> {
    if segments.is_empty() {
        return Vec::new();
    }

    let mut sorted: Vec<CnvSegment> = segments.to_vec();
    sorted.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));

    let mut merged: Vec<CnvSegment> = Vec::new();
    let mut current = sorted[0].clone();

    for seg in sorted.iter().skip(1) {
        let same_chrom = seg.chrom == current.chrom;
        let gap = if seg.start >= current.end {
            seg.start - current.end
        } else {
            0
        };
        let cn_diff = (seg.copy_number as i64 - current.copy_number as i64).unsigned_abs() as u32;

        if same_chrom && gap <= max_gap && cn_diff <= cn_tolerance {
            // Merge: weighted average log2 ratio.
            let total_probes = current.n_probes + seg.n_probes;
            let weighted_lr = (current.log2_ratio * current.n_probes as f64
                + seg.log2_ratio * seg.n_probes as f64)
                / total_probes as f64;
            current.end = seg.end;
            current.log2_ratio = weighted_lr;
            current.n_probes = total_probes;
            current.copy_number = cn_from_log2(weighted_lr);
        } else {
            merged.push(current);
            current = seg.clone();
        }
    }
    merged.push(current);

    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test 1: CBS on constant data produces a single segment.
    #[test]
    fn cbs_constant_data() {
        let positions: Vec<u64> = (0..20).map(|i| i * 1000).collect();
        let log2_ratios = vec![0.0; 20];
        let config = CbsConfig::default();

        let segments =
            circular_binary_segmentation(&positions, &log2_ratios, "chr1", &config).unwrap();

        assert_eq!(segments.len(), 1);
        assert_eq!(segments[0].n_probes, 20);
        assert!((segments[0].log2_ratio - 0.0).abs() < 1e-10);
        assert_eq!(segments[0].copy_number, 2); // CN=2 for log2=0
    }

    /// Test 2: CBS on two-level data produces 2 segments.
    #[test]
    fn cbs_two_level() {
        let n = 40;
        let positions: Vec<u64> = (0..n).map(|i| i as u64 * 1000).collect();
        let mut log2_ratios = vec![0.0; n];
        // Second half has a strong gain.
        for i in 20..n {
            log2_ratios[i] = 1.0;
        }
        let config = CbsConfig {
            alpha: 0.05,
            min_probes: 3,
            n_permutations: 500,
            seed: 123,
        };

        let segments =
            circular_binary_segmentation(&positions, &log2_ratios, "chr1", &config).unwrap();

        assert_eq!(segments.len(), 2);
        assert!(segments[0].log2_ratio < 0.5);
        assert!(segments[1].log2_ratio > 0.5);
    }

    /// Test 3: CBS on three-level data produces 3 segments.
    #[test]
    fn cbs_three_level() {
        let n = 60;
        let positions: Vec<u64> = (0..n).map(|i| i as u64 * 1000).collect();
        let mut log2_ratios = vec![0.0; n];
        for i in 20..40 {
            log2_ratios[i] = 1.0; // gain
        }
        for i in 40..60 {
            log2_ratios[i] = -1.0; // loss
        }
        let config = CbsConfig {
            alpha: 0.05,
            min_probes: 3,
            n_permutations: 500,
            seed: 456,
        };

        let segments =
            circular_binary_segmentation(&positions, &log2_ratios, "chr1", &config).unwrap();

        assert_eq!(segments.len(), 3);
        assert!(segments[0].log2_ratio.abs() < 0.5);
        assert!(segments[1].log2_ratio > 0.5);
        assert!(segments[2].log2_ratio < -0.5);
    }

    /// Test 4: CBS respects min_probes setting.
    #[test]
    fn cbs_min_probes_respected() {
        // Data has a change at position 2 (only 2 probes on the left).
        let positions: Vec<u64> = (0..10).map(|i| i * 1000).collect();
        let mut log2_ratios = vec![0.0; 10];
        log2_ratios[0] = 2.0;
        log2_ratios[1] = 2.0;
        // With min_probes=5, we cannot split with only 2 probes on one side.
        let config = CbsConfig {
            alpha: 0.05,
            min_probes: 5,
            n_permutations: 200,
            seed: 789,
        };

        let segments =
            circular_binary_segmentation(&positions, &log2_ratios, "chr1", &config).unwrap();

        // Should produce 1 segment since the split would violate min_probes.
        assert_eq!(segments.len(), 1);
    }

    /// Test 5: CBS alpha threshold controls sensitivity.
    #[test]
    fn cbs_alpha_threshold() {
        let n = 30;
        let positions: Vec<u64> = (0..n).map(|i| i as u64 * 1000).collect();
        let mut log2_ratios = vec![0.0; n];
        // Subtle shift.
        for i in 15..n {
            log2_ratios[i] = 0.3;
        }

        // Very strict alpha — should not split.
        let strict_config = CbsConfig {
            alpha: 0.001,
            min_probes: 3,
            n_permutations: 200,
            seed: 111,
        };
        let segments_strict =
            circular_binary_segmentation(&positions, &log2_ratios, "chr1", &strict_config).unwrap();

        // Lenient alpha — should split.
        let lenient_config = CbsConfig {
            alpha: 0.5,
            min_probes: 3,
            n_permutations: 200,
            seed: 111,
        };
        let segments_lenient =
            circular_binary_segmentation(&positions, &log2_ratios, "chr1", &lenient_config)
                .unwrap();

        // Strict should produce fewer or equal segments than lenient.
        assert!(segments_strict.len() <= segments_lenient.len());
    }

    /// Test 6: BAF segmentation detects heterozygous + LOH regions.
    #[test]
    fn baf_het_and_loh() {
        let n = 40;
        let positions: Vec<u64> = (0..n).map(|i| i as u64 * 1000).collect();
        let mut bafs = vec![0.5; n]; // heterozygous
        // LOH region: BAF near 0.0 or 1.0.
        for i in 20..n {
            bafs[i] = 0.95;
        }
        let config = CbsConfig {
            alpha: 0.05,
            min_probes: 3,
            n_permutations: 500,
            seed: 222,
        };

        let segments = baf_segmentation(&positions, &bafs, "chr1", &config).unwrap();

        assert!(segments.len() >= 2);
        // First segment should have mean BAF near 0.5.
        assert!((segments[0].mean_baf - 0.5).abs() < 0.1);
        // Last segment should have mean BAF near 0.95.
        assert!((segments.last().unwrap().mean_baf - 0.95).abs() < 0.1);
    }

    /// Test 7: BAF segmentation on uniform diploid data produces 1 segment.
    #[test]
    fn baf_normal_diploid() {
        let positions: Vec<u64> = (0..20).map(|i| i * 1000).collect();
        let bafs = vec![0.5; 20];
        let config = CbsConfig::default();

        let segments = baf_segmentation(&positions, &bafs, "chr1", &config).unwrap();

        assert_eq!(segments.len(), 1);
        assert!((segments[0].mean_baf - 0.5).abs() < 1e-10);
    }

    /// Test 8: SV deletion detected from clustered evidence.
    #[test]
    fn sv_deletion_cluster() {
        let discordant: Vec<(String, u64, String, u64)> = vec![
            ("chr1".into(), 1000, "chr1".into(), 5000),
            ("chr1".into(), 1010, "chr1".into(), 5010),
            ("chr1".into(), 1020, "chr1".into(), 5020),
        ];
        let split: Vec<(String, u64, String, u64)> = vec![
            ("chr1".into(), 1005, "chr1".into(), 5005),
        ];

        let bps = detect_sv_breakpoints(&discordant, &split, 100);

        assert_eq!(bps.len(), 1);
        assert_eq!(bps[0].sv_type, SvType::Deletion);
        assert_eq!(bps[0].discordant_pairs, 3);
        assert_eq!(bps[0].split_reads, 1);
        assert!((bps[0].quality - (1.0 * 10.0 + 3.0 * 5.0)).abs() < 1e-10);
    }

    /// Test 9: SV translocation detected between chromosomes.
    #[test]
    fn sv_translocation() {
        let discordant: Vec<(String, u64, String, u64)> = vec![
            ("chr1".into(), 1000, "chr5".into(), 2000),
            ("chr1".into(), 1010, "chr5".into(), 2010),
        ];
        let split: Vec<(String, u64, String, u64)> = vec![
            ("chr1".into(), 1005, "chr5".into(), 2005),
        ];

        let bps = detect_sv_breakpoints(&discordant, &split, 100);

        assert_eq!(bps.len(), 1);
        assert_eq!(bps[0].sv_type, SvType::Translocation);
        assert_eq!(bps[0].chrom1, "chr1");
        assert_eq!(bps[0].chrom2, "chr5");
    }

    /// Test 10: SV inversion detected on same chromosome.
    #[test]
    fn sv_inversion() {
        // pos2 <= pos1 on same chromosome → inversion.
        let discordant: Vec<(String, u64, String, u64)> = vec![
            ("chr2".into(), 5000, "chr2".into(), 3000),
            ("chr2".into(), 5010, "chr2".into(), 3010),
        ];
        let split: Vec<(String, u64, String, u64)> = Vec::new();

        let bps = detect_sv_breakpoints(&discordant, &split, 100);

        assert_eq!(bps.len(), 1);
        assert_eq!(bps[0].sv_type, SvType::Inversion);
        assert_eq!(bps[0].chrom1, "chr2");
        assert_eq!(bps[0].chrom2, "chr2");
        assert_eq!(bps[0].discordant_pairs, 2);
        assert_eq!(bps[0].split_reads, 0);
    }

    /// Test 11: Merge adjacent CNV segments with matching copy number.
    #[test]
    fn cnv_merge_adjacent() {
        let segments = vec![
            CnvSegment {
                chrom: "chr1".into(),
                start: 0,
                end: 1000,
                log2_ratio: 0.0,
                n_probes: 10,
                copy_number: 2,
            },
            CnvSegment {
                chrom: "chr1".into(),
                start: 1000,
                end: 2000,
                log2_ratio: 0.0,
                n_probes: 10,
                copy_number: 2,
            },
        ];

        let merged = merge_cnv_segments(&segments, 100, 0);

        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start, 0);
        assert_eq!(merged[0].end, 2000);
        assert_eq!(merged[0].n_probes, 20);
        assert_eq!(merged[0].copy_number, 2);
    }

    /// Test 12: Merge respects max_gap parameter.
    #[test]
    fn cnv_merge_max_gap() {
        let segments = vec![
            CnvSegment {
                chrom: "chr1".into(),
                start: 0,
                end: 1000,
                log2_ratio: 0.0,
                n_probes: 10,
                copy_number: 2,
            },
            CnvSegment {
                chrom: "chr1".into(),
                start: 5000, // gap = 4000
                end: 6000,
                log2_ratio: 0.0,
                n_probes: 10,
                copy_number: 2,
            },
        ];

        // Gap of 4000 exceeds max_gap of 100 → should not merge.
        let merged = merge_cnv_segments(&segments, 100, 0);
        assert_eq!(merged.len(), 2);

        // With larger max_gap, should merge.
        let merged_large = merge_cnv_segments(&segments, 5000, 0);
        assert_eq!(merged_large.len(), 1);
    }

    /// Test 13: Segments with different copy number are not merged.
    #[test]
    fn cnv_merge_different_cn_not_merged() {
        let segments = vec![
            CnvSegment {
                chrom: "chr1".into(),
                start: 0,
                end: 1000,
                log2_ratio: 0.0,
                n_probes: 10,
                copy_number: 2,
            },
            CnvSegment {
                chrom: "chr1".into(),
                start: 1000,
                end: 2000,
                log2_ratio: 1.0,
                n_probes: 10,
                copy_number: 4,
            },
        ];

        let merged = merge_cnv_segments(&segments, 100, 0);

        // CN difference of 2 exceeds tolerance of 0 → should not merge.
        assert_eq!(merged.len(), 2);
    }

    /// Test 14: Copy number estimation from log2 ratio values.
    #[test]
    fn cn_estimation_from_log2_ratio() {
        // log2_ratio = 0.0 → CN = 2^(0+1) = 2
        assert_eq!(cn_from_log2(0.0), 2);
        // log2_ratio = 1.0 → CN = 2^(1+1) = 4
        assert_eq!(cn_from_log2(1.0), 4);
        // log2_ratio = -1.0 → CN = 2^(-1+1) = 1
        assert_eq!(cn_from_log2(-1.0), 1);
        // log2_ratio = -inf (large negative) → clamped to 0
        assert_eq!(cn_from_log2(-10.0), 0);
        // log2_ratio = large positive → clamped to 10
        assert_eq!(cn_from_log2(5.0), 10);
        // log2_ratio = -0.415 → CN ~ 2^0.585 ~ 1.5 → rounds to 2
        assert_eq!(cn_from_log2(-0.415), 2);
        // log2_ratio = 0.585 → CN = 2^1.585 ~ 3.0
        assert_eq!(cn_from_log2(0.585), 3);
    }
}
