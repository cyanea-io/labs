//! Methylation analysis — CpG sites, differentially methylated regions, and CpG islands.
//!
//! This module provides tools for bisulfite sequencing analysis:
//!
//! - [`CpgSite`] — a single CpG dinucleotide with methylation counts
//! - [`call_methylation`] — identify CpG sites from bisulfite-seq read counts
//! - [`find_dmrs`] — detect differentially methylated regions between sample groups
//! - [`find_cpg_islands`] — locate CpG islands in a reference sequence
//! - [`bisulfite_convert`] — simulate in-silico bisulfite conversion

use std::collections::HashMap;

use cyanea_core::{CyaneaError, Result};

use crate::genomic::Strand;

/// A single CpG dinucleotide with methylation read counts.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CpgSite {
    /// Chromosome name.
    pub chrom: String,
    /// 0-based position of the cytosine in the CpG dinucleotide.
    pub position: u64,
    /// Strand on which the CpG was observed.
    pub strand: Strand,
    /// Number of reads showing methylation (C) at this position.
    pub methylated_reads: u32,
    /// Total number of reads covering this position (methylated + unmethylated).
    pub total_reads: u32,
}

impl CpgSite {
    /// Returns the beta value (methylation fraction) for this CpG site.
    ///
    /// Beta = methylated_reads / total_reads. Returns 0.0 if total_reads is 0.
    pub fn beta(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            self.methylated_reads as f64 / self.total_reads as f64
        }
    }
}

/// A differentially methylated region (DMR) between two sample groups.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DmRegion {
    /// Chromosome name.
    pub chrom: String,
    /// Start position (0-based, inclusive).
    pub start: u64,
    /// End position (0-based, exclusive).
    pub end: u64,
    /// Mean delta-beta (group1 mean beta - group2 mean beta) across CpGs in the region.
    pub mean_delta_beta: f64,
    /// Number of CpG sites in this region.
    pub n_cpgs: usize,
    /// P-value from Welch's t-test on beta values.
    pub p_value: f64,
}

/// A CpG island — a genomic region enriched for CpG dinucleotides.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CpgIsland {
    /// Chromosome name.
    pub chrom: String,
    /// Start position (0-based, inclusive).
    pub start: u64,
    /// End position (0-based, exclusive).
    pub end: u64,
    /// Number of CpG dinucleotides in the island.
    pub cpg_count: usize,
    /// Observed/expected CpG ratio.
    pub obs_exp_ratio: f64,
    /// GC content as a fraction (0.0–1.0).
    pub gc_content: f64,
}

/// Configuration for differentially methylated region (DMR) detection.
pub struct DmrConfig {
    /// Minimum absolute delta-beta to consider a site significant.
    pub min_delta_beta: f64,
    /// Maximum gap (in bp) between adjacent significant CpGs to merge into one DMR.
    pub max_gap: u64,
    /// Minimum number of CpG sites required in a DMR.
    pub min_cpgs: usize,
    /// Minimum total read coverage required at a CpG in both groups.
    pub min_coverage: u32,
}

impl Default for DmrConfig {
    fn default() -> Self {
        Self {
            min_delta_beta: 0.2,
            max_gap: 500,
            min_cpgs: 3,
            min_coverage: 5,
        }
    }
}

/// Call methylation at CpG sites from bisulfite sequencing counts.
///
/// For each position in `positions`, checks whether `reference[position]` and
/// `reference[position + 1]` form a CpG dinucleotide (C followed by G on the
/// forward strand). If so, creates a [`CpgSite`] with the corresponding counts.
/// Positions where the reference does not contain a CpG are skipped.
///
/// # Arguments
///
/// * `positions` — 0-based positions in the reference to check
/// * `c_counts` — number of methylated (C) reads at each position
/// * `t_counts` — number of unmethylated (T) reads at each position
/// * `chrom` — chromosome name
/// * `reference` — the reference sequence as bytes (e.g., `b"ACGTACGT"`)
pub fn call_methylation(
    positions: &[u64],
    c_counts: &[u32],
    t_counts: &[u32],
    chrom: &str,
    reference: &[u8],
) -> Vec<CpgSite> {
    let n = positions.len().min(c_counts.len()).min(t_counts.len());
    let ref_len = reference.len();
    let mut sites = Vec::new();

    for i in 0..n {
        let pos = positions[i] as usize;
        if pos + 1 >= ref_len {
            continue;
        }

        let base = reference[pos].to_ascii_uppercase();
        let next = reference[pos + 1].to_ascii_uppercase();

        if base == b'C' && next == b'G' {
            // Forward-strand CpG
            sites.push(CpgSite {
                chrom: chrom.to_string(),
                position: positions[i],
                strand: Strand::Forward,
                methylated_reads: c_counts[i],
                total_reads: c_counts[i] + t_counts[i],
            });
        } else if base == b'G' && pos > 0 && reference[pos - 1].to_ascii_uppercase() == b'C' {
            // Reverse-strand CpG (G preceded by C)
            // Skip: this would be captured when the C position itself is queried.
            continue;
        }
    }

    sites
}

/// Detect differentially methylated regions (DMRs) between two groups of samples.
///
/// For each CpG position present in both groups, computes mean beta values and
/// their difference (delta-beta). Sites with |delta_beta| >= `config.min_delta_beta`
/// and sufficient coverage are merged into regions when they lie within
/// `config.max_gap` base pairs of each other. Regions with fewer than
/// `config.min_cpgs` sites are discarded.
///
/// P-values are computed using Welch's t-test on beta values.
pub fn find_dmrs(
    group1: &[&[CpgSite]],
    group2: &[&[CpgSite]],
    config: &DmrConfig,
) -> Result<Vec<DmRegion>> {
    if group1.is_empty() || group2.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "Both groups must have at least one sample".to_string(),
        ));
    }

    // Collect beta values per (chrom, position) for each group.
    let collect_group = |group: &[&[CpgSite]]| -> HashMap<(String, u64), Vec<(f64, u32)>> {
        let mut map: HashMap<(String, u64), Vec<(f64, u32)>> = HashMap::new();
        for sample in group {
            for site in *sample {
                map.entry((site.chrom.clone(), site.position))
                    .or_default()
                    .push((site.beta(), site.total_reads));
            }
        }
        map
    };

    let g1_map = collect_group(group1);
    let g2_map = collect_group(group2);

    // Find positions present in both groups with sufficient coverage and delta-beta.
    struct SigSite {
        chrom: String,
        position: u64,
        g1_betas: Vec<f64>,
        g2_betas: Vec<f64>,
        delta_beta: f64,
    }

    let mut sig_sites: Vec<SigSite> = Vec::new();

    for (key, g1_entries) in &g1_map {
        if let Some(g2_entries) = g2_map.get(key) {
            // Check min coverage in both groups.
            let g1_pass = g1_entries.iter().all(|(_, cov)| *cov >= config.min_coverage);
            let g2_pass = g2_entries.iter().all(|(_, cov)| *cov >= config.min_coverage);
            if !g1_pass || !g2_pass {
                continue;
            }

            let g1_betas: Vec<f64> = g1_entries.iter().map(|(b, _)| *b).collect();
            let g2_betas: Vec<f64> = g2_entries.iter().map(|(b, _)| *b).collect();

            let mean1 = g1_betas.iter().sum::<f64>() / g1_betas.len() as f64;
            let mean2 = g2_betas.iter().sum::<f64>() / g2_betas.len() as f64;
            let delta = mean1 - mean2;

            if delta.abs() >= config.min_delta_beta {
                sig_sites.push(SigSite {
                    chrom: key.0.clone(),
                    position: key.1,
                    g1_betas,
                    g2_betas,
                    delta_beta: delta,
                });
            }
        }
    }

    // Sort by (chrom, position).
    sig_sites.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.position.cmp(&b.position)));

    // Merge adjacent sites within max_gap into regions.
    let mut regions: Vec<DmRegion> = Vec::new();
    let mut i = 0;

    while i < sig_sites.len() {
        let chrom = sig_sites[i].chrom.clone();
        let start = sig_sites[i].position;
        let mut end = sig_sites[i].position + 1;
        let mut all_g1_betas: Vec<f64> = sig_sites[i].g1_betas.clone();
        let mut all_g2_betas: Vec<f64> = sig_sites[i].g2_betas.clone();
        let mut delta_sum = sig_sites[i].delta_beta;
        let mut n_cpgs = 1;

        let mut j = i + 1;
        while j < sig_sites.len()
            && sig_sites[j].chrom == chrom
            && sig_sites[j].position <= end + config.max_gap - 1
        {
            end = sig_sites[j].position + 1;
            all_g1_betas.extend_from_slice(&sig_sites[j].g1_betas);
            all_g2_betas.extend_from_slice(&sig_sites[j].g2_betas);
            delta_sum += sig_sites[j].delta_beta;
            n_cpgs += 1;
            j += 1;
        }

        if n_cpgs >= config.min_cpgs {
            let mean_delta = delta_sum / n_cpgs as f64;
            let p_value = welch_t_test(&all_g1_betas, &all_g2_betas);

            regions.push(DmRegion {
                chrom,
                start,
                end,
                mean_delta_beta: mean_delta,
                n_cpgs,
                p_value,
            });
        }

        i = j;
    }

    Ok(regions)
}

/// Approximate the error function using the Abramowitz & Stegun formula (7.1.26).
fn erf_approx(x: f64) -> f64 {
    let sign = if x >= 0.0 { 1.0 } else { -1.0 };
    let x = x.abs();

    let p = 0.3275911;
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;

    let t = 1.0 / (1.0 + p * x);
    let poly = t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))));
    sign * (1.0 - poly * (-x * x).exp())
}

/// Approximate erfc(x) = 1 - erf(x).
fn erfc_approx(x: f64) -> f64 {
    1.0 - erf_approx(x)
}

/// Compute a two-sided p-value from Welch's t-test using normal approximation.
///
/// t = (mean1 - mean2) / sqrt(var1/n1 + var2/n2)
/// p = erfc(|t| / sqrt(2))
fn welch_t_test(group1: &[f64], group2: &[f64]) -> f64 {
    let n1 = group1.len() as f64;
    let n2 = group2.len() as f64;

    if n1 < 2.0 || n2 < 2.0 {
        return 1.0;
    }

    let mean1 = group1.iter().sum::<f64>() / n1;
    let mean2 = group2.iter().sum::<f64>() / n2;

    let var1 = group1.iter().map(|x| (x - mean1).powi(2)).sum::<f64>() / (n1 - 1.0);
    let var2 = group2.iter().map(|x| (x - mean2).powi(2)).sum::<f64>() / (n2 - 1.0);

    let se = (var1 / n1 + var2 / n2).sqrt();
    if se == 0.0 {
        return 1.0;
    }

    let t = (mean1 - mean2) / se;
    erfc_approx(t.abs() / std::f64::consts::SQRT_2)
}

/// Find CpG islands in a DNA sequence using a sliding window approach.
///
/// Scans the sequence with a 200 bp window and identifies regions where:
/// - GC content >= `min_gc`
/// - Observed/expected CpG ratio >= `min_obs_exp`
/// - Region length >= `min_length`
///
/// The observed/expected ratio is `(n_cpg * window_length) / (n_c * n_g)`.
/// Overlapping qualifying windows are merged into contiguous islands.
pub fn find_cpg_islands(
    sequence: &[u8],
    chrom: &str,
    min_length: u64,
    min_gc: f64,
    min_obs_exp: f64,
) -> Vec<CpgIsland> {
    let seq_len = sequence.len();
    let window_size: usize = 200;

    if seq_len < window_size {
        return Vec::new();
    }

    // Identify qualifying windows.
    let mut qualifying: Vec<(usize, usize)> = Vec::new(); // (start, end)

    for start in 0..=(seq_len - window_size) {
        let window = &sequence[start..start + window_size];
        let wlen = window_size as f64;

        let mut n_c: usize = 0;
        let mut n_g: usize = 0;
        let mut n_cpg: usize = 0;

        for j in 0..window.len() {
            let b = window[j].to_ascii_uppercase();
            if b == b'C' {
                n_c += 1;
                if j + 1 < window.len() && window[j + 1].to_ascii_uppercase() == b'G' {
                    n_cpg += 1;
                }
            } else if b == b'G' {
                n_g += 1;
            }
        }

        let gc_content = (n_c + n_g) as f64 / wlen;
        if gc_content < min_gc {
            continue;
        }

        let obs_exp = if n_c > 0 && n_g > 0 {
            (n_cpg as f64 * wlen) / (n_c as f64 * n_g as f64)
        } else {
            0.0
        };

        if obs_exp < min_obs_exp {
            continue;
        }

        qualifying.push((start, start + window_size));
    }

    if qualifying.is_empty() {
        return Vec::new();
    }

    // Merge overlapping windows.
    let mut merged: Vec<(usize, usize)> = Vec::new();
    let mut cur_start = qualifying[0].0;
    let mut cur_end = qualifying[0].1;

    for &(s, e) in &qualifying[1..] {
        if s <= cur_end {
            cur_end = cur_end.max(e);
        } else {
            merged.push((cur_start, cur_end));
            cur_start = s;
            cur_end = e;
        }
    }
    merged.push((cur_start, cur_end));

    // Build CpG islands from merged windows.
    let mut islands = Vec::new();
    for (start, end) in merged {
        let length = (end - start) as u64;
        if length < min_length {
            continue;
        }

        let region = &sequence[start..end];
        let mut n_c: usize = 0;
        let mut n_g: usize = 0;
        let mut n_cpg: usize = 0;

        for j in 0..region.len() {
            let b = region[j].to_ascii_uppercase();
            if b == b'C' {
                n_c += 1;
                if j + 1 < region.len() && region[j + 1].to_ascii_uppercase() == b'G' {
                    n_cpg += 1;
                }
            } else if b == b'G' {
                n_g += 1;
            }
        }

        let gc_content = (n_c + n_g) as f64 / region.len() as f64;
        let obs_exp = if n_c > 0 && n_g > 0 {
            (n_cpg as f64 * region.len() as f64) / (n_c as f64 * n_g as f64)
        } else {
            0.0
        };

        islands.push(CpgIsland {
            chrom: chrom.to_string(),
            start: start as u64,
            end: end as u64,
            cpg_count: n_cpg,
            obs_exp_ratio: obs_exp,
            gc_content,
        });
    }

    islands
}

/// Simulate bisulfite conversion of a DNA sequence.
///
/// In bisulfite sequencing, unmethylated cytosines are converted to uracil
/// (read as thymine), while methylated cytosines are protected. This function
/// converts all C bases to T, except at positions listed in `methylated_positions`,
/// which remain as C.
///
/// Other bases (G, A, T) are unaffected.
pub fn bisulfite_convert(sequence: &[u8], methylated_positions: &[u64]) -> Vec<u8> {
    let meth_set: std::collections::HashSet<u64> =
        methylated_positions.iter().copied().collect();

    sequence
        .iter()
        .enumerate()
        .map(|(i, &base)| {
            let upper = base.to_ascii_uppercase();
            if upper == b'C' && !meth_set.contains(&(i as u64)) {
                if base.is_ascii_uppercase() {
                    b'T'
                } else {
                    b't'
                }
            } else {
                base
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_beta_value_calculation() {
        let site = CpgSite {
            chrom: "chr1".to_string(),
            position: 100,
            strand: Strand::Forward,
            methylated_reads: 8,
            total_reads: 10,
        };
        assert!((site.beta() - 0.8).abs() < 1e-10);
    }

    #[test]
    fn test_beta_zero_coverage() {
        let site = CpgSite {
            chrom: "chr1".to_string(),
            position: 100,
            strand: Strand::Forward,
            methylated_reads: 0,
            total_reads: 0,
        };
        assert_eq!(site.beta(), 0.0);
    }

    #[test]
    fn test_cpg_identification_from_reference() {
        // Reference has CpG at position 2 (CG dinucleotide).
        let reference = b"AACGTTCGAA";
        let positions = vec![2, 6];
        let c_counts = vec![5, 3];
        let t_counts = vec![5, 7];

        let sites = call_methylation(&positions, &c_counts, &t_counts, "chr1", reference);
        assert_eq!(sites.len(), 2);

        assert_eq!(sites[0].position, 2);
        assert_eq!(sites[0].strand, Strand::Forward);
        assert_eq!(sites[0].methylated_reads, 5);
        assert_eq!(sites[0].total_reads, 10);

        assert_eq!(sites[1].position, 6);
        assert_eq!(sites[1].methylated_reads, 3);
        assert_eq!(sites[1].total_reads, 10);
    }

    #[test]
    fn test_non_cpg_position_skipped() {
        // Position 1 is 'A' followed by 'C', not CpG.
        let reference = b"AACGTT";
        let positions = vec![1]; // 'A' at pos 1
        let c_counts = vec![5];
        let t_counts = vec![5];

        let sites = call_methylation(&positions, &c_counts, &t_counts, "chr1", reference);
        assert!(sites.is_empty());
    }

    #[test]
    fn test_bisulfite_unmethylated_c_to_t() {
        let seq = b"ACGT";
        let methylated: &[u64] = &[];
        let converted = bisulfite_convert(seq, methylated);
        // C at position 1 should become T.
        assert_eq!(converted, b"ATGT");
    }

    #[test]
    fn test_bisulfite_methylated_c_stays_c() {
        let seq = b"ACGT";
        let methylated = &[1]; // position 1 is methylated
        let converted = bisulfite_convert(seq, methylated);
        // Methylated C stays as C.
        assert_eq!(converted, b"ACGT");
    }

    #[test]
    fn test_bisulfite_g_a_t_unaffected() {
        let seq = b"GATT";
        let methylated: &[u64] = &[];
        let converted = bisulfite_convert(seq, methylated);
        // No C in the sequence, everything stays the same.
        assert_eq!(converted, b"GATT");
    }

    #[test]
    fn test_cpg_island_found_on_cg_rich_sequence() {
        // Build a 250 bp CG-rich sequence: repeating "CG" pattern.
        let mut seq = Vec::new();
        for _ in 0..125 {
            seq.push(b'C');
            seq.push(b'G');
        }
        assert_eq!(seq.len(), 250);

        let islands = find_cpg_islands(&seq, "chr1", 200, 0.5, 0.6);
        assert!(!islands.is_empty());
        assert!(islands[0].gc_content >= 0.5);
        assert!(islands[0].obs_exp_ratio >= 0.6);
        assert!(islands[0].cpg_count > 0);
    }

    #[test]
    fn test_no_cpg_island_on_at_rich_sequence() {
        // 300 bp AT-rich sequence.
        let mut seq = Vec::new();
        for _ in 0..150 {
            seq.push(b'A');
            seq.push(b'T');
        }
        assert_eq!(seq.len(), 300);

        let islands = find_cpg_islands(&seq, "chr1", 200, 0.5, 0.6);
        assert!(islands.is_empty());
    }

    #[test]
    fn test_dmr_detected_between_differentially_methylated_groups() {
        // Group 1: high methylation at 5 CpG sites.
        let g1_sample: Vec<CpgSite> = (0..5)
            .map(|i| CpgSite {
                chrom: "chr1".to_string(),
                position: i * 100,
                strand: Strand::Forward,
                methylated_reads: 9,
                total_reads: 10,
            })
            .collect();

        // Group 2: low methylation at the same 5 sites.
        let g2_sample: Vec<CpgSite> = (0..5)
            .map(|i| CpgSite {
                chrom: "chr1".to_string(),
                position: i * 100,
                strand: Strand::Forward,
                methylated_reads: 1,
                total_reads: 10,
            })
            .collect();

        let g1_refs: Vec<&[CpgSite]> = vec![&g1_sample, &g1_sample];
        let g2_refs: Vec<&[CpgSite]> = vec![&g2_sample, &g2_sample];

        let config = DmrConfig {
            min_delta_beta: 0.2,
            max_gap: 500,
            min_cpgs: 3,
            min_coverage: 5,
        };

        let dmrs = find_dmrs(&g1_refs, &g2_refs, &config).unwrap();
        assert!(!dmrs.is_empty());
        assert!(dmrs[0].mean_delta_beta > 0.0);
        assert!(dmrs[0].n_cpgs >= 3);
        assert!(dmrs[0].p_value < 0.05);
    }

    #[test]
    fn test_dmr_none_when_groups_identical() {
        let sample: Vec<CpgSite> = (0..5)
            .map(|i| CpgSite {
                chrom: "chr1".to_string(),
                position: i * 100,
                strand: Strand::Forward,
                methylated_reads: 5,
                total_reads: 10,
            })
            .collect();

        let g1_refs: Vec<&[CpgSite]> = vec![&sample];
        let g2_refs: Vec<&[CpgSite]> = vec![&sample];

        let config = DmrConfig::default();
        let dmrs = find_dmrs(&g1_refs, &g2_refs, &config).unwrap();
        assert!(dmrs.is_empty());
    }

    #[test]
    fn test_dmr_min_cpgs_filter() {
        // Only 2 CpG sites — should not pass the min_cpgs=3 filter.
        let g1_sample: Vec<CpgSite> = (0..2)
            .map(|i| CpgSite {
                chrom: "chr1".to_string(),
                position: i * 100,
                strand: Strand::Forward,
                methylated_reads: 9,
                total_reads: 10,
            })
            .collect();

        let g2_sample: Vec<CpgSite> = (0..2)
            .map(|i| CpgSite {
                chrom: "chr1".to_string(),
                position: i * 100,
                strand: Strand::Forward,
                methylated_reads: 1,
                total_reads: 10,
            })
            .collect();

        let g1_refs: Vec<&[CpgSite]> = vec![&g1_sample, &g1_sample];
        let g2_refs: Vec<&[CpgSite]> = vec![&g2_sample, &g2_sample];

        let config = DmrConfig {
            min_delta_beta: 0.2,
            max_gap: 500,
            min_cpgs: 3,
            min_coverage: 5,
        };

        let dmrs = find_dmrs(&g1_refs, &g2_refs, &config).unwrap();
        assert!(dmrs.is_empty());
    }
}
