//! Haplotype analysis — EM phasing, haplotype blocks, and haplotype diversity.
//!
//! Phase genotypes into haplotypes using expectation-maximization,
//! detect haplotype blocks from LD structure, and compute diversity.

use cyanea_core::{CyaneaError, Result};

/// A haplotype: a sequence of alleles (0 or 1 for biallelic markers).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Haplotype {
    /// Allele at each SNP position (0 = reference, 1 = alternate).
    pub alleles: Vec<u8>,
}

/// Result of EM-based haplotype phasing.
#[derive(Debug, Clone)]
pub struct PhasedGenotypes {
    /// Phased haplotype pairs for each sample.
    pub haplotypes: Vec<(Haplotype, Haplotype)>,
    /// Population haplotype frequencies.
    pub frequencies: Vec<(Haplotype, f64)>,
    /// Log-likelihood of the final solution.
    pub log_likelihood: f64,
}

/// A contiguous block of SNPs in high LD.
#[derive(Debug, Clone)]
pub struct HaplotypeBlock {
    /// Start SNP index (0-based).
    pub start: usize,
    /// End SNP index (inclusive).
    pub end: usize,
    /// Number of SNPs in the block.
    pub n_snps: usize,
    /// Haplotype diversity within the block.
    pub diversity: f64,
}

/// Phase genotypes using the Expectation-Maximization algorithm.
///
/// Input: `genotypes[sample][snp]` with values 0 (hom-ref), 1 (het), 2 (hom-alt).
///
/// # Errors
///
/// Returns an error if genotypes is empty, samples have inconsistent SNP counts,
/// or genotype values are not in {0, 1, 2}.
pub fn phase_em(genotypes: &[Vec<u8>], max_iter: usize) -> Result<PhasedGenotypes> {
    if genotypes.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "at least one sample is required".into(),
        ));
    }
    let n_snps = genotypes[0].len();
    if n_snps == 0 {
        return Err(CyaneaError::InvalidInput(
            "at least one SNP is required".into(),
        ));
    }
    for (i, g) in genotypes.iter().enumerate() {
        if g.len() != n_snps {
            return Err(CyaneaError::InvalidInput(format!(
                "sample {} has {} SNPs, expected {}",
                i,
                g.len(),
                n_snps
            )));
        }
        for &v in g {
            if v > 2 {
                return Err(CyaneaError::InvalidInput(format!(
                    "genotype values must be 0, 1, or 2; got {}",
                    v
                )));
            }
        }
    }

    let n_samples = genotypes.len();

    // Enumerate all possible haplotypes (feasible for small n_snps).
    let n_haplotypes = 1usize << n_snps;
    let all_haps: Vec<Haplotype> = (0..n_haplotypes)
        .map(|i| {
            let alleles = (0..n_snps)
                .map(|bit| ((i >> bit) & 1) as u8)
                .collect();
            Haplotype { alleles }
        })
        .collect();

    // Initialize haplotype frequencies uniformly.
    let mut freqs = vec![1.0 / n_haplotypes as f64; n_haplotypes];

    // For each sample, enumerate compatible haplotype pairs.
    let compatible: Vec<Vec<(usize, usize)>> = genotypes
        .iter()
        .map(|g| {
            let mut pairs = Vec::new();
            for (hi, h1) in all_haps.iter().enumerate() {
                for (hj, h2) in all_haps.iter().enumerate() {
                    if hi > hj {
                        continue;
                    }
                    let compatible = h1
                        .alleles
                        .iter()
                        .zip(h2.alleles.iter())
                        .zip(g.iter())
                        .all(|((&a1, &a2), &gt)| a1 + a2 == gt);
                    if compatible {
                        pairs.push((hi, hj));
                    }
                }
            }
            pairs
        })
        .collect();

    // EM iterations.
    for _ in 0..max_iter {
        let mut new_freqs = vec![0.0f64; n_haplotypes];

        for (_si, pairs) in compatible.iter().enumerate() {
            // E-step: weight each compatible pair by product of frequencies.
            let mut weights: Vec<f64> = pairs
                .iter()
                .map(|&(hi, hj)| {
                    let f = freqs[hi] * freqs[hj];
                    if hi != hj { 2.0 * f } else { f }
                })
                .collect();

            let sum_w: f64 = weights.iter().sum();
            if sum_w > 0.0 {
                for w in &mut weights {
                    *w /= sum_w;
                }
            }

            // M-step: accumulate weighted haplotype counts.
            for (pi, &(hi, hj)) in pairs.iter().enumerate() {
                new_freqs[hi] += weights[pi];
                new_freqs[hj] += weights[pi];
            }
        }

        // Normalize frequencies.
        let total: f64 = new_freqs.iter().sum();
        if total > 0.0 {
            for f in &mut new_freqs {
                *f /= total;
            }
        }

        // Check convergence.
        let max_delta: f64 = freqs
            .iter()
            .zip(new_freqs.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        freqs = new_freqs;
        if max_delta < 1e-8 {
            break;
        }
    }

    // Assign best haplotype pair for each sample.
    let mut phased = Vec::with_capacity(n_samples);
    for pairs in &compatible {
        let best = pairs
            .iter()
            .max_by(|&&(hi1, hj1), &&(hi2, hj2)| {
                let f1 = freqs[hi1] * freqs[hj1];
                let f2 = freqs[hi2] * freqs[hj2];
                f1.partial_cmp(&f2).unwrap()
            })
            .unwrap_or(&(0, 0));
        phased.push((all_haps[best.0].clone(), all_haps[best.1].clone()));
    }

    // Compute log-likelihood.
    let ll: f64 = compatible
        .iter()
        .map(|pairs| {
            let p: f64 = pairs
                .iter()
                .map(|&(hi, hj)| {
                    let f = freqs[hi] * freqs[hj];
                    if hi != hj { 2.0 * f } else { f }
                })
                .sum();
            if p > 0.0 { p.ln() } else { 0.0 }
        })
        .sum();

    // Collect non-zero frequency haplotypes.
    let frequencies: Vec<(Haplotype, f64)> = all_haps
        .into_iter()
        .zip(freqs.iter())
        .filter(|(_, &f)| f > 1e-10)
        .map(|(h, &f)| (h, f))
        .collect();

    Ok(PhasedGenotypes {
        haplotypes: phased,
        frequencies,
        log_likelihood: ll,
    })
}

/// Detect haplotype blocks from genotype data based on LD structure.
///
/// Uses a sliding window approach: a block boundary is placed where
/// the average pairwise LD (r²) between adjacent SNPs drops below `threshold`.
///
/// Input: `genotypes[sample][snp]` with values 0, 1, 2.
///
/// # Errors
///
/// Returns an error if genotypes is empty or threshold is not in [0, 1].
pub fn haplotype_blocks(
    genotypes: &[Vec<u8>],
    threshold: f64,
) -> Result<Vec<HaplotypeBlock>> {
    if genotypes.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "at least one sample is required".into(),
        ));
    }
    if !(0.0..=1.0).contains(&threshold) {
        return Err(CyaneaError::InvalidInput(format!(
            "threshold must be in [0, 1], got {}",
            threshold
        )));
    }

    let n_snps = genotypes[0].len();
    if n_snps <= 1 {
        if n_snps == 1 {
            return Ok(vec![HaplotypeBlock {
                start: 0,
                end: 0,
                n_snps: 1,
                diversity: 0.0,
            }]);
        }
        return Ok(Vec::new());
    }

    // Compute pairwise r² between adjacent SNPs.
    let mut r2_adjacent = Vec::with_capacity(n_snps - 1);
    for j in 0..n_snps - 1 {
        let r2 = compute_r2(genotypes, j, j + 1);
        r2_adjacent.push(r2);
    }

    // Find block boundaries where r² drops below threshold.
    let mut blocks = Vec::new();
    let mut block_start = 0;

    for j in 0..r2_adjacent.len() {
        if r2_adjacent[j] < threshold {
            // End current block.
            let block = make_block(genotypes, block_start, j);
            blocks.push(block);
            block_start = j + 1;
        }
    }
    // Final block.
    let block = make_block(genotypes, block_start, n_snps - 1);
    blocks.push(block);

    Ok(blocks)
}

/// Haplotype diversity: Hd = (n/(n-1)) * (1 - Σ p_i²).
///
/// Analogous to expected heterozygosity applied to haplotypes.
pub fn haplotype_diversity(haplotypes: &[&Haplotype]) -> f64 {
    let n = haplotypes.len();
    if n <= 1 {
        return 0.0;
    }

    // Count frequency of each distinct haplotype.
    let mut freq_map: Vec<(&Haplotype, usize)> = Vec::new();
    for &h in haplotypes {
        if let Some(entry) = freq_map.iter_mut().find(|(hap, _)| hap.alleles == h.alleles) {
            entry.1 += 1;
        } else {
            freq_map.push((h, 1));
        }
    }

    let n_f = n as f64;
    let sum_p2: f64 = freq_map
        .iter()
        .map(|(_, count)| {
            let p = *count as f64 / n_f;
            p * p
        })
        .sum();

    (n_f / (n_f - 1.0)) * (1.0 - sum_p2)
}

/// Compute r² between two SNP positions.
fn compute_r2(genotypes: &[Vec<u8>], snp_a: usize, snp_b: usize) -> f64 {
    let n = genotypes.len() as f64;
    if n == 0.0 {
        return 0.0;
    }

    let mut sum_a = 0.0;
    let mut sum_b = 0.0;
    let mut sum_ab = 0.0;
    let mut sum_a2 = 0.0;
    let mut sum_b2 = 0.0;
    let mut count = 0.0;

    for g in genotypes {
        let a = g[snp_a] as f64;
        let b = g[snp_b] as f64;
        sum_a += a;
        sum_b += b;
        sum_ab += a * b;
        sum_a2 += a * a;
        sum_b2 += b * b;
        count += 1.0;
    }

    if count == 0.0 {
        return 0.0;
    }

    let mean_a = sum_a / count;
    let mean_b = sum_b / count;
    let var_a = sum_a2 / count - mean_a * mean_a;
    let var_b = sum_b2 / count - mean_b * mean_b;
    let cov = sum_ab / count - mean_a * mean_b;

    if var_a <= 0.0 || var_b <= 0.0 {
        return 0.0;
    }

    let r = cov / (var_a * var_b).sqrt();
    r * r
}

fn make_block(genotypes: &[Vec<u8>], start: usize, end: usize) -> HaplotypeBlock {
    let n_snps = end - start + 1;

    // Extract haplotype-like signatures from genotypes within the block.
    let haps: Vec<Haplotype> = genotypes
        .iter()
        .map(|g| Haplotype {
            alleles: g[start..=end].to_vec(),
        })
        .collect();
    let hap_refs: Vec<&Haplotype> = haps.iter().collect();
    let diversity = haplotype_diversity(&hap_refs);

    HaplotypeBlock {
        start,
        end,
        n_snps,
        diversity,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phase_em_homozygous() {
        // All homozygous: trivial phasing.
        let genotypes = vec![vec![0, 0], vec![2, 2]];
        let result = phase_em(&genotypes, 100).unwrap();
        // Sample 0: both haplotypes should be [0, 0].
        assert_eq!(result.haplotypes[0].0.alleles, vec![0, 0]);
        assert_eq!(result.haplotypes[0].1.alleles, vec![0, 0]);
        // Sample 1: both haplotypes should be [1, 1].
        assert_eq!(result.haplotypes[1].0.alleles, vec![1, 1]);
        assert_eq!(result.haplotypes[1].1.alleles, vec![1, 1]);
    }

    #[test]
    fn phase_em_simple_het() {
        // Single heterozygous individual with 2 SNPs, both het.
        let genotypes = vec![vec![1, 1]];
        let result = phase_em(&genotypes, 100).unwrap();
        // Should produce two haplotypes that sum to [1, 1].
        let h1 = &result.haplotypes[0].0;
        let h2 = &result.haplotypes[0].1;
        for i in 0..2 {
            assert_eq!(h1.alleles[i] + h2.alleles[i], 1);
        }
    }

    #[test]
    fn haplotype_blocks_high_ld() {
        // All SNPs in perfect LD → single block.
        let genotypes = vec![
            vec![0, 0, 0],
            vec![1, 1, 1],
            vec![2, 2, 2],
            vec![0, 0, 0],
            vec![2, 2, 2],
        ];
        let blocks = haplotype_blocks(&genotypes, 0.5).unwrap();
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].start, 0);
        assert_eq!(blocks[0].end, 2);
        assert_eq!(blocks[0].n_snps, 3);
    }

    #[test]
    fn haplotype_diversity_single() {
        // All identical haplotypes → diversity = 0.
        let h = Haplotype { alleles: vec![0, 1] };
        let haps = vec![&h, &h, &h];
        assert!((haplotype_diversity(&haps) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn haplotype_diversity_maximum() {
        // All unique haplotypes → maximum diversity.
        let h1 = Haplotype { alleles: vec![0] };
        let h2 = Haplotype { alleles: vec![1] };
        let haps = vec![&h1, &h2];
        let hd = haplotype_diversity(&haps);
        // Hd = (2/1) * (1 - 2*(0.5^2)) = 2 * 0.5 = 1.0
        assert!((hd - 1.0).abs() < 1e-10);
    }
}
