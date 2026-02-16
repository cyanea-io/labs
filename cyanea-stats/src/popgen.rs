//! Population genetics — allele frequencies, HWE, Fst, diversity, Tajima's D, LD, genotype PCA.
//!
//! All functions operate on 012-encoded genotype data (PLINK convention):
//! - `0` = homozygous reference
//! - `1` = heterozygous
//! - `2` = homozygous alternate
//! - `None` (or `255` for u8 variants) = missing
//!
//! No VCF dependency — genotypes are plain slices.

// Genetics nomenclature uses mixed-case allele names (A/a, B/b) — suppress warnings.
#![allow(non_snake_case)]

use cyanea_core::{CyaneaError, Result};

use crate::distribution::{ChiSquared, Distribution};

// ── Result types ─────────────────────────────────────────────────────────

/// Allele frequency summary for a single site.
#[derive(Debug, Clone)]
pub struct AlleleFrequencies {
    /// Reference allele frequency.
    pub freq_ref: f64,
    /// Alternate allele frequency.
    pub freq_alt: f64,
    /// Total non-missing allele count (2 × non-missing individuals).
    pub allele_count: usize,
    /// Number of missing genotypes.
    pub missing_count: usize,
    /// Observed heterozygosity (proportion of hets among non-missing).
    pub observed_het: f64,
    /// Expected heterozygosity under HWE (2pq).
    pub expected_het: f64,
}

/// Hardy-Weinberg equilibrium test result.
#[derive(Debug, Clone)]
pub struct HweResult {
    /// Chi-squared statistic (1 df).
    pub chi_squared: f64,
    /// P-value from chi-squared distribution.
    pub p_value: f64,
    /// Observed genotype counts \[hom-ref, het, hom-alt\].
    pub observed: [usize; 3],
    /// Expected genotype counts under HWE.
    pub expected: [f64; 3],
    /// Inbreeding coefficient F = 1 - Ho/He.
    pub inbreeding_coeff: f64,
}

/// Method used for Fst estimation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FstMethod {
    /// Weir & Cockerham (1984) ANOVA-based estimator.
    WeirCockerham,
    /// Hudson (1992) estimator.
    Hudson,
}

/// Fixation index (Fst) result.
#[derive(Debug, Clone)]
pub struct FstResult {
    /// Estimated Fst value.
    pub fst: f64,
    /// Method used.
    pub method: FstMethod,
}

/// Nucleotide diversity statistics.
#[derive(Debug, Clone)]
pub struct DiversityStats {
    /// Nucleotide diversity (pi) — average pairwise differences per site.
    pub pi: f64,
    /// Watterson's theta (theta_w).
    pub theta_w: f64,
    /// Number of segregating sites.
    pub segregating_sites: usize,
    /// Number of sequences.
    pub n_sequences: usize,
}

/// Tajima's D statistic.
#[derive(Debug, Clone)]
pub struct TajimaD {
    /// Tajima's D value.
    pub d: f64,
    /// Nucleotide diversity (pi).
    pub pi: f64,
    /// Watterson's theta.
    pub theta_w: f64,
    /// Number of segregating sites.
    pub segregating_sites: usize,
    /// Number of sequences.
    pub n_sequences: usize,
}

/// Linkage disequilibrium result between two sites.
#[derive(Debug, Clone)]
pub struct LdResult {
    /// Squared correlation coefficient.
    pub r_squared: f64,
    /// Normalized LD coefficient (|D| / D_max).
    pub d_prime: f64,
    /// Raw LD coefficient D.
    pub d: f64,
    /// Number of haplotypes used (2 × non-missing diploid individuals).
    pub n_haplotypes: usize,
}

// ── Internal helpers ─────────────────────────────────────────────────────

/// Count genotypes: returns (n_0, n_1, n_2, n_missing).
fn count_genotypes(genotypes: &[Option<u8>]) -> Result<(usize, usize, usize, usize)> {
    let mut n0 = 0usize;
    let mut n1 = 0usize;
    let mut n2 = 0usize;
    let mut n_miss = 0usize;
    for (i, g) in genotypes.iter().enumerate() {
        match g {
            Some(0) => n0 += 1,
            Some(1) => n1 += 1,
            Some(2) => n2 += 1,
            None => n_miss += 1,
            Some(v) => {
                return Err(CyaneaError::InvalidInput(format!(
                    "invalid genotype {} at index {} (expected 0, 1, or 2)",
                    v, i
                )));
            }
        }
    }
    Ok((n0, n1, n2, n_miss))
}

/// Allele frequencies (p=ref, q=alt) from genotype counts.
fn allele_freq_from_counts(n0: usize, n1: usize, n2: usize) -> (f64, f64) {
    let total = 2 * (n0 + n1 + n2);
    if total == 0 {
        return (0.0, 0.0);
    }
    let p = (2 * n0 + n1) as f64 / total as f64;
    let q = (2 * n2 + n1) as f64 / total as f64;
    (p, q)
}

/// Harmonic number H_n = sum_{i=1}^{n} 1/i.
fn harmonic(n: usize) -> f64 {
    (1..=n).map(|i| 1.0 / i as f64).sum()
}

/// Sum of squared reciprocals: sum_{i=1}^{n} 1/i^2.
fn harmonic_sq(n: usize) -> f64 {
    (1..=n).map(|i| 1.0 / (i as f64 * i as f64)).sum()
}

/// EM algorithm for haplotype frequency estimation from unphased diploid data.
///
/// Returns (p_AB, p_Ab, p_aB, p_ab) haplotype frequencies.
fn em_haplotype_freqs(
    site_a: &[Option<u8>],
    site_b: &[Option<u8>],
    max_iter: usize,
    tol: f64,
) -> Result<(f64, f64, f64, f64)> {
    // Count joint genotypes (only non-missing at both sites)
    let mut counts = [[0usize; 3]; 3]; // counts[ga][gb]
    let mut n = 0usize;
    for (ga, gb) in site_a.iter().zip(site_b.iter()) {
        if let (Some(a), Some(b)) = (ga, gb) {
            if *a > 2 || *b > 2 {
                return Err(CyaneaError::InvalidInput(
                    "invalid genotype value (expected 0, 1, or 2)".into(),
                ));
            }
            counts[*a as usize][*b as usize] += 1;
            n += 1;
        }
    }

    if n == 0 {
        return Err(CyaneaError::InvalidInput(
            "no non-missing genotype pairs".into(),
        ));
    }

    // Compute marginal allele frequencies for initialization
    let mut a_alt = 0usize;
    let mut b_alt = 0usize;
    for ga in 0..3 {
        for gb in 0..3 {
            a_alt += ga * counts[ga][gb];
            b_alt += gb * counts[ga][gb];
        }
    }
    let pa = a_alt as f64 / (2 * n) as f64; // freq of 'a' (alt at site A)
    let pb = b_alt as f64 / (2 * n) as f64; // freq of 'b' (alt at site B)

    // Initialize haplotype freqs assuming linkage equilibrium
    let mut p_ab = pa * pb; // ab = both alt
    let mut p_aB = pa * (1.0 - pb);
    let mut p_Ab = (1.0 - pa) * pb;
    let mut p_AB = (1.0 - pa) * (1.0 - pb);

    for _ in 0..max_iter {
        // E-step: estimate haplotype counts
        let mut h_AB = 0.0;
        let mut h_Ab = 0.0;
        let mut h_aB = 0.0;
        let mut h_ab = 0.0;

        for ga in 0..3usize {
            for gb in 0..3usize {
                let c = counts[ga][gb] as f64;
                if c == 0.0 {
                    continue;
                }
                // For each diploid genotype pair, distribute haplotypes
                match (ga, gb) {
                    (0, 0) => h_AB += 2.0 * c, // AA,BB → 2×AB
                    (0, 1) => {
                        h_AB += c;
                        h_Ab += c;
                    } // AA,Bb → AB + Ab
                    (0, 2) => h_Ab += 2.0 * c, // AA,bb → 2×Ab
                    (1, 0) => {
                        h_AB += c;
                        h_aB += c;
                    } // Aa,BB → AB + aB
                    (1, 1) => {
                        // Aa,Bb → ambiguous: AB+ab or Ab+aB
                        let denom = p_AB * p_ab + p_Ab * p_aB;
                        if denom > 0.0 {
                            let w = p_AB * p_ab / denom;
                            h_AB += c * w;
                            h_ab += c * w;
                            h_Ab += c * (1.0 - w);
                            h_aB += c * (1.0 - w);
                        } else {
                            // Ambiguous — split evenly
                            h_AB += 0.5 * c;
                            h_ab += 0.5 * c;
                            h_Ab += 0.5 * c;
                            h_aB += 0.5 * c;
                        }
                    }
                    (1, 2) => {
                        h_Ab += c;
                        h_ab += c;
                    } // Aa,bb → Ab + ab
                    (2, 0) => h_aB += 2.0 * c, // aa,BB → 2×aB
                    (2, 1) => {
                        h_aB += c;
                        h_ab += c;
                    } // aa,Bb → aB + ab
                    (2, 2) => h_ab += 2.0 * c, // aa,bb → 2×ab
                    _ => {}
                }
            }
        }

        // M-step: update frequencies
        let total = h_AB + h_Ab + h_aB + h_ab;
        let new_AB = h_AB / total;
        let new_Ab = h_Ab / total;
        let new_aB = h_aB / total;
        let new_ab = h_ab / total;

        let diff = (new_AB - p_AB).abs()
            + (new_Ab - p_Ab).abs()
            + (new_aB - p_aB).abs()
            + (new_ab - p_ab).abs();

        p_AB = new_AB;
        p_Ab = new_Ab;
        p_aB = new_aB;
        p_ab = new_ab;

        if diff < tol {
            break;
        }
    }

    Ok((p_AB, p_Ab, p_aB, p_ab))
}

// ── Public API ───────────────────────────────────────────────────────────

/// Compute allele frequencies and heterozygosity from 012-encoded genotypes.
///
/// # Errors
///
/// Returns an error if genotypes is empty, all genotypes are missing,
/// or any value is outside {0, 1, 2, None}.
pub fn allele_frequencies(genotypes: &[Option<u8>]) -> Result<AlleleFrequencies> {
    if genotypes.is_empty() {
        return Err(CyaneaError::InvalidInput("empty genotype array".into()));
    }
    let (n0, n1, n2, n_miss) = count_genotypes(genotypes)?;
    let n_called = n0 + n1 + n2;
    if n_called == 0 {
        return Err(CyaneaError::InvalidInput(
            "all genotypes are missing".into(),
        ));
    }
    let (p, q) = allele_freq_from_counts(n0, n1, n2);
    let obs_het = n1 as f64 / n_called as f64;
    let exp_het = 2.0 * p * q;

    Ok(AlleleFrequencies {
        freq_ref: p,
        freq_alt: q,
        allele_count: 2 * n_called,
        missing_count: n_miss,
        observed_het: obs_het,
        expected_het: exp_het,
    })
}

/// Convenience wrapper: compute allele frequencies from u8 genotypes
/// where 255 encodes missing data.
pub fn allele_frequencies_u8(genotypes: &[u8]) -> Result<AlleleFrequencies> {
    let opts: Vec<Option<u8>> = genotypes
        .iter()
        .map(|&g| if g == 255 { None } else { Some(g) })
        .collect();
    allele_frequencies(&opts)
}

/// Hardy-Weinberg equilibrium chi-squared test (1 df).
///
/// Tests whether observed genotype frequencies deviate from HWE expectations.
///
/// # Errors
///
/// Returns an error if the site is monomorphic (freq_alt == 0 or freq_ref == 0),
/// genotypes are empty, or all are missing.
pub fn hwe_test(genotypes: &[Option<u8>]) -> Result<HweResult> {
    if genotypes.is_empty() {
        return Err(CyaneaError::InvalidInput("empty genotype array".into()));
    }
    let (n0, n1, n2, _) = count_genotypes(genotypes)?;
    let n = n0 + n1 + n2;
    if n == 0 {
        return Err(CyaneaError::InvalidInput(
            "all genotypes are missing".into(),
        ));
    }
    let (p, q) = allele_freq_from_counts(n0, n1, n2);
    if p == 0.0 || q == 0.0 {
        return Err(CyaneaError::InvalidInput(
            "monomorphic site: HWE test requires polymorphism".into(),
        ));
    }

    let n_f = n as f64;
    let exp = [p * p * n_f, 2.0 * p * q * n_f, q * q * n_f];
    let obs = [n0, n1, n2];

    let chi_sq: f64 = obs
        .iter()
        .zip(exp.iter())
        .map(|(&o, &e)| {
            if e > 0.0 {
                (o as f64 - e).powi(2) / e
            } else {
                0.0
            }
        })
        .sum();

    let chi2_dist = ChiSquared::new(1.0)?;
    let p_value = 1.0 - chi2_dist.cdf(chi_sq);

    let obs_het = n1 as f64 / n_f;
    let exp_het = 2.0 * p * q;
    let f = if exp_het > 0.0 {
        1.0 - obs_het / exp_het
    } else {
        0.0
    };

    Ok(HweResult {
        chi_squared: chi_sq,
        p_value,
        observed: obs,
        expected: exp,
        inbreeding_coeff: f,
    })
}

/// Weir & Cockerham (1984) Fst estimator across multiple loci.
///
/// Each element of `pop1` and `pop2` is the genotype vector for one locus.
/// Uses ratio-of-averages: Fst = sum(a) / sum(a + b + c).
///
/// # Errors
///
/// Returns an error if no loci are provided or populations are empty.
pub fn fst_weir_cockerham(
    pop1: &[&[Option<u8>]],
    pop2: &[&[Option<u8>]],
) -> Result<FstResult> {
    if pop1.is_empty() || pop2.is_empty() {
        return Err(CyaneaError::InvalidInput("empty loci or populations".into()));
    }

    let mut sum_a = 0.0;
    let mut sum_abc = 0.0;

    for (locus1, locus2) in pop1.iter().zip(pop2.iter()) {
        let (n0_1, n1_1, n2_1, _) = count_genotypes(locus1)?;
        let (n0_2, n1_2, n2_2, _) = count_genotypes(locus2)?;
        let ni1 = (n0_1 + n1_1 + n2_1) as f64;
        let ni2 = (n0_2 + n1_2 + n2_2) as f64;
        if ni1 == 0.0 || ni2 == 0.0 {
            continue;
        }

        let r = 2.0; // number of populations
        let n_bar = (ni1 + ni2) / r;
        let n_total = ni1 + ni2;

        let (_, q1) = allele_freq_from_counts(n0_1, n1_1, n2_1);
        let (_, q2) = allele_freq_from_counts(n0_2, n1_2, n2_2);

        let p_bar = (ni1 * q1 + ni2 * q2) / n_total;

        let s_sq = (ni1 * (q1 - p_bar).powi(2) + ni2 * (q2 - p_bar).powi(2))
            / ((r - 1.0) * n_bar);

        let h_bar = (ni1 * n1_1 as f64 / ni1 + ni2 * n1_2 as f64 / ni2) / n_total;

        let nc = n_total - (ni1 * ni1 + ni2 * ni2) / n_total;
        let nc = nc / (r - 1.0);

        // Weir-Cockerham variance components
        let a = (n_bar / nc)
            * (s_sq - (1.0 / (n_bar - 1.0)) * (p_bar * (1.0 - p_bar) - (r - 1.0) / r * s_sq - 0.25 * h_bar));

        let b = (n_bar / (n_bar - 1.0))
            * (p_bar * (1.0 - p_bar) - (r - 1.0) / r * s_sq
                - (2.0 * n_bar - 1.0) / (4.0 * n_bar) * h_bar);

        let c = 0.5 * h_bar;

        sum_a += a;
        sum_abc += a + b + c;
    }

    let fst = if sum_abc > 0.0 {
        (sum_a / sum_abc).clamp(0.0, 1.0)
    } else {
        0.0
    };

    Ok(FstResult {
        fst,
        method: FstMethod::WeirCockerham,
    })
}

/// Hudson (1992) Fst estimator across multiple loci.
///
/// Fst = 1 - Hw / Hb, averaged across loci (ratio of averages).
///
/// # Errors
///
/// Returns an error if no loci are provided or populations are empty.
pub fn fst_hudson(
    pop1: &[&[Option<u8>]],
    pop2: &[&[Option<u8>]],
) -> Result<FstResult> {
    if pop1.is_empty() || pop2.is_empty() {
        return Err(CyaneaError::InvalidInput("empty loci or populations".into()));
    }

    let mut sum_num = 0.0;
    let mut sum_den = 0.0;

    for (locus1, locus2) in pop1.iter().zip(pop2.iter()) {
        let (n0_1, n1_1, n2_1, _) = count_genotypes(locus1)?;
        let (n0_2, n1_2, n2_2, _) = count_genotypes(locus2)?;
        let n1 = (n0_1 + n1_1 + n2_1) as f64;
        let n2_pop = (n0_2 + n1_2 + n2_2) as f64;
        if n1 == 0.0 || n2_pop == 0.0 {
            continue;
        }

        let (_, q1) = allele_freq_from_counts(n0_1, n1_1, n2_1);
        let (_, q2) = allele_freq_from_counts(n0_2, n1_2, n2_2);

        // Within-population heterozygosity (bias-corrected)
        let hw1 = 2.0 * q1 * (1.0 - q1) * n1 / (n1 - 1.0);
        let hw2 = 2.0 * q2 * (1.0 - q2) * n2_pop / (n2_pop - 1.0);
        let hw = (hw1 + hw2) / 2.0;

        // Between-population heterozygosity
        let hb = q1 * (1.0 - q2) + q2 * (1.0 - q1);

        // Hudson's numerator/denominator
        sum_num += hb - hw;
        sum_den += hb;
    }

    let fst = if sum_den > 0.0 {
        (sum_num / sum_den).clamp(0.0, 1.0)
    } else {
        0.0
    };

    Ok(FstResult {
        fst,
        method: FstMethod::Hudson,
    })
}

/// Multi-population Fst using the Weir-Cockerham estimator.
///
/// `populations[i][j]` is the genotype vector for population i, locus j.
/// All populations must have the same number of loci.
///
/// # Errors
///
/// Returns an error if fewer than 2 populations or no loci.
pub fn fst_multi_population(populations: &[&[&[Option<u8>]]]) -> Result<FstResult> {
    let r = populations.len();
    if r < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 populations for Fst".into(),
        ));
    }
    let n_loci = populations[0].len();
    if n_loci == 0 {
        return Err(CyaneaError::InvalidInput("no loci provided".into()));
    }
    for pop in populations.iter() {
        if pop.len() != n_loci {
            return Err(CyaneaError::InvalidInput(
                "all populations must have the same number of loci".into(),
            ));
        }
    }

    let mut sum_a = 0.0;
    let mut sum_abc = 0.0;
    let r_f = r as f64;

    for locus_idx in 0..n_loci {
        let mut ni = Vec::with_capacity(r);
        let mut qi = Vec::with_capacity(r);
        let mut hi = Vec::with_capacity(r);

        for pop in populations.iter() {
            let (n0, n1, n2, _) = count_genotypes(pop[locus_idx])?;
            let n = (n0 + n1 + n2) as f64;
            if n == 0.0 {
                ni.push(0.0);
                qi.push(0.0);
                hi.push(0.0);
                continue;
            }
            let (_, q) = allele_freq_from_counts(n0, n1, n2);
            ni.push(n);
            qi.push(q);
            hi.push(n1 as f64 / n);
        }

        let n_total: f64 = ni.iter().sum();
        if n_total == 0.0 {
            continue;
        }

        let n_bar = n_total / r_f;
        let p_bar: f64 = ni.iter().zip(qi.iter()).map(|(n, q)| n * q).sum::<f64>() / n_total;
        let s_sq: f64 = ni
            .iter()
            .zip(qi.iter())
            .map(|(n, q)| n * (q - p_bar).powi(2))
            .sum::<f64>()
            / ((r_f - 1.0) * n_bar);
        let h_bar: f64 = ni.iter().zip(hi.iter()).map(|(n, h)| n * h).sum::<f64>() / n_total;

        let nc = (n_total - ni.iter().map(|n| n * n).sum::<f64>() / n_total) / (r_f - 1.0);

        let a = (n_bar / nc)
            * (s_sq
                - (1.0 / (n_bar - 1.0))
                    * (p_bar * (1.0 - p_bar) - (r_f - 1.0) / r_f * s_sq - 0.25 * h_bar));

        let b = (n_bar / (n_bar - 1.0))
            * (p_bar * (1.0 - p_bar) - (r_f - 1.0) / r_f * s_sq
                - (2.0 * n_bar - 1.0) / (4.0 * n_bar) * h_bar);

        let c = 0.5 * h_bar;

        sum_a += a;
        sum_abc += a + b + c;
    }

    let fst = if sum_abc > 0.0 {
        (sum_a / sum_abc).clamp(0.0, 1.0)
    } else {
        0.0
    };

    Ok(FstResult {
        fst,
        method: FstMethod::WeirCockerham,
    })
}

/// Compute nucleotide diversity (pi) and Watterson's theta from a genotype matrix.
///
/// `genotype_matrix[i]` is the genotype vector for locus i.
/// `n_sequences` is the number of haploid sequences (typically 2× diploid individuals).
///
/// Pi is estimated as the average number of pairwise differences per site.
///
/// # Errors
///
/// Returns an error if `n_sequences` < 2 or no loci are provided.
pub fn diversity(
    genotype_matrix: &[&[Option<u8>]],
    n_sequences: usize,
) -> Result<DiversityStats> {
    if n_sequences < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 sequences for diversity".into(),
        ));
    }
    if genotype_matrix.is_empty() {
        return Err(CyaneaError::InvalidInput("no loci provided".into()));
    }

    let mut seg_sites = 0usize;
    let mut total_pi = 0.0;

    for locus in genotype_matrix {
        let (n0, n1, n2, _) = count_genotypes(locus)?;
        let n_called = n0 + n1 + n2;
        if n_called == 0 {
            continue;
        }
        let (_, q) = allele_freq_from_counts(n0, n1, n2);
        if q > 0.0 && q < 1.0 {
            seg_sites += 1;
        }
        // Per-site pi: for diploid 012 data, estimate from allele frequency
        let n_hap = (2 * n_called) as f64;
        let pi_site = 2.0 * q * (1.0 - q) * n_hap / (n_hap - 1.0);
        total_pi += pi_site;
    }

    let n_loci = genotype_matrix.len() as f64;
    let pi = total_pi / n_loci;

    let a1 = harmonic(n_sequences - 1);
    let theta_w = if a1 > 0.0 {
        (seg_sites as f64 / n_loci) / a1
    } else {
        0.0
    };

    Ok(DiversityStats {
        pi,
        theta_w,
        segregating_sites: seg_sites,
        n_sequences,
    })
}

/// Compute diversity statistics from pre-computed summary values.
///
/// Useful when you already have segregating site counts and average pairwise
/// differences (e.g., from a VCF summary).
///
/// # Errors
///
/// Returns an error if `n_sequences` < 2.
pub fn diversity_from_counts(
    segregating_sites: usize,
    n_sequences: usize,
    avg_pairwise_diff: f64,
) -> Result<DiversityStats> {
    if n_sequences < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 sequences for diversity".into(),
        ));
    }

    let a1 = harmonic(n_sequences - 1);
    let theta_w = if a1 > 0.0 {
        segregating_sites as f64 / a1
    } else {
        0.0
    };

    Ok(DiversityStats {
        pi: avg_pairwise_diff,
        theta_w,
        segregating_sites,
        n_sequences,
    })
}

/// Compute Tajima's D statistic.
///
/// Compares pi (average pairwise differences) to theta_w (from segregating sites).
/// Under neutrality, D ≈ 0. Positive D suggests balancing selection or population
/// contraction; negative D suggests purifying selection or expansion.
///
/// No p-value is provided (the null distribution requires coalescent simulation).
///
/// # Errors
///
/// Returns an error if `n_sequences` < 4 or `segregating_sites` == 0.
pub fn tajimas_d(
    segregating_sites: usize,
    n_sequences: usize,
    avg_pairwise_diff: f64,
) -> Result<TajimaD> {
    if n_sequences < 4 {
        return Err(CyaneaError::InvalidInput(
            "Tajima's D requires at least 4 sequences".into(),
        ));
    }
    if segregating_sites == 0 {
        return Err(CyaneaError::InvalidInput(
            "Tajima's D undefined with 0 segregating sites".into(),
        ));
    }

    let n = n_sequences as f64;
    let s = segregating_sites as f64;

    let a1 = harmonic(n_sequences - 1);
    let a2 = harmonic_sq(n_sequences - 1);

    let theta_w = s / a1;

    // Tajima's D variance components
    let b1 = (n + 1.0) / (3.0 * (n - 1.0));
    let b2 = 2.0 * (n * n + n + 3.0) / (9.0 * n * (n - 1.0));

    let c1 = b1 - 1.0 / a1;
    let c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1);

    let e1 = c1 / a1;
    let e2 = c2 / (a1 * a1 + a2);

    let var_d = e1 * s + e2 * s * (s - 1.0);

    let d = if var_d > 0.0 {
        (avg_pairwise_diff - theta_w) / var_d.sqrt()
    } else {
        0.0
    };

    Ok(TajimaD {
        d,
        pi: avg_pairwise_diff,
        theta_w,
        segregating_sites,
        n_sequences,
    })
}

/// Compute linkage disequilibrium (r², D', D) between two sites.
///
/// Uses the EM algorithm to estimate haplotype frequencies from unphased
/// diploid genotype data.
///
/// # Errors
///
/// Returns an error if sites have different lengths, either site is monomorphic,
/// or all genotypes are missing.
pub fn linkage_disequilibrium(
    site_a: &[Option<u8>],
    site_b: &[Option<u8>],
) -> Result<LdResult> {
    if site_a.len() != site_b.len() {
        return Err(CyaneaError::InvalidInput(
            "site vectors must have the same length".into(),
        ));
    }
    if site_a.is_empty() {
        return Err(CyaneaError::InvalidInput("empty site vectors".into()));
    }

    // Count non-missing pairs
    let mut n_pairs = 0usize;
    for (a, b) in site_a.iter().zip(site_b.iter()) {
        if a.is_some() && b.is_some() {
            n_pairs += 1;
        }
    }
    if n_pairs == 0 {
        return Err(CyaneaError::InvalidInput(
            "no non-missing genotype pairs".into(),
        ));
    }

    // Check that both sites are polymorphic
    let af_a = allele_frequencies(site_a)?;
    let af_b = allele_frequencies(site_b)?;
    if af_a.freq_alt == 0.0 || af_a.freq_alt == 1.0 {
        return Err(CyaneaError::InvalidInput(
            "site A is monomorphic: LD undefined".into(),
        ));
    }
    if af_b.freq_alt == 0.0 || af_b.freq_alt == 1.0 {
        return Err(CyaneaError::InvalidInput(
            "site B is monomorphic: LD undefined".into(),
        ));
    }

    let (p_AB, _p_Ab, _p_aB, _p_ab) = em_haplotype_freqs(site_a, site_b, 1000, 1e-10)?;

    // Marginal frequencies (A = ref at site A, a = alt at site A)
    let p_a = af_a.freq_alt; // freq of alt allele at site A
    let p_A = af_a.freq_ref;
    let p_b = af_b.freq_alt; // freq of alt allele at site B
    let p_B = af_b.freq_ref;

    // D = freq(AB) - freq(A)*freq(B)
    let d = p_AB - p_A * p_B;

    // D'
    let d_max = if d >= 0.0 {
        (p_A * p_b).min(p_a * p_B)
    } else {
        (p_A * p_B).min(p_a * p_b)
    };
    let d_prime = if d_max.abs() > 0.0 {
        (d / d_max).abs().min(1.0)
    } else {
        0.0
    };

    // r²
    let denom = p_A * p_a * p_B * p_b;
    let r_squared = if denom > 0.0 {
        (d * d / denom).min(1.0)
    } else {
        0.0
    };

    Ok(LdResult {
        r_squared,
        d_prime,
        d,
        n_haplotypes: 2 * n_pairs,
    })
}

/// PCA on a genotype matrix with missing-value imputation.
///
/// `genotype_matrix[i]` is the genotype vector for sample i (each sample across all loci).
/// Missing values are imputed with the per-site mean before delegating to
/// [`crate::reduction::pca`].
///
/// # Errors
///
/// Returns an error if the matrix is empty, samples have inconsistent lengths,
/// or `n_components` is invalid.
pub fn genotype_pca(
    genotype_matrix: &[&[Option<u8>]],
    n_components: usize,
) -> Result<crate::reduction::PcaResult> {
    if genotype_matrix.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "empty genotype matrix".into(),
        ));
    }
    let n_features = genotype_matrix[0].len();
    if n_features == 0 {
        return Err(CyaneaError::InvalidInput("zero loci".into()));
    }
    for (i, row) in genotype_matrix.iter().enumerate() {
        if row.len() != n_features {
            return Err(CyaneaError::InvalidInput(format!(
                "sample {} has {} loci, expected {}",
                i,
                row.len(),
                n_features
            )));
        }
    }

    let n_samples = genotype_matrix.len();

    // Compute per-site mean (for imputation)
    let mut site_sum = vec![0.0; n_features];
    let mut site_count = vec![0usize; n_features];
    for row in genotype_matrix {
        for (j, g) in row.iter().enumerate() {
            if let Some(v) = g {
                site_sum[j] += *v as f64;
                site_count[j] += 1;
            }
        }
    }
    let site_mean: Vec<f64> = (0..n_features)
        .map(|j| {
            if site_count[j] > 0 {
                site_sum[j] / site_count[j] as f64
            } else {
                0.0
            }
        })
        .collect();

    // Build imputed f64 matrix
    let mut data: Vec<Vec<f64>> = Vec::with_capacity(n_samples);
    for row in genotype_matrix {
        let imputed: Vec<f64> = row
            .iter()
            .enumerate()
            .map(|(j, g)| match g {
                Some(v) => *v as f64,
                None => site_mean[j],
            })
            .collect();
        data.push(imputed);
    }

    let refs: Vec<&[f64]> = data.iter().map(|v| v.as_slice()).collect();
    crate::reduction::pca(&refs, n_components)
}

// ── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // ── Allele frequencies ───────────────────────────────────────────

    #[test]
    fn af_basic() {
        // 2 hom-ref, 2 het, 1 hom-alt → p=0.6, q=0.4
        let g = [Some(0), Some(0), Some(1), Some(1), Some(2)];
        let af = allele_frequencies(&g).unwrap();
        assert!((af.freq_ref - 0.6).abs() < TOL);
        assert!((af.freq_alt - 0.4).abs() < TOL);
        assert_eq!(af.allele_count, 10);
        assert_eq!(af.missing_count, 0);
        assert!((af.observed_het - 0.4).abs() < TOL);
    }

    #[test]
    fn af_all_hom_ref() {
        let g = [Some(0), Some(0), Some(0)];
        let af = allele_frequencies(&g).unwrap();
        assert!((af.freq_ref - 1.0).abs() < TOL);
        assert!((af.freq_alt - 0.0).abs() < TOL);
        assert!((af.observed_het - 0.0).abs() < TOL);
    }

    #[test]
    fn af_all_het() {
        let g = [Some(1), Some(1), Some(1), Some(1)];
        let af = allele_frequencies(&g).unwrap();
        assert!((af.freq_ref - 0.5).abs() < TOL);
        assert!((af.freq_alt - 0.5).abs() < TOL);
        assert!((af.observed_het - 1.0).abs() < TOL);
    }

    #[test]
    fn af_with_missing() {
        let g = [Some(0), None, Some(1), Some(2), None];
        let af = allele_frequencies(&g).unwrap();
        assert_eq!(af.missing_count, 2);
        assert_eq!(af.allele_count, 6);
    }

    #[test]
    fn af_empty_error() {
        let g: [Option<u8>; 0] = [];
        assert!(allele_frequencies(&g).is_err());
    }

    #[test]
    fn af_invalid_genotype() {
        let g = [Some(0), Some(3), Some(1)];
        assert!(allele_frequencies(&g).is_err());
    }

    // ── HWE ──────────────────────────────────────────────────────────

    #[test]
    fn hwe_perfect_equilibrium() {
        // p=0.5 → expected 25% AA, 50% Aa, 25% aa
        // 25 hom-ref, 50 het, 25 hom-alt
        let mut g = Vec::new();
        for _ in 0..25 {
            g.push(Some(0));
        }
        for _ in 0..50 {
            g.push(Some(1));
        }
        for _ in 0..25 {
            g.push(Some(2));
        }
        let result = hwe_test(&g).unwrap();
        assert!(result.chi_squared < 0.01, "chi_sq={}", result.chi_squared);
        assert!(result.p_value > 0.9);
        assert!(result.inbreeding_coeff.abs() < 0.01);
    }

    #[test]
    fn hwe_excess_heterozygosity() {
        // All het → strong deviation
        let g: Vec<Option<u8>> = vec![Some(1); 100];
        let result = hwe_test(&g).unwrap();
        assert!(result.chi_squared > 10.0);
        assert!(result.p_value < 0.01);
        assert!(result.inbreeding_coeff < -0.5);
    }

    #[test]
    fn hwe_deficit_het() {
        // 50 hom-ref, 0 het, 50 hom-alt → extreme het deficit
        let mut g = Vec::new();
        for _ in 0..50 {
            g.push(Some(0));
        }
        for _ in 0..50 {
            g.push(Some(2));
        }
        let result = hwe_test(&g).unwrap();
        assert!(result.chi_squared > 50.0);
        assert!(result.inbreeding_coeff > 0.9);
    }

    #[test]
    fn hwe_textbook() {
        // 298 AA, 489 Aa, 213 aa (n=1000)
        // Classic textbook example: close to HWE
        let mut g = Vec::new();
        for _ in 0..298 {
            g.push(Some(0));
        }
        for _ in 0..489 {
            g.push(Some(1));
        }
        for _ in 0..213 {
            g.push(Some(2));
        }
        let result = hwe_test(&g).unwrap();
        // These genotype frequencies are close to HWE for p≈0.5425
        // so chi-squared should be small
        assert!(result.chi_squared < 10.0);
    }

    #[test]
    fn hwe_monomorphic_error() {
        let g = [Some(0), Some(0), Some(0)];
        assert!(hwe_test(&g).is_err());
    }

    #[test]
    fn hwe_with_missing() {
        let mut g: Vec<Option<u8>> = Vec::new();
        for _ in 0..20 {
            g.push(Some(0));
        }
        for _ in 0..40 {
            g.push(Some(1));
        }
        for _ in 0..20 {
            g.push(Some(2));
        }
        for _ in 0..10 {
            g.push(None);
        }
        let result = hwe_test(&g).unwrap();
        assert!(result.p_value > 0.5);
    }

    // ── Fst ──────────────────────────────────────────────────────────

    #[test]
    fn fst_identical_populations() {
        let locus: Vec<Option<u8>> = vec![Some(0), Some(1), Some(1), Some(2)];
        let pop1: Vec<&[Option<u8>]> = vec![locus.as_slice()];
        let pop2: Vec<&[Option<u8>]> = vec![locus.as_slice()];
        let result = fst_weir_cockerham(&pop1, &pop2).unwrap();
        assert!(result.fst < 0.05, "fst={}", result.fst);
    }

    #[test]
    fn fst_fully_differentiated() {
        let pop1_locus: Vec<Option<u8>> = vec![Some(0); 50];
        let pop2_locus: Vec<Option<u8>> = vec![Some(2); 50];
        let pop1: Vec<&[Option<u8>]> = vec![pop1_locus.as_slice()];
        let pop2: Vec<&[Option<u8>]> = vec![pop2_locus.as_slice()];
        let result = fst_weir_cockerham(&pop1, &pop2).unwrap();
        assert!(result.fst > 0.8, "fst={}", result.fst);
    }

    #[test]
    fn fst_wc_known() {
        // Two pops with moderate differentiation
        let pop1_locus: Vec<Option<u8>> = vec![Some(0), Some(0), Some(0), Some(1), Some(1)];
        let pop2_locus: Vec<Option<u8>> = vec![Some(1), Some(1), Some(2), Some(2), Some(2)];
        let pop1: Vec<&[Option<u8>]> = vec![pop1_locus.as_slice()];
        let pop2: Vec<&[Option<u8>]> = vec![pop2_locus.as_slice()];
        let result = fst_weir_cockerham(&pop1, &pop2).unwrap();
        assert!(result.fst > 0.1 && result.fst < 0.9, "fst={}", result.fst);
        assert_eq!(result.method, FstMethod::WeirCockerham);
    }

    #[test]
    fn fst_hudson_known() {
        let pop1_locus: Vec<Option<u8>> = vec![Some(0), Some(0), Some(0), Some(1), Some(1)];
        let pop2_locus: Vec<Option<u8>> = vec![Some(1), Some(1), Some(2), Some(2), Some(2)];
        let pop1: Vec<&[Option<u8>]> = vec![pop1_locus.as_slice()];
        let pop2: Vec<&[Option<u8>]> = vec![pop2_locus.as_slice()];
        let result = fst_hudson(&pop1, &pop2).unwrap();
        assert!(result.fst > 0.1 && result.fst < 0.9, "fst={}", result.fst);
        assert_eq!(result.method, FstMethod::Hudson);
    }

    #[test]
    fn fst_multi_locus() {
        // Multiple loci should give a more robust estimate
        let p1_l1: Vec<Option<u8>> = vec![Some(0), Some(0), Some(1), Some(1)];
        let p1_l2: Vec<Option<u8>> = vec![Some(0), Some(1), Some(1), Some(2)];
        let p2_l1: Vec<Option<u8>> = vec![Some(1), Some(2), Some(2), Some(2)];
        let p2_l2: Vec<Option<u8>> = vec![Some(1), Some(1), Some(2), Some(2)];
        let pop1: Vec<&[Option<u8>]> = vec![p1_l1.as_slice(), p1_l2.as_slice()];
        let pop2: Vec<&[Option<u8>]> = vec![p2_l1.as_slice(), p2_l2.as_slice()];
        let result = fst_weir_cockerham(&pop1, &pop2).unwrap();
        assert!(result.fst >= 0.0 && result.fst <= 1.0);
    }

    #[test]
    fn fst_three_populations() {
        let l1: Vec<Option<u8>> = vec![Some(0), Some(0), Some(1)];
        let l2: Vec<Option<u8>> = vec![Some(1), Some(2), Some(2)];
        let l3: Vec<Option<u8>> = vec![Some(0), Some(1), Some(2)];
        let pop1: Vec<&[Option<u8>]> = vec![l1.as_slice()];
        let pop2: Vec<&[Option<u8>]> = vec![l2.as_slice()];
        let pop3: Vec<&[Option<u8>]> = vec![l3.as_slice()];
        let pops: Vec<&[&[Option<u8>]]> = vec![pop1.as_slice(), pop2.as_slice(), pop3.as_slice()];
        let result = super::fst_multi_population(&pops).unwrap();
        assert!(result.fst >= 0.0 && result.fst <= 1.0);
    }

    #[test]
    fn fst_unequal_sizes() {
        let pop1_locus: Vec<Option<u8>> = vec![Some(0), Some(1)];
        let pop2_locus: Vec<Option<u8>> = vec![Some(1), Some(2), Some(2), Some(2), Some(2)];
        let pop1: Vec<&[Option<u8>]> = vec![pop1_locus.as_slice()];
        let pop2: Vec<&[Option<u8>]> = vec![pop2_locus.as_slice()];
        let result = fst_hudson(&pop1, &pop2).unwrap();
        assert!(result.fst >= 0.0 && result.fst <= 1.0);
    }

    #[test]
    fn fst_error_empty() {
        let pop1: Vec<&[Option<u8>]> = vec![];
        let pop2: Vec<&[Option<u8>]> = vec![];
        assert!(fst_weir_cockerham(&pop1, &pop2).is_err());
    }

    // ── Diversity ────────────────────────────────────────────────────

    #[test]
    fn diversity_known_pi() {
        // 4 loci, 2 segregating. With known allele frequencies we can verify pi.
        let l1: Vec<Option<u8>> = vec![Some(0), Some(0), Some(1), Some(2)]; // q=0.375
        let l2: Vec<Option<u8>> = vec![Some(0), Some(0), Some(0), Some(0)]; // monomorphic
        let l3: Vec<Option<u8>> = vec![Some(1), Some(1), Some(1), Some(1)]; // q=0.5
        let l4: Vec<Option<u8>> = vec![Some(0), Some(0), Some(0), Some(0)]; // monomorphic
        let matrix: Vec<&[Option<u8>]> =
            vec![l1.as_slice(), l2.as_slice(), l3.as_slice(), l4.as_slice()];
        let result = diversity(&matrix, 8).unwrap();
        assert!(result.pi > 0.0);
        assert_eq!(result.segregating_sites, 2);
        assert!(result.theta_w > 0.0);
    }

    #[test]
    fn diversity_known_theta_w() {
        // S=10, n=20 → theta_w = S / H_{19}
        let result = diversity_from_counts(10, 20, 5.0).unwrap();
        let h19 = harmonic(19);
        assert!((result.theta_w - 10.0 / h19).abs() < TOL);
    }

    #[test]
    fn diversity_zero_segregating() {
        let l1: Vec<Option<u8>> = vec![Some(0), Some(0), Some(0)];
        let matrix: Vec<&[Option<u8>]> = vec![l1.as_slice()];
        let result = diversity(&matrix, 6).unwrap();
        assert_eq!(result.segregating_sites, 0);
        assert!((result.pi - 0.0).abs() < TOL);
    }

    #[test]
    fn diversity_error_n_lt_2() {
        assert!(diversity_from_counts(5, 1, 2.0).is_err());
    }

    // ── Tajima's D ───────────────────────────────────────────────────

    #[test]
    fn tajimas_d_neutral() {
        // When pi ≈ theta_w, D ≈ 0
        let n = 20;
        let s = 10;
        let a1 = harmonic(n - 1);
        let theta_w = s as f64 / a1;
        let result = tajimas_d(s, n, theta_w).unwrap();
        assert!(result.d.abs() < 0.01, "d={}", result.d);
    }

    #[test]
    fn tajimas_d_positive() {
        // pi >> theta_w → positive D (balancing selection / contraction)
        let result = tajimas_d(5, 20, 10.0).unwrap();
        assert!(result.d > 0.0, "d={}", result.d);
    }

    #[test]
    fn tajimas_d_negative() {
        // pi << theta_w → negative D (purifying selection / expansion)
        let result = tajimas_d(50, 20, 1.0).unwrap();
        assert!(result.d < 0.0, "d={}", result.d);
    }

    #[test]
    fn tajimas_d_textbook() {
        // Textbook example: n=10, S=12, pi=3.5
        let result = tajimas_d(12, 10, 3.5).unwrap();
        // D should be negative since pi (3.5) < S/a1 (12/2.829 ≈ 4.24)
        assert!(result.d < 0.0, "d={}", result.d);
        assert_eq!(result.segregating_sites, 12);
        assert_eq!(result.n_sequences, 10);
    }

    #[test]
    fn tajimas_d_n_lt_4_error() {
        assert!(tajimas_d(5, 3, 2.0).is_err());
    }

    #[test]
    fn tajimas_d_zero_s_error() {
        assert!(tajimas_d(0, 20, 0.0).is_err());
    }

    // ── Linkage disequilibrium ───────────────────────────────────────

    #[test]
    fn ld_perfect() {
        // Perfect LD: sites perfectly correlated
        let site_a: Vec<Option<u8>> = vec![Some(0), Some(0), Some(2), Some(2)];
        let site_b: Vec<Option<u8>> = vec![Some(0), Some(0), Some(2), Some(2)];
        let result = linkage_disequilibrium(&site_a, &site_b).unwrap();
        assert!(result.r_squared > 0.9, "r²={}", result.r_squared);
        assert!(result.d_prime > 0.9, "D'={}", result.d_prime);
    }

    #[test]
    fn ld_no_ld() {
        // Sites with no LD (independent)
        // Create a larger sample to get stable LD estimates
        let mut site_a = Vec::new();
        let mut site_b = Vec::new();
        // AB, Ab, aB, ab equally frequent → no LD
        for _ in 0..25 {
            site_a.push(Some(0));
            site_b.push(Some(0));
        }
        for _ in 0..25 {
            site_a.push(Some(0));
            site_b.push(Some(2));
        }
        for _ in 0..25 {
            site_a.push(Some(2));
            site_b.push(Some(0));
        }
        for _ in 0..25 {
            site_a.push(Some(2));
            site_b.push(Some(2));
        }
        let result = linkage_disequilibrium(&site_a, &site_b).unwrap();
        assert!(result.r_squared < 0.05, "r²={}", result.r_squared);
    }

    #[test]
    fn ld_intermediate() {
        // Some LD but not perfect
        let site_a: Vec<Option<u8>> =
            vec![Some(0), Some(0), Some(0), Some(2), Some(2), Some(1)];
        let site_b: Vec<Option<u8>> =
            vec![Some(0), Some(0), Some(2), Some(2), Some(2), Some(1)];
        let result = linkage_disequilibrium(&site_a, &site_b).unwrap();
        assert!(result.r_squared >= 0.0 && result.r_squared <= 1.0);
        assert!(result.d_prime >= 0.0 && result.d_prime <= 1.0);
    }

    #[test]
    fn ld_d_prime_vs_r_squared() {
        // D' and r² measure different things: D' can be 1 when r² < 1
        let site_a: Vec<Option<u8>> =
            vec![Some(0), Some(0), Some(0), Some(0), Some(1), Some(1), Some(2)];
        let site_b: Vec<Option<u8>> =
            vec![Some(0), Some(0), Some(0), Some(1), Some(1), Some(2), Some(2)];
        let result = linkage_disequilibrium(&site_a, &site_b).unwrap();
        // Both should be in valid range
        assert!(result.r_squared >= 0.0 && result.r_squared <= 1.0);
        assert!(result.d_prime >= 0.0 && result.d_prime <= 1.0);
    }

    #[test]
    fn ld_with_missing() {
        let site_a: Vec<Option<u8>> = vec![Some(0), Some(0), None, Some(2), Some(2)];
        let site_b: Vec<Option<u8>> = vec![Some(0), None, Some(0), Some(2), Some(2)];
        let result = linkage_disequilibrium(&site_a, &site_b).unwrap();
        assert!(result.n_haplotypes > 0);
        assert!(result.r_squared >= 0.0);
    }

    #[test]
    fn ld_monomorphic_error() {
        let site_a: Vec<Option<u8>> = vec![Some(0), Some(0), Some(0)];
        let site_b: Vec<Option<u8>> = vec![Some(0), Some(1), Some(2)];
        assert!(linkage_disequilibrium(&site_a, &site_b).is_err());
    }

    #[test]
    fn ld_length_mismatch() {
        let site_a: Vec<Option<u8>> = vec![Some(0), Some(1)];
        let site_b: Vec<Option<u8>> = vec![Some(0), Some(1), Some(2)];
        assert!(linkage_disequilibrium(&site_a, &site_b).is_err());
    }

    // ── Genotype PCA ─────────────────────────────────────────────────

    #[test]
    fn pca_basic() {
        let s1: Vec<Option<u8>> = vec![Some(0), Some(0), Some(1), Some(2)];
        let s2: Vec<Option<u8>> = vec![Some(0), Some(1), Some(1), Some(2)];
        let s3: Vec<Option<u8>> = vec![Some(2), Some(2), Some(1), Some(0)];
        let s4: Vec<Option<u8>> = vec![Some(2), Some(1), Some(0), Some(0)];
        let matrix: Vec<&[Option<u8>]> =
            vec![s1.as_slice(), s2.as_slice(), s3.as_slice(), s4.as_slice()];
        let result = genotype_pca(&matrix, 2).unwrap();
        assert_eq!(result.components.len(), 2);
        assert_eq!(result.transformed.len(), 4);
        assert_eq!(result.transformed[0].len(), 2);
    }

    #[test]
    fn pca_missing_imputation() {
        // Missing values should be imputed without error
        let s1: Vec<Option<u8>> = vec![Some(0), None, Some(1)];
        let s2: Vec<Option<u8>> = vec![Some(1), Some(1), Some(2)];
        let s3: Vec<Option<u8>> = vec![Some(2), Some(2), None];
        let matrix: Vec<&[Option<u8>]> = vec![s1.as_slice(), s2.as_slice(), s3.as_slice()];
        let result = genotype_pca(&matrix, 2).unwrap();
        assert_eq!(result.transformed.len(), 3);
    }

    #[test]
    fn pca_variance() {
        // Variance ratios should be non-negative and sum ≤ 1
        let s1: Vec<Option<u8>> = vec![Some(0), Some(0), Some(2), Some(2)];
        let s2: Vec<Option<u8>> = vec![Some(0), Some(1), Some(1), Some(2)];
        let s3: Vec<Option<u8>> = vec![Some(2), Some(2), Some(0), Some(0)];
        let s4: Vec<Option<u8>> = vec![Some(1), Some(1), Some(1), Some(1)];
        let s5: Vec<Option<u8>> = vec![Some(0), Some(2), Some(0), Some(2)];
        let matrix: Vec<&[Option<u8>]> = vec![
            s1.as_slice(),
            s2.as_slice(),
            s3.as_slice(),
            s4.as_slice(),
            s5.as_slice(),
        ];
        let result = genotype_pca(&matrix, 3).unwrap();
        for r in &result.explained_variance_ratio {
            assert!(*r >= 0.0);
        }
        let total: f64 = result.explained_variance_ratio.iter().sum();
        assert!(total <= 1.0 + 0.01);
    }

    #[test]
    fn pca_population_structure() {
        // Two distinct populations should separate along PC1
        let mut matrix: Vec<Vec<Option<u8>>> = Vec::new();
        // Pop A: mostly 0s
        for _ in 0..10 {
            matrix.push(vec![Some(0), Some(0), Some(0), Some(0), Some(1)]);
        }
        // Pop B: mostly 2s
        for _ in 0..10 {
            matrix.push(vec![Some(2), Some(2), Some(2), Some(2), Some(1)]);
        }
        let refs: Vec<&[Option<u8>]> = matrix.iter().map(|v| v.as_slice()).collect();
        let result = genotype_pca(&refs, 2).unwrap();
        // PC1 should explain most variance
        assert!(result.explained_variance_ratio[0] > 0.5);
        // Pop A and Pop B should be on opposite sides of PC1
        let mean_a: f64 = result.transformed[..10].iter().map(|t| t[0]).sum::<f64>() / 10.0;
        let mean_b: f64 = result.transformed[10..].iter().map(|t| t[0]).sum::<f64>() / 10.0;
        assert!((mean_a - mean_b).abs() > 0.1);
    }

    #[test]
    fn pca_empty_error() {
        let matrix: Vec<&[Option<u8>]> = vec![];
        assert!(genotype_pca(&matrix, 1).is_err());
    }

    // ── Edge cases ───────────────────────────────────────────────────

    #[test]
    fn u8_convenience() {
        let g: Vec<u8> = vec![0, 1, 2, 255, 0];
        let af = allele_frequencies_u8(&g).unwrap();
        assert_eq!(af.missing_count, 1);
        assert_eq!(af.allele_count, 8);
    }

    #[test]
    fn all_missing_error() {
        let g = [None, None, None];
        assert!(allele_frequencies(&g).is_err());
    }

    #[test]
    fn af_consistency() {
        // freq_ref + freq_alt should always sum to 1
        let g = [Some(0), Some(0), Some(1), Some(2), Some(2), Some(2)];
        let af = allele_frequencies(&g).unwrap();
        assert!((af.freq_ref + af.freq_alt - 1.0).abs() < TOL);
    }
}
