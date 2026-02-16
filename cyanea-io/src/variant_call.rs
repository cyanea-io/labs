//! Bayesian variant caller with genotype likelihoods, indel support,
//! strand bias detection, and VCF output.
//!
//! Transforms [`Pileup`] data into genotyped variant calls using a diploid
//! genotype model with configurable priors, quality filters, and strand bias
//! detection via Fisher's exact test.

use cyanea_omics::variant::{Variant, VariantFilter};
use cyanea_stats::testing::fisher_exact;

use crate::pileup::{base_index, Pileup, PileupColumn};

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Configuration for the Bayesian variant caller.
#[derive(Debug, Clone)]
pub struct VariantCallConfig {
    /// Minimum depth to consider a site (default: 5).
    pub min_depth: u32,
    /// Minimum base quality to include a read (default: 20).
    pub min_base_quality: u8,
    /// Minimum mapping quality to include a read (default: 20).
    pub min_mapq: u8,
    /// Minimum alternate allele read count (default: 2).
    pub min_alt_count: u32,
    /// Minimum QUAL to emit a variant (default: 10.0).
    pub min_qual: f64,
    /// Minimum genotype quality (default: 20.0).
    pub min_gq: f64,
    /// Strand bias Fisher p-value threshold (default: 0.001).
    pub strand_bias_threshold: f64,
    /// Prior probability of a heterozygous genotype (default: 0.001).
    pub prior_het: f64,
    /// Prior probability of a homozygous alternate genotype (default: 0.0005).
    pub prior_hom_alt: f64,
}

impl Default for VariantCallConfig {
    fn default() -> Self {
        Self {
            min_depth: 5,
            min_base_quality: 20,
            min_mapq: 20,
            min_alt_count: 2,
            min_qual: 10.0,
            min_gq: 20.0,
            strand_bias_threshold: 0.001,
            prior_het: 0.001,
            prior_hom_alt: 0.0005,
        }
    }
}

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

/// Diploid genotype call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Genotype {
    /// Homozygous reference (0/0).
    HomRef,
    /// Heterozygous (0/1).
    Het,
    /// Homozygous alternate (1/1).
    HomAlt,
}

/// A called variant with genotype information.
#[derive(Debug, Clone)]
pub struct CalledVariant {
    /// Underlying variant record (chrom, pos, ref, alt, quality, filter).
    pub variant: Variant,
    /// Called genotype.
    pub genotype: Genotype,
    /// Phred-scaled genotype quality.
    pub genotype_quality: f64,
    /// Total depth at this site.
    pub depth: u32,
    /// Allele depth: [ref_count, alt_count].
    pub allele_depth: [u32; 2],
    /// Forward-strand allele depth: [ref_fwd, alt_fwd].
    pub allele_depth_fwd: [u32; 2],
    /// Reverse-strand allele depth: [ref_rev, alt_rev].
    pub allele_depth_rev: [u32; 2],
    /// Strand bias Fisher exact test p-value.
    pub strand_bias_p: f64,
    /// Phred-scaled genotype likelihoods (PL): [hom-ref, het, hom-alt].
    /// Normalized so the best genotype has PL = 0.
    pub genotype_likelihoods: [u32; 3],
}

/// Summary statistics for a set of called variants.
#[derive(Debug, Clone)]
pub struct VariantCallStats {
    /// Total number of called variants.
    pub total: usize,
    /// Number of single nucleotide variants.
    pub snvs: usize,
    /// Number of insertions.
    pub insertions: usize,
    /// Number of deletions.
    pub deletions: usize,
    /// Number of heterozygous calls.
    pub het_count: usize,
    /// Number of homozygous alternate calls.
    pub hom_alt_count: usize,
    /// Transition/transversion ratio (NaN if no transversions).
    pub ti_tv_ratio: f64,
    /// Number of variants passing all filters.
    pub pass_count: usize,
}

// ---------------------------------------------------------------------------
// Log-space math helpers
// ---------------------------------------------------------------------------

/// Numerically stable log-sum-exp: ln(exp(a) + exp(b)).
fn ln_sum_exp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let max = a.max(b);
    max + ((a - max).exp() + (b - max).exp()).ln()
}

/// Numerically stable log-sum-exp for three values.
fn ln_sum_exp3(a: f64, b: f64, c: f64) -> f64 {
    ln_sum_exp(ln_sum_exp(a, b), c)
}

// ---------------------------------------------------------------------------
// Genotype likelihood computation
// ---------------------------------------------------------------------------

/// Compute log-likelihoods for three diploid genotypes at an SNP site.
///
/// Returns `[ln P(data|RR), ln P(data|RA), ln P(data|AA)]` where
/// R = reference allele, A = alternate allele.
fn genotype_likelihoods_snp(
    bases: &[u8],
    qualities: &[u8],
    mapqs: &[u8],
    ref_base: u8,
    alt_base: u8,
    config: &VariantCallConfig,
) -> [f64; 3] {
    let ref_upper = ref_base.to_ascii_uppercase();
    let alt_upper = alt_base.to_ascii_uppercase();

    let mut ll_rr = 0.0_f64; // ln P(data | hom-ref)
    let mut ll_ra = 0.0_f64; // ln P(data | het)
    let mut ll_aa = 0.0_f64; // ln P(data | hom-alt)

    for i in 0..bases.len() {
        let q = qualities[i];
        let mq = mapqs[i];

        if q < config.min_base_quality || mq < config.min_mapq {
            continue;
        }

        let base_upper = bases[i].to_ascii_uppercase();
        let e = 10.0_f64.powf(-(q as f64) / 10.0); // error probability

        // P(base | genotype) for each genotype
        let p_rr = if base_upper == ref_upper {
            1.0 - e
        } else {
            e / 3.0
        };
        let p_aa = if base_upper == alt_upper {
            1.0 - e
        } else {
            e / 3.0
        };
        let p_ra = 0.5 * p_rr_for_base(base_upper, ref_upper, e)
            + 0.5 * p_rr_for_base(base_upper, alt_upper, e);

        ll_rr += p_rr.ln();
        ll_ra += p_ra.ln();
        ll_aa += p_aa.ln();
    }

    [ll_rr, ll_ra, ll_aa]
}

/// Probability of observing `base` given the true allele is `allele`, with error rate `e`.
fn p_rr_for_base(base: u8, allele: u8, e: f64) -> f64 {
    if base == allele {
        1.0 - e
    } else {
        e / 3.0
    }
}

/// Given log-likelihoods and priors, compute QUAL, GQ, genotype, and PL.
///
/// Returns `(qual, gq, genotype, pl)`.
fn call_genotype(likelihoods: [f64; 3], config: &VariantCallConfig) -> (f64, f64, Genotype, [u32; 3]) {
    let prior_rr = (1.0 - config.prior_het - config.prior_hom_alt).ln();
    let prior_ra = config.prior_het.ln();
    let prior_aa = config.prior_hom_alt.ln();

    // Log-posteriors (unnormalized)
    let post_rr = likelihoods[0] + prior_rr;
    let post_ra = likelihoods[1] + prior_ra;
    let post_aa = likelihoods[2] + prior_aa;

    // Normalize
    let log_total = ln_sum_exp3(post_rr, post_ra, post_aa);
    let norm_rr = post_rr - log_total;
    let norm_ra = post_ra - log_total;
    let norm_aa = post_aa - log_total;

    // QUAL = -10 * log10(P(RR | data))
    let qual = -10.0 * norm_rr * std::f64::consts::LOG10_E;
    let qual = qual.max(0.0);

    // Phred-scaled likelihoods (PL) — from raw log-likelihoods, not posteriors
    let max_ll = likelihoods[0].max(likelihoods[1]).max(likelihoods[2]);
    let pl: [u32; 3] = [
        ((-10.0 * (likelihoods[0] - max_ll) * std::f64::consts::LOG10_E).round() as u32),
        ((-10.0 * (likelihoods[1] - max_ll) * std::f64::consts::LOG10_E).round() as u32),
        ((-10.0 * (likelihoods[2] - max_ll) * std::f64::consts::LOG10_E).round() as u32),
    ];

    // Genotype = best posterior
    let genotype = if norm_ra >= norm_rr && norm_ra >= norm_aa {
        Genotype::Het
    } else if norm_aa >= norm_rr {
        Genotype::HomAlt
    } else {
        Genotype::HomRef
    };

    // GQ = second-best PL - best PL = second-best PL (since best = 0)
    let mut sorted_pl = [pl[0], pl[1], pl[2]];
    sorted_pl.sort_unstable();
    let gq = sorted_pl[1] as f64; // second smallest (best is 0)

    (qual, gq, genotype, pl)
}

// ---------------------------------------------------------------------------
// Strand bias
// ---------------------------------------------------------------------------

/// Compute strand bias p-value using Fisher's exact test.
///
/// Constructs a 2x2 contingency table:
/// ```text
///           Forward  Reverse
/// Ref          a        b
/// Alt          c        d
/// ```
fn strand_bias(ref_fwd: usize, ref_rev: usize, alt_fwd: usize, alt_rev: usize) -> f64 {
    let table = [[ref_fwd, ref_rev], [alt_fwd, alt_rev]];
    fisher_exact(&table).map(|r| r.p_value).unwrap_or(1.0)
}

// ---------------------------------------------------------------------------
// Filters
// ---------------------------------------------------------------------------

/// Apply quality filters to a variant call.
fn compute_filters(
    qual: f64,
    gq: f64,
    depth: u32,
    sb_p: f64,
    config: &VariantCallConfig,
) -> VariantFilter {
    let mut reasons = Vec::new();

    if depth < config.min_depth {
        reasons.push("LowDepth".to_string());
    }
    if qual < config.min_qual {
        reasons.push("LowQual".to_string());
    }
    if gq < config.min_gq {
        reasons.push("LowGQ".to_string());
    }
    if sb_p < config.strand_bias_threshold {
        reasons.push("StrandBias".to_string());
    }

    if reasons.is_empty() {
        VariantFilter::Pass
    } else {
        VariantFilter::Fail(reasons)
    }
}

// ---------------------------------------------------------------------------
// SNP calling
// ---------------------------------------------------------------------------

/// Count strand-specific alleles from a pileup column.
///
/// Returns (ref_fwd, ref_rev, alt_fwd, alt_rev) for a given alt base.
fn count_stranded_alleles(col: &PileupColumn, ref_base: u8, alt_base: u8) -> (u32, u32, u32, u32) {
    let ref_upper = ref_base.to_ascii_uppercase();
    let alt_upper = alt_base.to_ascii_uppercase();
    let ref_lower = ref_upper.to_ascii_lowercase();
    let alt_lower = alt_upper.to_ascii_lowercase();

    let mut ref_fwd = 0u32;
    let mut ref_rev = 0u32;
    let mut alt_fwd = 0u32;
    let mut alt_rev = 0u32;

    for &b in &col.bases {
        if b == ref_upper {
            ref_fwd += 1;
        } else if b == ref_lower {
            ref_rev += 1;
        } else if b == alt_upper {
            alt_fwd += 1;
        } else if b == alt_lower {
            alt_rev += 1;
        }
    }

    (ref_fwd, ref_rev, alt_fwd, alt_rev)
}

/// Find the best alternate allele at a position (highest count, excluding ref and deletions).
fn best_alt_allele(col: &PileupColumn, ref_base: u8) -> Option<(u8, u32)> {
    let ref_idx = base_index(ref_base.to_ascii_uppercase());

    let mut best_idx = usize::MAX;
    let mut best_count = 0u32;

    for i in 0..5 {
        // A=0, C=1, G=2, T=3, N=4 (skip del=5)
        if i == ref_idx {
            continue;
        }
        if col.base_counts[i] > best_count {
            best_count = col.base_counts[i];
            best_idx = i;
        }
    }

    if best_idx == usize::MAX || best_count == 0 {
        return None;
    }

    let alt_base = match best_idx {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    };

    Some((alt_base, best_count))
}

/// Call SNPs from a single pileup column.
fn call_snp_at_position(col: &PileupColumn, config: &VariantCallConfig) -> Option<CalledVariant> {
    if col.depth < config.min_depth {
        return None;
    }

    let ref_base = col.ref_base.to_ascii_uppercase();
    if ref_base == b'N' {
        return None;
    }

    let (alt_base, alt_count) = best_alt_allele(col, ref_base)?;
    if alt_count < config.min_alt_count {
        return None;
    }

    // Compute genotype likelihoods
    let lls = genotype_likelihoods_snp(
        &col.bases,
        &col.qualities,
        &col.mapping_qualities,
        ref_base,
        alt_base,
        config,
    );

    let (qual, gq, genotype, pl) = call_genotype(lls, config);

    // Don't emit hom-ref calls
    if genotype == Genotype::HomRef {
        return None;
    }

    // Strand-specific counts
    let (ref_fwd, ref_rev, alt_fwd, alt_rev) = count_stranded_alleles(col, ref_base, alt_base);
    let sb_p = strand_bias(ref_fwd as usize, ref_rev as usize, alt_fwd as usize, alt_rev as usize);

    let filter = compute_filters(qual, gq, col.depth, sb_p, config);

    let mut variant = Variant::new(
        &col.rname,
        col.pos + 1, // VCF is 1-based
        vec![ref_base],
        vec![vec![alt_base]],
    )
    .ok()?;
    variant.quality = Some(qual);
    variant.filter = filter;

    Some(CalledVariant {
        variant,
        genotype,
        genotype_quality: gq,
        depth: col.depth,
        allele_depth: [ref_fwd + ref_rev, alt_fwd + alt_rev],
        allele_depth_fwd: [ref_fwd, alt_fwd],
        allele_depth_rev: [ref_rev, alt_rev],
        strand_bias_p: sb_p,
        genotype_likelihoods: pl,
    })
}

// ---------------------------------------------------------------------------
// Deletion calling
// ---------------------------------------------------------------------------

/// Scan for deletion runs across consecutive pileup columns.
fn call_deletions(pileup: &Pileup, config: &VariantCallConfig) -> Vec<CalledVariant> {
    let mut variants = Vec::new();
    let cols = &pileup.columns;
    let mut i = 0;

    while i < cols.len() {
        let del_count = cols[i].base_counts[5]; // deletion marker index
        if del_count < config.min_alt_count || cols[i].depth < config.min_depth {
            i += 1;
            continue;
        }

        // Find the extent of consecutive deletion positions
        let start = i;
        let mut end = i;
        while end + 1 < cols.len()
            && cols[end + 1].base_counts[5] >= config.min_alt_count
            && cols[end + 1].pos == cols[end].pos + 1
        {
            end += 1;
        }

        // Need an anchor base before the deletion
        let anchor_idx = if start > 0 && cols[start - 1].pos == cols[start].pos - 1 {
            Some(start - 1)
        } else {
            None
        };

        if let Some(anchor_idx) = anchor_idx {
            let anchor_col = &cols[anchor_idx];
            let anchor_base = anchor_col.ref_base.to_ascii_uppercase();
            if anchor_base != b'N' {
                // Build REF = anchor + deleted bases, ALT = anchor
                let mut ref_allele = vec![anchor_base];
                for col in cols.iter().take(end + 1).skip(start) {
                    let del_base = col.ref_base.to_ascii_uppercase();
                    ref_allele.push(if del_base == b'N' { b'N' } else { del_base });
                }
                let alt_allele = vec![anchor_base];

                // Use deletion count from first deletion position for likelihood
                let total_depth = anchor_col.depth;
                let ref_count = total_depth.saturating_sub(del_count);

                let (lls, ref_fwd, ref_rev, alt_fwd, alt_rev) =
                    deletion_likelihoods(&cols[start..=end], anchor_col, config);

                let (qual, gq, genotype, pl) = call_genotype(lls, config);

                if genotype == Genotype::HomRef {
                    i = end + 1;
                    continue;
                }

                let sb_p = strand_bias(
                    ref_fwd as usize,
                    ref_rev as usize,
                    alt_fwd as usize,
                    alt_rev as usize,
                );
                let filter = compute_filters(qual, gq, total_depth, sb_p, config);

                if let Ok(mut variant) = Variant::new(
                    &anchor_col.rname,
                    anchor_col.pos + 1,
                    ref_allele,
                    vec![alt_allele],
                ) {
                    variant.quality = Some(qual);
                    variant.filter = filter;

                    variants.push(CalledVariant {
                        variant,
                        genotype,
                        genotype_quality: gq,
                        depth: total_depth,
                        allele_depth: [ref_count, del_count],
                        allele_depth_fwd: [ref_fwd, alt_fwd],
                        allele_depth_rev: [ref_rev, alt_rev],
                        strand_bias_p: sb_p,
                        genotype_likelihoods: pl,
                    });
                }

                // Skip past the deletion region
                i = end + 1;
                continue;
            }
        }

        i = end + 1;
    }

    variants
}

/// Compute genotype likelihoods for a deletion.
///
/// Returns (log-likelihoods, ref_fwd, ref_rev, del_fwd, del_rev).
fn deletion_likelihoods(
    del_cols: &[PileupColumn],
    anchor_col: &PileupColumn,
    _config: &VariantCallConfig,
) -> ([f64; 3], u32, u32, u32, u32) {
    // Count strand-specific deletion evidence from the first deletion column
    let first = &del_cols[0];
    let mut del_fwd = 0u32;
    let mut del_rev = 0u32;
    for &b in &first.bases {
        if b == b'*' {
            del_fwd += 1;
        } else if b == b'#' {
            del_rev += 1;
        }
    }

    // Count non-deletion reads at anchor
    let ref_upper = anchor_col.ref_base.to_ascii_uppercase();
    let mut ref_fwd = 0u32;
    let mut ref_rev = 0u32;
    for &b in &anchor_col.bases {
        let bu = b.to_ascii_uppercase();
        if bu == ref_upper {
            if b.is_ascii_uppercase() {
                ref_fwd += 1;
            } else {
                ref_rev += 1;
            }
        }
    }

    let del_total = del_fwd + del_rev;
    let ref_total = ref_fwd + ref_rev;

    // Simple likelihood model: each read either supports del or ref
    let avg_qual = 30.0_f64; // assume moderate quality for deletions
    let e = 10.0_f64.powf(-avg_qual / 10.0);

    let mut ll_rr = 0.0_f64;
    let mut ll_ra = 0.0_f64;
    let mut ll_aa = 0.0_f64;

    for _ in 0..ref_total {
        ll_rr += (1.0 - e).ln();
        ll_ra += (0.5 * (1.0 - e) + 0.5 * e).ln();
        ll_aa += e.ln();
    }
    for _ in 0..del_total {
        ll_rr += e.ln();
        ll_ra += (0.5 * e + 0.5 * (1.0 - e)).ln();
        ll_aa += (1.0 - e).ln();
    }

    ([ll_rr, ll_ra, ll_aa], ref_fwd, ref_rev, del_fwd, del_rev)
}

// ---------------------------------------------------------------------------
// Insertion calling
// ---------------------------------------------------------------------------

/// Call insertions from insertion evidence in the pileup.
fn call_insertions(pileup: &Pileup, config: &VariantCallConfig) -> Vec<CalledVariant> {
    let mut variants = Vec::new();

    for col in &pileup.columns {
        if col.insertions.is_empty() || col.depth < config.min_depth {
            continue;
        }

        let anchor_base = col.ref_base.to_ascii_uppercase();
        if anchor_base == b'N' {
            continue;
        }

        // Group identical insertion sequences (case-insensitive) and count
        let mut ins_groups: std::collections::HashMap<Vec<u8>, (u32, u32, u32)> =
            std::collections::HashMap::new();
        for ev in &col.insertions {
            let key: Vec<u8> = ev.bases.iter().map(|b| b.to_ascii_uppercase()).collect();
            let entry = ins_groups.entry(key).or_insert((0, 0, 0));
            entry.0 += 1; // total count
            // Check strand from first base case
            if !ev.bases.is_empty() && ev.bases[0].is_ascii_uppercase() {
                entry.1 += 1; // fwd
            } else {
                entry.2 += 1; // rev
            }
        }

        // Process each insertion allele
        for (ins_seq, (count, ins_fwd, ins_rev)) in &ins_groups {
            if *count < config.min_alt_count {
                continue;
            }

            // REF = anchor base, ALT = anchor + inserted sequence
            let mut alt_allele = vec![anchor_base];
            alt_allele.extend_from_slice(ins_seq);

            // Count ref (non-insertion) reads at this position
            let ref_total = col.depth - col.insertions.len().min(col.depth as usize) as u32;
            let ref_fwd = ref_total / 2; // approximate
            let ref_rev = ref_total - ref_fwd;

            // Simple genotype likelihood for insertion
            let avg_qual = if !col.insertions.is_empty() {
                let sum: u32 = col
                    .insertions
                    .iter()
                    .map(|e| e.mean_quality as u32)
                    .sum();
                (sum / col.insertions.len() as u32).max(1) as f64
            } else {
                30.0
            };

            let e = 10.0_f64.powf(-avg_qual / 10.0);
            let mut ll_rr = 0.0_f64;
            let mut ll_ra = 0.0_f64;
            let mut ll_aa = 0.0_f64;

            // Reads supporting reference (no insertion)
            let no_ins = col.depth.saturating_sub(*count);
            for _ in 0..no_ins {
                ll_rr += (1.0 - e).ln();
                ll_ra += (0.5 * (1.0 - e) + 0.5 * e).ln();
                ll_aa += e.ln();
            }
            // Reads supporting insertion
            for _ in 0..*count {
                ll_rr += e.ln();
                ll_ra += (0.5 * e + 0.5 * (1.0 - e)).ln();
                ll_aa += (1.0 - e).ln();
            }

            let (qual, gq, genotype, pl) = call_genotype([ll_rr, ll_ra, ll_aa], config);

            if genotype == Genotype::HomRef {
                continue;
            }

            let sb_p = strand_bias(
                ref_fwd as usize,
                ref_rev as usize,
                *ins_fwd as usize,
                *ins_rev as usize,
            );
            let filter = compute_filters(qual, gq, col.depth, sb_p, config);

            if let Ok(mut variant) = Variant::new(
                &col.rname,
                col.pos + 1,
                vec![anchor_base],
                vec![alt_allele],
            ) {
                variant.quality = Some(qual);
                variant.filter = filter;

                variants.push(CalledVariant {
                    variant,
                    genotype,
                    genotype_quality: gq,
                    depth: col.depth,
                    allele_depth: [no_ins, *count],
                    allele_depth_fwd: [ref_fwd, *ins_fwd],
                    allele_depth_rev: [ref_rev, *ins_rev],
                    strand_bias_p: sb_p,
                    genotype_likelihoods: pl,
                });
            }
        }
    }

    variants
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Call variants from a single pileup using the Bayesian genotype model.
///
/// Calls SNPs, deletions, and insertions. Results are sorted by position.
pub fn call_variants(pileup: &Pileup, config: &VariantCallConfig) -> Vec<CalledVariant> {
    let mut variants = Vec::new();

    // SNPs
    for col in &pileup.columns {
        if let Some(v) = call_snp_at_position(col, config) {
            variants.push(v);
        }
    }

    // Deletions
    variants.extend(call_deletions(pileup, config));

    // Insertions
    variants.extend(call_insertions(pileup, config));

    // Sort by position
    variants.sort_by_key(|v| v.variant.position);
    variants
}

/// Call variants from multiple pileups (one per reference sequence).
pub fn call_variants_all(pileups: &[Pileup], config: &VariantCallConfig) -> Vec<CalledVariant> {
    let mut variants = Vec::new();
    for p in pileups {
        variants.extend(call_variants(p, config));
    }
    variants
}

/// Compute summary statistics for a set of called variants.
pub fn variant_call_stats(variants: &[CalledVariant]) -> VariantCallStats {
    let mut snvs = 0usize;
    let mut insertions = 0usize;
    let mut deletions = 0usize;
    let mut het_count = 0usize;
    let mut hom_alt_count = 0usize;
    let mut transitions = 0usize;
    let mut transversions = 0usize;
    let mut pass_count = 0usize;

    for cv in variants {
        let v = &cv.variant;

        if v.is_snv() {
            snvs += 1;
            if v.is_transition() {
                transitions += 1;
            } else {
                transversions += 1;
            }
        } else if v.ref_allele.len() < v.alt_alleles[0].len() {
            insertions += 1;
        } else if v.ref_allele.len() > v.alt_alleles[0].len() {
            deletions += 1;
        }

        match cv.genotype {
            Genotype::Het => het_count += 1,
            Genotype::HomAlt => hom_alt_count += 1,
            Genotype::HomRef => {}
        }

        if cv.variant.filter == VariantFilter::Pass {
            pass_count += 1;
        }
    }

    let ti_tv_ratio = if transversions > 0 {
        transitions as f64 / transversions as f64
    } else {
        f64::NAN
    };

    VariantCallStats {
        total: variants.len(),
        snvs,
        insertions,
        deletions,
        het_count,
        hom_alt_count,
        ti_tv_ratio,
        pass_count,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::pileup;
    use crate::sam::SamRecord;
    use std::collections::HashMap;

    fn make_record(
        qname: &str,
        flag: u16,
        rname: &str,
        pos: u64,
        mapq: u8,
        cigar: &str,
        seq: &str,
        qual: &str,
    ) -> SamRecord {
        SamRecord {
            qname: qname.to_string(),
            flag,
            rname: rname.to_string(),
            pos,
            mapq,
            cigar: cigar.to_string(),
            rnext: "*".to_string(),
            pnext: 0,
            tlen: 0,
            sequence: seq.to_string(),
            quality: qual.to_string(),
        }
    }

    fn lenient_config() -> VariantCallConfig {
        VariantCallConfig {
            min_depth: 1,
            min_base_quality: 0,
            min_mapq: 0,
            min_alt_count: 1,
            min_qual: 0.0,
            min_gq: 0.0,
            strand_bias_threshold: 0.0, // never filter on SB
            prior_het: 0.001,
            prior_hom_alt: 0.0005,
        }
    }

    // -- Config defaults --

    #[test]
    fn test_config_defaults() {
        let config = VariantCallConfig::default();
        assert_eq!(config.min_depth, 5);
        assert_eq!(config.min_base_quality, 20);
        assert_eq!(config.min_mapq, 20);
        assert_eq!(config.min_alt_count, 2);
        assert!((config.min_qual - 10.0).abs() < f64::EPSILON);
        assert!((config.min_gq - 20.0).abs() < f64::EPSILON);
        assert!((config.strand_bias_threshold - 0.001).abs() < f64::EPSILON);
        assert!((config.prior_het - 0.001).abs() < f64::EPSILON);
        assert!((config.prior_hom_alt - 0.0005).abs() < f64::EPSILON);
    }

    // -- Obvious hom-alt --

    #[test]
    fn test_obvious_hom_alt() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAAA".to_vec());

        // 20 reads all with T at position 0 (ref A)
        let records: Vec<SamRecord> = (0..20)
            .map(|i| make_record(&format!("r{i}"), 0, "chr1", 1, 60, "1M", "T", "I"))
            .collect();

        let pileups = pileup(&records, Some(&reference)).unwrap();
        let config = lenient_config();
        let variants = call_variants(&pileups[0], &config);

        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].genotype, Genotype::HomAlt);
        assert!(variants[0].variant.quality.unwrap() > 30.0);
        assert_eq!(variants[0].allele_depth, [0, 20]);
    }

    // -- Obvious het --

    #[test]
    fn test_obvious_het() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAAA".to_vec());

        // 10 ref + 10 alt reads
        let mut records: Vec<SamRecord> = (0..10)
            .map(|i| make_record(&format!("rA{i}"), 0, "chr1", 1, 60, "1M", "A", "I"))
            .collect();
        records.extend((0..10).map(|i| {
            make_record(&format!("rT{i}"), 0, "chr1", 1, 60, "1M", "T", "I")
        }));

        let pileups = pileup(&records, Some(&reference)).unwrap();
        let config = lenient_config();
        let variants = call_variants(&pileups[0], &config);

        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].genotype, Genotype::Het);
        assert_eq!(variants[0].allele_depth, [10, 10]);
    }

    // -- Hom-ref emits nothing --

    #[test]
    fn test_hom_ref_no_variant() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAAA".to_vec());

        let records: Vec<SamRecord> = (0..20)
            .map(|i| make_record(&format!("r{i}"), 0, "chr1", 1, 60, "1M", "A", "I"))
            .collect();

        let pileups = pileup(&records, Some(&reference)).unwrap();
        let config = lenient_config();
        let variants = call_variants(&pileups[0], &config);

        assert!(variants.is_empty());
    }

    // -- Low depth filter --

    #[test]
    fn test_low_depth_filtered() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAAA".to_vec());

        // Only 2 reads, but config requires min_depth = 5
        let records = vec![
            make_record("r1", 0, "chr1", 1, 60, "1M", "T", "I"),
            make_record("r2", 0, "chr1", 1, 60, "1M", "T", "I"),
        ];

        let pileups = pileup(&records, Some(&reference)).unwrap();
        let config = VariantCallConfig::default(); // min_depth = 5
        let variants = call_variants(&pileups[0], &config);

        assert!(variants.is_empty());
    }

    // -- Strand bias --

    #[test]
    fn test_strand_bias_detected() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAAA".to_vec());

        // All alt reads on forward strand, all ref on reverse
        let mut records: Vec<SamRecord> = (0..10)
            .map(|i| make_record(&format!("rF{i}"), 0, "chr1", 1, 60, "1M", "T", "I"))
            .collect();
        records.extend((0..10).map(|i| {
            make_record(&format!("rR{i}"), 16, "chr1", 1, 60, "1M", "A", "I")
        }));

        let pileups = pileup(&records, Some(&reference)).unwrap();
        let mut config = lenient_config();
        config.strand_bias_threshold = 0.05; // reasonable threshold
        let variants = call_variants(&pileups[0], &config);

        assert_eq!(variants.len(), 1);
        // All alt forward, all ref reverse → should trigger strand bias
        assert!(variants[0].strand_bias_p < 0.05);
        assert!(matches!(variants[0].variant.filter, VariantFilter::Fail(_)));
    }

    // -- No strand bias → PASS --

    #[test]
    fn test_no_strand_bias_passes() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAAA".to_vec());

        // Mix of forward and reverse for both ref and alt
        let mut records: Vec<SamRecord> = Vec::new();
        for i in 0..5 {
            records.push(make_record(&format!("rAf{i}"), 0, "chr1", 1, 60, "1M", "A", "I"));
            records.push(make_record(&format!("rAr{i}"), 16, "chr1", 1, 60, "1M", "A", "I"));
            records.push(make_record(&format!("rTf{i}"), 0, "chr1", 1, 60, "1M", "T", "I"));
            records.push(make_record(&format!("rTr{i}"), 16, "chr1", 1, 60, "1M", "T", "I"));
        }

        let pileups = pileup(&records, Some(&reference)).unwrap();
        let config = lenient_config();
        let variants = call_variants(&pileups[0], &config);

        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].variant.filter, VariantFilter::Pass);
    }

    // -- Quality varies with base quality --

    #[test]
    fn test_quality_varies_with_base_quality() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAAA".to_vec());

        // High-quality alt reads (Phred ~40: 'I')
        let records_hq: Vec<SamRecord> = (0..20)
            .map(|i| make_record(&format!("r{i}"), 0, "chr1", 1, 60, "1M", "T", "I"))
            .collect();

        // Low-quality alt reads (Phred ~2: '#')
        let records_lq: Vec<SamRecord> = (0..20)
            .map(|i| make_record(&format!("r{i}"), 0, "chr1", 1, 60, "1M", "T", "#"))
            .collect();

        let pileups_hq = pileup(&records_hq, Some(&reference)).unwrap();
        let pileups_lq = pileup(&records_lq, Some(&reference)).unwrap();
        let config = lenient_config();

        let vars_hq = call_variants(&pileups_hq[0], &config);
        let vars_lq = call_variants(&pileups_lq[0], &config);

        assert!(!vars_hq.is_empty());
        assert!(!vars_lq.is_empty());
        assert!(vars_hq[0].variant.quality.unwrap() > vars_lq[0].variant.quality.unwrap());
    }

    // -- Deletion calling --

    #[test]
    fn test_deletion_called() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"ACGTACGT".to_vec());

        // Multiple reads with a deletion: 2M2D2M
        let records: Vec<SamRecord> = (0..10)
            .map(|i| {
                make_record(
                    &format!("r{i}"),
                    0,
                    "chr1",
                    1,
                    60,
                    "2M2D2M",
                    "ACGT",
                    "IIII",
                )
            })
            .collect();

        let pileups = pileup(&records, Some(&reference)).unwrap();
        let config = lenient_config();
        let variants = call_variants(&pileups[0], &config);

        // Should find the deletion
        let dels: Vec<_> = variants.iter().filter(|v| v.variant.is_indel()).collect();
        assert!(!dels.is_empty());
        // The deletion REF should be longer than ALT
        assert!(dels[0].variant.ref_allele.len() > dels[0].variant.alt_alleles[0].len());
    }

    // -- Insertion calling --

    #[test]
    fn test_insertion_called() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"ACGTACGT".to_vec());

        // Multiple reads with an insertion: 3M2I3M
        let records: Vec<SamRecord> = (0..10)
            .map(|i| {
                make_record(
                    &format!("r{i}"),
                    0,
                    "chr1",
                    1,
                    60,
                    "3M2I3M",
                    "ACGTTAAA",
                    "IIIIIIII",
                )
            })
            .collect();

        let pileups = pileup(&records, Some(&reference)).unwrap();
        let mut config = lenient_config();
        config.min_alt_count = 2;
        let variants = call_variants(&pileups[0], &config);

        // Should find the insertion
        let ins: Vec<_> = variants
            .iter()
            .filter(|v| v.variant.ref_allele.len() < v.variant.alt_alleles[0].len())
            .collect();
        assert!(!ins.is_empty());
        // ALT should be anchor + inserted bases
        assert!(ins[0].variant.alt_alleles[0].len() > ins[0].variant.ref_allele.len());
    }

    // -- Genotype likelihood math --

    #[test]
    fn test_genotype_likelihoods_all_ref() {
        let config = lenient_config();
        let bases = vec![b'A'; 10];
        let quals = vec![30u8; 10];
        let mapqs = vec![60u8; 10];

        let lls = genotype_likelihoods_snp(&bases, &quals, &mapqs, b'A', b'T', &config);

        // RR should have the highest likelihood
        assert!(lls[0] > lls[1]);
        assert!(lls[0] > lls[2]);
    }

    #[test]
    fn test_genotype_likelihoods_all_alt() {
        let config = lenient_config();
        let bases = vec![b'T'; 10];
        let quals = vec![30u8; 10];
        let mapqs = vec![60u8; 10];

        let lls = genotype_likelihoods_snp(&bases, &quals, &mapqs, b'A', b'T', &config);

        // AA should have the highest likelihood
        assert!(lls[2] > lls[1]);
        assert!(lls[2] > lls[0]);
    }

    #[test]
    fn test_genotype_likelihoods_het() {
        let config = lenient_config();
        let mut bases = vec![b'A'; 10];
        bases.extend(vec![b'T'; 10]);
        let quals = vec![30u8; 20];
        let mapqs = vec![60u8; 20];

        let lls = genotype_likelihoods_snp(&bases, &quals, &mapqs, b'A', b'T', &config);

        // RA should have the highest likelihood for 50/50 mix
        assert!(lls[1] > lls[0]);
        assert!(lls[1] > lls[2]);
    }

    // -- QUAL/GQ/PL correctness --

    #[test]
    fn test_qual_gq_pl() {
        let config = lenient_config();
        let bases = vec![b'T'; 20];
        let quals = vec![30u8; 20];
        let mapqs = vec![60u8; 20];

        let lls = genotype_likelihoods_snp(&bases, &quals, &mapqs, b'A', b'T', &config);
        let (qual, gq, genotype, pl) = call_genotype(lls, &config);

        assert_eq!(genotype, Genotype::HomAlt);
        assert!(qual > 50.0); // Strong evidence against RR
        assert!(gq > 10.0); // Good genotype quality
        assert_eq!(pl[2], 0); // Best genotype (HomAlt) should have PL = 0
        assert!(pl[0] > 0); // HomRef should have positive PL
    }

    // -- variant_call_stats --

    #[test]
    fn test_variant_call_stats_ti_tv() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AACCCC".to_vec());

        // Transition: A→G
        let mut records: Vec<SamRecord> = (0..20)
            .map(|i| make_record(&format!("r{i}"), 0, "chr1", 1, 60, "1M", "G", "I"))
            .collect();
        // Transversion: C→A at pos 3
        records.extend((0..20).map(|i| {
            make_record(&format!("s{i}"), 0, "chr1", 3, 60, "1M", "A", "I")
        }));

        let pileups = pileup(&records, Some(&reference)).unwrap();
        let config = lenient_config();
        let variants = call_variants(&pileups[0], &config);

        let stats = variant_call_stats(&variants);
        assert_eq!(stats.snvs, 2);
        assert!((stats.ti_tv_ratio - 1.0).abs() < f64::EPSILON); // 1 Ti, 1 Tv
    }

    #[test]
    fn test_variant_call_stats_counts() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAAA".to_vec());

        // het variant
        let mut records: Vec<SamRecord> = (0..10)
            .map(|i| make_record(&format!("rA{i}"), 0, "chr1", 1, 60, "1M", "A", "I"))
            .collect();
        records.extend((0..10).map(|i| {
            make_record(&format!("rT{i}"), 0, "chr1", 1, 60, "1M", "T", "I")
        }));

        let pileups = pileup(&records, Some(&reference)).unwrap();
        let config = lenient_config();
        let variants = call_variants(&pileups[0], &config);

        let stats = variant_call_stats(&variants);
        assert_eq!(stats.total, 1);
        assert_eq!(stats.het_count, 1);
        assert_eq!(stats.hom_alt_count, 0);
    }

    // -- Empty pileup --

    #[test]
    fn test_empty_pileup() {
        let p = Pileup {
            rname: "chr1".to_string(),
            columns: vec![],
        };
        let config = VariantCallConfig::default();
        let variants = call_variants(&p, &config);
        assert!(variants.is_empty());
    }

    // -- Multiple references --

    #[test]
    fn test_call_variants_all() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAAA".to_vec());
        reference.insert("chr2".to_string(), b"CCCC".to_vec());

        let records = vec![
            make_record("r1", 0, "chr1", 1, 60, "1M", "T", "I"),
            make_record("r2", 0, "chr1", 1, 60, "1M", "T", "I"),
            make_record("r3", 0, "chr1", 1, 60, "1M", "T", "I"),
            make_record("r4", 0, "chr1", 1, 60, "1M", "T", "I"),
            make_record("r5", 0, "chr1", 1, 60, "1M", "T", "I"),
            make_record("r6", 0, "chr2", 1, 60, "1M", "A", "I"),
            make_record("r7", 0, "chr2", 1, 60, "1M", "A", "I"),
            make_record("r8", 0, "chr2", 1, 60, "1M", "A", "I"),
            make_record("r9", 0, "chr2", 1, 60, "1M", "A", "I"),
            make_record("r10", 0, "chr2", 1, 60, "1M", "A", "I"),
        ];

        let pileups = pileup(&records, Some(&reference)).unwrap();
        let config = lenient_config();
        let variants = call_variants_all(&pileups, &config);

        assert!(variants.len() >= 2);
        // Should have variants from both chromosomes
        let chroms: Vec<_> = variants.iter().map(|v| v.variant.chrom.as_str()).collect();
        assert!(chroms.contains(&"chr1"));
        assert!(chroms.contains(&"chr2"));
    }

    // -- Multi-allelic picks best alt --

    #[test]
    fn test_multi_allelic_picks_best() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAAA".to_vec());

        // 15 T reads, 5 G reads (T should be picked as best alt)
        let mut records: Vec<SamRecord> = (0..15)
            .map(|i| make_record(&format!("rT{i}"), 0, "chr1", 1, 60, "1M", "T", "I"))
            .collect();
        records.extend((0..5).map(|i| {
            make_record(&format!("rG{i}"), 0, "chr1", 1, 60, "1M", "G", "I")
        }));

        let pileups = pileup(&records, Some(&reference)).unwrap();
        let config = lenient_config();
        let variants = call_variants(&pileups[0], &config);

        assert!(!variants.is_empty());
        assert_eq!(variants[0].variant.alt_alleles[0], vec![b'T']);
    }

    // -- ln_sum_exp --

    #[test]
    fn test_ln_sum_exp() {
        // ln(e^0 + e^0) = ln(2)
        let result = ln_sum_exp(0.0, 0.0);
        assert!((result - 2.0_f64.ln()).abs() < 1e-10);

        // When one is -inf, returns the other
        assert!((ln_sum_exp(f64::NEG_INFINITY, 5.0) - 5.0).abs() < 1e-10);
        assert!((ln_sum_exp(5.0, f64::NEG_INFINITY) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_ln_sum_exp3() {
        // ln(e^0 + e^0 + e^0) = ln(3)
        let result = ln_sum_exp3(0.0, 0.0, 0.0);
        assert!((result - 3.0_f64.ln()).abs() < 1e-10);
    }
}
