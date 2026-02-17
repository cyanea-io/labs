//! Gene set enrichment and over-representation analysis.
//!
//! Provides two standard methods for determining whether predefined sets of
//! genes show statistically significant enrichment:
//!
//! - **Over-representation analysis** ([`ora`]) — tests whether significant
//!   genes are enriched in gene sets using the hypergeometric distribution.
//! - **GSEA preranked** ([`gsea_preranked`]) — walks a pre-ranked gene list and
//!   computes enrichment scores with permutation-based significance
//!   (Subramanian et al. 2005).
//!
//! Both methods apply Benjamini-Hochberg correction via
//! [`crate::correction::benjamini_hochberg`].

use std::collections::HashSet;

use cyanea_core::{CyaneaError, Result};

use crate::correction;
use crate::testing::ln_choose;

// ── Gene set representation ─────────────────────────────────────────────────

/// A named set of gene indices.
#[derive(Debug, Clone)]
pub struct GeneSet {
    /// Name of the gene set (e.g., pathway name).
    pub name: String,
    /// Gene indices (0-based) belonging to this set.
    pub genes: Vec<usize>,
}

// ── ORA types ───────────────────────────────────────────────────────────────

/// Result of over-representation analysis for one gene set.
#[derive(Debug, Clone)]
pub struct OraResult {
    /// Name of the gene set.
    pub gene_set: String,
    /// Number of significant genes found in this set (k).
    pub overlap: usize,
    /// Expected overlap under the null hypothesis.
    pub expected: f64,
    /// Effective gene set size within the universe (K).
    pub gene_set_size: usize,
    /// Hypergeometric upper-tail p-value: P(X >= k).
    pub p_value: f64,
    /// Benjamini-Hochberg adjusted p-value.
    pub p_adjusted: f64,
}

/// Over-representation analysis using the hypergeometric test.
///
/// For each gene set, tests whether `significant` genes are over-represented
/// using P(X >= k) where X ~ Hypergeometric(N, K, n).
///
/// - `significant`: indices of significant genes.
/// - `gene_sets`: named gene sets to test.
/// - `n_total`: size of the gene universe (N). All gene indices must be in
///   `[0, n_total)`.
///
/// Returns results sorted by ascending p-value, with BH-corrected p-values.
///
/// # Example
///
/// ```
/// use cyanea_stats::enrichment::{GeneSet, ora};
///
/// let significant = vec![0, 1, 2, 3, 4];
/// let gene_sets = vec![
///     GeneSet { name: "pathway_A".into(), genes: vec![0, 1, 2, 10, 11] },
///     GeneSet { name: "pathway_B".into(), genes: vec![50, 51, 52, 53, 54] },
/// ];
/// let results = ora(&significant, &gene_sets, 100).unwrap();
/// assert!(results[0].gene_set == "pathway_A");
/// assert!(results[0].overlap == 3);
/// ```
pub fn ora(
    significant: &[usize],
    gene_sets: &[GeneSet],
    n_total: usize,
) -> Result<Vec<OraResult>> {
    // Validate inputs
    if n_total == 0 {
        return Err(CyaneaError::InvalidInput(
            "ora: n_total must be > 0".into(),
        ));
    }
    if gene_sets.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "ora: gene_sets must be non-empty".into(),
        ));
    }
    for &idx in significant {
        if idx >= n_total {
            return Err(CyaneaError::InvalidInput(format!(
                "ora: significant gene index {} >= n_total {}",
                idx, n_total,
            )));
        }
    }
    for gs in gene_sets {
        for &idx in &gs.genes {
            if idx >= n_total {
                return Err(CyaneaError::InvalidInput(format!(
                    "ora: gene set '{}' contains index {} >= n_total {}",
                    gs.name, idx, n_total,
                )));
            }
        }
    }

    let sig_set: HashSet<usize> = significant.iter().copied().collect();
    let n = sig_set.len(); // sample size (number of significant genes)

    let mut results = Vec::with_capacity(gene_sets.len());

    for gs in gene_sets {
        // Effective gene set size within the universe
        let gs_unique: HashSet<usize> = gs.genes.iter().copied().collect();
        let big_k = gs_unique.len(); // success population size

        // Overlap: significant genes in this set
        let k = sig_set.intersection(&gs_unique).count();

        // Expected overlap: E[X] = n * K / N
        let expected = n as f64 * big_k as f64 / n_total as f64;

        // Hypergeometric upper-tail: P(X >= k) = Σ_{i=k}^{min(n,K)} P(X = i)
        let p_value = hypergeometric_upper_tail(k, n, big_k, n_total);

        results.push(OraResult {
            gene_set: gs.name.clone(),
            overlap: k,
            expected,
            gene_set_size: big_k,
            p_value,
            p_adjusted: 1.0, // filled in below
        });
    }

    // BH correction
    let raw_p: Vec<f64> = results.iter().map(|r| r.p_value).collect();
    let adj_p = correction::benjamini_hochberg(&raw_p)?;
    for (r, &padj) in results.iter_mut().zip(adj_p.iter()) {
        r.p_adjusted = padj;
    }

    // Sort by p-value ascending
    results.sort_by(|a, b| a.p_value.total_cmp(&b.p_value));

    Ok(results)
}

/// Hypergeometric upper-tail probability: P(X >= k).
///
/// X ~ Hypergeometric(N, K, n) where:
/// - N = total population
/// - K = success states in population
/// - n = draws (sample size)
/// - k = observed successes
fn hypergeometric_upper_tail(k: usize, n: usize, big_k: usize, big_n: usize) -> f64 {
    if k == 0 {
        return 1.0;
    }
    let max_i = n.min(big_k);
    if k > max_i {
        return 0.0;
    }

    // Sum PMF from k to min(n, K) in log-space for stability
    let mut sum = 0.0_f64;
    let log_denom = ln_choose(big_n, n);
    for i in k..=max_i {
        // Check that n - i <= N - K (otherwise PMF is 0)
        if n < i || big_n - big_k < n - i {
            continue;
        }
        let log_p = ln_choose(big_k, i) + ln_choose(big_n - big_k, n - i) - log_denom;
        sum += log_p.exp();
    }
    sum.min(1.0)
}

// ── GSEA types ──────────────────────────────────────────────────────────────

/// Result of GSEA for one gene set.
#[derive(Debug, Clone)]
pub struct GseaResult {
    /// Name of the gene set.
    pub gene_set: String,
    /// Enrichment score (ES): maximum deviation of the running sum.
    pub enrichment_score: f64,
    /// Normalized enrichment score (NES): ES / mean(|ES_null|) for matching sign.
    pub normalized_es: f64,
    /// Permutation-based p-value.
    pub p_value: f64,
    /// Benjamini-Hochberg adjusted p-value.
    pub p_adjusted: f64,
    /// Number of genes in the leading edge.
    pub leading_edge_size: usize,
    /// Effective gene set size (genes present in the ranked list).
    pub gene_set_size: usize,
}

/// GSEA preranked: enrichment analysis on a pre-ranked gene list.
///
/// Implements the Subramanian et al. (2005) algorithm:
/// 1. Walk down the ranked list, incrementing the running sum at set members
///    (weighted by |score|^weight) and decrementing at non-members.
/// 2. ES = maximum absolute deviation of the running sum.
/// 3. Significance via gene-label permutation.
/// 4. NES = ES / mean(|ES_null|) for the same sign.
/// 5. BH correction across gene sets.
///
/// - `genes`: gene indices in rank order (best first).
/// - `scores`: ranking metric per gene (e.g., log2FC, signed statistic).
///   Must have the same length as `genes`.
/// - `gene_sets`: named gene sets to test.
/// - `weight`: exponent for score weighting (1.0 = classic GSEA, 0.0 = unweighted Kolmogorov-Smirnov).
/// - `n_permutations`: number of permutations for p-value estimation (e.g., 1000).
///
/// Returns results sorted by ascending p-value, with BH-corrected p-values.
///
/// # Example
///
/// ```
/// use cyanea_stats::enrichment::{GeneSet, gsea_preranked};
///
/// // Genes 0-4 at top of list, 5-9 at bottom
/// let genes: Vec<usize> = (0..20).collect();
/// let scores: Vec<f64> = (0..20).rev().map(|i| i as f64).collect();
/// let gene_sets = vec![
///     GeneSet { name: "top_genes".into(), genes: vec![0, 1, 2, 3, 4] },
/// ];
/// let results = gsea_preranked(&genes, &scores, &gene_sets, 1.0, 100).unwrap();
/// assert!(results[0].enrichment_score > 0.0);
/// ```
pub fn gsea_preranked(
    genes: &[usize],
    scores: &[f64],
    gene_sets: &[GeneSet],
    weight: f64,
    n_permutations: usize,
) -> Result<Vec<GseaResult>> {
    // Validate inputs
    if genes.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "gsea_preranked: genes must be non-empty".into(),
        ));
    }
    if genes.len() != scores.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "gsea_preranked: genes length ({}) != scores length ({})",
            genes.len(),
            scores.len(),
        )));
    }
    if gene_sets.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "gsea_preranked: gene_sets must be non-empty".into(),
        ));
    }
    if n_permutations == 0 {
        return Err(CyaneaError::InvalidInput(
            "gsea_preranked: n_permutations must be > 0".into(),
        ));
    }

    let gene_set_in_list: HashSet<usize> = genes.iter().copied().collect();

    let mut results = Vec::with_capacity(gene_sets.len());

    for gs in gene_sets {
        // Effective set: gene set members present in the ranked list
        let set_members: HashSet<usize> = gs.genes.iter().copied()
            .filter(|g| gene_set_in_list.contains(g))
            .collect();
        let n_h = set_members.len();

        if n_h == 0 {
            // No overlap with ranked list — no enrichment
            results.push(GseaResult {
                gene_set: gs.name.clone(),
                enrichment_score: 0.0,
                normalized_es: 0.0,
                p_value: 1.0,
                p_adjusted: 1.0,
                leading_edge_size: 0,
                gene_set_size: 0,
            });
            continue;
        }

        // Compute observed ES
        let (es, leading_edge) = compute_es(genes, scores, &set_members, n_h, weight);

        // Permutation null distribution
        let mut rng = Xorshift64::new(42);
        let mut null_es = Vec::with_capacity(n_permutations);

        // We permute score assignments: shuffle scores and recompute ES
        let mut perm_scores: Vec<f64> = scores.to_vec();
        for _ in 0..n_permutations {
            fisher_yates_shuffle(&mut perm_scores, &mut rng);
            let (perm_es, _) = compute_es(genes, &perm_scores, &set_members, n_h, weight);
            null_es.push(perm_es);
        }

        // p-value: fraction of null ES as or more extreme than observed
        // Separate positive and negative
        let p_value = if es >= 0.0 {
            let count = null_es.iter().filter(|&&x| x >= es).count();
            (count as f64 / n_permutations as f64).max(1.0 / n_permutations as f64)
        } else {
            let count = null_es.iter().filter(|&&x| x <= es).count();
            (count as f64 / n_permutations as f64).max(1.0 / n_permutations as f64)
        };

        // NES: normalize by mean of |ES_null| for same sign
        let normalized_es = if es >= 0.0 {
            let pos_null: Vec<f64> = null_es.iter().copied().filter(|&x| x >= 0.0).collect();
            if pos_null.is_empty() {
                0.0
            } else {
                let mean_pos: f64 = pos_null.iter().sum::<f64>() / pos_null.len() as f64;
                if mean_pos > 1e-15 { es / mean_pos } else { 0.0 }
            }
        } else {
            let neg_null: Vec<f64> = null_es.iter().copied().filter(|&x| x < 0.0).collect();
            if neg_null.is_empty() {
                0.0
            } else {
                let mean_neg: f64 = neg_null.iter().map(|x| x.abs()).sum::<f64>() / neg_null.len() as f64;
                if mean_neg > 1e-15 { -(es.abs() / mean_neg) } else { 0.0 }
            }
        };

        results.push(GseaResult {
            gene_set: gs.name.clone(),
            enrichment_score: es,
            normalized_es,
            p_value,
            p_adjusted: 1.0, // filled in below
            leading_edge_size: leading_edge,
            gene_set_size: n_h,
        });
    }

    // BH correction
    let raw_p: Vec<f64> = results.iter().map(|r| r.p_value).collect();
    let adj_p = correction::benjamini_hochberg(&raw_p)?;
    for (r, &padj) in results.iter_mut().zip(adj_p.iter()) {
        r.p_adjusted = padj;
    }

    // Sort by p-value ascending
    results.sort_by(|a, b| a.p_value.total_cmp(&b.p_value));

    Ok(results)
}

/// Compute the enrichment score for one gene set.
///
/// Returns (ES, leading_edge_size).
fn compute_es(
    genes: &[usize],
    scores: &[f64],
    set_members: &HashSet<usize>,
    n_h: usize,
    weight: f64,
) -> (f64, usize) {
    let n = genes.len();

    // N_R: sum of |score|^weight for set members
    let n_r: f64 = genes.iter().zip(scores.iter())
        .filter(|(&g, _)| set_members.contains(&g))
        .map(|(_, &s)| s.abs().powf(weight))
        .sum();

    // Miss penalty per non-member step
    let n_miss = n - n_h;
    let miss_penalty = if n_miss > 0 { 1.0 / n_miss as f64 } else { 0.0 };

    let mut running_sum = 0.0_f64;
    let mut max_dev = 0.0_f64;
    let mut max_dev_pos = 0usize; // position of max |deviation|
    let mut max_dev_signed = 0.0_f64;

    for (i, (&g, &s)) in genes.iter().zip(scores.iter()).enumerate() {
        if set_members.contains(&g) {
            let hit_score = if n_r > 1e-15 {
                s.abs().powf(weight) / n_r
            } else {
                1.0 / n_h as f64
            };
            running_sum += hit_score;
        } else {
            running_sum -= miss_penalty;
        }

        if running_sum.abs() > max_dev {
            max_dev = running_sum.abs();
            max_dev_signed = running_sum;
            max_dev_pos = i;
        }
    }

    // Leading edge: set members before (and including) the peak
    let leading_edge = genes[..=max_dev_pos].iter()
        .filter(|g| set_members.contains(g))
        .count();

    (max_dev_signed, leading_edge)
}

// ── Simple xorshift64 PRNG ──────────────────────────────────────────────────

/// Minimal xorshift64 PRNG for deterministic permutations.
struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        // Ensure non-zero state
        Self { state: if seed == 0 { 1 } else { seed } }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    /// Generate a random index in [0, n).
    fn next_usize(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

/// Fisher-Yates shuffle using our xorshift PRNG.
fn fisher_yates_shuffle(slice: &mut [f64], rng: &mut Xorshift64) {
    let n = slice.len();
    for i in (1..n).rev() {
        let j = rng.next_usize(i + 1);
        slice.swap(i, j);
    }
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ORA tests ──────────────────────────────────────────────────────

    #[test]
    fn ora_basic_enrichment() {
        // 5 significant genes out of 100, gene set has 3 of them
        let significant = vec![0, 1, 2, 3, 4];
        let gene_sets = vec![
            GeneSet { name: "enriched".into(), genes: vec![0, 1, 2, 50, 51] },
        ];
        let results = ora(&significant, &gene_sets, 100).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].overlap, 3);
        assert_eq!(results[0].gene_set_size, 5);
        // 3 out of 5 significant genes in a set of 5/100 is very enriched
        assert!(results[0].p_value < 0.01, "p={}", results[0].p_value);
    }

    #[test]
    fn ora_no_overlap() {
        let significant = vec![0, 1, 2];
        let gene_sets = vec![
            GeneSet { name: "disjoint".into(), genes: vec![50, 51, 52] },
        ];
        let results = ora(&significant, &gene_sets, 100).unwrap();
        assert_eq!(results[0].overlap, 0);
        assert!((results[0].p_value - 1.0).abs() < 1e-10, "p={}", results[0].p_value);
    }

    #[test]
    fn ora_complete_overlap() {
        // All significant genes are in the set
        let significant = vec![0, 1, 2, 3, 4];
        let gene_sets = vec![
            GeneSet { name: "complete".into(), genes: vec![0, 1, 2, 3, 4] },
        ];
        let results = ora(&significant, &gene_sets, 1000).unwrap();
        assert_eq!(results[0].overlap, 5);
        assert!(results[0].p_value < 0.001, "p={}", results[0].p_value);
    }

    #[test]
    fn ora_empty_significant() {
        let significant: Vec<usize> = vec![];
        let gene_sets = vec![
            GeneSet { name: "any".into(), genes: vec![0, 1, 2] },
        ];
        let results = ora(&significant, &gene_sets, 100).unwrap();
        assert_eq!(results[0].overlap, 0);
        assert!((results[0].p_value - 1.0).abs() < 1e-10);
    }

    #[test]
    fn ora_expected_overlap() {
        // E[k] = n * K / N
        let significant: Vec<usize> = (0..10).collect();
        let gene_sets = vec![
            GeneSet { name: "gs".into(), genes: (0..20).collect() },
        ];
        let results = ora(&significant, &gene_sets, 100).unwrap();
        // Expected = 10 * 20 / 100 = 2.0
        assert!((results[0].expected - 2.0).abs() < 1e-10, "expected={}", results[0].expected);
    }

    #[test]
    fn ora_bh_correction_applied() {
        // Multiple gene sets: check that p_adjusted >= p_value
        let significant: Vec<usize> = (0..5).collect();
        let gene_sets = vec![
            GeneSet { name: "gs1".into(), genes: vec![0, 1, 2, 50, 51] },
            GeneSet { name: "gs2".into(), genes: vec![60, 61, 62, 63, 64] },
            GeneSet { name: "gs3".into(), genes: vec![0, 1, 70, 71, 72] },
        ];
        let results = ora(&significant, &gene_sets, 100).unwrap();
        for r in &results {
            assert!(
                r.p_adjusted >= r.p_value - 1e-15,
                "{}: padj={} < p={}",
                r.gene_set, r.p_adjusted, r.p_value,
            );
        }
    }

    #[test]
    fn ora_sorted_by_pvalue() {
        let significant: Vec<usize> = (0..5).collect();
        let gene_sets = vec![
            GeneSet { name: "gs1".into(), genes: vec![0, 1, 2, 50, 51] },
            GeneSet { name: "gs2".into(), genes: vec![60, 61, 62, 63, 64] },
            GeneSet { name: "gs3".into(), genes: vec![0, 1, 2, 3, 4] },
        ];
        let results = ora(&significant, &gene_sets, 100).unwrap();
        for w in results.windows(2) {
            assert!(w[0].p_value <= w[1].p_value + 1e-15);
        }
    }

    #[test]
    fn ora_error_index_too_large() {
        let significant = vec![100]; // >= n_total
        let gene_sets = vec![
            GeneSet { name: "gs".into(), genes: vec![0] },
        ];
        assert!(ora(&significant, &gene_sets, 100).is_err());
    }

    #[test]
    fn ora_error_gene_set_index_too_large() {
        let significant = vec![0];
        let gene_sets = vec![
            GeneSet { name: "gs".into(), genes: vec![200] },
        ];
        assert!(ora(&significant, &gene_sets, 100).is_err());
    }

    #[test]
    fn ora_error_empty_gene_sets() {
        assert!(ora(&[0], &[], 100).is_err());
    }

    #[test]
    fn ora_error_zero_universe() {
        let gene_sets = vec![GeneSet { name: "gs".into(), genes: vec![] }];
        assert!(ora(&[], &gene_sets, 0).is_err());
    }

    #[test]
    fn ora_duplicate_genes_deduplicated() {
        // Duplicate indices in significant and gene sets should be deduplicated
        let significant = vec![0, 0, 1, 1, 2];
        let gene_sets = vec![
            GeneSet { name: "gs".into(), genes: vec![0, 0, 1, 50, 50] },
        ];
        let results = ora(&significant, &gene_sets, 100).unwrap();
        // Unique significant: {0, 1, 2}, unique gs: {0, 1, 50}
        assert_eq!(results[0].overlap, 2); // {0, 1}
        assert_eq!(results[0].gene_set_size, 3); // {0, 1, 50}
    }

    // ── GSEA tests ─────────────────────────────────────────────────────

    #[test]
    fn gsea_genes_at_top() {
        // Gene set members concentrated at top of ranked list → positive ES
        let genes: Vec<usize> = (0..50).collect();
        let scores: Vec<f64> = (0..50).rev().map(|i| i as f64).collect();
        let gene_sets = vec![
            GeneSet { name: "top".into(), genes: vec![0, 1, 2, 3, 4] },
        ];
        let results = gsea_preranked(&genes, &scores, &gene_sets, 1.0, 200).unwrap();
        assert!(results[0].enrichment_score > 0.0, "ES={}", results[0].enrichment_score);
    }

    #[test]
    fn gsea_genes_at_bottom() {
        // Gene set members concentrated at bottom → negative ES
        let genes: Vec<usize> = (0..50).collect();
        let scores: Vec<f64> = (0..50).rev().map(|i| i as f64).collect();
        let gene_sets = vec![
            GeneSet { name: "bottom".into(), genes: vec![45, 46, 47, 48, 49] },
        ];
        let results = gsea_preranked(&genes, &scores, &gene_sets, 1.0, 200).unwrap();
        assert!(results[0].enrichment_score < 0.0, "ES={}", results[0].enrichment_score);
    }

    #[test]
    fn gsea_uniform_distribution() {
        // Gene set members evenly distributed → ES near 0, not significant
        let genes: Vec<usize> = (0..100).collect();
        let scores: Vec<f64> = (0..100).rev().map(|i| i as f64).collect();
        let gene_sets = vec![
            GeneSet { name: "uniform".into(), genes: vec![0, 20, 40, 60, 80] },
        ];
        let results = gsea_preranked(&genes, &scores, &gene_sets, 1.0, 200).unwrap();
        // ES should be relatively small
        assert!(results[0].enrichment_score.abs() < 0.5, "ES={}", results[0].enrichment_score);
    }

    #[test]
    fn gsea_unweighted() {
        // weight=0 gives unweighted KS-like statistic
        let genes: Vec<usize> = (0..20).collect();
        let scores: Vec<f64> = (0..20).rev().map(|i| i as f64).collect();
        let gene_sets = vec![
            GeneSet { name: "top".into(), genes: vec![0, 1, 2, 3, 4] },
        ];
        let results = gsea_preranked(&genes, &scores, &gene_sets, 0.0, 100).unwrap();
        assert!(results[0].enrichment_score > 0.0);
    }

    #[test]
    fn gsea_weight_one_classic() {
        // weight=1.0 is classic GSEA
        let genes: Vec<usize> = (0..20).collect();
        let scores: Vec<f64> = (0..20).rev().map(|i| i as f64 + 1.0).collect();
        let gene_sets = vec![
            GeneSet { name: "top".into(), genes: vec![0, 1, 2, 3, 4] },
        ];
        let results = gsea_preranked(&genes, &scores, &gene_sets, 1.0, 100).unwrap();
        assert!(results[0].enrichment_score > 0.0);
    }

    #[test]
    fn gsea_nes_sign_matches_es() {
        let genes: Vec<usize> = (0..50).collect();
        let scores: Vec<f64> = (0..50).rev().map(|i| i as f64).collect();
        let gene_sets = vec![
            GeneSet { name: "top".into(), genes: vec![0, 1, 2, 3, 4] },
            GeneSet { name: "bottom".into(), genes: vec![45, 46, 47, 48, 49] },
        ];
        let results = gsea_preranked(&genes, &scores, &gene_sets, 1.0, 200).unwrap();
        for r in &results {
            if r.enrichment_score > 0.0 {
                assert!(r.normalized_es >= 0.0, "NES={} but ES={}", r.normalized_es, r.enrichment_score);
            } else if r.enrichment_score < 0.0 {
                assert!(r.normalized_es <= 0.0, "NES={} but ES={}", r.normalized_es, r.enrichment_score);
            }
        }
    }

    #[test]
    fn gsea_leading_edge_size() {
        let genes: Vec<usize> = (0..20).collect();
        let scores: Vec<f64> = (0..20).rev().map(|i| i as f64).collect();
        let gene_sets = vec![
            GeneSet { name: "top".into(), genes: vec![0, 1, 2, 3, 4] },
        ];
        let results = gsea_preranked(&genes, &scores, &gene_sets, 1.0, 100).unwrap();
        // Leading edge should be between 1 and gene_set_size
        assert!(results[0].leading_edge_size >= 1);
        assert!(results[0].leading_edge_size <= results[0].gene_set_size);
    }

    #[test]
    fn gsea_pvalue_in_range() {
        let genes: Vec<usize> = (0..30).collect();
        let scores: Vec<f64> = (0..30).rev().map(|i| i as f64).collect();
        let gene_sets = vec![
            GeneSet { name: "gs".into(), genes: vec![0, 1, 2, 15, 16] },
        ];
        let results = gsea_preranked(&genes, &scores, &gene_sets, 1.0, 100).unwrap();
        assert!(results[0].p_value >= 0.0 && results[0].p_value <= 1.0,
            "p={}", results[0].p_value);
        assert!(results[0].p_adjusted >= 0.0 && results[0].p_adjusted <= 1.0,
            "padj={}", results[0].p_adjusted);
    }

    #[test]
    fn gsea_deterministic_with_fixed_seed() {
        let genes: Vec<usize> = (0..30).collect();
        let scores: Vec<f64> = (0..30).rev().map(|i| i as f64).collect();
        let gene_sets = vec![
            GeneSet { name: "gs".into(), genes: vec![0, 1, 2, 3, 4] },
        ];
        let r1 = gsea_preranked(&genes, &scores, &gene_sets, 1.0, 100).unwrap();
        let r2 = gsea_preranked(&genes, &scores, &gene_sets, 1.0, 100).unwrap();
        assert!((r1[0].enrichment_score - r2[0].enrichment_score).abs() < 1e-15);
        assert!((r1[0].p_value - r2[0].p_value).abs() < 1e-15);
    }

    #[test]
    fn gsea_no_overlap_with_list() {
        // Gene set has no members in the ranked list
        let genes: Vec<usize> = (0..10).collect();
        let scores: Vec<f64> = (0..10).rev().map(|i| i as f64).collect();
        let gene_sets = vec![
            GeneSet { name: "outside".into(), genes: vec![100, 101, 102] },
        ];
        let results = gsea_preranked(&genes, &scores, &gene_sets, 1.0, 100).unwrap();
        assert_eq!(results[0].gene_set_size, 0);
        assert!((results[0].enrichment_score - 0.0).abs() < 1e-15);
        assert!((results[0].p_value - 1.0).abs() < 1e-15);
    }

    #[test]
    fn gsea_single_gene_set_single_gene() {
        let genes: Vec<usize> = (0..10).collect();
        let scores: Vec<f64> = (0..10).rev().map(|i| i as f64 + 1.0).collect();
        let gene_sets = vec![
            GeneSet { name: "single".into(), genes: vec![0] },
        ];
        let results = gsea_preranked(&genes, &scores, &gene_sets, 1.0, 100).unwrap();
        assert_eq!(results[0].gene_set_size, 1);
        assert!(results[0].enrichment_score > 0.0);
    }

    #[test]
    fn gsea_error_empty_genes() {
        let gene_sets = vec![GeneSet { name: "gs".into(), genes: vec![0] }];
        assert!(gsea_preranked(&[], &[], &gene_sets, 1.0, 100).is_err());
    }

    #[test]
    fn gsea_error_mismatched_lengths() {
        let gene_sets = vec![GeneSet { name: "gs".into(), genes: vec![0] }];
        assert!(gsea_preranked(&[0, 1], &[1.0], &gene_sets, 1.0, 100).is_err());
    }

    #[test]
    fn gsea_error_empty_gene_sets() {
        assert!(gsea_preranked(&[0], &[1.0], &[], 1.0, 100).is_err());
    }

    #[test]
    fn gsea_error_zero_permutations() {
        let gene_sets = vec![GeneSet { name: "gs".into(), genes: vec![0] }];
        assert!(gsea_preranked(&[0], &[1.0], &gene_sets, 1.0, 0).is_err());
    }

    #[test]
    fn gsea_multiple_sets_bh_correction() {
        let genes: Vec<usize> = (0..50).collect();
        let scores: Vec<f64> = (0..50).rev().map(|i| i as f64).collect();
        let gene_sets = vec![
            GeneSet { name: "gs1".into(), genes: vec![0, 1, 2, 3, 4] },
            GeneSet { name: "gs2".into(), genes: vec![25, 26, 27, 28, 29] },
            GeneSet { name: "gs3".into(), genes: vec![45, 46, 47, 48, 49] },
        ];
        let results = gsea_preranked(&genes, &scores, &gene_sets, 1.0, 200).unwrap();
        // p_adjusted >= p_value for all
        for r in &results {
            assert!(
                r.p_adjusted >= r.p_value - 1e-15,
                "{}: padj={} < p={}",
                r.gene_set, r.p_adjusted, r.p_value,
            );
        }
        // Sorted by p-value
        for w in results.windows(2) {
            assert!(w[0].p_value <= w[1].p_value + 1e-15);
        }
    }

    // ── Xorshift PRNG tests ────────────────────────────────────────────

    #[test]
    fn xorshift_deterministic() {
        let mut rng1 = Xorshift64::new(42);
        let mut rng2 = Xorshift64::new(42);
        for _ in 0..100 {
            assert_eq!(rng1.next_u64(), rng2.next_u64());
        }
    }

    #[test]
    fn xorshift_different_seeds() {
        let mut rng1 = Xorshift64::new(1);
        let mut rng2 = Xorshift64::new(2);
        // At least one of the first 10 values should differ
        let differs = (0..10).any(|_| rng1.next_u64() != rng2.next_u64());
        assert!(differs);
    }
}
