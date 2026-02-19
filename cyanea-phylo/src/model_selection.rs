//! Information-theoretic model selection: AIC, BIC, AICc, LRT, ModelFinder.
//!
//! Supports comparing substitution models via likelihood-based criteria.

use crate::generic_likelihood::generic_tree_likelihood;
use crate::models::GammaRates;
use crate::subst_model::SubstitutionModel;
use crate::tree::PhyloTree;
use cyanea_core::Result;

/// Results for a single model evaluation.
#[derive(Debug, Clone)]
pub struct ModelResult {
    pub name: String,
    pub log_likelihood: f64,
    pub n_params: usize,
    pub n_sites: usize,
    pub aic: f64,
    pub bic: f64,
    pub aicc: f64,
}

/// Results from comparing multiple models.
#[derive(Debug, Clone)]
pub struct ModelSelectionResult {
    pub models: Vec<ModelResult>,
    pub best_aic: String,
    pub best_bic: String,
}

/// Result of a likelihood ratio test.
#[derive(Debug, Clone)]
pub struct LrtResult {
    pub stat: f64,
    pub df: usize,
    pub p_value: f64,
    pub reject: bool,
}

/// Akaike Information Criterion: AIC = -2 ln(L) + 2k
pub fn aic(log_likelihood: f64, n_params: usize) -> f64 {
    -2.0 * log_likelihood + 2.0 * n_params as f64
}

/// Bayesian Information Criterion: BIC = -2 ln(L) + k ln(n)
pub fn bic(log_likelihood: f64, n_params: usize, n_sites: usize) -> f64 {
    -2.0 * log_likelihood + n_params as f64 * (n_sites as f64).ln()
}

/// Corrected AIC for small samples: AICc = AIC + 2k(k+1)/(n-k-1)
pub fn aicc(log_likelihood: f64, n_params: usize, n_sites: usize) -> f64 {
    let k = n_params as f64;
    let n = n_sites as f64;
    let base = aic(log_likelihood, n_params);
    if n > k + 1.0 {
        base + 2.0 * k * (k + 1.0) / (n - k - 1.0)
    } else {
        f64::INFINITY
    }
}

/// Likelihood ratio test between nested models.
///
/// Tests whether the more complex model (with higher likelihood) is
/// significantly better than the simpler model.
///
/// `ll_null`: log-likelihood of the null (simpler) model.
/// `ll_alt`: log-likelihood of the alternative (more complex) model.
/// `df`: difference in number of free parameters.
pub fn lrt(ll_null: f64, ll_alt: f64, df: usize) -> LrtResult {
    let stat = 2.0 * (ll_alt - ll_null);
    let stat = stat.max(0.0); // Guard against numerical issues

    // Chi-squared p-value: P(X > stat) where X ~ chi-squared(df)
    // Using regularized gamma: P(chi2 > x) = 1 - P(df/2, x/2)
    let p_value = 1.0 - chi2_cdf(stat, df);

    LrtResult {
        stat,
        df,
        p_value,
        reject: p_value < 0.05,
    }
}

/// Compare multiple substitution models and rank by AIC/BIC.
///
/// Each candidate model is evaluated on the given tree and sequences.
/// Models are provided as `(name, model)` pairs.
pub fn model_finder(
    tree: &PhyloTree,
    sequences: &[&[u8]],
    candidate_models: &[(&str, &dyn SubstitutionModel)],
    gamma: Option<&GammaRates>,
) -> Result<ModelSelectionResult> {
    let n_taxa = tree.leaf_count();
    let n_branches = if n_taxa > 2 { 2 * n_taxa - 3 } else { 1 };
    let gamma_params = if gamma.is_some() { 1 } else { 0 };

    let seq_len = if sequences.is_empty() {
        0
    } else {
        sequences[0].len()
    };

    let mut results = Vec::new();

    for &(name, model) in candidate_models {
        let ll = generic_tree_likelihood(tree, sequences, model, gamma)?;
        let k = model.n_free_params() + n_branches + gamma_params;

        results.push(ModelResult {
            name: name.to_string(),
            log_likelihood: ll,
            n_params: k,
            n_sites: seq_len,
            aic: aic(ll, k),
            bic: bic(ll, k, seq_len),
            aicc: aicc(ll, k, seq_len),
        });
    }

    // Sort by AIC.
    results.sort_by(|a, b| a.aic.partial_cmp(&b.aic).unwrap_or(std::cmp::Ordering::Equal));

    let best_aic = results
        .first()
        .map(|r| r.name.clone())
        .unwrap_or_default();

    // Find best BIC.
    let best_bic = results
        .iter()
        .min_by(|a, b| a.bic.partial_cmp(&b.bic).unwrap_or(std::cmp::Ordering::Equal))
        .map(|r| r.name.clone())
        .unwrap_or_default();

    Ok(ModelSelectionResult {
        models: results,
        best_aic,
        best_bic,
    })
}

/// Chi-squared CDF: P(X <= x) where X ~ chi-squared(df).
fn chi2_cdf(x: f64, df: usize) -> f64 {
    if df == 0 || x <= 0.0 {
        return 0.0;
    }
    let a = df as f64 / 2.0;
    let x_half = x / 2.0;
    let p = crate::models::gamma_regularized(a, x_half);
    // Clamp to [0, 1] for numerical robustness.
    p.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::subst_model::{Hky85Model, Jc69Model};

    #[test]
    fn aic_known_values() {
        // AIC = -2*(-100) + 2*5 = 210
        assert!((aic(-100.0, 5) - 210.0).abs() < 1e-10);
        // AIC = -2*(-50) + 2*2 = 104
        assert!((aic(-50.0, 2) - 104.0).abs() < 1e-10);
    }

    #[test]
    fn bic_known_values() {
        // BIC = -2*(-100) + 5*ln(100) = 200 + 5*4.605... ≈ 223.026
        let val = bic(-100.0, 5, 100);
        assert!((val - 223.0259).abs() < 0.01);
    }

    #[test]
    fn aicc_converges_to_aic_for_large_n() {
        let aic_val = aic(-100.0, 3);
        let aicc_val = aicc(-100.0, 3, 10000);
        assert!(
            (aic_val - aicc_val).abs() < 0.01,
            "AICc should approach AIC for large n: {} vs {}",
            aicc_val, aic_val
        );
    }

    #[test]
    fn lrt_nested_models() {
        // Simple nested model test: JC69 (0 params) vs HKY85 (4 params)
        let result = lrt(-120.0, -110.0, 4);
        assert!(result.stat > 0.0, "stat should be positive: {}", result.stat);
        // stat = 2*((-110) - (-120)) = 20, should be significant
        assert!((result.stat - 20.0).abs() < 1e-10, "stat should be 20: {}", result.stat);
        assert!((0.0..=1.0).contains(&result.p_value), "p_value {} out of [0,1]", result.p_value);
    }

    #[test]
    fn lrt_p_value_bounds() {
        let result = lrt(-100.0, -100.0, 1);
        assert!(result.p_value >= 0.0 && result.p_value <= 1.0);
        // stat = 0, p-value should be 1.0
        assert!(result.p_value > 0.9, "p = {} should be ~1.0", result.p_value);
    }

    #[test]
    fn model_finder_ranks_jc69_on_jc69_data() {
        let tree = PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let seqs: Vec<Vec<u8>> = vec![
            b"ACGTACGTACGT".to_vec(),
            b"ACGTACGTACGT".to_vec(),
            b"TGCATGCATGCA".to_vec(),
            b"TGCATGCATGCA".to_vec(),
        ];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        let jc = Jc69Model::new();
        let hky = Hky85Model::new(2.0, [0.25; 4]).unwrap();

        let candidates: Vec<(&str, &dyn SubstitutionModel)> =
            vec![("JC69", &jc), ("HKY85", &hky)];

        let result = model_finder(&tree, &refs, &candidates, None).unwrap();
        assert_eq!(result.models.len(), 2);
        // JC69 should have better (lower) BIC since data is uniform and JC69 has fewer params
        assert_eq!(
            result.best_bic, "JC69",
            "JC69 should be best by BIC on uniform-frequency data"
        );
    }

    #[test]
    fn model_finder_on_protein_data() {
        let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.1,(C:0.3,D:0.1):0.1);").unwrap();
        let seqs: Vec<Vec<u8>> = vec![
            b"ACDEFGHIKL".to_vec(),
            b"ACDEFGHIKL".to_vec(),
            b"MNPQRSTVWY".to_vec(),
            b"MNPQRSTVWY".to_vec(),
        ];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        let lg = crate::protein_models::LgModel::new();
        let wag = crate::protein_models::WagModel::new();

        let candidates: Vec<(&str, &dyn SubstitutionModel)> =
            vec![("LG", &lg), ("WAG", &wag)];

        let result = model_finder(&tree, &refs, &candidates, None).unwrap();
        assert_eq!(result.models.len(), 2);
        // Both should have finite likelihoods
        for m in &result.models {
            assert!(m.log_likelihood.is_finite());
            assert!(m.aic.is_finite());
        }
    }

    #[test]
    fn lrt_large_stat_is_significant() {
        let result = lrt(-200.0, -100.0, 1);
        assert!(result.reject, "large LRT stat should reject null");
        assert!(result.p_value < 0.001);
    }
}
