//! Spatial transcriptomics deconvolution — estimating cell type proportions
//! within each spot from bulk-like expression profiles.
//!
//! Implements NNLS-based and regression-based deconvolution, plus enrichment
//! scoring for cell type signatures.

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Cell type signature: a set of marker genes with expected expression levels.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CellTypeSignature {
    /// Cell type name.
    pub cell_type: String,
    /// Gene names in the signature.
    pub genes: Vec<String>,
    /// Expected expression values (reference profile).
    pub weights: Vec<f64>,
}

/// Deconvolution result for a single spot.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SpotDeconvolution {
    /// Spot index.
    pub spot_idx: usize,
    /// Estimated cell type proportions (sums to ~1).
    pub proportions: Vec<f64>,
    /// Residual error (reconstruction quality).
    pub residual: f64,
}

/// Full deconvolution result.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DeconvolutionResult {
    /// Cell type names.
    pub cell_types: Vec<String>,
    /// Per-spot deconvolution results.
    pub spots: Vec<SpotDeconvolution>,
    /// Mean residual across all spots.
    pub mean_residual: f64,
}

/// Enrichment score for a cell type at a spot.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct EnrichmentScore {
    /// Spot index.
    pub spot_idx: usize,
    /// Cell type name.
    pub cell_type: String,
    /// Enrichment score (higher = more enriched).
    pub score: f64,
    /// p-value from permutation test.
    pub p_value: f64,
}

// ---------------------------------------------------------------------------
// NNLS deconvolution
// ---------------------------------------------------------------------------

/// Deconvolve spatial spots using non-negative least squares (NNLS).
///
/// Given spot expression profiles and cell type reference signatures,
/// estimates the proportion of each cell type in each spot.
///
/// `expression` is spot-major: `expression[spot][gene]`.
/// `gene_names` are the gene identifiers for columns.
/// `signatures` defines reference profiles for each cell type.
///
/// Uses iterative NNLS: starts with uniform proportions and iteratively
/// adjusts to minimize reconstruction error while keeping proportions >= 0.
pub fn nnls_deconvolve(
    expression: &[Vec<f64>],
    gene_names: &[String],
    signatures: &[CellTypeSignature],
) -> Result<DeconvolutionResult> {
    let n_spots = expression.len();
    if n_spots == 0 {
        return Err(CyaneaError::InvalidInput("empty expression".into()));
    }
    let n_genes = gene_names.len();
    if expression[0].len() != n_genes {
        return Err(CyaneaError::InvalidInput(
            "gene count mismatch".into(),
        ));
    }
    if signatures.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "need at least 1 cell type signature".into(),
        ));
    }

    let n_types = signatures.len();

    // Build reference matrix: reference[type][gene] (only shared genes)
    // Find shared genes
    let mut shared_genes: Vec<usize> = Vec::new(); // indices into gene_names
    let mut sig_gene_indices: Vec<Vec<Option<usize>>> = vec![Vec::new(); n_types];

    for (gi, gname) in gene_names.iter().enumerate() {
        let mut in_all = true;
        for (ti, sig) in signatures.iter().enumerate() {
            if let Some(pos) = sig.genes.iter().position(|g| g == gname) {
                if sig_gene_indices[ti].len() <= gi {
                    sig_gene_indices[ti].resize(gi + 1, None);
                }
                sig_gene_indices[ti][gi] = Some(pos);
            } else {
                in_all = false;
            }
        }
        if in_all {
            shared_genes.push(gi);
        }
    }

    // Use any gene found in at least one signature
    let mut used_genes: Vec<usize> = Vec::new();
    for (gi, _) in gene_names.iter().enumerate() {
        let in_any = signatures.iter().any(|sig| sig.genes.iter().any(|g| g == &gene_names[gi]));
        if in_any {
            used_genes.push(gi);
        }
    }

    if used_genes.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "no shared genes between expression data and signatures".into(),
        ));
    }

    // Build reference matrix for used genes
    let n_used = used_genes.len();
    let mut reference: Vec<Vec<f64>> = vec![vec![0.0; n_used]; n_types];
    for (ti, sig) in signatures.iter().enumerate() {
        let sig_idx: std::collections::HashMap<&str, usize> = sig
            .genes
            .iter()
            .enumerate()
            .map(|(i, g)| (g.as_str(), i))
            .collect();
        for (ui, &gi) in used_genes.iter().enumerate() {
            if let Some(&pos) = sig_idx.get(gene_names[gi].as_str()) {
                reference[ti][ui] = sig.weights[pos];
            }
        }
    }

    // NNLS for each spot
    let mut spots = Vec::with_capacity(n_spots);
    let mut total_residual = 0.0;

    for si in 0..n_spots {
        // Extract observed values for used genes
        let observed: Vec<f64> = used_genes.iter().map(|&gi| expression[si][gi]).collect();

        // Iterative NNLS (coordinate descent)
        let mut proportions = vec![1.0 / n_types as f64; n_types];
        let max_iter = 200;

        for _ in 0..max_iter {
            let mut changed = false;
            for t in 0..n_types {
                // Compute residual without type t's contribution
                let mut residual_without: Vec<f64> = observed.clone();
                for t2 in 0..n_types {
                    if t2 == t {
                        continue;
                    }
                    for g in 0..n_used {
                        residual_without[g] -= proportions[t2] * reference[t2][g];
                    }
                }

                // Optimal proportion for type t: dot(residual, reference[t]) / dot(reference[t], reference[t])
                let dot_rr: f64 = reference[t].iter().map(|x| x * x).sum();
                if dot_rr < 1e-30 {
                    proportions[t] = 0.0;
                    continue;
                }
                let dot_res: f64 = residual_without
                    .iter()
                    .zip(reference[t].iter())
                    .map(|(r, f)| r * f)
                    .sum();
                let new_p = (dot_res / dot_rr).max(0.0);

                if (new_p - proportions[t]).abs() > 1e-8 {
                    changed = true;
                }
                proportions[t] = new_p;
            }
            if !changed {
                break;
            }
        }

        // Normalize to sum to 1
        let sum: f64 = proportions.iter().sum();
        if sum > 0.0 {
            for p in &mut proportions {
                *p /= sum;
            }
        }

        // Compute residual
        let mut residual = 0.0;
        for g in 0..n_used {
            let predicted: f64 = (0..n_types)
                .map(|t| proportions[t] * reference[t][g])
                .sum();
            residual += (observed[g] - predicted).powi(2);
        }
        let residual = (residual / n_used as f64).sqrt();
        total_residual += residual;

        spots.push(SpotDeconvolution {
            spot_idx: si,
            proportions,
            residual,
        });
    }

    let cell_types = signatures.iter().map(|s| s.cell_type.clone()).collect();

    Ok(DeconvolutionResult {
        cell_types,
        spots,
        mean_residual: total_residual / n_spots as f64,
    })
}

// ---------------------------------------------------------------------------
// Enrichment scoring
// ---------------------------------------------------------------------------

/// Score enrichment of cell type signatures at each spot.
///
/// For each spot and cell type, computes a score based on the mean z-scored
/// expression of the signature genes. Higher scores indicate higher enrichment.
///
/// `expression` is spot-major: `expression[spot][gene]`.
/// `gene_names` are column labels.
/// `signatures` defines marker genes per cell type.
/// `n_permutations` controls permutation-based p-value estimation.
pub fn score_enrichment(
    expression: &[Vec<f64>],
    gene_names: &[String],
    signatures: &[CellTypeSignature],
    n_permutations: usize,
    seed: u64,
) -> Result<Vec<EnrichmentScore>> {
    let n_spots = expression.len();
    if n_spots == 0 {
        return Err(CyaneaError::InvalidInput("empty expression".into()));
    }
    let n_genes = gene_names.len();

    // Z-score expression per gene
    let mut z_expr: Vec<Vec<f64>> = expression.to_vec();
    for g in 0..n_genes {
        let mean = expression.iter().map(|r| r[g]).sum::<f64>() / n_spots as f64;
        let var = expression.iter().map(|r| (r[g] - mean).powi(2)).sum::<f64>() / n_spots as f64;
        let std = var.sqrt().max(1e-10);
        for s in 0..n_spots {
            z_expr[s][g] = (expression[s][g] - mean) / std;
        }
    }

    // Gene index map
    let gene_idx: std::collections::HashMap<&str, usize> = gene_names
        .iter()
        .enumerate()
        .map(|(i, g)| (g.as_str(), i))
        .collect();

    let mut rng_state = if seed == 0 { 0x5EED_CAFE_BABE } else { seed };
    let xorshift = |state: &mut u64| -> u64 {
        let mut x = *state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        *state = x;
        x
    };

    let mut results = Vec::new();

    for sig in signatures {
        // Find signature gene indices
        let sig_indices: Vec<usize> = sig
            .genes
            .iter()
            .filter_map(|g| gene_idx.get(g.as_str()).copied())
            .collect();

        if sig_indices.is_empty() {
            continue;
        }

        let n_sig = sig_indices.len();

        for si in 0..n_spots {
            // Score = mean z-score of signature genes
            let score: f64 = sig_indices.iter().map(|&gi| z_expr[si][gi]).sum::<f64>() / n_sig as f64;

            // Permutation: sample random gene sets of same size
            let mut perm_ge = 0usize;
            let all_indices: Vec<usize> = (0..n_genes).collect();
            for _ in 0..n_permutations {
                // Random sample of n_sig genes
                let mut sampled = Vec::with_capacity(n_sig);
                let mut remaining: Vec<usize> = all_indices.clone();
                for _ in 0..n_sig {
                    if remaining.is_empty() {
                        break;
                    }
                    let idx = (xorshift(&mut rng_state) as usize) % remaining.len();
                    sampled.push(remaining.swap_remove(idx));
                }
                let perm_score: f64 =
                    sampled.iter().map(|&gi| z_expr[si][gi]).sum::<f64>() / n_sig as f64;
                if perm_score >= score {
                    perm_ge += 1;
                }
            }

            let p_value = if n_permutations > 0 {
                (perm_ge as f64 + 1.0) / (n_permutations as f64 + 1.0)
            } else {
                1.0
            };

            results.push(EnrichmentScore {
                spot_idx: si,
                cell_type: sig.cell_type.clone(),
                score,
                p_value,
            });
        }
    }

    Ok(results)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_signatures() -> Vec<CellTypeSignature> {
        vec![
            CellTypeSignature {
                cell_type: "T_cell".into(),
                genes: vec!["CD3D".into(), "CD3E".into(), "CD8A".into()],
                weights: vec![10.0, 8.0, 6.0],
            },
            CellTypeSignature {
                cell_type: "Macrophage".into(),
                genes: vec!["CD68".into(), "CD163".into(), "CSF1R".into()],
                weights: vec![9.0, 7.0, 5.0],
            },
        ]
    }

    #[test]
    fn test_nnls_pure_spot() {
        let gene_names = vec![
            "CD3D".into(),
            "CD3E".into(),
            "CD8A".into(),
            "CD68".into(),
            "CD163".into(),
            "CSF1R".into(),
        ];
        // Spot with pure T cell expression
        let expression = vec![vec![10.0, 8.0, 6.0, 0.0, 0.0, 0.0]];
        let signatures = make_signatures();

        let result = nnls_deconvolve(&expression, &gene_names, &signatures).unwrap();
        assert_eq!(result.cell_types.len(), 2);
        assert_eq!(result.spots.len(), 1);

        let props = &result.spots[0].proportions;
        // T cell proportion should dominate
        assert!(props[0] > 0.8, "T cell proportion should be high: {}", props[0]);
    }

    #[test]
    fn test_nnls_mixed_spot() {
        let gene_names = vec![
            "CD3D".into(),
            "CD3E".into(),
            "CD8A".into(),
            "CD68".into(),
            "CD163".into(),
            "CSF1R".into(),
        ];
        // Mixed spot: half T cell, half macrophage
        let expression = vec![vec![5.0, 4.0, 3.0, 4.5, 3.5, 2.5]];
        let signatures = make_signatures();

        let result = nnls_deconvolve(&expression, &gene_names, &signatures).unwrap();
        let props = &result.spots[0].proportions;
        // Both should be substantial
        assert!(props[0] > 0.2 && props[0] < 0.8);
        assert!(props[1] > 0.2 && props[1] < 0.8);
        // Should sum to ~1
        let sum: f64 = props.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_nnls_multiple_spots() {
        let gene_names = vec![
            "CD3D".into(),
            "CD3E".into(),
            "CD8A".into(),
            "CD68".into(),
            "CD163".into(),
            "CSF1R".into(),
        ];
        let expression = vec![
            vec![10.0, 8.0, 6.0, 0.0, 0.0, 0.0], // T cell
            vec![0.0, 0.0, 0.0, 9.0, 7.0, 5.0],   // Macrophage
            vec![5.0, 4.0, 3.0, 4.5, 3.5, 2.5],   // Mixed
        ];
        let signatures = make_signatures();

        let result = nnls_deconvolve(&expression, &gene_names, &signatures).unwrap();
        assert_eq!(result.spots.len(), 3);
        assert!(result.mean_residual >= 0.0);
    }

    #[test]
    fn test_nnls_no_shared_genes() {
        let gene_names = vec!["UNRELATED1".into(), "UNRELATED2".into()];
        let expression = vec![vec![1.0, 2.0]];
        let signatures = make_signatures();

        let result = nnls_deconvolve(&expression, &gene_names, &signatures);
        assert!(result.is_err());
    }

    #[test]
    fn test_enrichment_scoring() {
        let gene_names = vec![
            "CD3D".into(),
            "CD3E".into(),
            "CD8A".into(),
            "CD68".into(),
            "CD163".into(),
            "CSF1R".into(),
        ];
        // Spot 0: T cell enriched, Spot 1: Macrophage enriched
        let expression = vec![
            vec![20.0, 15.0, 12.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 18.0, 14.0, 10.0],
        ];
        let signatures = make_signatures();

        let results = score_enrichment(&expression, &gene_names, &signatures, 99, 42).unwrap();
        // 2 spots × 2 cell types = 4 scores
        assert_eq!(results.len(), 4);

        // T cell enrichment should be higher at spot 0
        let t_spot0 = results
            .iter()
            .find(|r| r.spot_idx == 0 && r.cell_type == "T_cell")
            .unwrap();
        let t_spot1 = results
            .iter()
            .find(|r| r.spot_idx == 1 && r.cell_type == "T_cell")
            .unwrap();
        assert!(t_spot0.score > t_spot1.score);
    }

    #[test]
    fn test_enrichment_p_values() {
        let gene_names: Vec<String> = (0..20).map(|i| format!("Gene{}", i)).collect();
        let expression = vec![vec![1.0; 20]; 5];
        let signatures = vec![CellTypeSignature {
            cell_type: "Test".into(),
            genes: vec!["Gene0".into(), "Gene1".into()],
            weights: vec![1.0, 1.0],
        }];

        let results = score_enrichment(&expression, &gene_names, &signatures, 99, 1).unwrap();
        for r in &results {
            assert!(r.p_value > 0.0 && r.p_value <= 1.0);
        }
    }

    #[test]
    fn test_proportions_sum_to_one() {
        let gene_names = vec!["G1".into(), "G2".into(), "G3".into(), "G4".into()];
        let expression = vec![vec![5.0, 3.0, 7.0, 2.0]];
        let signatures = vec![
            CellTypeSignature {
                cell_type: "TypeA".into(),
                genes: vec!["G1".into(), "G2".into()],
                weights: vec![10.0, 8.0],
            },
            CellTypeSignature {
                cell_type: "TypeB".into(),
                genes: vec!["G3".into(), "G4".into()],
                weights: vec![9.0, 6.0],
            },
        ];

        let result = nnls_deconvolve(&expression, &gene_names, &signatures).unwrap();
        let sum: f64 = result.spots[0].proportions.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-10,
            "proportions should sum to 1, got {}",
            sum
        );
    }
}
