//! CellChat-style cell-cell communication analysis for spatial transcriptomics.
//!
//! Extends the basic ligand-receptor scoring in `spatial.rs` with:
//! - A curated L-R interaction database
//! - Multi-subunit complex support
//! - Pathway-level aggregation
//! - Communication probability based on mass-action kinetics
//! - Spatial distance weighting

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// A ligand-receptor interaction pair (possibly multi-subunit).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LrPair {
    /// Interaction name (e.g., "CXCL12_CXCR4").
    pub name: String,
    /// Ligand gene(s) — multiple for complexes.
    pub ligand_genes: Vec<String>,
    /// Receptor gene(s) — multiple for complexes.
    pub receptor_genes: Vec<String>,
    /// Signaling pathway (e.g., "CXCL", "WNT", "NOTCH").
    pub pathway: String,
    /// Interaction type (e.g., "Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact").
    pub interaction_type: String,
}

/// Result of a cell-cell communication analysis.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CommunicationResult {
    /// L-R pair name.
    pub lr_pair: String,
    /// Source cell type / domain.
    pub source: String,
    /// Target cell type / domain.
    pub target: String,
    /// Communication probability (0-1).
    pub probability: f64,
    /// p-value from permutation test.
    pub p_value: f64,
    /// Pathway name.
    pub pathway: String,
}

/// Aggregated pathway-level communication.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PathwayCommunication {
    /// Pathway name.
    pub pathway: String,
    /// Source cell type.
    pub source: String,
    /// Target cell type.
    pub target: String,
    /// Aggregated communication strength (sum of probabilities).
    pub strength: f64,
    /// Number of significant L-R pairs.
    pub n_significant: usize,
}

/// Parameters for communication analysis.
#[derive(Debug, Clone)]
pub struct CommParams {
    /// Distance decay factor for spatial weighting (σ in Gaussian kernel).
    pub distance_sigma: f64,
    /// Number of permutations for significance testing.
    pub n_permutations: usize,
    /// p-value threshold for significance.
    pub p_threshold: f64,
    /// Minimum expression fraction in a cell type to consider.
    pub min_pct: f64,
    /// Random seed.
    pub seed: u64,
}

impl Default for CommParams {
    fn default() -> Self {
        Self {
            distance_sigma: 100.0,
            n_permutations: 100,
            p_threshold: 0.05,
            min_pct: 0.1,
            seed: 42,
        }
    }
}

// ---------------------------------------------------------------------------
// L-R database
// ---------------------------------------------------------------------------

/// Return a curated set of ligand-receptor interactions for demo/analysis.
///
/// Includes major signaling families: chemokines, growth factors, Notch,
/// WNT, integrins, and immune checkpoints.
pub fn demo_lr_database() -> Vec<LrPair> {
    vec![
        LrPair {
            name: "CXCL12_CXCR4".into(),
            ligand_genes: vec!["CXCL12".into()],
            receptor_genes: vec!["CXCR4".into()],
            pathway: "CXCL".into(),
            interaction_type: "Secreted Signaling".into(),
        },
        LrPair {
            name: "CCL2_CCR2".into(),
            ligand_genes: vec!["CCL2".into()],
            receptor_genes: vec!["CCR2".into()],
            pathway: "CCL".into(),
            interaction_type: "Secreted Signaling".into(),
        },
        LrPair {
            name: "DLL1_NOTCH1".into(),
            ligand_genes: vec!["DLL1".into()],
            receptor_genes: vec!["NOTCH1".into()],
            pathway: "NOTCH".into(),
            interaction_type: "Cell-Cell Contact".into(),
        },
        LrPair {
            name: "JAG1_NOTCH2".into(),
            ligand_genes: vec!["JAG1".into()],
            receptor_genes: vec!["NOTCH2".into()],
            pathway: "NOTCH".into(),
            interaction_type: "Cell-Cell Contact".into(),
        },
        LrPair {
            name: "WNT3A_FZD1_LRP6".into(),
            ligand_genes: vec!["WNT3A".into()],
            receptor_genes: vec!["FZD1".into(), "LRP6".into()],
            pathway: "WNT".into(),
            interaction_type: "Secreted Signaling".into(),
        },
        LrPair {
            name: "VEGFA_FLT1".into(),
            ligand_genes: vec!["VEGFA".into()],
            receptor_genes: vec!["FLT1".into()],
            pathway: "VEGF".into(),
            interaction_type: "Secreted Signaling".into(),
        },
        LrPair {
            name: "TGFB1_TGFBR1_TGFBR2".into(),
            ligand_genes: vec!["TGFB1".into()],
            receptor_genes: vec!["TGFBR1".into(), "TGFBR2".into()],
            pathway: "TGFb".into(),
            interaction_type: "Secreted Signaling".into(),
        },
        LrPair {
            name: "COL1A1_ITGA1_ITGB1".into(),
            ligand_genes: vec!["COL1A1".into()],
            receptor_genes: vec!["ITGA1".into(), "ITGB1".into()],
            pathway: "COLLAGEN".into(),
            interaction_type: "ECM-Receptor".into(),
        },
        LrPair {
            name: "FN1_ITGAV_ITGB3".into(),
            ligand_genes: vec!["FN1".into()],
            receptor_genes: vec!["ITGAV".into(), "ITGB3".into()],
            pathway: "FN1".into(),
            interaction_type: "ECM-Receptor".into(),
        },
        LrPair {
            name: "CD274_PDCD1".into(),
            ligand_genes: vec!["CD274".into()],
            receptor_genes: vec!["PDCD1".into()],
            pathway: "PD-L1".into(),
            interaction_type: "Cell-Cell Contact".into(),
        },
        LrPair {
            name: "EGF_EGFR".into(),
            ligand_genes: vec!["EGF".into()],
            receptor_genes: vec!["EGFR".into()],
            pathway: "EGF".into(),
            interaction_type: "Secreted Signaling".into(),
        },
        LrPair {
            name: "HGF_MET".into(),
            ligand_genes: vec!["HGF".into()],
            receptor_genes: vec!["MET".into()],
            pathway: "HGF".into(),
            interaction_type: "Secreted Signaling".into(),
        },
    ]
}

// ---------------------------------------------------------------------------
// Communication analysis
// ---------------------------------------------------------------------------

/// Compute cell-cell communication probabilities using spatial expression data.
///
/// `expression` is cell-major: `expression[cell][gene]`.
/// `gene_names` maps column index to gene name.
/// `cell_types` assigns each cell to a cell type string.
/// `coords` are `(x, y)` per cell.
/// `lr_pairs` is the set of L-R interactions to test.
///
/// The communication probability for an L-R pair between source type A and
/// target type B is computed as:
///
/// `P = mean(L_i * R_j * w_ij)` for cells i in A, j in B
///
/// where `L_i` is the geometric mean of ligand subunit expression at cell i,
/// `R_j` is the geometric mean of receptor subunit expression at cell j,
/// and `w_ij = exp(-d²/(2σ²))` is a Gaussian spatial weight.
pub fn analyze_communication(
    expression: &[Vec<f64>],
    gene_names: &[String],
    cell_types: &[String],
    coords: &[(f64, f64)],
    lr_pairs: &[LrPair],
    params: &CommParams,
) -> Result<Vec<CommunicationResult>> {
    let n_cells = expression.len();
    if n_cells == 0 {
        return Err(CyaneaError::InvalidInput("empty expression".into()));
    }
    if gene_names.len() != expression[0].len()
        || cell_types.len() != n_cells
        || coords.len() != n_cells
    {
        return Err(CyaneaError::InvalidInput(
            "dimension mismatch".into(),
        ));
    }

    // Build gene name → index map
    let gene_idx: std::collections::HashMap<&str, usize> = gene_names
        .iter()
        .enumerate()
        .map(|(i, g)| (g.as_str(), i))
        .collect();

    // Unique cell types
    let mut type_set: Vec<String> = cell_types.to_vec();
    type_set.sort();
    type_set.dedup();

    // Cells per type
    let cells_by_type: std::collections::HashMap<&str, Vec<usize>> = type_set
        .iter()
        .map(|t| {
            let members: Vec<usize> = (0..n_cells)
                .filter(|&i| cell_types[i] == *t)
                .collect();
            (t.as_str(), members)
        })
        .collect();

    // Spatial weight
    let sigma2 = 2.0 * params.distance_sigma * params.distance_sigma;

    let mut results = Vec::new();
    let mut rng_state = if params.seed == 0 {
        0x5EED_DEAD_BEEF_CAFE
    } else {
        params.seed
    };

    for lr in lr_pairs {
        // Find gene indices for ligand and receptor subunits
        let lig_indices: Vec<usize> = lr
            .ligand_genes
            .iter()
            .filter_map(|g| gene_idx.get(g.as_str()).copied())
            .collect();
        let rec_indices: Vec<usize> = lr
            .receptor_genes
            .iter()
            .filter_map(|g| gene_idx.get(g.as_str()).copied())
            .collect();

        // Skip if genes not found
        if lig_indices.len() != lr.ligand_genes.len()
            || rec_indices.len() != lr.receptor_genes.len()
        {
            continue;
        }

        // Compute ligand/receptor scores per cell (geometric mean of subunits)
        let lig_scores: Vec<f64> = (0..n_cells)
            .map(|i| {
                let prod: f64 = lig_indices.iter().map(|&g| expression[i][g].max(0.0)).product();
                prod.powf(1.0 / lig_indices.len() as f64)
            })
            .collect();
        let rec_scores: Vec<f64> = (0..n_cells)
            .map(|i| {
                let prod: f64 = rec_indices.iter().map(|&g| expression[i][g].max(0.0)).product();
                prod.powf(1.0 / rec_indices.len() as f64)
            })
            .collect();

        // Test all source → target type pairs
        for src_type in &type_set {
            for tgt_type in &type_set {
                let src_cells = &cells_by_type[src_type.as_str()];
                let tgt_cells = &cells_by_type[tgt_type.as_str()];

                if src_cells.is_empty() || tgt_cells.is_empty() {
                    continue;
                }

                // Check minimum expression fraction
                let src_pct = src_cells.iter().filter(|&&i| lig_scores[i] > 0.0).count() as f64
                    / src_cells.len() as f64;
                let tgt_pct = tgt_cells.iter().filter(|&&i| rec_scores[i] > 0.0).count() as f64
                    / tgt_cells.len() as f64;

                if src_pct < params.min_pct || tgt_pct < params.min_pct {
                    continue;
                }

                // Compute communication probability
                let prob = compute_comm_prob(
                    &lig_scores,
                    &rec_scores,
                    src_cells,
                    tgt_cells,
                    coords,
                    sigma2,
                );

                // Permutation test
                let mut perm_ge = 0usize;
                let mut shuffled_types = cell_types.to_vec();
                for _ in 0..params.n_permutations {
                    // Shuffle cell type labels
                    for i in (1..n_cells).rev() {
                        rng_state ^= rng_state << 13;
                        rng_state ^= rng_state >> 7;
                        rng_state ^= rng_state << 17;
                        let j = (rng_state as usize) % (i + 1);
                        shuffled_types.swap(i, j);
                    }
                    let perm_src: Vec<usize> = (0..n_cells)
                        .filter(|&i| shuffled_types[i] == *src_type)
                        .collect();
                    let perm_tgt: Vec<usize> = (0..n_cells)
                        .filter(|&i| shuffled_types[i] == *tgt_type)
                        .collect();
                    if perm_src.is_empty() || perm_tgt.is_empty() {
                        continue;
                    }
                    let perm_prob = compute_comm_prob(
                        &lig_scores,
                        &rec_scores,
                        &perm_src,
                        &perm_tgt,
                        coords,
                        sigma2,
                    );
                    if perm_prob >= prob {
                        perm_ge += 1;
                    }
                }

                let p_value = (perm_ge as f64 + 1.0) / (params.n_permutations as f64 + 1.0);

                results.push(CommunicationResult {
                    lr_pair: lr.name.clone(),
                    source: src_type.clone(),
                    target: tgt_type.clone(),
                    probability: prob,
                    p_value,
                    pathway: lr.pathway.clone(),
                });
            }
        }
    }

    Ok(results)
}

/// Compute mean L*R*w across source-target cell pairs.
fn compute_comm_prob(
    lig_scores: &[f64],
    rec_scores: &[f64],
    src_cells: &[usize],
    tgt_cells: &[usize],
    coords: &[(f64, f64)],
    sigma2: f64,
) -> f64 {
    let mut total = 0.0_f64;
    let mut count = 0usize;

    for &i in src_cells {
        for &j in tgt_cells {
            if i == j {
                continue;
            }
            let d2 = (coords[i].0 - coords[j].0).powi(2) + (coords[i].1 - coords[j].1).powi(2);
            let w = (-d2 / sigma2).exp();
            total += lig_scores[i] * rec_scores[j] * w;
            count += 1;
        }
    }

    if count > 0 {
        total / count as f64
    } else {
        0.0
    }
}

/// Aggregate communication results to pathway level.
pub fn aggregate_pathways(
    results: &[CommunicationResult],
    p_threshold: f64,
) -> Vec<PathwayCommunication> {
    let mut pathway_map: std::collections::HashMap<(String, String, String), (f64, usize)> =
        std::collections::HashMap::new();

    for r in results {
        let key = (r.pathway.clone(), r.source.clone(), r.target.clone());
        let entry = pathway_map.entry(key).or_insert((0.0, 0));
        entry.0 += r.probability;
        if r.p_value < p_threshold {
            entry.1 += 1;
        }
    }

    let mut pathways: Vec<PathwayCommunication> = pathway_map
        .into_iter()
        .map(|((pathway, source, target), (strength, n_sig))| PathwayCommunication {
            pathway,
            source,
            target,
            strength,
            n_significant: n_sig,
        })
        .collect();

    pathways.sort_by(|a, b| {
        b.strength
            .partial_cmp(&a.strength)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    pathways
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_demo_lr_database() {
        let db = demo_lr_database();
        assert_eq!(db.len(), 12);
        // Check multi-subunit receptor
        let wnt = db.iter().find(|p| p.name == "WNT3A_FZD1_LRP6").unwrap();
        assert_eq!(wnt.receptor_genes.len(), 2);
        assert_eq!(wnt.pathway, "WNT");
    }

    #[test]
    fn test_analyze_communication_basic() {
        // Simple scenario: 2 cell types, 1 L-R pair
        let gene_names = vec!["LIG1".into(), "REC1".into()];
        let expression = vec![
            vec![10.0, 0.0], // type A, high ligand
            vec![10.0, 0.0], // type A, high ligand
            vec![0.0, 10.0], // type B, high receptor
            vec![0.0, 10.0], // type B, high receptor
        ];
        let cell_types = vec!["A".into(), "A".into(), "B".into(), "B".into()];
        let coords = vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)];
        let lr_pairs = vec![LrPair {
            name: "LIG1_REC1".into(),
            ligand_genes: vec!["LIG1".into()],
            receptor_genes: vec!["REC1".into()],
            pathway: "TEST".into(),
            interaction_type: "Secreted Signaling".into(),
        }];

        let params = CommParams {
            distance_sigma: 10.0,
            n_permutations: 50,
            p_threshold: 0.05,
            min_pct: 0.1,
            seed: 42,
        };

        let results = analyze_communication(
            &expression,
            &gene_names,
            &cell_types,
            &coords,
            &lr_pairs,
            &params,
        )
        .unwrap();

        // Should have A→B with high probability (ligand in A, receptor in B)
        let a_to_b = results
            .iter()
            .find(|r| r.source == "A" && r.target == "B")
            .unwrap();
        assert!(a_to_b.probability > 0.0);
    }

    #[test]
    fn test_communication_no_expression() {
        // Genes not in dataset → no results
        let gene_names = vec!["GENE1".into(), "GENE2".into()];
        let expression = vec![vec![1.0, 2.0]; 4];
        let cell_types = vec!["A".into(), "A".into(), "B".into(), "B".into()];
        let coords = vec![(0.0, 0.0); 4];
        let lr_pairs = vec![LrPair {
            name: "MISSING_L_R".into(),
            ligand_genes: vec!["NOTHERE".into()],
            receptor_genes: vec!["ALSO_NOT".into()],
            pathway: "X".into(),
            interaction_type: "Secreted Signaling".into(),
        }];

        let params = CommParams::default();
        let results = analyze_communication(
            &expression,
            &gene_names,
            &cell_types,
            &coords,
            &lr_pairs,
            &params,
        )
        .unwrap();
        assert!(results.is_empty());
    }

    #[test]
    fn test_multi_subunit() {
        // Multi-subunit receptor: both subunits must be expressed
        let gene_names = vec!["LIG".into(), "REC_A".into(), "REC_B".into()];
        let expression = vec![
            vec![10.0, 0.0, 0.0], // type A
            vec![0.0, 10.0, 10.0], // type B: both receptor subunits
            vec![0.0, 10.0, 0.0], // type C: only one subunit
        ];
        let cell_types = vec!["A".into(), "B".into(), "C".into()];
        let coords = vec![(0.0, 0.0), (5.0, 0.0), (10.0, 0.0)];
        let lr_pairs = vec![LrPair {
            name: "LIG_RECA_RECB".into(),
            ligand_genes: vec!["LIG".into()],
            receptor_genes: vec!["REC_A".into(), "REC_B".into()],
            pathway: "TEST".into(),
            interaction_type: "Secreted Signaling".into(),
        }];

        let params = CommParams {
            distance_sigma: 20.0,
            n_permutations: 0,
            min_pct: 0.0,
            ..Default::default()
        };

        let results = analyze_communication(
            &expression,
            &gene_names,
            &cell_types,
            &coords,
            &lr_pairs,
            &params,
        )
        .unwrap();

        // A→B should have positive probability (both receptor subunits)
        let a_b = results.iter().find(|r| r.source == "A" && r.target == "B");
        // A→C should have zero probability (only one receptor subunit, geometric mean = 0)
        let a_c = results.iter().find(|r| r.source == "A" && r.target == "C");

        if let Some(ab) = a_b {
            assert!(ab.probability > 0.0);
        }
        if let Some(ac) = a_c {
            assert!(ac.probability < 0.01);
        }
    }

    #[test]
    fn test_aggregate_pathways() {
        let results = vec![
            CommunicationResult {
                lr_pair: "L1_R1".into(),
                source: "A".into(),
                target: "B".into(),
                probability: 0.5,
                p_value: 0.01,
                pathway: "WNT".into(),
            },
            CommunicationResult {
                lr_pair: "L2_R2".into(),
                source: "A".into(),
                target: "B".into(),
                probability: 0.3,
                p_value: 0.1,
                pathway: "WNT".into(),
            },
        ];

        let pathways = aggregate_pathways(&results, 0.05);
        assert_eq!(pathways.len(), 1);
        assert!((pathways[0].strength - 0.8).abs() < 1e-10);
        assert_eq!(pathways[0].n_significant, 1); // only first has p < 0.05
    }

    #[test]
    fn test_spatial_distance_weighting() {
        // Closer cells should have higher communication probability
        let gene_names = vec!["L".into(), "R".into()];
        let expression = vec![
            vec![10.0, 0.0], // source
            vec![0.0, 10.0], // close target
            vec![0.0, 10.0], // far target
        ];
        let cell_types = vec!["src".into(), "close".into(), "far".into()];
        let coords = vec![(0.0, 0.0), (1.0, 0.0), (100.0, 0.0)];
        let lr_pairs = vec![LrPair {
            name: "L_R".into(),
            ligand_genes: vec!["L".into()],
            receptor_genes: vec!["R".into()],
            pathway: "T".into(),
            interaction_type: "Secreted Signaling".into(),
        }];
        let params = CommParams {
            distance_sigma: 5.0,
            n_permutations: 0,
            min_pct: 0.0,
            ..Default::default()
        };

        let results = analyze_communication(
            &expression,
            &gene_names,
            &cell_types,
            &coords,
            &lr_pairs,
            &params,
        )
        .unwrap();

        let to_close = results
            .iter()
            .find(|r| r.source == "src" && r.target == "close");
        let to_far = results
            .iter()
            .find(|r| r.source == "src" && r.target == "far");

        if let (Some(close), Some(far)) = (to_close, to_far) {
            assert!(close.probability > far.probability);
        }
    }
}
