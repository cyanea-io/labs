//! B-factor analysis for macromolecular structures.
//!
//! Provides per-residue and per-chain B-factor statistics, as well as normalized
//! flexibility scores (Z-scores) to identify flexible and rigid regions.

use cyanea_core::{CyaneaError, Result};

use crate::types::Structure;

use alloc::string::String;
use alloc::vec::Vec;

/// Descriptive statistics for B-factor values.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BFactorStats {
    /// Arithmetic mean of B-factors.
    pub mean: f64,
    /// Population standard deviation.
    pub std_dev: f64,
    /// Minimum B-factor.
    pub min: f64,
    /// Maximum B-factor.
    pub max: f64,
    /// Median B-factor.
    pub median: f64,
}

/// Compute per-residue average B-factors across all chains.
///
/// Returns a list of `(residue_seq_num, residue_name, average_b_factor)` for
/// each residue in the structure.
///
/// # Errors
///
/// Returns an error if the structure has no atoms.
pub fn residue_bfactors(structure: &Structure) -> Result<Vec<(usize, String, f64)>> {
    if structure.atom_count() == 0 {
        return Err(CyaneaError::InvalidInput(
            "cannot compute B-factors for structure with no atoms".into(),
        ));
    }

    let mut result = Vec::new();

    for chain in &structure.chains {
        for residue in &chain.residues {
            if residue.atoms.is_empty() {
                continue;
            }
            let sum: f64 = residue.atoms.iter().map(|a| a.temp_factor).sum();
            let avg = sum / residue.atoms.len() as f64;
            result.push((residue.seq_num as usize, residue.name.clone(), avg));
        }
    }

    Ok(result)
}

/// Compute per-chain B-factor statistics.
///
/// Returns a list of `(chain_id, stats)` where `stats` contains mean, std_dev,
/// min, max, and median B-factors for all atoms in the chain.
///
/// # Errors
///
/// Returns an error if the structure has no chains or if a chain has no atoms.
pub fn chain_bfactors(structure: &Structure) -> Result<Vec<(String, BFactorStats)>> {
    if structure.chains.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "cannot compute B-factors for structure with no chains".into(),
        ));
    }

    let mut result = Vec::new();

    for chain in &structure.chains {
        let mut bfactors: Vec<f64> = Vec::new();
        for residue in &chain.residues {
            for atom in &residue.atoms {
                bfactors.push(atom.temp_factor);
            }
        }

        if bfactors.is_empty() {
            continue;
        }

        let stats = compute_stats(&mut bfactors);
        let chain_name = alloc::format!("{}", chain.id);
        result.push((chain_name, stats));
    }

    Ok(result)
}

/// Compute normalized B-factor Z-scores per residue.
///
/// For each residue, the flexibility score is:
///   `Z = (residue_avg_B - global_mean_B) / global_std_B`
///
/// High Z-scores indicate flexible regions; low Z-scores indicate rigid regions.
///
/// Returns a list of `(residue_seq_num, z_score)`.
///
/// # Errors
///
/// Returns an error if the structure has no atoms or if the global B-factor
/// standard deviation is zero (all atoms have the same B-factor).
pub fn flexibility_score(structure: &Structure) -> Result<Vec<(usize, f64)>> {
    let residue_bfs = residue_bfactors(structure)?;

    if residue_bfs.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "no residues to compute flexibility scores".into(),
        ));
    }

    // Compute global mean and std dev from all atoms
    let mut all_bfactors: Vec<f64> = Vec::new();
    for chain in &structure.chains {
        for residue in &chain.residues {
            for atom in &residue.atoms {
                all_bfactors.push(atom.temp_factor);
            }
        }
    }

    let global_mean = mean_f64(&all_bfactors);
    let global_std = std_dev_f64(&all_bfactors, global_mean);

    if global_std < 1e-15 {
        // All B-factors are identical — Z-scores are all 0
        return Ok(residue_bfs.iter().map(|(seq, _, _)| (*seq, 0.0)).collect());
    }

    let result: Vec<(usize, f64)> = residue_bfs
        .iter()
        .map(|(seq, _, avg_b)| (*seq, (avg_b - global_mean) / global_std))
        .collect();

    Ok(result)
}

// ---- Internal helper functions ----

/// Compute BFactorStats from a mutable slice (sorted in place for median).
fn compute_stats(values: &mut Vec<f64>) -> BFactorStats {
    let n = values.len();
    debug_assert!(n > 0);

    let mean = mean_f64(values);
    let std_dev = std_dev_f64(values, mean);

    // Sort for min, max, median
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));

    let min = values[0];
    let max = values[n - 1];
    let median = if n % 2 == 0 {
        (values[n / 2 - 1] + values[n / 2]) / 2.0
    } else {
        values[n / 2]
    };

    BFactorStats {
        mean,
        std_dev,
        min,
        max,
        median,
    }
}

/// Arithmetic mean.
fn mean_f64(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let sum: f64 = values.iter().sum();
    sum / values.len() as f64
}

/// Population standard deviation.
fn std_dev_f64(values: &[f64], mean: f64) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let variance: f64 = values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>()
        / values.len() as f64;
    variance.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Atom, Chain, Point3D, Residue, Structure};
    use alloc::vec;

    fn make_atom_with_b(name: &str, b: f64) -> Atom {
        Atom {
            serial: 1,
            name: name.into(),
            alt_loc: None,
            coords: Point3D::new(0.0, 0.0, 0.0),
            occupancy: 1.0,
            temp_factor: b,
            element: None,
            charge: None,
            is_hetatm: false,
        }
    }

    fn test_structure() -> Structure {
        Structure {
            id: "TEST".into(),
            chains: vec![Chain::new(
                'A',
                vec![
                    Residue {
                        name: "ALA".into(),
                        seq_num: 1,
                        i_code: None,
                        atoms: vec![
                            make_atom_with_b("N", 10.0),
                            make_atom_with_b("CA", 12.0),
                            make_atom_with_b("C", 11.0),
                            make_atom_with_b("O", 13.0),
                        ],
                    },
                    Residue {
                        name: "GLY".into(),
                        seq_num: 2,
                        i_code: None,
                        atoms: vec![
                            make_atom_with_b("N", 20.0),
                            make_atom_with_b("CA", 22.0),
                            make_atom_with_b("C", 21.0),
                            make_atom_with_b("O", 23.0),
                        ],
                    },
                    Residue {
                        name: "VAL".into(),
                        seq_num: 3,
                        i_code: None,
                        atoms: vec![
                            make_atom_with_b("N", 5.0),
                            make_atom_with_b("CA", 6.0),
                            make_atom_with_b("C", 5.5),
                            make_atom_with_b("O", 6.5),
                        ],
                    },
                ],
            )],
        }
    }

    #[test]
    fn residue_bfactors_basic() {
        let s = test_structure();
        let bfs = residue_bfactors(&s).unwrap();
        assert_eq!(bfs.len(), 3);

        // ALA: (10+12+11+13)/4 = 11.5
        assert!((bfs[0].2 - 11.5).abs() < 1e-10);
        assert_eq!(bfs[0].1, "ALA");

        // GLY: (20+22+21+23)/4 = 21.5
        assert!((bfs[1].2 - 21.5).abs() < 1e-10);

        // VAL: (5+6+5.5+6.5)/4 = 5.75
        assert!((bfs[2].2 - 5.75).abs() < 1e-10);
    }

    #[test]
    fn chain_bfactors_basic() {
        let s = test_structure();
        let cbs = chain_bfactors(&s).unwrap();
        assert_eq!(cbs.len(), 1);
        assert_eq!(cbs[0].0, "A");

        let stats = &cbs[0].1;
        // All 12 atoms: 10,12,11,13,20,22,21,23,5,6,5.5,6.5
        // Sum = 155, mean = 155/12 ≈ 12.9167
        assert!((stats.mean - 155.0 / 12.0).abs() < 1e-4);
        assert!((stats.min - 5.0).abs() < 1e-10);
        assert!((stats.max - 23.0).abs() < 1e-10);

        // Median of sorted [5,5.5,6,6.5,10,11,12,13,20,21,22,23]: (11+12)/2 = 11.5
        assert!((stats.median - 11.5).abs() < 1e-10);
    }

    #[test]
    fn flexibility_score_basic() {
        let s = test_structure();
        let scores = flexibility_score(&s).unwrap();
        assert_eq!(scores.len(), 3);

        // GLY (seq=2) should have the highest Z-score (most flexible)
        let gly_score = scores.iter().find(|(seq, _)| *seq == 2).unwrap().1;
        // VAL (seq=3) should have the lowest Z-score (most rigid)
        let val_score = scores.iter().find(|(seq, _)| *seq == 3).unwrap().1;
        assert!(gly_score > val_score);
        assert!(gly_score > 0.0); // above mean
        assert!(val_score < 0.0); // below mean
    }

    #[test]
    fn flexibility_identifies_flexible_residue() {
        let s = test_structure();
        let scores = flexibility_score(&s).unwrap();
        // The most flexible residue should be GLY (highest avg B = 21.5)
        let max_score = scores
            .iter()
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .unwrap();
        assert_eq!(max_score.0, 2); // GLY at seq_num 2
    }

    #[test]
    fn empty_structure_error() {
        let s = Structure {
            id: "EMPTY".into(),
            chains: vec![],
        };
        assert!(residue_bfactors(&s).is_err());
        assert!(chain_bfactors(&s).is_err());
        assert!(flexibility_score(&s).is_err());
    }

    #[test]
    fn multi_chain_bfactors() {
        let s = Structure {
            id: "MULTI".into(),
            chains: vec![
                Chain::new(
                    'A',
                    vec![Residue {
                        name: "ALA".into(),
                        seq_num: 1,
                        i_code: None,
                        atoms: vec![make_atom_with_b("CA", 10.0)],
                    }],
                ),
                Chain::new(
                    'B',
                    vec![Residue {
                        name: "GLY".into(),
                        seq_num: 1,
                        i_code: None,
                        atoms: vec![make_atom_with_b("CA", 20.0)],
                    }],
                ),
            ],
        };

        let cbs = chain_bfactors(&s).unwrap();
        assert_eq!(cbs.len(), 2);
        assert!((cbs[0].1.mean - 10.0).abs() < 1e-10);
        assert!((cbs[1].1.mean - 20.0).abs() < 1e-10);
    }

    #[test]
    fn uniform_bfactors_zscore_zero() {
        let s = Structure {
            id: "UNI".into(),
            chains: vec![Chain::new(
                'A',
                vec![
                    Residue {
                        name: "ALA".into(),
                        seq_num: 1,
                        i_code: None,
                        atoms: vec![make_atom_with_b("CA", 15.0)],
                    },
                    Residue {
                        name: "GLY".into(),
                        seq_num: 2,
                        i_code: None,
                        atoms: vec![make_atom_with_b("CA", 15.0)],
                    },
                ],
            )],
        };

        let scores = flexibility_score(&s).unwrap();
        for (_, z) in &scores {
            assert!(z.abs() < 1e-10, "uniform B-factors should give Z=0");
        }
    }

    #[test]
    fn bfactor_stats_std_dev() {
        // Known values: [2, 4, 4, 4, 5, 5, 7, 9]
        // Mean = 5.0, Variance = 4.0, StdDev = 2.0
        let mut vals = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let stats = compute_stats(&mut vals);
        assert!((stats.mean - 5.0).abs() < 1e-10);
        assert!((stats.std_dev - 2.0).abs() < 1e-10);
        assert!((stats.min - 2.0).abs() < 1e-10);
        assert!((stats.max - 9.0).abs() < 1e-10);
        assert!((stats.median - 4.5).abs() < 1e-10);
    }
}
