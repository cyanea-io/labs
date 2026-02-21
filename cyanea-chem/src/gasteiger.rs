//! Gasteiger-Marsili partial charge calculation.
//!
//! Iterative partial equalization of orbital electronegativity for computing
//! partial atomic charges. Used by force field electrostatic terms.

use cyanea_core::Result;

use crate::molecule::Molecule;

/// Electronegativity coefficients (a, b, c) for the equation:
/// χ = a + b·q + c·q²
/// where q is the partial charge on the atom.
///
/// Source: Gasteiger & Marsili, Tetrahedron 36, 3219 (1980).
struct ElectroParams {
    a: f64,
    b: f64,
    c: f64,
}

/// Get Gasteiger electronegativity parameters for a given atomic number and
/// bond-order-sum based hybridization state.
fn electro_params(atomic_number: u8, bond_order_sum: f64) -> ElectroParams {
    match atomic_number {
        // Hydrogen
        1 => ElectroParams { a: 7.17, b: 6.24, c: -0.56 },
        // Carbon
        6 => {
            if bond_order_sum > 3.5 {
                // sp (triple bond)
                ElectroParams { a: 10.39, b: 9.45, c: 0.73 }
            } else if bond_order_sum > 2.5 {
                // sp2 (double bond or aromatic)
                ElectroParams { a: 8.79, b: 9.32, c: 1.51 }
            } else {
                // sp3
                ElectroParams { a: 7.98, b: 9.18, c: 1.88 }
            }
        }
        // Nitrogen
        7 => {
            if bond_order_sum > 3.5 {
                ElectroParams { a: 15.68, b: 11.70, c: -0.27 }
            } else if bond_order_sum > 2.5 {
                ElectroParams { a: 12.87, b: 11.15, c: 0.85 }
            } else {
                ElectroParams { a: 11.54, b: 10.82, c: 1.36 }
            }
        }
        // Oxygen
        8 => {
            if bond_order_sum > 1.5 {
                // sp2 (C=O)
                ElectroParams { a: 17.07, b: 13.79, c: 0.47 }
            } else {
                // sp3 (C-O-H, C-O-C)
                ElectroParams { a: 14.18, b: 12.92, c: 1.39 }
            }
        }
        // Fluorine
        9 => ElectroParams { a: 14.66, b: 13.85, c: 2.31 },
        // Silicon
        14 => ElectroParams { a: 5.60, b: 6.00, c: 1.20 },
        // Phosphorus
        15 => ElectroParams { a: 8.90, b: 8.24, c: 0.96 },
        // Sulfur
        16 => {
            if bond_order_sum > 2.5 {
                ElectroParams { a: 12.00, b: 9.88, c: 1.58 }
            } else {
                ElectroParams { a: 10.14, b: 9.13, c: 1.38 }
            }
        }
        // Chlorine
        17 => ElectroParams { a: 11.00, b: 9.69, c: 1.35 },
        // Bromine
        35 => ElectroParams { a: 10.08, b: 8.47, c: 1.16 },
        // Iodine
        53 => ElectroParams { a: 9.90, b: 7.96, c: 0.96 },
        // Selenium
        34 => ElectroParams { a: 10.00, b: 8.80, c: 1.20 },
        // Default: use carbon sp3 as fallback
        _ => ElectroParams { a: 7.98, b: 9.18, c: 1.88 },
    }
}

/// Electronegativity at partial charge q.
fn electronegativity(params: &ElectroParams, q: f64) -> f64 {
    params.a + params.b * q + params.c * q * q
}

/// Compute Gasteiger-Marsili partial charges for all atoms in a molecule.
///
/// Uses 6 iterations of iterative partial equalization of orbital
/// electronegativity with a damping factor of 0.5^iteration.
///
/// # Example
///
/// ```
/// use cyanea_chem::{parse_smiles, gasteiger_charges};
///
/// let mol = parse_smiles("CCO").unwrap();
/// let charges = gasteiger_charges(&mol).unwrap();
/// assert_eq!(charges.len(), 3);
/// // Oxygen should be negative
/// assert!(charges[2] < 0.0);
/// ```
pub fn gasteiger_charges(mol: &Molecule) -> Result<Vec<f64>> {
    let n = mol.atom_count();
    if n == 0 {
        return Ok(Vec::new());
    }

    // Compute bond order sums for hybridization assignment
    let bond_order_sums = compute_bond_order_sums(mol);

    // Get electronegativity params per atom
    let params: Vec<ElectroParams> = (0..n)
        .map(|i| {
            let bos = bond_order_sums[i] + mol.atoms[i].implicit_hydrogens as f64;
            electro_params(mol.atoms[i].atomic_number, bos)
        })
        .collect();

    let mut charges = vec![0.0_f64; n];
    let n_iterations = 6;

    for iteration in 0..n_iterations {
        let damping = 0.5_f64.powi(iteration as i32 + 1);
        let mut delta = vec![0.0_f64; n];

        for bond in &mol.bonds {
            let a1 = bond.atom1;
            let a2 = bond.atom2;

            let chi1 = electronegativity(&params[a1], charges[a1]);
            let chi2 = electronegativity(&params[a2], charges[a2]);

            // Charge flows from less electronegative to more electronegative
            let diff = chi2 - chi1;

            // Scale by the more electronegative atom's parameters
            let scale = if diff > 0.0 {
                // a2 is more electronegative — charge flows a1 → a2
                params[a2].a + params[a2].b + params[a2].c
            } else {
                params[a1].a + params[a1].b + params[a1].c
            };

            if scale.abs() < 1e-12 {
                continue;
            }

            let transfer = damping * diff / scale;
            delta[a1] += transfer;
            delta[a2] -= transfer;
        }

        for i in 0..n {
            charges[i] += delta[i];
        }
    }

    // Handle formal charges — add them directly
    for i in 0..n {
        charges[i] += mol.atoms[i].formal_charge as f64;
    }

    Ok(charges)
}

/// Compute sum of bond orders for each atom (excluding implicit H).
fn compute_bond_order_sums(mol: &Molecule) -> Vec<f64> {
    let mut sums = vec![0.0_f64; mol.atom_count()];
    for bond in &mol.bonds {
        let order = bond.order.as_f64();
        sums[bond.atom1] += order;
        sums[bond.atom2] += order;
    }
    sums
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn charges_sum_near_zero_for_neutral() {
        let mol = parse_smiles("CCO").unwrap();
        let charges = gasteiger_charges(&mol).unwrap();
        let sum: f64 = charges.iter().sum();
        assert!(sum.abs() < 0.01, "charge sum = {sum}, expected ~0");
    }

    #[test]
    fn oxygen_is_negative() {
        let mol = parse_smiles("CCO").unwrap();
        let charges = gasteiger_charges(&mol).unwrap();
        // Oxygen is atom index 2
        assert!(charges[2] < 0.0, "O charge = {}, expected < 0", charges[2]);
    }

    #[test]
    fn methane_all_zero() {
        let mol = parse_smiles("C").unwrap();
        let charges = gasteiger_charges(&mol).unwrap();
        assert_eq!(charges.len(), 1);
        assert!(charges[0].abs() < 0.01, "C charge = {}", charges[0]);
    }

    #[test]
    fn acetic_acid_carbonyl_o_more_negative() {
        // CC(=O)O — atom 0=C, 1=C, 2=O(=), 3=O(-H)
        let mol = parse_smiles("CC(=O)O").unwrap();
        let charges = gasteiger_charges(&mol).unwrap();
        // Both oxygens should be negative
        assert!(charges[2] < 0.0, "carbonyl O = {}", charges[2]);
        assert!(charges[3] < 0.0, "hydroxyl O = {}", charges[3]);
    }

    #[test]
    fn empty_molecule() {
        use crate::molecule::Molecule;
        let mol = Molecule::new("empty".into(), vec![], vec![]);
        let charges = gasteiger_charges(&mol).unwrap();
        assert!(charges.is_empty());
    }
}
