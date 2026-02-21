//! Chemical reaction handling: SMIRKS application, enumeration, atom-atom mapping,
//! and retrosynthetic disconnection.
//!
//! Uses the SMARTS pattern matcher with atom-map extensions for reaction templates.

use cyanea_core::{CyaneaError, Result};

use crate::canon::canonical_smiles;
use crate::molecule::{Bond, BondOrder, BondStereo, MolAtom, Molecule};
use crate::smarts::{parse_smarts, smarts_find_all, SmartsPattern};
use crate::smiles::parse_smiles;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// A parsed chemical reaction (SMIRKS).
#[derive(Debug, Clone)]
pub struct Reaction {
    /// Reactant SMARTS pattern.
    pub reactant_pattern: SmartsPattern,
    /// Product SMARTS template.
    pub product_template: SmartsPattern,
    /// Atom map: (reactant_atom_idx, map_number) pairs.
    pub reactant_maps: Vec<(usize, u16)>,
    /// Atom map: (product_atom_idx, map_number) pairs.
    pub product_maps: Vec<(usize, u16)>,
}

/// A reaction product.
#[derive(Debug, Clone)]
pub struct ReactionProduct {
    /// The product molecule.
    pub molecule: Molecule,
    /// Canonical SMILES of the product.
    pub smiles: String,
}

/// Atom-atom mapping between reactant and product.
#[derive(Debug, Clone)]
pub struct AtomAtomMapping {
    /// Pairs of (reactant_atom_idx, product_atom_idx) for mapped atoms.
    pub mapping: Vec<(usize, usize)>,
    /// Reactant atom indices not mapped to any product atom.
    pub unmapped_reactant: Vec<usize>,
    /// Product atom indices not mapped to any reactant atom.
    pub unmapped_product: Vec<usize>,
}

/// A retrosynthetic disconnection.
#[derive(Debug, Clone)]
pub struct Disconnection {
    /// Name of the transform applied.
    pub transform_name: String,
    /// SMIRKS string (product>>reactant direction).
    pub smirks: String,
    /// Precursor molecules as canonical SMILES.
    pub precursors: Vec<String>,
}

// ---------------------------------------------------------------------------
// Reaction parsing
// ---------------------------------------------------------------------------

/// Parse a reaction SMIRKS string (`reactant>>product`).
///
/// # Example
///
/// ```
/// use cyanea_chem::reaction::parse_reaction;
///
/// let rxn = parse_reaction("[C:1][OH:2]>>[C:1][O:2]C").unwrap();
/// assert!(!rxn.reactant_maps.is_empty());
/// ```
pub fn parse_reaction(smirks: &str) -> Result<Reaction> {
    let parts: Vec<&str> = smirks.split(">>").collect();
    if parts.len() != 2 {
        return Err(CyaneaError::Parse(
            "SMIRKS must contain exactly one '>>' separator".into(),
        ));
    }

    let reactant_str = parts[0].trim();
    let product_str = parts[1].trim();

    if reactant_str.is_empty() || product_str.is_empty() {
        return Err(CyaneaError::Parse(
            "SMIRKS reactant and product sides must not be empty".into(),
        ));
    }

    let reactant_pattern = parse_smarts(reactant_str)?;
    let product_template = parse_smarts(product_str)?;

    let reactant_maps = reactant_pattern.atom_maps();
    let product_maps = product_template.atom_maps();

    Ok(Reaction {
        reactant_pattern,
        product_template,
        reactant_maps,
        product_maps,
    })
}

// ---------------------------------------------------------------------------
// Reaction application
// ---------------------------------------------------------------------------

/// Apply a reaction (SMIRKS) to a molecule.
///
/// Matches the reactant pattern against the molecule, then constructs the product
/// by modifying mapped atoms according to the product template.
///
/// # Example
///
/// ```
/// use cyanea_chem::reaction::{parse_reaction, apply_reaction};
/// use cyanea_chem::parse_smiles;
///
/// let mol = parse_smiles("CCO").unwrap();
/// let rxn = parse_reaction("[C:1][OH:2]>>[C:1][O:2]C").unwrap();
/// let products = apply_reaction(&mol, &rxn).unwrap();
/// assert!(!products.is_empty());
/// ```
pub fn apply_reaction(mol: &Molecule, reaction: &Reaction) -> Result<Vec<ReactionProduct>> {
    let matches = smarts_find_all(mol, &reaction.reactant_pattern);
    if matches.is_empty() {
        return Ok(Vec::new());
    }

    let mut products = Vec::new();
    let mut seen_smiles = std::collections::HashSet::new();

    for smatch in &matches {
        if let Some(product) = build_product(mol, reaction, &smatch.atom_mapping) {
            let smi = canonical_smiles(&product);
            if seen_smiles.insert(smi.clone()) {
                products.push(ReactionProduct {
                    molecule: product,
                    smiles: smi,
                });
            }
        }
    }

    Ok(products)
}

/// Build a product molecule from a reactant match.
fn build_product(
    mol: &Molecule,
    reaction: &Reaction,
    match_mapping: &[(usize, usize)],
) -> Option<Molecule> {
    // Map from pattern atom index to target atom index
    let pattern_to_target: std::collections::HashMap<usize, usize> =
        match_mapping.iter().copied().collect();

    // Map from atom-map number to target atom index
    let mut map_to_target = std::collections::HashMap::new();
    for &(pattern_idx, map_num) in &reaction.reactant_maps {
        if let Some(&target_idx) = pattern_to_target.get(&pattern_idx) {
            map_to_target.insert(map_num, target_idx);
        }
    }

    // Start with a copy of the molecule
    let mut atoms = mol.atoms.clone();
    let mut bonds = mol.bonds.clone();

    // Track which bonds to remove (by marking)
    let mut bonds_to_remove = std::collections::HashSet::new();

    // Process product template: modify mapped atoms
    for &(prod_idx, map_num) in &reaction.product_maps {
        if let Some(&target_idx) = map_to_target.get(&map_num) {
            let prod_atom = &reaction.product_template.atoms[prod_idx];
            // Apply product atom properties to the target atom
            apply_product_atom_props(&mut atoms[target_idx], prod_atom);
        }
    }

    // Process product bonds: check for new bonds between mapped atoms
    for bond in &reaction.product_template.bonds {
        let a1_map = reaction.product_maps.iter().find(|&&(idx, _)| idx == bond.atom1);
        let a2_map = reaction.product_maps.iter().find(|&&(idx, _)| idx == bond.atom2);

        if let (Some(&(_, map1)), Some(&(_, map2))) = (a1_map, a2_map) {
            if let (Some(&t1), Some(&t2)) = (map_to_target.get(&map1), map_to_target.get(&map2)) {
                // Check if this bond already exists
                let existing = mol.get_bond(t1, t2);
                let prod_order = smarts_bond_to_bond_order(&bond.expr);

                if existing.is_none() {
                    // Add new bond
                    bonds.push(Bond {
                        atom1: t1,
                        atom2: t2,
                        order: prod_order,
                        is_aromatic: prod_order == BondOrder::Aromatic,
                        stereo: BondStereo::None,
                    });
                } else if existing.map(|b| b.order) != Some(prod_order) {
                    // Modify existing bond order
                    if let Some(bi) = mol.adjacency[t1]
                        .iter()
                        .find(|&&(n, _)| n == t2)
                        .map(|&(_, bi)| bi)
                    {
                        bonds[bi] = Bond {
                            atom1: t1,
                            atom2: t2,
                            order: prod_order,
                            is_aromatic: prod_order == BondOrder::Aromatic,
                            stereo: BondStereo::None,
                        };
                    }
                }
            }
        }
    }

    // Check for bonds in reactant but not in product (bonds to remove)
    for bond in &reaction.reactant_pattern.bonds {
        let a1_map = reaction.reactant_maps.iter().find(|&&(idx, _)| idx == bond.atom1);
        let a2_map = reaction.reactant_maps.iter().find(|&&(idx, _)| idx == bond.atom2);

        if let (Some(&(_, map1)), Some(&(_, map2))) = (a1_map, a2_map) {
            // Check if this bond exists in product template
            let in_product = reaction.product_template.bonds.iter().any(|pb| {
                let pb1_map = reaction.product_maps.iter().find(|&&(idx, _)| idx == pb.atom1);
                let pb2_map = reaction.product_maps.iter().find(|&&(idx, _)| idx == pb.atom2);
                match (pb1_map, pb2_map) {
                    (Some(&(_, m1)), Some(&(_, m2))) => {
                        (m1 == map1 && m2 == map2) || (m1 == map2 && m2 == map1)
                    }
                    _ => false,
                }
            });

            if !in_product {
                // Remove this bond
                if let (Some(&t1), Some(&t2)) = (map_to_target.get(&map1), map_to_target.get(&map2)) {
                    for (bi, b) in bonds.iter().enumerate() {
                        if (b.atom1 == t1 && b.atom2 == t2) || (b.atom1 == t2 && b.atom2 == t1) {
                            bonds_to_remove.insert(bi);
                        }
                    }
                }
            }
        }
    }

    // Remove marked bonds
    let bonds: Vec<Bond> = bonds
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !bonds_to_remove.contains(i))
        .map(|(_, b)| b)
        .collect();

    // Recalculate implicit hydrogens
    let mut product = Molecule::new("product".into(), atoms, bonds);
    recalculate_implicit_hydrogens(&mut product);

    Some(product)
}

fn apply_product_atom_props(atom: &mut MolAtom, smarts_atom: &crate::smarts::SmartsAtom) {
    // Extract atomic number from the SMARTS expression if it specifies one
    if let Some(num) = extract_atomic_num(&smarts_atom.expr) {
        atom.atomic_number = num;
    }
}

fn extract_atomic_num(expr: &crate::smarts::AtomExpr) -> Option<u8> {
    use crate::smarts::{AtomExpr, AtomPrimitive};
    match expr {
        AtomExpr::Prim(AtomPrimitive::AtomicNum(n)) => Some(*n),
        AtomExpr::And(terms) => terms.iter().find_map(|t| extract_atomic_num(t)),
        AtomExpr::Or(terms) => terms.first().and_then(|t| extract_atomic_num(t)),
        _ => None,
    }
}

fn smarts_bond_to_bond_order(expr: &crate::smarts::BondExpr) -> BondOrder {
    use crate::smarts::BondExpr;
    match expr {
        BondExpr::Single => BondOrder::Single,
        BondExpr::Double => BondOrder::Double,
        BondExpr::Triple => BondOrder::Triple,
        BondExpr::Aromatic => BondOrder::Aromatic,
        BondExpr::Any => BondOrder::Single,
        _ => BondOrder::Single,
    }
}

fn recalculate_implicit_hydrogens(mol: &mut Molecule) {
    for i in 0..mol.atoms.len() {
        let atom = &mol.atoms[i];
        let valence = crate::element::element_by_number(atom.atomic_number)
            .map(|e| e.valence)
            .unwrap_or(0);

        let bond_sum: f64 = mol.adjacency[i]
            .iter()
            .map(|&(_, bi)| mol.bonds[bi].order.as_f64())
            .sum();

        let charge_adj = atom.formal_charge.abs() as f64;
        let target_valence = valence as f64 - charge_adj;
        let implicit_h = (target_valence - bond_sum).max(0.0) as u8;
        mol.atoms[i].implicit_hydrogens = implicit_h;
    }
}

// ---------------------------------------------------------------------------
// Reaction enumeration
// ---------------------------------------------------------------------------

/// Enumerate reaction products from multiple reactant sets.
///
/// For a multi-component reaction (reactant SMIRKS with `.` separator),
/// generates the Cartesian product of reactant sets and applies the reaction
/// to each combination. Results are deduplicated by canonical SMILES.
pub fn enumerate_reactions(
    reactant_sets: &[Vec<&str>],
    smirks: &str,
) -> Result<Vec<ReactionProduct>> {
    let reaction = parse_reaction(smirks)?;
    let mut all_products = Vec::new();
    let mut seen = std::collections::HashSet::new();

    if reactant_sets.len() == 1 {
        // Single-component reaction
        for smi in &reactant_sets[0] {
            let mol = parse_smiles(smi)?;
            let products = apply_reaction(&mol, &reaction)?;
            for p in products {
                if seen.insert(p.smiles.clone()) {
                    all_products.push(p);
                }
            }
        }
    } else if reactant_sets.len() == 2 {
        // Two-component: combine each pair
        for smi1 in &reactant_sets[0] {
            for smi2 in &reactant_sets[1] {
                // Combine as `smi1.smi2`
                let combined = format!("{}.{}", smi1, smi2);
                if let Ok(mol) = parse_smiles(&combined) {
                    let products = apply_reaction(&mol, &reaction)?;
                    for p in products {
                        if seen.insert(p.smiles.clone()) {
                            all_products.push(p);
                        }
                    }
                }
            }
        }
    } else {
        // General multi-component: flatten into pairwise
        let mut current: Vec<String> = reactant_sets[0].iter().map(|s| s.to_string()).collect();
        for set in &reactant_sets[1..] {
            let mut next = Vec::new();
            for existing in &current {
                for smi in set.iter() {
                    next.push(format!("{}.{}", existing, smi));
                }
            }
            current = next;
        }
        for combined in &current {
            if let Ok(mol) = parse_smiles(combined) {
                let products = apply_reaction(&mol, &reaction)?;
                for p in products {
                    if seen.insert(p.smiles.clone()) {
                        all_products.push(p);
                    }
                }
            }
        }
    }

    Ok(all_products)
}

// ---------------------------------------------------------------------------
// Atom-atom mapping
// ---------------------------------------------------------------------------

/// Compute atom-atom mapping between a reactant and product molecule.
///
/// Uses maximum common substructure (MCS) to find corresponding atoms,
/// then heuristically assigns remaining unmapped atoms.
pub fn atom_atom_map(reactant: &Molecule, product: &Molecule) -> AtomAtomMapping {
    let mcs = crate::scaffold::maximum_common_substructure(reactant, product, 1000)
        .unwrap_or_else(|_| crate::scaffold::McsResult {
            mcs_atoms: 0,
            mcs_bonds: 0,
            mapping_a: Vec::new(),
            mapping_b: Vec::new(),
        });

    let mut mapping: Vec<(usize, usize)> = Vec::new();
    let mut mapped_reactant = std::collections::HashSet::new();
    let mut mapped_product = std::collections::HashSet::new();

    // mapping_a: Vec<(mol_a_atom, node_idx)>, mapping_b: Vec<(mol_b_atom, node_idx)>
    // Match by node_idx to get reactant↔product atom pairs
    for &(a_atom, node_idx) in &mcs.mapping_a {
        if let Some(&(b_atom, _)) = mcs.mapping_b.iter().find(|&&(_, ni)| ni == node_idx) {
            mapping.push((a_atom, b_atom));
            mapped_reactant.insert(a_atom);
            mapped_product.insert(b_atom);
        }
    }

    // Heuristic: map remaining atoms by element
    let unmapped_r: Vec<usize> = (0..reactant.atom_count())
        .filter(|i| !mapped_reactant.contains(i))
        .collect();
    let unmapped_p: Vec<usize> = (0..product.atom_count())
        .filter(|i| !mapped_product.contains(i))
        .collect();

    let mut still_unmapped_r = Vec::new();
    let mut still_unmapped_p: Vec<usize> = unmapped_p.clone();

    for &ri in &unmapped_r {
        let r_elem = reactant.atoms[ri].atomic_number;
        if let Some(pos) = still_unmapped_p.iter().position(|&pi| {
            product.atoms[pi].atomic_number == r_elem
        }) {
            let pi = still_unmapped_p.remove(pos);
            mapping.push((ri, pi));
        } else {
            still_unmapped_r.push(ri);
        }
    }

    AtomAtomMapping {
        mapping,
        unmapped_reactant: still_unmapped_r,
        unmapped_product: still_unmapped_p,
    }
}

// ---------------------------------------------------------------------------
// Retrosynthetic disconnection
// ---------------------------------------------------------------------------

/// Common retrosynthetic transforms as (name, forward SMIRKS).
/// For retrosynthesis, we apply the reverse (product → reactants).
const RETRO_TRANSFORMS: &[(&str, &str)] = &[
    ("Amide bond formation", "[C:1](=O)[NH:2]>>[C:1](=O)O.[NH2:2]"),
    ("Ester hydrolysis", "[C:1](=O)[O:2]>>[C:1](=O)O.[OH:2]"),
    ("Suzuki coupling", "[c:1][c:2]>>[c:1]B(O)O.[c:2]Br"),
    ("Reductive amination", "[C:1][NH:2]>>[C:1]=O.[NH2:2]"),
    ("Williamson ether", "[C:1][O:2][C:3]>>[C:1]Br.[OH:2][C:3]"),
    ("Fischer esterification", "[C:1](=O)[O:2][C:3]>>[C:1](=O)O.[OH:2][C:3]"),
    ("N-alkylation", "[N:1][C:2]>>[NH:1].[C:2]Br"),
    ("Grignard addition", "[C:1]([OH:2])[C:3]>>[C:1]=O.[C:3]MgBr"),
    ("Michael addition", "[C:1][C:2][C:3](=O)>>[C:1]=[C:2].[C:3](=O)"),
    ("Aldol condensation", "[C:1]([OH:2])[C:3][C:4](=O)>>[C:1](=O).[C:3][C:4](=O)"),
    ("Wittig olefination", "[C:1]=[C:2]>>[C:1]=O.[C:2]"),
    ("Diels-Alder", "[C:1]1[C:2]=[C:3][C:4][C:5]=[C:6]1>>[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]"),
    ("Heck coupling", "[c:1][C:2]=[C:3]>>[c:1]Br.[C:2]=[C:3]"),
    ("Buchwald-Hartwig", "[c:1][N:2]>>[c:1]Br.[NH:2]"),
    ("Sonogashira coupling", "[c:1][C:2]#[C:3]>>[c:1]Br.[C:2]#[C:3]"),
];

/// Find retrosynthetic disconnections for a target molecule.
///
/// Applies common transform rules in reverse to identify potential precursors.
///
/// # Example
///
/// ```
/// use cyanea_chem::reaction::retrosynthetic_disconnections;
/// use cyanea_chem::parse_smiles;
///
/// let mol = parse_smiles("CC(=O)NC").unwrap(); // N-methylacetamide
/// let disconnections = retrosynthetic_disconnections(&mol);
/// // Should find at least the amide disconnection
/// ```
pub fn retrosynthetic_disconnections(mol: &Molecule) -> Vec<Disconnection> {
    let mut results = Vec::new();

    for &(name, smirks) in RETRO_TRANSFORMS {
        // Parse the forward reaction and apply the reactant pattern (which represents
        // the substructure we're looking for in the target)
        let parts: Vec<&str> = smirks.split(">>").collect();
        if parts.len() != 2 {
            continue;
        }

        // For retrosynthesis: we match the reactant side (left of >>) against our molecule,
        // and the products (right of >>) are the precursors
        let reactant_smarts = parts[0];
        let precursor_smiles_template = parts[1];

        let pattern = match parse_smarts(reactant_smarts) {
            Ok(p) => p,
            Err(_) => continue,
        };

        let matches = smarts_find_all(mol, &pattern);
        if matches.is_empty() {
            continue;
        }

        // Collect precursor SMILES from the template
        let precursors: Vec<String> = precursor_smiles_template
            .split('.')
            .filter(|s| !s.is_empty())
            .map(|s| {
                // Strip atom map numbers for cleaner output
                strip_atom_maps(s)
            })
            .collect();

        if !precursors.is_empty() {
            results.push(Disconnection {
                transform_name: name.to_string(),
                smirks: smirks.to_string(),
                precursors,
            });
        }
    }

    results
}

/// Strip atom map numbers (`:N`) from SMILES-like strings for cleaner output.
fn strip_atom_maps(s: &str) -> String {
    let mut result = String::with_capacity(s.len());
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b':' && i > 0 {
            // Skip `:digits`
            i += 1;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
        } else {
            result.push(bytes[i] as char);
            i += 1;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn parse_reaction_basic() {
        let rxn = parse_reaction("[C:1][O:2]>>[C:1]=[O:2]").unwrap();
        assert!(!rxn.reactant_maps.is_empty());
        assert!(!rxn.product_maps.is_empty());
    }

    #[test]
    fn parse_reaction_invalid() {
        assert!(parse_reaction("").is_err());
        assert!(parse_reaction("CC").is_err()); // no >>
        assert!(parse_reaction(">>CC").is_err()); // empty reactant
        assert!(parse_reaction("CC>>").is_err()); // empty product
    }

    #[test]
    fn apply_reaction_basic() {
        // Oxidation: C-OH → C=O
        let mol = parse_smiles("CCO").unwrap();
        let rxn = parse_reaction("[C:1][OH:2]>>[C:1]=[O:2]").unwrap();
        let products = apply_reaction(&mol, &rxn).unwrap();
        // Should produce at least one product
        assert!(!products.is_empty(), "expected products");
    }

    #[test]
    fn apply_reaction_no_match() {
        let mol = parse_smiles("CCC").unwrap(); // no OH group
        let rxn = parse_reaction("[C:1][OH:2]>>[C:1]=[O:2]").unwrap();
        let products = apply_reaction(&mol, &rxn).unwrap();
        assert!(products.is_empty());
    }

    #[test]
    fn apply_reaction_dedup() {
        // Symmetric molecule — same product from multiple matches
        let mol = parse_smiles("OCC(O)CCO").unwrap();
        let rxn = parse_reaction("[C:1][OH:2]>>[C:1]=[O:2]").unwrap();
        let products = apply_reaction(&mol, &rxn).unwrap();
        // Products should be deduplicated by canonical SMILES
        let unique: std::collections::HashSet<_> = products.iter().map(|p| &p.smiles).collect();
        assert_eq!(unique.len(), products.len());
    }

    #[test]
    fn enumerate_single_component() {
        let reactants = vec![vec!["CCO", "CCCO"]];
        let smirks = "[C:1][OH:2]>>[C:1]=[O:2]";
        let products = enumerate_reactions(&reactants, smirks).unwrap();
        // Should get products from both reactants
        assert!(!products.is_empty());
    }

    #[test]
    fn atom_atom_map_basic() {
        let reactant = parse_smiles("CCO").unwrap();
        let product = parse_smiles("CC=O").unwrap();
        let mapping = atom_atom_map(&reactant, &product);
        assert!(!mapping.mapping.is_empty());
    }

    #[test]
    fn atom_atom_map_identical() {
        let mol = parse_smiles("CCO").unwrap();
        let mapping = atom_atom_map(&mol, &mol);
        assert_eq!(mapping.mapping.len(), 3);
        assert!(mapping.unmapped_reactant.is_empty());
        assert!(mapping.unmapped_product.is_empty());
    }

    #[test]
    fn retrosynthetic_amide() {
        let mol = parse_smiles("CC(=O)NC").unwrap();
        let disconnections = retrosynthetic_disconnections(&mol);
        let names: Vec<&str> = disconnections.iter().map(|d| d.transform_name.as_str()).collect();
        assert!(
            names.contains(&"Amide bond formation"),
            "expected amide disconnection, got {:?}",
            names
        );
    }

    #[test]
    fn retrosynthetic_ester() {
        let mol = parse_smiles("CC(=O)OC").unwrap();
        let disconnections = retrosynthetic_disconnections(&mol);
        let names: Vec<&str> = disconnections.iter().map(|d| d.transform_name.as_str()).collect();
        // Should find ester hydrolysis or Fischer esterification
        assert!(
            names.iter().any(|n| n.contains("ster")),
            "expected ester disconnection, got {:?}",
            names
        );
    }

    #[test]
    fn retrosynthetic_no_match() {
        let mol = parse_smiles("C").unwrap(); // methane — no disconnections
        let disconnections = retrosynthetic_disconnections(&mol);
        assert!(disconnections.is_empty());
    }

    #[test]
    fn retrosynthetic_has_precursors() {
        let mol = parse_smiles("CC(=O)NC").unwrap();
        let disconnections = retrosynthetic_disconnections(&mol);
        for d in &disconnections {
            assert!(!d.precursors.is_empty(), "disconnection {} has no precursors", d.transform_name);
        }
    }

    #[test]
    fn strip_atom_maps_basic() {
        assert_eq!(strip_atom_maps("[C:1](=O)O"), "[C](=O)O");
        assert_eq!(strip_atom_maps("[NH2:2]"), "[NH2]");
        assert_eq!(strip_atom_maps("CC"), "CC");
    }

    #[test]
    fn reaction_product_has_smiles() {
        let mol = parse_smiles("CCO").unwrap();
        let rxn = parse_reaction("[C:1][OH:2]>>[C:1]=[O:2]").unwrap();
        let products = apply_reaction(&mol, &rxn).unwrap();
        for p in &products {
            assert!(!p.smiles.is_empty());
        }
    }

    #[test]
    fn parse_reaction_maps() {
        let rxn = parse_reaction("[C:1][N:2]>>[C:1]=[N:2]").unwrap();
        assert_eq!(rxn.reactant_maps.len(), 2);
        assert_eq!(rxn.product_maps.len(), 2);
        // Map numbers should match
        let r_maps: Vec<u16> = rxn.reactant_maps.iter().map(|&(_, m)| m).collect();
        let p_maps: Vec<u16> = rxn.product_maps.iter().map(|&(_, m)| m).collect();
        assert_eq!(r_maps, p_maps);
    }

    #[test]
    fn retrosynthetic_ether() {
        let mol = parse_smiles("COCC").unwrap();
        let disconnections = retrosynthetic_disconnections(&mol);
        // Should find Williamson ether synthesis
        let names: Vec<&str> = disconnections.iter().map(|d| d.transform_name.as_str()).collect();
        assert!(
            names.iter().any(|n| n.contains("ether") || n.contains("Williamson")),
            "expected ether disconnection, got {:?}",
            names
        );
    }

    #[test]
    fn retrosynthetic_secondary_amine() {
        let mol = parse_smiles("CNCC").unwrap();
        let disconnections = retrosynthetic_disconnections(&mol);
        // Should find reductive amination or N-alkylation
        assert!(!disconnections.is_empty(), "expected at least one disconnection");
    }
}
