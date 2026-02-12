//! MACCS 166-key structural fingerprints.
//!
//! MACCS (Molecular ACCess System) keys are one of the most widely used structural
//! fingerprint types. Each of the 166 bit positions corresponds to the presence or
//! absence of a specific structural feature: element presence, ring topology, bond
//! types, functional groups, and atom/ring count thresholds.
//!
//! This implementation evaluates keys directly from the molecular graph using
//! element counts, ring detection, bond iteration, and local neighborhood checks.
//! Complex SMARTS-based keys that would require full substructure matching are
//! approximated where possible or left unset.

use crate::fingerprint::Fingerprint;
use crate::molecule::{BondOrder, Molecule};
use crate::ring;

// Atomic numbers for convenience.
const H: u8 = 1;
const B: u8 = 5;
const C: u8 = 6;
const N: u8 = 7;
const O: u8 = 8;
const F: u8 = 9;
const SI: u8 = 14;
const P: u8 = 15;
const S: u8 = 16;
const CL: u8 = 17;
const K: u8 = 19;
const CA: u8 = 20;
const AS: u8 = 33;
const SE: u8 = 34;
const BR: u8 = 35;
const LI: u8 = 3;
const NA: u8 = 11;
const AL: u8 = 13;
const I: u8 = 53;
const FE: u8 = 26;
const ZN: u8 = 30;
const MG: u8 = 12;
const CU: u8 = 29;
const MN: u8 = 25;
const CO: u8 = 27;
const GE: u8 = 32;
const SN: u8 = 50;
const TE: u8 = 52;

/// Compute a MACCS-like 166-key structural fingerprint.
///
/// Each bit position corresponds to a specific structural feature from the MACCS
/// key definitions. Keys cover element presence, ring topology, bond types,
/// functional groups, and count-based thresholds. This is a practical implementation
/// that evaluates ~80+ of the 166 keys using the molecular graph directly.
///
/// # Example
///
/// ```
/// use cyanea_chem::{parse_smiles, maccs_fingerprint, tanimoto_similarity};
///
/// let aspirin = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
/// let fp = maccs_fingerprint(&aspirin);
/// assert_eq!(fp.nbits(), 166);
/// assert!(fp.count_ones() > 10); // aspirin has many features
/// ```
pub fn maccs_fingerprint(mol: &Molecule) -> Fingerprint {
    let mut fp = Fingerprint::new(166);

    if mol.atom_count() == 0 {
        return fp;
    }

    let rings = ring::find_sssr(mol);

    // Element counts indexed by atomic number (0..120).
    let mut element_counts = [0u32; 120];
    for atom in &mol.atoms {
        if (atom.atomic_number as usize) < 120 {
            element_counts[atom.atomic_number as usize] += 1;
        }
    }

    set_element_keys(&mut fp, &element_counts, mol);
    set_ring_keys(&mut fp, &rings, mol);
    set_bond_keys(&mut fp, mol);
    set_functional_group_keys(&mut fp, mol);
    set_count_keys(&mut fp, &element_counts, &rings, mol);

    fp
}

// ---------------------------------------------------------------------------
// Element presence keys
// ---------------------------------------------------------------------------

/// Set fingerprint bits for element-presence MACCS keys.
fn set_element_keys(fp: &mut Fingerprint, counts: &[u32; 120], mol: &Molecule) {
    // Key 103: Group IIIA atom (B, Al, Ga, In, Tl) present
    if counts[B as usize] > 0 || counts[AL as usize] > 0 {
        fp.set_bit(103);
    }

    // Key 110: Isotope present
    if mol.atoms.iter().any(|a| a.isotope.is_some()) {
        fp.set_bit(110);
    }

    // Key 111: >= 4 heavy atoms
    if mol.heavy_atom_count() >= 4 {
        fp.set_bit(111);
    }

    // Key 112: >= 8 heavy atoms
    if mol.heavy_atom_count() >= 8 {
        fp.set_bit(112);
    }

    // Key 113: >= 16 heavy atoms
    if mol.heavy_atom_count() >= 16 {
        fp.set_bit(113);
    }

    // Key 114: >= 32 heavy atoms
    if mol.heavy_atom_count() >= 32 {
        fp.set_bit(114);
    }

    // Key 115: Phosphorus present
    if counts[P as usize] > 0 {
        fp.set_bit(115);
    }

    // Key 116: Sulfur present
    if counts[S as usize] > 0 {
        fp.set_bit(116);
    }

    // Key 117: Silicon present
    if counts[SI as usize] > 0 {
        fp.set_bit(117);
    }

    // Key 118: Oxygen present
    if counts[O as usize] > 0 {
        fp.set_bit(118);
    }

    // Key 119: Nitrogen present
    if counts[N as usize] > 0 {
        fp.set_bit(119);
    }

    // Key 120: Halogen present (F, Cl, Br, I)
    if counts[F as usize] > 0
        || counts[CL as usize] > 0
        || counts[BR as usize] > 0
        || counts[I as usize] > 0
    {
        fp.set_bit(120);
    }

    // Key 122: At least 2 different elements (beyond H)
    let distinct_heavy: u32 = counts[2..120]
        .iter()
        .filter(|&&c| c > 0)
        .count() as u32;
    if distinct_heavy >= 2 {
        fp.set_bit(122);
    }

    // Key 123: At least 3 different elements (beyond H)
    if distinct_heavy >= 3 {
        fp.set_bit(123);
    }

    // Key 124: At least 4 different elements (beyond H)
    if distinct_heavy >= 4 {
        fp.set_bit(124);
    }

    // Key 126: Carbon present
    if counts[C as usize] > 0 {
        fp.set_bit(126);
    }

    // Key 127: Sodium, potassium, or other alkali metal present
    if counts[LI as usize] > 0 || counts[NA as usize] > 0 || counts[K as usize] > 0 {
        fp.set_bit(127);
    }

    // Key 128: Alkaline earth present (Mg, Ca)
    if counts[MG as usize] > 0 || counts[CA as usize] > 0 {
        fp.set_bit(128);
    }

    // Key 129: Transition metal present (Fe, Cu, Zn, Mn, Co)
    if counts[FE as usize] > 0
        || counts[CU as usize] > 0
        || counts[ZN as usize] > 0
        || counts[MN as usize] > 0
        || counts[CO as usize] > 0
    {
        fp.set_bit(129);
    }

    // Key 130: Germanium or Tin present
    if counts[GE as usize] > 0 || counts[SN as usize] > 0 {
        fp.set_bit(130);
    }

    // Key 131: Selenium or Tellurium present
    if counts[SE as usize] > 0 || counts[TE as usize] > 0 {
        fp.set_bit(131);
    }

    // Key 132: Arsenic present
    if counts[AS as usize] > 0 {
        fp.set_bit(132);
    }

    // Key 133: >= 1 hydrogen (implicit or explicit)
    if counts[H as usize] > 0 || mol.atoms.iter().any(|a| a.implicit_hydrogens > 0) {
        fp.set_bit(133);
    }

    // Key 134: >= 2 carbons
    if counts[C as usize] >= 2 {
        fp.set_bit(134);
    }

    // Key 135: >= 4 carbons
    if counts[C as usize] >= 4 {
        fp.set_bit(135);
    }

    // Key 136: >= 8 carbons
    if counts[C as usize] >= 8 {
        fp.set_bit(136);
    }

    // Key 137: >= 16 carbons
    if counts[C as usize] >= 16 {
        fp.set_bit(137);
    }

    // Key 138: >= 2 heteroatoms (non-C, non-H heavy atoms)
    let hetero_count: u32 = (2..120)
        .filter(|&i| i != C as usize)
        .map(|i| counts[i])
        .sum();
    if hetero_count >= 2 {
        fp.set_bit(138);
    }

    // Key 139: Chlorine present
    if counts[CL as usize] > 0 {
        fp.set_bit(139);
    }

    // Key 140: Bromine present
    if counts[BR as usize] > 0 {
        fp.set_bit(140);
    }

    // Key 141: Iodine present
    if counts[I as usize] > 0 {
        fp.set_bit(141);
    }

    // Key 142: Fluorine present
    if counts[F as usize] > 0 {
        fp.set_bit(142);
    }
}

// ---------------------------------------------------------------------------
// Ring keys
// ---------------------------------------------------------------------------

/// Set fingerprint bits for ring-topology MACCS keys.
fn set_ring_keys(fp: &mut Fingerprint, rings: &[Vec<usize>], mol: &Molecule) {
    // Key 125: Any ring present
    if !rings.is_empty() {
        fp.set_bit(125);
    }

    for ring in rings {
        let size = ring.len();
        match size {
            3 => fp.set_bit(145), // Key 145: 3-membered ring
            4 => fp.set_bit(146), // Key 146: 4-membered ring
            5 => fp.set_bit(147), // Key 147: 5-membered ring
            6 => fp.set_bit(148), // Key 148: 6-membered ring
            _ if size >= 7 => fp.set_bit(149), // Key 149: 7+ membered ring
            _ => {}
        }

        // Check if ring is aromatic (all atoms in ring are aromatic)
        let all_aromatic = ring.iter().all(|&idx| mol.atoms[idx].is_aromatic);
        if all_aromatic {
            fp.set_bit(162); // Key 162: Aromatic ring present
        }

        // Key 165: Ring contains nitrogen
        if ring.iter().any(|&idx| mol.atoms[idx].atomic_number == N) {
            fp.set_bit(165);
        }

        // Check for heteroatom in ring
        if ring
            .iter()
            .any(|&idx| mol.atoms[idx].atomic_number != C && mol.atoms[idx].atomic_number != H)
        {
            fp.set_bit(102); // Key 102: Heterocyclic ring
        }

        // Check for oxygen in ring
        if ring.iter().any(|&idx| mol.atoms[idx].atomic_number == O) {
            fp.set_bit(101); // Key 101: Ring with oxygen
        }

        // Check for sulfur in ring
        if ring.iter().any(|&idx| mol.atoms[idx].atomic_number == S) {
            fp.set_bit(100); // Key 100: Ring with sulfur
        }
    }
}

// ---------------------------------------------------------------------------
// Bond type keys
// ---------------------------------------------------------------------------

/// Set fingerprint bits for bond-type MACCS keys.
fn set_bond_keys(fp: &mut Fingerprint, mol: &Molecule) {
    let mut has_double = false;
    let mut has_triple = false;
    let mut has_aromatic_bond = false;
    let mut has_c_double_c = false;
    let mut has_c_single_n = false;
    let mut has_c_single_o = false;
    let mut has_c_double_o = false;
    let mut has_c_double_n = false;
    let mut has_n_single_o = false;
    let mut has_n_double_o = false;
    let mut has_s_double_o = false;
    let mut has_c_single_s = false;
    let mut has_c_single_c = false;
    let mut has_c_triple_c = false;
    let mut has_c_triple_n = false;
    let mut has_n_single_n = false;
    let mut has_o_single_o = false;
    let mut has_s_single_s = false;
    let mut has_c_single_halogen = false;
    let mut has_c_single_f = false;
    let mut has_c_single_cl = false;
    let mut has_c_single_br = false;
    let mut has_c_single_i = false;
    let mut has_p_single_o = false;
    let mut has_p_double_o = false;
    let mut has_c_single_p = false;
    let mut has_n_single_s = false;
    let mut has_c_single_si = false;

    for bond in &mol.bonds {
        let a1 = mol.atoms[bond.atom1].atomic_number;
        let a2 = mol.atoms[bond.atom2].atomic_number;
        let (lo, hi) = if a1 <= a2 { (a1, a2) } else { (a2, a1) };

        match bond.order {
            BondOrder::Double => {
                has_double = true;
                if lo == C && hi == C {
                    has_c_double_c = true;
                }
                if (lo == C && hi == O) || (lo == O && hi == C) {
                    has_c_double_o = true;
                }
                if (lo == C && hi == N) || (lo == N && hi == C) {
                    has_c_double_n = true;
                }
                if (lo == N && hi == O) || (lo == O && hi == N) {
                    has_n_double_o = true;
                }
                if (lo == O && hi == S) || (lo == S && hi == O) {
                    has_s_double_o = true;
                }
                if (lo == O && hi == P) || (lo == P && hi == O) {
                    has_p_double_o = true;
                }
            }
            BondOrder::Triple => {
                has_triple = true;
                if lo == C && hi == C {
                    has_c_triple_c = true;
                }
                if (lo == C && hi == N) || (lo == N && hi == C) {
                    has_c_triple_n = true;
                }
            }
            BondOrder::Aromatic => {
                has_aromatic_bond = true;
            }
            BondOrder::Single => {
                if lo == C && hi == C {
                    has_c_single_c = true;
                }
                if (lo == C && hi == N) || (lo == N && hi == C) {
                    has_c_single_n = true;
                }
                if (lo == C && hi == O) || (lo == O && hi == C) {
                    has_c_single_o = true;
                }
                if (lo == N && hi == O) || (lo == O && hi == N) {
                    has_n_single_o = true;
                }
                if (lo == C && hi == S) || (lo == S && hi == C) {
                    has_c_single_s = true;
                }
                if lo == N && hi == N {
                    has_n_single_n = true;
                }
                if lo == O && hi == O {
                    has_o_single_o = true;
                }
                if lo == S && hi == S {
                    has_s_single_s = true;
                }
                if (lo == C && hi == F) || (lo == F && hi == C) {
                    has_c_single_f = true;
                    has_c_single_halogen = true;
                }
                if (lo == C && hi == CL) || (lo == CL && hi == C) {
                    has_c_single_cl = true;
                    has_c_single_halogen = true;
                }
                if (lo == C && hi == BR) || (lo == BR && hi == C) {
                    has_c_single_br = true;
                    has_c_single_halogen = true;
                }
                if (lo == C && hi == I) || (lo == I && hi == C) {
                    has_c_single_i = true;
                    has_c_single_halogen = true;
                }
                if (lo == O && hi == P) || (lo == P && hi == O) {
                    has_p_single_o = true;
                }
                if (lo == C && hi == P) || (lo == P && hi == C) {
                    has_c_single_p = true;
                }
                if (lo == N && hi == S) || (lo == S && hi == N) {
                    has_n_single_s = true;
                }
                if (lo == C && hi == SI) || (lo == SI && hi == C) {
                    has_c_single_si = true;
                }
            }
        }
    }

    // Key 143: Double bond present
    if has_double {
        fp.set_bit(143);
    }
    // Key 144: Triple bond present
    if has_triple {
        fp.set_bit(144);
    }
    // Key 163: Aromatic bond present
    if has_aromatic_bond {
        fp.set_bit(163);
    }

    // Bond-pair keys
    // Key 153: C=O present
    if has_c_double_o {
        fp.set_bit(153);
    }
    // Key 154: C-N present
    if has_c_single_n {
        fp.set_bit(154);
    }
    // Key 158: N-O present (single bond)
    if has_n_single_o {
        fp.set_bit(158);
    }

    // Additional bond keys
    // Key 1: C=C
    if has_c_double_c {
        fp.set_bit(1);
    }
    // Key 2: C#C
    if has_c_triple_c {
        fp.set_bit(2);
    }
    // Key 3: C#N
    if has_c_triple_n {
        fp.set_bit(3);
    }
    // Key 4: C=N
    if has_c_double_n {
        fp.set_bit(4);
    }
    // Key 5: N=O
    if has_n_double_o {
        fp.set_bit(5);
    }
    // Key 6: S=O
    if has_s_double_o {
        fp.set_bit(6);
    }
    // Key 7: C-O
    if has_c_single_o {
        fp.set_bit(7);
    }
    // Key 8: C-S
    if has_c_single_s {
        fp.set_bit(8);
    }
    // Key 9: C-C
    if has_c_single_c {
        fp.set_bit(9);
    }
    // Key 10: N-N
    if has_n_single_n {
        fp.set_bit(10);
    }
    // Key 11: O-O
    if has_o_single_o {
        fp.set_bit(11);
    }
    // Key 12: S-S
    if has_s_single_s {
        fp.set_bit(12);
    }
    // Key 13: C-halogen
    if has_c_single_halogen {
        fp.set_bit(13);
    }
    // Key 14: C-F
    if has_c_single_f {
        fp.set_bit(14);
    }
    // Key 15: C-Cl
    if has_c_single_cl {
        fp.set_bit(15);
    }
    // Key 16: C-Br
    if has_c_single_br {
        fp.set_bit(16);
    }
    // Key 17: C-I
    if has_c_single_i {
        fp.set_bit(17);
    }
    // Key 18: P-O
    if has_p_single_o || has_p_double_o {
        fp.set_bit(18);
    }
    // Key 19: C-P
    if has_c_single_p {
        fp.set_bit(19);
    }
    // Key 20: N-S
    if has_n_single_s {
        fp.set_bit(20);
    }
    // Key 21: C-Si
    if has_c_single_si {
        fp.set_bit(21);
    }
}

// ---------------------------------------------------------------------------
// Functional group keys (local neighborhood patterns)
// ---------------------------------------------------------------------------

/// Set fingerprint bits for functional-group MACCS keys, using local
/// atom-neighborhood inspection rather than full substructure search.
fn set_functional_group_keys(fp: &mut Fingerprint, mol: &Molecule) {
    let mut has_carboxyl = false; // C(=O)O pattern (acid, ester)
    let mut has_amide = false; // C(=O)N pattern
    let mut has_sulfonyl = false; // S(=O) pattern
    let mut has_hydroxyl = false; // O-H (O with implicit H and single bond to heavy atom)
    let mut has_amine_nh = false; // N-H (N with implicit H)
    let mut has_aldehyde = false; // H-C(=O) (C with implicit H and C=O)
    let mut has_ester = false; // C(=O)O-C
    let mut has_nitro = false; // N(=O)O or N(=O)(=O)
    let mut has_thiol = false; // S-H
    let mut has_phosphate = false; // P(=O)(O)(O)
    let mut has_c_o_c = false; // ether: C-O-C
    let mut has_c_n_c = false; // tertiary amine-like: C-N-C
    let mut has_charged_n = false; // positively charged nitrogen
    let mut has_charged_o = false; // negatively charged oxygen

    for (i, atom) in mol.atoms.iter().enumerate() {
        let an = atom.atomic_number;

        // Check for charged atoms
        if an == N && atom.formal_charge > 0 {
            has_charged_n = true;
        }
        if an == O && atom.formal_charge < 0 {
            has_charged_o = true;
        }

        // Hydroxyl: oxygen with at least 1 implicit H, bonded to non-H atom
        if an == O && atom.implicit_hydrogens >= 1 {
            has_hydroxyl = true;
        }

        // Amine N-H: nitrogen with at least 1 implicit H
        if an == N && atom.implicit_hydrogens >= 1 {
            has_amine_nh = true;
        }

        // Thiol: sulfur with at least 1 implicit H
        if an == S && atom.implicit_hydrogens >= 1 {
            has_thiol = true;
        }

        // Patterns centered on carbon: C(=O)O, C(=O)N, aldehyde
        if an == C {
            let neighbors = &mol.adjacency[i];
            let mut has_double_o = false;
            let mut has_single_o = false;
            let mut has_single_n = false;
            let mut single_o_idx: Option<usize> = None;

            for &(nbr, bond_idx) in neighbors {
                let nbr_an = mol.atoms[nbr].atomic_number;
                let bond = &mol.bonds[bond_idx];
                match (nbr_an, bond.order) {
                    (O, BondOrder::Double) => has_double_o = true,
                    (O, BondOrder::Single) => {
                        has_single_o = true;
                        single_o_idx = Some(nbr);
                    }
                    (N, BondOrder::Single) => has_single_n = true,
                    _ => {}
                }
            }

            if has_double_o && has_single_o {
                has_carboxyl = true;
                // Check if the single O is connected to C (ester) vs H (acid)
                if let Some(o_idx) = single_o_idx {
                    for &(o_nbr, _) in &mol.adjacency[o_idx] {
                        if o_nbr != i && mol.atoms[o_nbr].atomic_number == C {
                            has_ester = true;
                        }
                    }
                }
            }
            if has_double_o && has_single_n {
                has_amide = true;
            }
            // Aldehyde: C(=O) with implicit H on the carbon
            if has_double_o && atom.implicit_hydrogens >= 1 {
                has_aldehyde = true;
            }
        }

        // S(=O) pattern
        if an == S {
            for &(_, bond_idx) in &mol.adjacency[i] {
                let nbr_an = mol.atoms[mol.bonds[bond_idx].atom1].atomic_number;
                let nbr_an2 = mol.atoms[mol.bonds[bond_idx].atom2].atomic_number;
                let other_an = if mol.bonds[bond_idx].atom1 == i {
                    nbr_an2
                } else {
                    nbr_an
                };
                if other_an == O && mol.bonds[bond_idx].order == BondOrder::Double {
                    has_sulfonyl = true;
                    break;
                }
            }
        }

        // Nitro group: N bonded to two oxygens via double/single bonds, or N(=O)=O
        if an == N {
            let mut o_bond_count = 0u32;
            for &(nbr, bond_idx) in &mol.adjacency[i] {
                if mol.atoms[nbr].atomic_number == O {
                    let order = mol.bonds[bond_idx].order;
                    if order == BondOrder::Double || order == BondOrder::Single {
                        o_bond_count += 1;
                    }
                }
            }
            if o_bond_count >= 2 {
                has_nitro = true;
            }
        }

        // Phosphate: P with double-bond O and at least one single-bond O
        if an == P {
            let mut p_double_o = false;
            let mut p_single_o_count = 0u32;
            for &(nbr, bond_idx) in &mol.adjacency[i] {
                if mol.atoms[nbr].atomic_number == O {
                    match mol.bonds[bond_idx].order {
                        BondOrder::Double => p_double_o = true,
                        BondOrder::Single => p_single_o_count += 1,
                        _ => {}
                    }
                }
            }
            if p_double_o && p_single_o_count >= 1 {
                has_phosphate = true;
            }
        }

        // Ether: oxygen bonded to two carbons (C-O-C)
        if an == O && atom.implicit_hydrogens == 0 && atom.formal_charge == 0 {
            let c_neighbors: Vec<_> = mol.adjacency[i]
                .iter()
                .filter(|&&(nbr, bond_idx)| {
                    mol.atoms[nbr].atomic_number == C
                        && mol.bonds[bond_idx].order == BondOrder::Single
                })
                .collect();
            if c_neighbors.len() >= 2 {
                has_c_o_c = true;
            }
        }

        // C-N-C: nitrogen bonded to at least two carbons via single bonds
        if an == N {
            let c_single_count = mol.adjacency[i]
                .iter()
                .filter(|&&(nbr, bond_idx)| {
                    mol.atoms[nbr].atomic_number == C
                        && mol.bonds[bond_idx].order == BondOrder::Single
                })
                .count();
            if c_single_count >= 2 {
                has_c_n_c = true;
            }
        }
    }

    // Key 150: C(=O)O (carboxyl/ester)
    if has_carboxyl {
        fp.set_bit(150);
    }
    // Key 151: C(=O)N (amide)
    if has_amide {
        fp.set_bit(151);
    }
    // Key 152: S(=O) (sulfonyl/sulfoxide)
    if has_sulfonyl {
        fp.set_bit(152);
    }
    // Key 159: O-H (hydroxyl)
    if has_hydroxyl {
        fp.set_bit(159);
    }
    // Key 161: N-H (amine)
    if has_amine_nh {
        fp.set_bit(161);
    }

    // Additional functional group keys
    // Key 22: Aldehyde
    if has_aldehyde {
        fp.set_bit(22);
    }
    // Key 23: Ester
    if has_ester {
        fp.set_bit(23);
    }
    // Key 24: Nitro group
    if has_nitro {
        fp.set_bit(24);
    }
    // Key 25: Thiol
    if has_thiol {
        fp.set_bit(25);
    }
    // Key 26: Phosphate
    if has_phosphate {
        fp.set_bit(26);
    }
    // Key 27: Ether (C-O-C)
    if has_c_o_c {
        fp.set_bit(27);
    }
    // Key 28: Tertiary amine (C-N-C)
    if has_c_n_c {
        fp.set_bit(28);
    }
    // Key 29: Charged nitrogen
    if has_charged_n {
        fp.set_bit(29);
    }
    // Key 30: Charged oxygen
    if has_charged_o {
        fp.set_bit(30);
    }
}

// ---------------------------------------------------------------------------
// Count-based keys
// ---------------------------------------------------------------------------

/// Set fingerprint bits for count-threshold MACCS keys.
fn set_count_keys(
    fp: &mut Fingerprint,
    counts: &[u32; 120],
    rings: &[Vec<usize>],
    mol: &Molecule,
) {
    // Key 155: >= 2 oxygen atoms
    if counts[O as usize] >= 2 {
        fp.set_bit(155);
    }
    // Key 156: >= 2 nitrogen atoms
    if counts[N as usize] >= 2 {
        fp.set_bit(156);
    }
    // Key 157: >= 2 sulfur atoms
    if counts[S as usize] >= 2 {
        fp.set_bit(157);
    }

    // Count ring atoms
    let mut ring_atom_flags = vec![false; mol.atom_count()];
    for ring in rings {
        for &idx in ring {
            ring_atom_flags[idx] = true;
        }
    }
    let ring_atom_count = ring_atom_flags.iter().filter(|&&f| f).count();

    // Key 160: >= 4 ring atoms
    if ring_atom_count >= 4 {
        fp.set_bit(160);
    }

    // Key 164: >= 2 rings
    if rings.len() >= 2 {
        fp.set_bit(164);
    }

    // Key 31: >= 3 rings
    if rings.len() >= 3 {
        fp.set_bit(31);
    }

    // Key 32: >= 6 ring atoms
    if ring_atom_count >= 6 {
        fp.set_bit(32);
    }

    // Key 33: >= 10 ring atoms
    if ring_atom_count >= 10 {
        fp.set_bit(33);
    }

    // Key 34: >= 2 aromatic rings
    let aromatic_ring_count = rings
        .iter()
        .filter(|ring| ring.iter().all(|&idx| mol.atoms[idx].is_aromatic))
        .count();
    if aromatic_ring_count >= 2 {
        fp.set_bit(34);
    }

    // Key 35: >= 3 oxygen atoms
    if counts[O as usize] >= 3 {
        fp.set_bit(35);
    }

    // Key 36: >= 4 oxygen atoms
    if counts[O as usize] >= 4 {
        fp.set_bit(36);
    }

    // Key 37: >= 3 nitrogen atoms
    if counts[N as usize] >= 3 {
        fp.set_bit(37);
    }

    // Key 38: >= 2 halogens
    let halogen_count =
        counts[F as usize] + counts[CL as usize] + counts[BR as usize] + counts[I as usize];
    if halogen_count >= 2 {
        fp.set_bit(38);
    }

    // Key 39: >= 3 halogens
    if halogen_count >= 3 {
        fp.set_bit(39);
    }

    // Key 40: Bond count thresholds
    if mol.bond_count() >= 4 {
        fp.set_bit(40); // >= 4 bonds
    }
    if mol.bond_count() >= 8 {
        fp.set_bit(41); // Key 41: >= 8 bonds
    }
    if mol.bond_count() >= 16 {
        fp.set_bit(42); // Key 42: >= 16 bonds
    }
    if mol.bond_count() >= 32 {
        fp.set_bit(43); // Key 43: >= 32 bonds
    }

    // Key 44: >= 2 double bonds
    let double_bond_count = mol
        .bonds
        .iter()
        .filter(|b| b.order == BondOrder::Double)
        .count();
    if double_bond_count >= 2 {
        fp.set_bit(44);
    }

    // Key 45: >= 4 double bonds
    if double_bond_count >= 4 {
        fp.set_bit(45);
    }

    // Key 46: Degree >= 3 atom present (branching point)
    if mol.atoms.iter().enumerate().any(|(i, _)| mol.degree(i) >= 3) {
        fp.set_bit(46);
    }

    // Key 47: Degree >= 4 atom present
    if mol.atoms.iter().enumerate().any(|(i, _)| mol.degree(i) >= 4) {
        fp.set_bit(47);
    }

    // Key 48: >= 2 heteroatoms in rings
    let hetero_in_ring = ring_atom_flags
        .iter()
        .enumerate()
        .filter(|&(i, &in_ring)| {
            in_ring && mol.atoms[i].atomic_number != C && mol.atoms[i].atomic_number != H
        })
        .count();
    if hetero_in_ring >= 2 {
        fp.set_bit(48);
    }

    // Key 49: Formal charge present
    if mol.atoms.iter().any(|a| a.formal_charge != 0) {
        fp.set_bit(49);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fingerprint::tanimoto_similarity;
    use crate::smiles::parse_smiles;

    #[test]
    fn benzene_aromatic_ring_key() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let fp = maccs_fingerprint(&mol);

        // Aromatic ring present (key 162)
        assert!(fp.get_bit(162), "benzene should have aromatic ring key");
        // Ring present (key 125)
        assert!(fp.get_bit(125), "benzene should have ring key");
        // 6-membered ring (key 148)
        assert!(fp.get_bit(148), "benzene should have 6-membered ring key");
        // Aromatic bond present (key 163)
        assert!(fp.get_bit(163), "benzene should have aromatic bond key");
        // Carbon present (key 126)
        assert!(fp.get_bit(126), "benzene should have carbon key");
        // No oxygen (key 118)
        assert!(!fp.get_bit(118), "benzene should not have oxygen key");
        // No nitrogen (key 119)
        assert!(!fp.get_bit(119), "benzene should not have nitrogen key");
    }

    #[test]
    fn ethanol_oxygen_and_hydroxyl() {
        let mol = parse_smiles("CCO").unwrap();
        let fp = maccs_fingerprint(&mol);

        // Oxygen present (key 118)
        assert!(fp.get_bit(118), "ethanol should have oxygen key");
        // O-H hydroxyl (key 159)
        assert!(fp.get_bit(159), "ethanol should have hydroxyl key");
        // Carbon present (key 126)
        assert!(fp.get_bit(126), "ethanol should have carbon key");
        // No ring (key 125)
        assert!(!fp.get_bit(125), "ethanol should not have ring key");
        // C-O bond (key 7)
        assert!(fp.get_bit(7), "ethanol should have C-O bond key");
        // C-C bond (key 9)
        assert!(fp.get_bit(9), "ethanol should have C-C bond key");
    }

    #[test]
    fn chlorobenzene_halogen_and_aromatic() {
        let mol = parse_smiles("Clc1ccccc1").unwrap();
        let fp = maccs_fingerprint(&mol);

        // Chlorine present (key 139)
        assert!(fp.get_bit(139), "chlorobenzene should have Cl key");
        // Halogen present (key 120)
        assert!(fp.get_bit(120), "chlorobenzene should have halogen key");
        // Aromatic ring (key 162)
        assert!(fp.get_bit(162), "chlorobenzene should have aromatic ring key");
        // C-Cl bond (key 15)
        assert!(fp.get_bit(15), "chlorobenzene should have C-Cl bond key");
    }

    #[test]
    fn aspirin_many_keys() {
        // Aspirin: CC(=O)Oc1ccccc1C(=O)O
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let fp = maccs_fingerprint(&mol);

        // C=O present (key 153)
        assert!(fp.get_bit(153), "aspirin should have C=O key");
        // Aromatic ring (key 162)
        assert!(fp.get_bit(162), "aspirin should have aromatic ring key");
        // Oxygen present (key 118)
        assert!(fp.get_bit(118), "aspirin should have oxygen key");
        // Carbon present (key 126)
        assert!(fp.get_bit(126), "aspirin should have carbon key");
        // >= 2 oxygens (key 155)
        assert!(fp.get_bit(155), "aspirin should have >= 2 O key");
        // C(=O)O carboxyl (key 150)
        assert!(fp.get_bit(150), "aspirin should have carboxyl key");
        // Ester (key 23) — the acetyl O-C(aromatic) linkage
        assert!(fp.get_bit(23), "aspirin should have ester key");
        // Ring present (key 125)
        assert!(fp.get_bit(125), "aspirin should have ring key");
        // Many bits set overall
        assert!(
            fp.count_ones() > 10,
            "aspirin should have >10 bits set, got {}",
            fp.count_ones()
        );
    }

    #[test]
    fn deterministic_fingerprint() {
        let mol = parse_smiles("c1ccc(O)cc1").unwrap();
        let fp1 = maccs_fingerprint(&mol);
        let fp2 = maccs_fingerprint(&mol);
        assert_eq!(fp1.nbits(), fp2.nbits());
        // Check all bits are identical
        for i in 0..166 {
            assert_eq!(
                fp1.get_bit(i),
                fp2.get_bit(i),
                "bit {} differs between runs",
                i
            );
        }
    }

    #[test]
    fn empty_molecule_no_keys() {
        let mol = Molecule::new(String::new(), vec![], vec![]);
        let fp = maccs_fingerprint(&mol);
        assert_eq!(fp.count_ones(), 0, "empty molecule should have no bits set");
    }

    #[test]
    fn methane_minimal_keys() {
        let mol = parse_smiles("C").unwrap();
        let fp = maccs_fingerprint(&mol);

        // Carbon present (key 126)
        assert!(fp.get_bit(126), "methane should have carbon key");
        // Hydrogen present (key 133)
        assert!(fp.get_bit(133), "methane should have hydrogen key");
        // No ring
        assert!(!fp.get_bit(125), "methane should not have ring key");
        // No oxygen
        assert!(!fp.get_bit(118), "methane should not have oxygen key");
        // Few bits overall
        assert!(
            fp.count_ones() < 10,
            "methane should have few bits, got {}",
            fp.count_ones()
        );
    }

    #[test]
    fn tanimoto_similar_maccs_reasonable() {
        // Phenol and cresol should be quite similar
        let phenol = parse_smiles("Oc1ccccc1").unwrap();
        let cresol = parse_smiles("Cc1ccc(O)cc1").unwrap();
        let fp1 = maccs_fingerprint(&phenol);
        let fp2 = maccs_fingerprint(&cresol);
        let sim = tanimoto_similarity(&fp1, &fp2);
        assert!(
            sim > 0.5,
            "phenol/cresol MACCS Tanimoto should be >0.5, got {sim}"
        );
        assert!(
            sim < 1.0,
            "phenol/cresol MACCS Tanimoto should be <1.0, got {sim}"
        );
    }
}
