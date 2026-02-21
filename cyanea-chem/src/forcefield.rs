//! Force field energy calculations and minimization.
//!
//! Implements UFF (Universal Force Field) and MMFF94 energy terms:
//! bond stretching, angle bending, torsion, van der Waals, and electrostatic.
//! Provides steepest descent and conjugate gradient minimization.

use cyanea_core::{CyaneaError, Result};

use crate::conformer::Conformer;
use crate::gasteiger::gasteiger_charges;
use crate::molecule::{BondOrder, Molecule};
use crate::ring;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Energy contributions from each force field term.
#[derive(Debug, Clone)]
pub struct EnergyComponents {
    pub bond_stretch: f64,
    pub angle_bend: f64,
    pub torsion: f64,
    pub van_der_waals: f64,
    pub electrostatic: f64,
    pub out_of_plane: f64,
    pub total: f64,
}

/// Result of energy minimization.
#[derive(Debug, Clone)]
pub struct MinimizeResult {
    pub conformer: Conformer,
    pub initial_energy: f64,
    pub final_energy: f64,
    pub n_steps: usize,
    pub converged: bool,
    pub energy_components: EnergyComponents,
}

/// Minimization method.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MinimizeMethod {
    SteepestDescent,
    ConjugateGradient,
}

/// Minimization configuration.
#[derive(Debug, Clone)]
pub struct MinimizeConfig {
    pub max_steps: usize,
    pub gradient_threshold: f64,
    pub method: MinimizeMethod,
}

impl Default for MinimizeConfig {
    fn default() -> Self {
        MinimizeConfig {
            max_steps: 500,
            gradient_threshold: 0.1,
            method: MinimizeMethod::SteepestDescent,
        }
    }
}

// ---------------------------------------------------------------------------
// UFF atom types
// ---------------------------------------------------------------------------

/// UFF atom type determined by element and hybridization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UffAtomType {
    H,
    C3, C2, CR, C1,
    N3, N2, NR, N1,
    O3, O2, OR,
    S3, S2,
    P3,
    F, Cl, Br, I,
    Si3, Se3,
}

/// UFF parameters for an atom type.
struct UffParams {
    r_cov: f64,      // covalent radius (Å)
    theta0: f64,     // ideal bond angle (radians)
    x_vdw: f64,      // vdW distance (Å)
    d_vdw: f64,      // vdW well depth (kcal/mol)
    chi: f64,        // electronegativity (Rappé)
}

fn uff_params(at: UffAtomType) -> UffParams {
    use std::f64::consts::PI;
    match at {
        UffAtomType::H    => UffParams { r_cov: 0.354, theta0: 180.0*PI/180.0, x_vdw: 2.886, d_vdw: 0.044, chi: 2.20 },
        UffAtomType::C3   => UffParams { r_cov: 0.757, theta0: 109.47*PI/180.0, x_vdw: 3.851, d_vdw: 0.105, chi: 2.55 },
        UffAtomType::C2   => UffParams { r_cov: 0.732, theta0: 120.0*PI/180.0, x_vdw: 3.851, d_vdw: 0.105, chi: 2.55 },
        UffAtomType::CR   => UffParams { r_cov: 0.729, theta0: 120.0*PI/180.0, x_vdw: 3.851, d_vdw: 0.105, chi: 2.55 },
        UffAtomType::C1   => UffParams { r_cov: 0.706, theta0: 180.0*PI/180.0, x_vdw: 3.851, d_vdw: 0.105, chi: 2.55 },
        UffAtomType::N3   => UffParams { r_cov: 0.700, theta0: 106.7*PI/180.0, x_vdw: 3.660, d_vdw: 0.069, chi: 3.04 },
        UffAtomType::N2   => UffParams { r_cov: 0.685, theta0: 120.0*PI/180.0, x_vdw: 3.660, d_vdw: 0.069, chi: 3.04 },
        UffAtomType::NR   => UffParams { r_cov: 0.683, theta0: 120.0*PI/180.0, x_vdw: 3.660, d_vdw: 0.069, chi: 3.04 },
        UffAtomType::N1   => UffParams { r_cov: 0.656, theta0: 180.0*PI/180.0, x_vdw: 3.660, d_vdw: 0.069, chi: 3.04 },
        UffAtomType::O3   => UffParams { r_cov: 0.658, theta0: 104.51*PI/180.0, x_vdw: 3.500, d_vdw: 0.060, chi: 3.44 },
        UffAtomType::O2   => UffParams { r_cov: 0.634, theta0: 120.0*PI/180.0, x_vdw: 3.500, d_vdw: 0.060, chi: 3.44 },
        UffAtomType::OR   => UffParams { r_cov: 0.639, theta0: 120.0*PI/180.0, x_vdw: 3.500, d_vdw: 0.060, chi: 3.44 },
        UffAtomType::S3   => UffParams { r_cov: 1.016, theta0: 92.2*PI/180.0, x_vdw: 4.035, d_vdw: 0.274, chi: 2.58 },
        UffAtomType::S2   => UffParams { r_cov: 0.992, theta0: 120.0*PI/180.0, x_vdw: 4.035, d_vdw: 0.274, chi: 2.58 },
        UffAtomType::P3   => UffParams { r_cov: 1.018, theta0: 93.8*PI/180.0, x_vdw: 4.147, d_vdw: 0.305, chi: 2.19 },
        UffAtomType::F    => UffParams { r_cov: 0.668, theta0: 180.0*PI/180.0, x_vdw: 3.364, d_vdw: 0.050, chi: 3.98 },
        UffAtomType::Cl   => UffParams { r_cov: 1.033, theta0: 180.0*PI/180.0, x_vdw: 3.947, d_vdw: 0.227, chi: 3.16 },
        UffAtomType::Br   => UffParams { r_cov: 1.176, theta0: 180.0*PI/180.0, x_vdw: 4.189, d_vdw: 0.251, chi: 2.96 },
        UffAtomType::I    => UffParams { r_cov: 1.333, theta0: 180.0*PI/180.0, x_vdw: 4.500, d_vdw: 0.339, chi: 2.66 },
        UffAtomType::Si3  => UffParams { r_cov: 1.116, theta0: 109.47*PI/180.0, x_vdw: 4.295, d_vdw: 0.402, chi: 1.90 },
        UffAtomType::Se3  => UffParams { r_cov: 1.170, theta0: 90.6*PI/180.0, x_vdw: 4.205, d_vdw: 0.291, chi: 2.55 },
    }
}

/// Assign UFF atom types to each atom in the molecule.
pub fn assign_uff_types(mol: &Molecule) -> Result<Vec<UffAtomType>> {
    let mut types = Vec::with_capacity(mol.atom_count());
    let bond_order_sums = compute_bond_order_sums(mol);

    for i in 0..mol.atom_count() {
        let atom = &mol.atoms[i];
        let bos = bond_order_sums[i];
        let at = match atom.atomic_number {
            1 => UffAtomType::H,
            6 => {
                if atom.is_aromatic { UffAtomType::CR }
                else if bos > 3.5 { UffAtomType::C1 }
                else if bos > 2.5 { UffAtomType::C2 }
                else { UffAtomType::C3 }
            }
            7 => {
                if atom.is_aromatic { UffAtomType::NR }
                else if bos > 3.5 { UffAtomType::N1 }
                else if bos > 2.5 { UffAtomType::N2 }
                else { UffAtomType::N3 }
            }
            8 => {
                if atom.is_aromatic { UffAtomType::OR }
                else if bos > 1.5 { UffAtomType::O2 }
                else { UffAtomType::O3 }
            }
            16 => {
                if bos > 2.5 { UffAtomType::S2 } else { UffAtomType::S3 }
            }
            15 => UffAtomType::P3,
            9 => UffAtomType::F,
            17 => UffAtomType::Cl,
            35 => UffAtomType::Br,
            53 => UffAtomType::I,
            14 => UffAtomType::Si3,
            34 => UffAtomType::Se3,
            _ => UffAtomType::C3, // fallback
        };
        types.push(at);
    }
    Ok(types)
}

// ---------------------------------------------------------------------------
// MMFF94 atom types
// ---------------------------------------------------------------------------

/// Simplified MMFF94 atom type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mmff94AtomType {
    CR,   // Alkyl carbon sp3
    C2,   // Vinylic carbon sp2
    C3,   // Carbonyl carbon sp2
    CAR,  // Aromatic carbon
    C1,   // Acetylenic carbon sp
    NR,   // Amine nitrogen sp3
    N2,   // Imine nitrogen sp2
    NAR,  // Aromatic nitrogen
    N1,   // Nitrile nitrogen sp
    OR,   // Ether/alcohol oxygen sp3
    O2,   // Carbonyl oxygen sp2
    OAR,  // Aromatic oxygen (furan)
    SR,   // Thiol/thioether sulfur
    S2,   // Thione sulfur
    SAR,  // Aromatic sulfur (thiophene)
    PR,   // Phosphorus sp3
    F,
    Cl,
    Br,
    I,
    H,
    HO,   // Hydrogen on oxygen
    HN,   // Hydrogen on nitrogen
    Si,
    Se,
    // Generic fallback
    DU,
}

/// MMFF94 vdW parameters.
#[allow(dead_code)]
struct Mmff94VdwParams {
    alpha: f64,  // polarizability
    n_eff: f64,  // effective electrons
    a_type: f64, // scale factor
    g: f64,      // well-depth factor
    da: u8,      // donor/acceptor: 0=neither, 1=donor, 2=acceptor
}

fn mmff94_vdw_params(at: Mmff94AtomType) -> Mmff94VdwParams {
    match at {
        Mmff94AtomType::CR  => Mmff94VdwParams { alpha: 1.050, n_eff: 2.490, a_type: 3.890, g: 1.282, da: 0 },
        Mmff94AtomType::C2  => Mmff94VdwParams { alpha: 1.350, n_eff: 2.490, a_type: 3.890, g: 1.282, da: 0 },
        Mmff94AtomType::C3  => Mmff94VdwParams { alpha: 1.100, n_eff: 2.490, a_type: 3.890, g: 1.282, da: 0 },
        Mmff94AtomType::CAR => Mmff94VdwParams { alpha: 1.350, n_eff: 2.490, a_type: 3.890, g: 1.282, da: 0 },
        Mmff94AtomType::C1  => Mmff94VdwParams { alpha: 1.300, n_eff: 2.490, a_type: 3.890, g: 1.282, da: 0 },
        Mmff94AtomType::NR  => Mmff94VdwParams { alpha: 1.000, n_eff: 2.820, a_type: 3.890, g: 1.282, da: 1 },
        Mmff94AtomType::N2  => Mmff94VdwParams { alpha: 1.100, n_eff: 2.820, a_type: 3.890, g: 1.282, da: 1 },
        Mmff94AtomType::NAR => Mmff94VdwParams { alpha: 1.100, n_eff: 2.820, a_type: 3.890, g: 1.282, da: 2 },
        Mmff94AtomType::N1  => Mmff94VdwParams { alpha: 1.000, n_eff: 2.820, a_type: 3.890, g: 1.282, da: 2 },
        Mmff94AtomType::OR  => Mmff94VdwParams { alpha: 0.700, n_eff: 3.150, a_type: 3.890, g: 1.282, da: 2 },
        Mmff94AtomType::O2  => Mmff94VdwParams { alpha: 0.700, n_eff: 3.150, a_type: 3.890, g: 1.282, da: 2 },
        Mmff94AtomType::OAR => Mmff94VdwParams { alpha: 0.700, n_eff: 3.150, a_type: 3.890, g: 1.282, da: 2 },
        Mmff94AtomType::SR  => Mmff94VdwParams { alpha: 3.000, n_eff: 3.480, a_type: 4.250, g: 1.345, da: 0 },
        Mmff94AtomType::S2  => Mmff94VdwParams { alpha: 3.000, n_eff: 3.480, a_type: 4.250, g: 1.345, da: 0 },
        Mmff94AtomType::SAR => Mmff94VdwParams { alpha: 3.000, n_eff: 3.480, a_type: 4.250, g: 1.345, da: 0 },
        Mmff94AtomType::PR  => Mmff94VdwParams { alpha: 1.600, n_eff: 3.480, a_type: 4.150, g: 1.345, da: 0 },
        Mmff94AtomType::F   => Mmff94VdwParams { alpha: 0.350, n_eff: 3.480, a_type: 3.480, g: 1.282, da: 2 },
        Mmff94AtomType::Cl  => Mmff94VdwParams { alpha: 2.315, n_eff: 3.480, a_type: 3.947, g: 1.345, da: 0 },
        Mmff94AtomType::Br  => Mmff94VdwParams { alpha: 3.400, n_eff: 3.480, a_type: 4.189, g: 1.359, da: 0 },
        Mmff94AtomType::I   => Mmff94VdwParams { alpha: 5.500, n_eff: 3.480, a_type: 4.500, g: 1.404, da: 0 },
        Mmff94AtomType::H   => Mmff94VdwParams { alpha: 0.250, n_eff: 0.800, a_type: 3.340, g: 1.112, da: 0 },
        Mmff94AtomType::HO  => Mmff94VdwParams { alpha: 0.250, n_eff: 0.800, a_type: 3.340, g: 1.112, da: 1 },
        Mmff94AtomType::HN  => Mmff94VdwParams { alpha: 0.250, n_eff: 0.800, a_type: 3.340, g: 1.112, da: 1 },
        Mmff94AtomType::Si  => Mmff94VdwParams { alpha: 1.600, n_eff: 2.490, a_type: 4.295, g: 1.345, da: 0 },
        Mmff94AtomType::Se  => Mmff94VdwParams { alpha: 3.400, n_eff: 3.480, a_type: 4.205, g: 1.359, da: 0 },
        Mmff94AtomType::DU  => Mmff94VdwParams { alpha: 1.000, n_eff: 2.490, a_type: 3.890, g: 1.282, da: 0 },
    }
}

/// Assign MMFF94 atom types to each atom in the molecule.
pub fn assign_mmff94_types(mol: &Molecule) -> Result<Vec<Mmff94AtomType>> {
    let mut types = Vec::with_capacity(mol.atom_count());
    let bond_order_sums = compute_bond_order_sums(mol);

    for i in 0..mol.atom_count() {
        let atom = &mol.atoms[i];
        let bos = bond_order_sums[i];
        let at = match atom.atomic_number {
            1 => {
                // Check what H is bonded to
                let mut bonded_to_o = false;
                let mut bonded_to_n = false;
                for &(nb, _) in &mol.adjacency[i] {
                    match mol.atoms[nb].atomic_number {
                        8 => bonded_to_o = true,
                        7 => bonded_to_n = true,
                        _ => {}
                    }
                }
                if bonded_to_o { Mmff94AtomType::HO }
                else if bonded_to_n { Mmff94AtomType::HN }
                else { Mmff94AtomType::H }
            }
            6 => {
                if atom.is_aromatic { Mmff94AtomType::CAR }
                else {
                    let has_triple = mol.adjacency[i].iter().any(|&(_, bi)| {
                        mol.bonds[bi].order == BondOrder::Triple
                    });
                    let has_double_o = mol.adjacency[i].iter().any(|&(nb, bi)| {
                        mol.atoms[nb].atomic_number == 8 && mol.bonds[bi].order == BondOrder::Double
                    });
                    if has_triple { Mmff94AtomType::C1 }
                    else if has_double_o { Mmff94AtomType::C3 }
                    else if bos > 2.5 { Mmff94AtomType::C2 }
                    else { Mmff94AtomType::CR }
                }
            }
            7 => {
                if atom.is_aromatic { Mmff94AtomType::NAR }
                else if bos > 3.5 { Mmff94AtomType::N1 }
                else if bos > 2.5 { Mmff94AtomType::N2 }
                else { Mmff94AtomType::NR }
            }
            8 => {
                if atom.is_aromatic { Mmff94AtomType::OAR }
                else if bos > 1.5 { Mmff94AtomType::O2 }
                else { Mmff94AtomType::OR }
            }
            16 => {
                if atom.is_aromatic { Mmff94AtomType::SAR }
                else if bos > 2.5 { Mmff94AtomType::S2 }
                else { Mmff94AtomType::SR }
            }
            15 => Mmff94AtomType::PR,
            9 => Mmff94AtomType::F,
            17 => Mmff94AtomType::Cl,
            35 => Mmff94AtomType::Br,
            53 => Mmff94AtomType::I,
            14 => Mmff94AtomType::Si,
            34 => Mmff94AtomType::Se,
            _ => Mmff94AtomType::DU,
        };
        types.push(at);
    }
    Ok(types)
}

// ---------------------------------------------------------------------------
// UFF energy
// ---------------------------------------------------------------------------

/// Compute UFF energy for a molecule with 3D coordinates.
pub fn uff_energy(mol: &Molecule, conf: &Conformer) -> Result<EnergyComponents> {
    if conf.len() != mol.atom_count() {
        return Err(CyaneaError::InvalidInput(
            "conformer atom count doesn't match molecule".into(),
        ));
    }
    let types = assign_uff_types(mol)?;
    let charges = gasteiger_charges(mol).unwrap_or_else(|_| vec![0.0; mol.atom_count()]);
    let rings = ring::find_sssr(mol);
    let ring_bonds = ring_bond_set(mol, &rings);

    let e_bond = uff_bond_stretch(mol, conf, &types);
    let e_angle = uff_angle_bend(mol, conf, &types);
    let e_torsion = uff_torsion(mol, conf, &types, &ring_bonds);
    let e_vdw = uff_vdw(mol, conf, &types);
    let e_elec = electrostatic_energy(mol, conf, &charges);
    let total = e_bond + e_angle + e_torsion + e_vdw + e_elec;

    Ok(EnergyComponents {
        bond_stretch: e_bond,
        angle_bend: e_angle,
        torsion: e_torsion,
        van_der_waals: e_vdw,
        electrostatic: e_elec,
        out_of_plane: 0.0,
        total,
    })
}

/// Compute MMFF94 energy for a molecule with 3D coordinates.
pub fn mmff94_energy(mol: &Molecule, conf: &Conformer) -> Result<EnergyComponents> {
    if conf.len() != mol.atom_count() {
        return Err(CyaneaError::InvalidInput(
            "conformer atom count doesn't match molecule".into(),
        ));
    }
    let types = assign_mmff94_types(mol)?;
    let uff_types = assign_uff_types(mol)?;
    let charges = gasteiger_charges(mol).unwrap_or_else(|_| vec![0.0; mol.atom_count()]);
    let rings = ring::find_sssr(mol);
    let ring_bonds = ring_bond_set(mol, &rings);

    // MMFF94 uses similar functional forms to UFF with different parameters
    let e_bond = mmff94_bond_stretch(mol, conf, &types);
    let e_angle = uff_angle_bend(mol, conf, &uff_types); // reuse UFF angle bend logic
    let e_torsion = uff_torsion(mol, conf, &uff_types, &ring_bonds);
    let e_vdw = mmff94_vdw(mol, conf, &types);
    let e_elec = electrostatic_energy(mol, conf, &charges);
    let e_oop = mmff94_oop(mol, conf);
    let total = e_bond + e_angle + e_torsion + e_vdw + e_elec + e_oop;

    Ok(EnergyComponents {
        bond_stretch: e_bond,
        angle_bend: e_angle,
        torsion: e_torsion,
        van_der_waals: e_vdw,
        electrostatic: e_elec,
        out_of_plane: e_oop,
        total,
    })
}

// ---------------------------------------------------------------------------
// UFF bond stretching: E = 0.5 * k * (r - r0)^2
// ---------------------------------------------------------------------------

fn uff_bond_stretch(mol: &Molecule, conf: &Conformer, types: &[UffAtomType]) -> f64 {
    let mut energy = 0.0;
    for bond in &mol.bonds {
        let p1 = uff_params(types[bond.atom1]);
        let p2 = uff_params(types[bond.atom2]);

        // Bond order correction
        let bo_corr = match bond.order {
            BondOrder::Single => 0.0,
            BondOrder::Aromatic => -0.0332,
            BondOrder::Double => -0.0668,
            BondOrder::Triple => -0.0997,
        };

        let r0 = p1.r_cov + p2.r_cov + bo_corr;
        let r = conf.distance(bond.atom1, bond.atom2);

        // Force constant: 664.12 * Z1*Z2 / r0^3 (simplified)
        let chi_prod = p1.chi * p2.chi;
        let k = 664.12 * chi_prod / (r0 * r0 * r0);
        let k = k.min(2000.0); // cap to avoid blow-up

        energy += 0.5 * k * (r - r0) * (r - r0);
    }
    energy
}

// ---------------------------------------------------------------------------
// UFF angle bending: E = 0.5 * k * (theta - theta0)^2
// ---------------------------------------------------------------------------

fn uff_angle_bend(mol: &Molecule, conf: &Conformer, types: &[UffAtomType]) -> f64 {
    let mut energy = 0.0;

    for j in 0..mol.atom_count() {
        let neighbors = &mol.adjacency[j];
        if neighbors.len() < 2 {
            continue;
        }
        let theta0 = uff_params(types[j]).theta0;
        let k_theta = 50.0; // kcal/mol/rad^2 (simplified)

        for a in 0..neighbors.len() {
            for b in (a + 1)..neighbors.len() {
                let i = neighbors[a].0;
                let k = neighbors[b].0;
                let theta = conf.angle(i, j, k);
                let diff = theta - theta0;
                energy += 0.5 * k_theta * diff * diff;
            }
        }
    }
    energy
}

// ---------------------------------------------------------------------------
// UFF torsion: E = 0.5 * V * (1 - cos(n * phi))
// ---------------------------------------------------------------------------

fn uff_torsion(
    mol: &Molecule,
    conf: &Conformer,
    types: &[UffAtomType],
    ring_bonds: &[usize],
) -> f64 {
    let mut energy = 0.0;

    for (bi, bond) in mol.bonds.iter().enumerate() {
        // Skip terminal and ring bonds for torsion
        if mol.degree(bond.atom1) < 2 || mol.degree(bond.atom2) < 2 {
            continue;
        }

        let j = bond.atom1;
        let k = bond.atom2;

        // Determine torsion barrier and periodicity
        let (v, n) = torsion_params(types[j], types[k], bond.order);
        if v.abs() < 1e-10 {
            continue;
        }

        // Reduced torsion barrier for ring bonds
        let v = if ring_bonds.contains(&bi) { v * 0.5 } else { v };

        for &(i, _) in &mol.adjacency[j] {
            if i == k {
                continue;
            }
            for &(l, _) in &mol.adjacency[k] {
                if l == j || l == i {
                    continue;
                }
                let phi = conf.dihedral(i, j, k, l);
                energy += 0.5 * v * (1.0 - (n as f64 * phi).cos());
            }
        }
    }
    energy
}

fn torsion_params(t1: UffAtomType, t2: UffAtomType, order: BondOrder) -> (f64, u32) {
    match order {
        BondOrder::Double => (45.0, 2),  // High barrier, period 2
        BondOrder::Triple => (0.0, 1),   // No torsion for triple bonds
        BondOrder::Aromatic => (3.0, 2), // Moderate for aromatic
        BondOrder::Single => {
            // sp3-sp3 → V=1, n=3; sp3-sp2 → V=1, n=6; sp2-sp2 → V=5, n=2
            let is_sp2_1 = matches!(t1, UffAtomType::C2 | UffAtomType::CR | UffAtomType::N2 | UffAtomType::NR | UffAtomType::O2 | UffAtomType::OR);
            let is_sp2_2 = matches!(t2, UffAtomType::C2 | UffAtomType::CR | UffAtomType::N2 | UffAtomType::NR | UffAtomType::O2 | UffAtomType::OR);
            match (is_sp2_1, is_sp2_2) {
                (true, true) => (5.0, 2),
                (true, false) | (false, true) => (1.0, 6),
                (false, false) => (1.0, 3),
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Van der Waals: Lennard-Jones 12-6
// ---------------------------------------------------------------------------

fn uff_vdw(mol: &Molecule, conf: &Conformer, types: &[UffAtomType]) -> f64 {
    let mut energy = 0.0;
    let n = mol.atom_count();

    // Build 1-2 and 1-3 exclusion set
    let excluded = build_exclusion_set(mol);

    for i in 0..n {
        for j in (i + 1)..n {
            if excluded.contains(&(i, j)) {
                continue;
            }
            let p1 = uff_params(types[i]);
            let p2 = uff_params(types[j]);

            let x_ij = (p1.x_vdw * p2.x_vdw).sqrt(); // geometric mean
            let d_ij = (p1.d_vdw * p2.d_vdw).sqrt();

            let r = conf.distance(i, j);
            if r < 0.5 || r > 10.0 {
                continue; // skip unreasonable distances
            }

            let ratio = x_ij / r;
            let r6 = ratio * ratio * ratio * ratio * ratio * ratio;
            let r12 = r6 * r6;

            energy += d_ij * (r12 - 2.0 * r6);
        }
    }
    energy
}

fn mmff94_vdw(mol: &Molecule, conf: &Conformer, types: &[Mmff94AtomType]) -> f64 {
    let mut energy = 0.0;
    let n = mol.atom_count();
    let excluded = build_exclusion_set(mol);

    for i in 0..n {
        for j in (i + 1)..n {
            if excluded.contains(&(i, j)) {
                continue;
            }
            let p1 = mmff94_vdw_params(types[i]);
            let p2 = mmff94_vdw_params(types[j]);

            // MMFF94 buffered 14-7 potential (simplified as 12-6 with MMFF params)
            let r_star = (p1.a_type + p2.a_type) / 2.0;
            let eps = (p1.g * p2.g).sqrt() * 0.05;

            let r = conf.distance(i, j);
            if r < 0.5 || r > 10.0 {
                continue;
            }

            let ratio = r_star / r;
            let r6 = ratio * ratio * ratio * ratio * ratio * ratio;
            let r12 = r6 * r6;

            energy += eps * (r12 - 2.0 * r6);
        }
    }
    energy
}

// ---------------------------------------------------------------------------
// MMFF94 bond stretching with cubic term
// ---------------------------------------------------------------------------

fn mmff94_bond_stretch(mol: &Molecule, conf: &Conformer, types: &[Mmff94AtomType]) -> f64 {
    let mut energy = 0.0;
    for bond in &mol.bonds {
        let r0 = mmff94_eq_bond_length(types[bond.atom1], types[bond.atom2], bond.order);
        let r = conf.distance(bond.atom1, bond.atom2);
        let k = mmff94_bond_k(types[bond.atom1], types[bond.atom2], bond.order);
        let dr = r - r0;
        // MMFF94 cubic stretch: E = 0.5*k*dr^2 * (1 + cs*dr)
        let cs = -2.0; // cubic stretch constant
        energy += 0.5 * k * dr * dr * (1.0 + cs * dr);
    }
    energy
}

fn mmff94_eq_bond_length(t1: Mmff94AtomType, t2: Mmff94AtomType, order: BondOrder) -> f64 {
    // Simplified: use covalent radii sum with bond order correction
    let r1 = mmff94_cov_radius(t1);
    let r2 = mmff94_cov_radius(t2);
    let bo_corr = match order {
        BondOrder::Single => 0.0,
        BondOrder::Aromatic => -0.04,
        BondOrder::Double => -0.08,
        BondOrder::Triple => -0.12,
    };
    r1 + r2 + bo_corr
}

fn mmff94_cov_radius(at: Mmff94AtomType) -> f64 {
    match at {
        Mmff94AtomType::H | Mmff94AtomType::HO | Mmff94AtomType::HN => 0.33,
        Mmff94AtomType::CR | Mmff94AtomType::C2 | Mmff94AtomType::C3 |
        Mmff94AtomType::CAR | Mmff94AtomType::C1 => 0.77,
        Mmff94AtomType::NR | Mmff94AtomType::N2 | Mmff94AtomType::NAR |
        Mmff94AtomType::N1 => 0.70,
        Mmff94AtomType::OR | Mmff94AtomType::O2 | Mmff94AtomType::OAR => 0.66,
        Mmff94AtomType::F => 0.64,
        Mmff94AtomType::SR | Mmff94AtomType::S2 | Mmff94AtomType::SAR => 1.04,
        Mmff94AtomType::Cl => 0.99,
        Mmff94AtomType::Br => 1.14,
        Mmff94AtomType::I => 1.33,
        Mmff94AtomType::PR => 1.10,
        Mmff94AtomType::Si => 1.17,
        Mmff94AtomType::Se => 1.17,
        Mmff94AtomType::DU => 0.77,
    }
}

fn mmff94_bond_k(t1: Mmff94AtomType, t2: Mmff94AtomType, order: BondOrder) -> f64 {
    let base = match order {
        BondOrder::Single => 300.0,
        BondOrder::Aromatic => 400.0,
        BondOrder::Double => 600.0,
        BondOrder::Triple => 800.0,
    };
    // Adjust by element — heavier atoms have softer bonds
    let scale1 = mmff94_bond_scale(t1);
    let scale2 = mmff94_bond_scale(t2);
    base * scale1 * scale2
}

fn mmff94_bond_scale(at: Mmff94AtomType) -> f64 {
    match at {
        Mmff94AtomType::H | Mmff94AtomType::HO | Mmff94AtomType::HN => 0.8,
        Mmff94AtomType::CR | Mmff94AtomType::C2 | Mmff94AtomType::C3 |
        Mmff94AtomType::CAR | Mmff94AtomType::C1 => 1.0,
        Mmff94AtomType::NR | Mmff94AtomType::N2 | Mmff94AtomType::NAR |
        Mmff94AtomType::N1 => 1.1,
        Mmff94AtomType::OR | Mmff94AtomType::O2 | Mmff94AtomType::OAR => 1.2,
        Mmff94AtomType::F => 1.3,
        _ => 0.9,
    }
}

// ---------------------------------------------------------------------------
// MMFF94 out-of-plane bending
// ---------------------------------------------------------------------------

fn mmff94_oop(mol: &Molecule, conf: &Conformer) -> f64 {
    let mut energy = 0.0;
    let k_oop = 15.0; // kcal/mol/rad^2

    for j in 0..mol.atom_count() {
        let atom = &mol.atoms[j];
        // Only sp2 atoms (aromatic or double-bonded) with exactly 3 neighbors
        let is_sp2 = atom.is_aromatic || mol.adjacency[j].iter().any(|&(_, bi)| {
            mol.bonds[bi].order == BondOrder::Double
        });
        if !is_sp2 || mol.adjacency[j].len() != 3 {
            continue;
        }

        let neighbors: Vec<usize> = mol.adjacency[j].iter().map(|&(n, _)| n).collect();
        // Wilson angle: angle of atom j out of the plane of its 3 neighbors
        let chi = wilson_angle(conf, j, neighbors[0], neighbors[1], neighbors[2]);
        energy += 0.5 * k_oop * chi * chi;
    }
    energy
}

fn wilson_angle(conf: &Conformer, center: usize, a: usize, b: usize, c: usize) -> f64 {
    let va = sub3(conf.coords[a], conf.coords[center]);
    let vb = sub3(conf.coords[b], conf.coords[center]);
    let vc = sub3(conf.coords[c], conf.coords[center]);

    // Normal to plane of a, b, c
    let n = cross3(vb, vc);
    let n_len = norm3(n);
    let va_len = norm3(va);

    if n_len < 1e-12 || va_len < 1e-12 {
        return 0.0;
    }

    let sin_chi = dot3v(va, n) / (va_len * n_len);
    sin_chi.clamp(-1.0, 1.0).asin()
}

// ---------------------------------------------------------------------------
// Electrostatic: Coulomb with distance-dependent dielectric
// ---------------------------------------------------------------------------

fn electrostatic_energy(mol: &Molecule, conf: &Conformer, charges: &[f64]) -> f64 {
    let mut energy = 0.0;
    let n = mol.atom_count();
    let excluded = build_exclusion_set(mol);
    let coulomb_const = 332.0637; // kcal*Å/(mol*e^2)

    for i in 0..n {
        for j in (i + 1)..n {
            if excluded.contains(&(i, j)) {
                continue;
            }
            let r = conf.distance(i, j);
            if r < 0.5 {
                continue;
            }
            // Distance-dependent dielectric: ε = r
            energy += coulomb_const * charges[i] * charges[j] / (r * r);
        }
    }
    energy
}

// ---------------------------------------------------------------------------
// Minimization
// ---------------------------------------------------------------------------

/// Minimize UFF energy.
pub fn uff_minimize(
    mol: &Molecule,
    conformer: &Conformer,
    config: &MinimizeConfig,
) -> Result<MinimizeResult> {
    minimize_impl(mol, conformer, config, ForceFieldKind::Uff)
}

/// Minimize MMFF94 energy.
pub fn mmff94_minimize(
    mol: &Molecule,
    conformer: &Conformer,
    config: &MinimizeConfig,
) -> Result<MinimizeResult> {
    minimize_impl(mol, conformer, config, ForceFieldKind::Mmff94)
}

#[derive(Clone, Copy)]
enum ForceFieldKind {
    Uff,
    Mmff94,
}

fn energy_fn(mol: &Molecule, conf: &Conformer, kind: ForceFieldKind) -> Result<EnergyComponents> {
    match kind {
        ForceFieldKind::Uff => uff_energy(mol, conf),
        ForceFieldKind::Mmff94 => mmff94_energy(mol, conf),
    }
}

fn minimize_impl(
    mol: &Molecule,
    conformer: &Conformer,
    config: &MinimizeConfig,
    kind: ForceFieldKind,
) -> Result<MinimizeResult> {
    let initial_components = energy_fn(mol, conformer, kind)?;
    let initial_energy = initial_components.total;

    let mut current = conformer.clone();
    let mut current_energy = initial_energy;
    let n = mol.atom_count();
    let dx = 0.001; // Å for numerical gradient

    let mut prev_gradient = vec![[0.0_f64; 3]; n];
    let mut direction = vec![[0.0_f64; 3]; n];
    let mut converged = false;
    let mut step_count = 0;

    for step in 0..config.max_steps {
        step_count = step + 1;

        // Compute gradient via central differences
        let mut gradient = vec![[0.0_f64; 3]; n];
        for i in 0..n {
            for dim in 0..3 {
                let mut plus = current.clone();
                let mut minus = current.clone();
                plus.coords[i][dim] += dx;
                minus.coords[i][dim] -= dx;

                let e_plus = energy_fn(mol, &plus, kind)?.total;
                let e_minus = energy_fn(mol, &minus, kind)?.total;
                gradient[i][dim] = (e_plus - e_minus) / (2.0 * dx);
            }
        }

        // Check convergence: gradient norm
        let grad_norm: f64 = gradient.iter().map(|g| g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sum::<f64>().sqrt();
        if grad_norm < config.gradient_threshold {
            converged = true;
            break;
        }

        // Compute search direction
        match config.method {
            MinimizeMethod::SteepestDescent => {
                for i in 0..n {
                    direction[i] = [-gradient[i][0], -gradient[i][1], -gradient[i][2]];
                }
            }
            MinimizeMethod::ConjugateGradient => {
                if step == 0 {
                    for i in 0..n {
                        direction[i] = [-gradient[i][0], -gradient[i][1], -gradient[i][2]];
                    }
                } else {
                    // Fletcher-Reeves
                    let num: f64 = gradient.iter().map(|g| g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sum();
                    let den: f64 = prev_gradient.iter().map(|g| g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sum();
                    let beta = if den > 1e-30 { num / den } else { 0.0 };
                    let beta = beta.min(2.0); // reset if β is too large

                    for i in 0..n {
                        direction[i] = [
                            -gradient[i][0] + beta * direction[i][0],
                            -gradient[i][1] + beta * direction[i][1],
                            -gradient[i][2] + beta * direction[i][2],
                        ];
                    }
                }
            }
        }

        // Line search: try step sizes 0.01, 0.005, 0.001
        let mut best_step = 0.0;
        let mut best_energy = current_energy;

        for &alpha in &[0.02, 0.01, 0.005, 0.001] {
            let mut trial = current.clone();
            for i in 0..n {
                trial.coords[i][0] += alpha * direction[i][0];
                trial.coords[i][1] += alpha * direction[i][1];
                trial.coords[i][2] += alpha * direction[i][2];
            }
            let e = energy_fn(mol, &trial, kind)?.total;
            if e < best_energy {
                best_energy = e;
                best_step = alpha;
            }
        }

        if best_step > 0.0 {
            for i in 0..n {
                current.coords[i][0] += best_step * direction[i][0];
                current.coords[i][1] += best_step * direction[i][1];
                current.coords[i][2] += best_step * direction[i][2];
            }
            current_energy = best_energy;
        } else {
            // No improvement — reduce step or stop
            break;
        }

        prev_gradient = gradient;
    }

    let final_components = energy_fn(mol, &current, kind)?;

    Ok(MinimizeResult {
        conformer: current,
        initial_energy,
        final_energy: final_components.total,
        n_steps: step_count,
        converged,
        energy_components: final_components,
    })
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn compute_bond_order_sums(mol: &Molecule) -> Vec<f64> {
    let mut sums = vec![0.0_f64; mol.atom_count()];
    for bond in &mol.bonds {
        let order = bond.order.as_f64();
        sums[bond.atom1] += order;
        sums[bond.atom2] += order;
    }
    sums
}

fn ring_bond_set(mol: &Molecule, rings: &[Vec<usize>]) -> Vec<usize> {
    let mut set = Vec::new();
    for ring in rings {
        for i in 0..ring.len() {
            let a1 = ring[i];
            let a2 = ring[(i + 1) % ring.len()];
            for &(neighbor, bond_idx) in &mol.adjacency[a1] {
                if neighbor == a2 && !set.contains(&bond_idx) {
                    set.push(bond_idx);
                }
            }
        }
    }
    set
}

/// Build exclusion set: 1-2 (bonded) and 1-3 (angle) pairs.
fn build_exclusion_set(mol: &Molecule) -> std::collections::HashSet<(usize, usize)> {
    let mut excluded = std::collections::HashSet::new();

    // 1-2 exclusions (bonded atoms)
    for bond in &mol.bonds {
        let (a, b) = if bond.atom1 < bond.atom2 {
            (bond.atom1, bond.atom2)
        } else {
            (bond.atom2, bond.atom1)
        };
        excluded.insert((a, b));
    }

    // 1-3 exclusions (atoms separated by 2 bonds)
    for j in 0..mol.atom_count() {
        let neighbors = &mol.adjacency[j];
        for a in 0..neighbors.len() {
            for b in (a + 1)..neighbors.len() {
                let i = neighbors[a].0;
                let k = neighbors[b].0;
                let (lo, hi) = if i < k { (i, k) } else { (k, i) };
                excluded.insert((lo, hi));
            }
        }
    }

    excluded
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot3v(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformer::Conformer;
    use crate::smiles::parse_smiles;

    fn ethane_with_coords() -> (Molecule, Conformer) {
        let mol = parse_smiles("CC").unwrap();
        let conf = Conformer::new(vec![
            [0.0, 0.0, 0.0],
            [1.54, 0.0, 0.0],
        ]);
        (mol, conf)
    }

    fn water_with_coords() -> (Molecule, Conformer) {
        let mol = parse_smiles("[OH2]").unwrap();
        let conf = Conformer::new(vec![[0.0, 0.0, 0.0]]);
        (mol, conf)
    }

    #[test]
    fn uff_type_assignment() {
        let mol = parse_smiles("CCO").unwrap();
        let types = assign_uff_types(&mol).unwrap();
        assert_eq!(types[0], UffAtomType::C3); // methyl C
        assert_eq!(types[1], UffAtomType::C3); // methylene C
        assert_eq!(types[2], UffAtomType::O3); // hydroxyl O
    }

    #[test]
    fn uff_type_aromatic() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let types = assign_uff_types(&mol).unwrap();
        for t in &types {
            assert_eq!(*t, UffAtomType::CR);
        }
    }

    #[test]
    fn mmff94_type_assignment() {
        let mol = parse_smiles("CC(=O)O").unwrap();
        let types = assign_mmff94_types(&mol).unwrap();
        assert_eq!(types[0], Mmff94AtomType::CR);  // methyl C
        assert_eq!(types[1], Mmff94AtomType::C3);  // carbonyl C
        assert_eq!(types[2], Mmff94AtomType::O2);  // C=O
        assert_eq!(types[3], Mmff94AtomType::OR);  // C-O-H
    }

    #[test]
    fn uff_energy_ethane() {
        let (mol, conf) = ethane_with_coords();
        let e = uff_energy(&mol, &conf).unwrap();
        // Should produce a finite energy
        assert!(e.total.is_finite(), "total = {}", e.total);
        assert!(e.bond_stretch >= 0.0);
    }

    #[test]
    fn mmff94_energy_ethane() {
        let (mol, conf) = ethane_with_coords();
        let e = mmff94_energy(&mol, &conf).unwrap();
        assert!(e.total.is_finite(), "total = {}", e.total);
    }

    #[test]
    fn energy_atom_count_mismatch() {
        let mol = parse_smiles("CCO").unwrap();
        let conf = Conformer::new(vec![[0.0; 3], [1.0, 0.0, 0.0]]); // 2 atoms, mol has 3
        assert!(uff_energy(&mol, &conf).is_err());
    }

    #[test]
    fn uff_energy_components_positive() {
        let mol = parse_smiles("CCCC").unwrap();
        let conf = Conformer::new(vec![
            [0.0, 0.0, 0.0],
            [1.54, 0.0, 0.0],
            [2.31, 1.26, 0.0],
            [3.85, 1.26, 0.0],
        ]);
        let e = uff_energy(&mol, &conf).unwrap();
        assert!(e.total.is_finite());
        assert!(e.bond_stretch.is_finite());
        assert!(e.angle_bend.is_finite());
    }

    #[test]
    fn uff_minimize_basic() {
        let mol = parse_smiles("CC").unwrap();
        // Start with slightly wrong C-C distance
        let conf = Conformer::new(vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0], // too far — should pull closer
        ]);
        let config = MinimizeConfig {
            max_steps: 50,
            gradient_threshold: 1.0,
            method: MinimizeMethod::SteepestDescent,
        };
        let result = uff_minimize(&mol, &conf, &config).unwrap();
        assert!(result.final_energy <= result.initial_energy,
            "energy should decrease: {} -> {}", result.initial_energy, result.final_energy);
    }

    #[test]
    fn mmff94_minimize_basic() {
        let mol = parse_smiles("CC").unwrap();
        let conf = Conformer::new(vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ]);
        let config = MinimizeConfig {
            max_steps: 50,
            gradient_threshold: 1.0,
            method: MinimizeMethod::SteepestDescent,
        };
        let result = mmff94_minimize(&mol, &conf, &config).unwrap();
        assert!(result.final_energy <= result.initial_energy);
    }

    #[test]
    fn conjugate_gradient_minimize() {
        let mol = parse_smiles("CC").unwrap();
        let conf = Conformer::new(vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ]);
        let config = MinimizeConfig {
            max_steps: 50,
            gradient_threshold: 1.0,
            method: MinimizeMethod::ConjugateGradient,
        };
        let result = uff_minimize(&mol, &conf, &config).unwrap();
        assert!(result.final_energy <= result.initial_energy);
    }

    #[test]
    fn exclusion_set_ethane() {
        let mol = parse_smiles("CC").unwrap();
        let excl = build_exclusion_set(&mol);
        // Bond: (0,1)
        assert!(excl.contains(&(0, 1)));
    }

    #[test]
    fn minimize_result_steps() {
        let mol = parse_smiles("CC").unwrap();
        let conf = Conformer::new(vec![
            [0.0, 0.0, 0.0],
            [1.54, 0.0, 0.0], // near equilibrium
        ]);
        let config = MinimizeConfig {
            max_steps: 100,
            gradient_threshold: 0.01,
            method: MinimizeMethod::SteepestDescent,
        };
        let result = uff_minimize(&mol, &conf, &config).unwrap();
        assert!(result.n_steps <= 100);
    }

    #[test]
    fn uff_energy_single_atom() {
        let (mol, conf) = water_with_coords();
        let e = uff_energy(&mol, &conf).unwrap();
        assert!(e.total.is_finite());
        // Single atom: no bonds, no angles, no torsions
        assert!((e.bond_stretch).abs() < 1e-10);
    }

    #[test]
    fn benzene_uff_types() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let types = assign_uff_types(&mol).unwrap();
        assert!(types.iter().all(|t| *t == UffAtomType::CR));
    }

    #[test]
    fn mmff94_aromatic_nitrogen() {
        let mol = parse_smiles("c1ccncc1").unwrap(); // pyridine
        let types = assign_mmff94_types(&mol).unwrap();
        let n_idx = mol.atoms.iter().position(|a| a.atomic_number == 7).unwrap();
        assert_eq!(types[n_idx], Mmff94AtomType::NAR);
    }
}
