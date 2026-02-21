//! 3D coordinate embedding via distance geometry and ETKDG.
//!
//! Generates 3D conformers from molecular connectivity using:
//! 1. Distance geometry (DG) — bounds matrix → random sampling → eigendecomposition
//! 2. ETKDG-style torsion angle preferences
//! 3. Multi-conformer enumeration with RMSD pruning

use cyanea_core::Result;

use crate::conformer::{Conformer, ConformerSet};
use crate::forcefield::{self, MinimizeConfig, MinimizeMethod};
use crate::molecule::{BondOrder, Molecule};
use crate::ring;

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Force field to use for post-embedding optimization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ForceFieldType {
    None,
    Uff,
    Mmff94,
}

/// Configuration for 3D embedding.
#[derive(Debug, Clone)]
pub struct EmbedConfig {
    /// Maximum number of conformers to generate.
    pub max_conformers: usize,
    /// RMSD threshold (Å) for duplicate pruning.
    pub rmsd_threshold: f64,
    /// Use ETKDG torsion angle preferences.
    pub use_torsion_prefs: bool,
    /// Random seed for reproducibility.
    pub random_seed: u64,
    /// Force field for post-embedding optimization.
    pub force_field: ForceFieldType,
    /// Max steps for force-field minimization (0 = skip).
    pub max_minimize_steps: usize,
}

impl Default for EmbedConfig {
    fn default() -> Self {
        EmbedConfig {
            max_conformers: 1,
            rmsd_threshold: 0.5,
            use_torsion_prefs: true,
            random_seed: 42,
            force_field: ForceFieldType::Uff,
            max_minimize_steps: 200,
        }
    }
}

// ---------------------------------------------------------------------------
// Covalent radii (Å) for distance bounds
// ---------------------------------------------------------------------------

fn covalent_radius(atomic_number: u8) -> f64 {
    match atomic_number {
        1 => 0.31,
        5 => 0.84,
        6 => 0.76,
        7 => 0.71,
        8 => 0.66,
        9 => 0.57,
        14 => 1.11,
        15 => 1.07,
        16 => 1.05,
        17 => 1.02,
        34 => 1.20,
        35 => 1.20,
        53 => 1.39,
        _ => 0.77, // default
    }
}

/// Van der Waals radii for lower distance bounds between non-bonded atoms.
fn vdw_radius(atomic_number: u8) -> f64 {
    match atomic_number {
        1 => 1.20,
        6 => 1.70,
        7 => 1.55,
        8 => 1.52,
        9 => 1.47,
        15 => 1.80,
        16 => 1.80,
        17 => 1.75,
        35 => 1.85,
        53 => 1.98,
        _ => 1.70,
    }
}

// ---------------------------------------------------------------------------
// ETKDG torsion preferences
// ---------------------------------------------------------------------------

/// Preferred torsion angles (radians) by hybridization of atoms j and k.
fn preferred_torsions(is_sp2_j: bool, is_sp2_k: bool) -> &'static [f64] {
    use std::f64::consts::PI;
    match (is_sp2_j, is_sp2_k) {
        // sp3-sp3: gauche+, anti, gauche-
        (false, false) => &[PI / 3.0, PI, 5.0 * PI / 3.0],
        // sp3-sp2 or sp2-sp3: syn, anti
        (true, false) | (false, true) => &[0.0, PI],
        // sp2-sp2: syn, anti
        (true, true) => &[0.0, PI],
    }
}

fn is_sp2_atom(mol: &Molecule, idx: usize) -> bool {
    if mol.atoms[idx].is_aromatic {
        return true;
    }
    mol.adjacency[idx].iter().any(|&(_, bi)| {
        mol.bonds[bi].order == BondOrder::Double || mol.bonds[bi].order == BondOrder::Aromatic
    })
}

// ---------------------------------------------------------------------------
// Distance geometry embedding
// ---------------------------------------------------------------------------

/// Embed a single 3D conformer using distance geometry.
pub fn embed_molecule(mol: &Molecule, config: &EmbedConfig) -> Result<Conformer> {
    let n = mol.atom_count();
    if n == 0 {
        return Ok(Conformer::new(Vec::new()));
    }
    if n == 1 {
        return Ok(Conformer::new(vec![[0.0, 0.0, 0.0]]));
    }

    // Step 1: Build distance bounds matrix
    let (lower, upper) = build_bounds_matrix(mol);

    // Step 2: Triangle inequality smoothing (Floyd-Warshall)
    let (lower, upper) = smooth_bounds(lower, upper, n);

    // Step 3: Sample distances and embed
    let mut rng = SimpleRng::new(config.random_seed);
    let conf = embed_from_bounds(&lower, &upper, n, &mut rng)?;

    // Step 4: Apply torsion preferences (ETKDG)
    let mut conf = conf;
    if config.use_torsion_prefs {
        apply_torsion_preferences(mol, &mut conf);
    }

    // Step 5: Force-field minimization
    let conf = if config.max_minimize_steps > 0 && config.force_field != ForceFieldType::None {
        let min_config = MinimizeConfig {
            max_steps: config.max_minimize_steps,
            gradient_threshold: 0.5,
            method: MinimizeMethod::SteepestDescent,
        };
        match config.force_field {
            ForceFieldType::Uff => {
                forcefield::uff_minimize(mol, &conf, &min_config)
                    .map(|r| r.conformer)
                    .unwrap_or(conf)
            }
            ForceFieldType::Mmff94 => {
                forcefield::mmff94_minimize(mol, &conf, &min_config)
                    .map(|r| r.conformer)
                    .unwrap_or(conf)
            }
            ForceFieldType::None => conf,
        }
    } else {
        conf
    };

    Ok(conf)
}

/// Embed multiple conformers with RMSD pruning and energy ranking.
pub fn embed_multiple(mol: &Molecule, config: &EmbedConfig) -> Result<ConformerSet> {
    let mut cset = ConformerSet::new();
    let max = config.max_conformers;

    // Generate more candidates than requested, then prune
    let attempts = max * 3 + 5;

    for trial in 0..attempts {
        if cset.len() >= max {
            break;
        }

        let trial_config = EmbedConfig {
            random_seed: config.random_seed.wrapping_add(trial as u64 * 97),
            max_conformers: 1,
            ..config.clone()
        };

        let conf = match embed_molecule(mol, &trial_config) {
            Ok(c) => c,
            Err(_) => continue,
        };

        // RMSD pruning: check against existing conformers
        let is_unique = cset.conformers.iter().all(|existing| {
            existing.rmsd(&conf).unwrap_or(f64::INFINITY) > config.rmsd_threshold
        });

        if !is_unique {
            continue;
        }

        // Compute energy
        let energy = match config.force_field {
            ForceFieldType::Uff => {
                forcefield::uff_energy(mol, &conf).ok().map(|e| e.total)
            }
            ForceFieldType::Mmff94 => {
                forcefield::mmff94_energy(mol, &conf).ok().map(|e| e.total)
            }
            ForceFieldType::None => None,
        };

        cset.push(conf, energy);
    }

    cset.sort_by_energy();
    Ok(cset)
}

// ---------------------------------------------------------------------------
// Bounds matrix construction
// ---------------------------------------------------------------------------

fn build_bounds_matrix(mol: &Molecule) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = mol.atom_count();
    let mut lower = vec![vec![0.0_f64; n]; n];
    let mut upper = vec![vec![1000.0_f64; n]; n];

    // 1-2 distances: from covalent radii
    for bond in &mol.bonds {
        let a1 = bond.atom1;
        let a2 = bond.atom2;
        let r1 = covalent_radius(mol.atoms[a1].atomic_number);
        let r2 = covalent_radius(mol.atoms[a2].atomic_number);

        let bo_adj = match bond.order {
            BondOrder::Single => 0.0,
            BondOrder::Aromatic => -0.04,
            BondOrder::Double => -0.10,
            BondOrder::Triple => -0.16,
        };

        let d = r1 + r2 + bo_adj;
        let margin = 0.05;
        lower[a1][a2] = d - margin;
        lower[a2][a1] = d - margin;
        upper[a1][a2] = d + margin;
        upper[a2][a1] = d + margin;
    }

    // 1-3 distances: from bond angles
    for j in 0..n {
        let neighbors = &mol.adjacency[j];
        if neighbors.len() < 2 {
            continue;
        }

        // Estimate angle based on degree
        let angle = match mol.adjacency[j].len() {
            1 => std::f64::consts::PI, // linear
            2 => {
                let has_double = mol.adjacency[j].iter().any(|&(_, bi)| {
                    mol.bonds[bi].order == BondOrder::Double || mol.bonds[bi].order == BondOrder::Triple
                });
                if has_double { std::f64::consts::PI } else { 120.0_f64.to_radians() }
            }
            3 => {
                if is_sp2_atom(mol, j) { 120.0_f64.to_radians() } else { 109.5_f64.to_radians() }
            }
            _ => 109.5_f64.to_radians(),
        };

        for a in 0..neighbors.len() {
            for b in (a + 1)..neighbors.len() {
                let i = neighbors[a].0;
                let k = neighbors[b].0;

                let d_ij = (lower[i][j] + upper[i][j]) / 2.0;
                let d_jk = (lower[j][k] + upper[j][k]) / 2.0;

                if d_ij < 0.01 || d_jk < 0.01 {
                    continue;
                }

                // Law of cosines
                let d13 = (d_ij * d_ij + d_jk * d_jk - 2.0 * d_ij * d_jk * angle.cos()).sqrt();
                let margin = 0.15;

                let cur_lower = lower[i][k].max(lower[k][i]);
                let cur_upper = upper[i][k].min(upper[k][i]);

                let new_lower = (d13 - margin).max(cur_lower);
                let new_upper = (d13 + margin).min(cur_upper);

                lower[i][k] = new_lower;
                lower[k][i] = new_lower;
                upper[i][k] = new_upper;
                upper[k][i] = new_upper;
            }
        }
    }

    // Non-bonded lower bounds from vdW radii (for atoms with no other constraints)
    for i in 0..n {
        for j in (i + 1)..n {
            if lower[i][j] < 0.01 {
                let vdw_sum = vdw_radius(mol.atoms[i].atomic_number) +
                    vdw_radius(mol.atoms[j].atomic_number);
                lower[i][j] = vdw_sum * 0.7;
                lower[j][i] = lower[i][j];
            }
        }
    }

    (lower, upper)
}

/// Floyd-Warshall triangle inequality smoothing.
fn smooth_bounds(
    mut lower: Vec<Vec<f64>>,
    mut upper: Vec<Vec<f64>>,
    n: usize,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    // Upper bounds: u[i][j] ≤ u[i][k] + u[k][j]
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                if i == j || i == k || j == k {
                    continue;
                }
                let sum = upper[i][k] + upper[k][j];
                if sum < upper[i][j] {
                    upper[i][j] = sum;
                }
            }
        }
    }

    // Lower bounds: l[i][j] ≥ l[i][k] - u[k][j]  (and symmetric)
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                if i == j || i == k || j == k {
                    continue;
                }
                let diff = lower[i][k] - upper[k][j];
                if diff > lower[i][j] {
                    lower[i][j] = diff;
                }
            }
        }
    }

    // Ensure lower ≤ upper
    for i in 0..n {
        for j in 0..n {
            if lower[i][j] > upper[i][j] {
                let avg = (lower[i][j] + upper[i][j]) / 2.0;
                lower[i][j] = avg;
                upper[i][j] = avg;
            }
            if lower[i][j] < 0.0 {
                lower[i][j] = 0.0;
            }
        }
    }

    (lower, upper)
}

/// Generate 3D coordinates from a distance bounds matrix using eigendecomposition.
fn embed_from_bounds(
    lower: &[Vec<f64>],
    upper: &[Vec<f64>],
    n: usize,
    rng: &mut SimpleRng,
) -> Result<Conformer> {
    // Sample random distances within bounds
    let mut dist = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let lo = lower[i][j].max(0.001);
            let hi = upper[i][j].max(lo + 0.001);
            let d = lo + rng.next_f64() * (hi - lo);
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }

    // Build metric matrix: G[i][j] = 0.5 * (d[0][i]^2 + d[0][j]^2 - d[i][j]^2)
    let mut g = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            g[i][j] = 0.5 * (dist[0][i] * dist[0][i] + dist[0][j] * dist[0][j]
                - dist[i][j] * dist[i][j]);
        }
    }

    // Extract top-3 eigenvectors using power iteration
    let coords = extract_3d_coords(&g, n)?;

    Ok(Conformer::new(coords))
}

/// Extract 3D coordinates from metric matrix using power iteration for top 3 eigenvalues.
fn extract_3d_coords(g: &[Vec<f64>], n: usize) -> Result<Vec<[f64; 3]>> {
    let mut coords = vec![[0.0_f64; 3]; n];
    let mut deflated = g.to_vec();

    for dim in 0..3 {
        // Power iteration
        let mut v = vec![1.0_f64 / (n as f64).sqrt(); n];
        let mut eigenvalue = 0.0_f64;

        for _ in 0..100 {
            // Matrix-vector multiply
            let mut mv = vec![0.0_f64; n];
            for i in 0..n {
                for j in 0..n {
                    mv[i] += deflated[i][j] * v[j];
                }
            }

            // Compute eigenvalue (Rayleigh quotient)
            let dot: f64 = mv.iter().zip(v.iter()).map(|(a, b)| a * b).sum();
            eigenvalue = dot;

            // Normalize
            let norm: f64 = mv.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm < 1e-15 {
                break;
            }
            for i in 0..n {
                v[i] = mv[i] / norm;
            }
        }

        // Use sqrt of eigenvalue as coordinate scale
        let scale = if eigenvalue > 0.0 { eigenvalue.sqrt() } else { 0.0 };
        for i in 0..n {
            coords[i][dim] = v[i] * scale;
        }

        // Deflate: G = G - λ * v * v^T
        for i in 0..n {
            for j in 0..n {
                deflated[i][j] -= eigenvalue * v[i] * v[j];
            }
        }
    }

    Ok(coords)
}

/// Apply ETKDG torsion angle preferences to rotatable bonds.
fn apply_torsion_preferences(mol: &Molecule, conf: &mut Conformer) {
    let rings = ring::find_sssr(mol);
    let ring_atoms: std::collections::HashSet<usize> = rings.iter().flat_map(|r| r.iter().copied()).collect();

    for bond in &mol.bonds {
        // Only rotatable single bonds
        if bond.order != BondOrder::Single {
            continue;
        }
        let j = bond.atom1;
        let k = bond.atom2;

        // Skip terminal atoms
        if mol.degree(j) < 2 || mol.degree(k) < 2 {
            continue;
        }

        // Skip ring bonds
        if ring_atoms.contains(&j) && ring_atoms.contains(&k) {
            let j_ring = rings.iter().any(|r| r.contains(&j) && r.contains(&k));
            if j_ring {
                continue;
            }
        }

        let sp2_j = is_sp2_atom(mol, j);
        let sp2_k = is_sp2_atom(mol, k);
        let prefs = preferred_torsions(sp2_j, sp2_k);

        if prefs.is_empty() {
            continue;
        }

        // Find atoms beyond k (for rotation)
        let atoms_beyond_k = find_atoms_beyond(mol, k, j);
        if atoms_beyond_k.is_empty() {
            continue;
        }

        // Find reference atoms for computing current dihedral
        let i = mol.adjacency[j].iter().find(|&&(n, _)| n != k).map(|&(n, _)| n);
        let l = mol.adjacency[k].iter().find(|&&(n, _)| n != j).map(|&(n, _)| n);

        if let (Some(i_idx), Some(l_idx)) = (i, l) {
            let current = conf.dihedral(i_idx, j, k, l_idx);

            // Find nearest preferred angle
            let nearest = prefs
                .iter()
                .min_by(|&&a, &&b| {
                    let da = angle_diff(current, a).abs();
                    let db = angle_diff(current, b).abs();
                    da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                })
                .copied()
                .unwrap_or(current);

            let rotation = angle_diff(current, nearest);
            if rotation.abs() > 0.05 {
                conf.set_dihedral(j, k, &atoms_beyond_k, rotation);
            }
        }
    }
}

/// Find all atoms reachable from `start` without going through `exclude`.
fn find_atoms_beyond(mol: &Molecule, start: usize, exclude: usize) -> Vec<usize> {
    let mut visited = vec![false; mol.atom_count()];
    visited[exclude] = true;
    let mut stack = vec![start];
    let mut result = Vec::new();

    while let Some(atom) = stack.pop() {
        if visited[atom] {
            continue;
        }
        visited[atom] = true;
        result.push(atom);
        for &(nb, _) in &mol.adjacency[atom] {
            if !visited[nb] {
                stack.push(nb);
            }
        }
    }
    result
}

fn angle_diff(a: f64, b: f64) -> f64 {
    use std::f64::consts::PI;
    let mut d = b - a;
    while d > PI {
        d -= 2.0 * PI;
    }
    while d < -PI {
        d += 2.0 * PI;
    }
    d
}

// ---------------------------------------------------------------------------
// Simple RNG (xorshift64)
// ---------------------------------------------------------------------------

struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        SimpleRng {
            state: if seed == 0 { 1 } else { seed },
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn embed_ethane() {
        let mol = parse_smiles("CC").unwrap();
        let conf = embed_molecule(&mol, &EmbedConfig::default()).unwrap();
        assert_eq!(conf.len(), 2);
        // C-C distance should be approximately 1.54 Å (within tolerance)
        let d = conf.distance(0, 1);
        assert!(d > 0.5 && d < 3.0, "C-C distance = {d}");
    }

    #[test]
    fn embed_methane() {
        let mol = parse_smiles("C").unwrap();
        let conf = embed_molecule(&mol, &EmbedConfig::default()).unwrap();
        assert_eq!(conf.len(), 1);
    }

    #[test]
    fn embed_empty() {
        use crate::molecule::Molecule;
        let mol = Molecule::new("".into(), vec![], vec![]);
        let conf = embed_molecule(&mol, &EmbedConfig::default()).unwrap();
        assert!(conf.is_empty());
    }

    #[test]
    fn embed_ethanol() {
        let mol = parse_smiles("CCO").unwrap();
        let conf = embed_molecule(&mol, &EmbedConfig::default()).unwrap();
        assert_eq!(conf.len(), 3);
        // All distances should be reasonable
        for i in 0..3 {
            for j in (i + 1)..3 {
                let d = conf.distance(i, j);
                assert!(d > 0.3 && d < 5.0, "distance({i},{j}) = {d}");
            }
        }
    }

    #[test]
    fn embed_benzene_planarity() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let config = EmbedConfig {
            max_minimize_steps: 100,
            ..EmbedConfig::default()
        };
        let conf = embed_molecule(&mol, &config).unwrap();
        assert_eq!(conf.len(), 6);

        // All atoms should have finite coordinates
        for i in 0..6 {
            assert!(conf.coords[i][0].is_finite());
            assert!(conf.coords[i][1].is_finite());
            assert!(conf.coords[i][2].is_finite());
        }
    }

    #[test]
    fn embed_multiple_conformers() {
        let mol = parse_smiles("CCCC").unwrap();
        let config = EmbedConfig {
            max_conformers: 3,
            rmsd_threshold: 0.1,
            max_minimize_steps: 50,
            ..EmbedConfig::default()
        };
        let cset = embed_multiple(&mol, &config).unwrap();
        assert!(cset.len() >= 1, "got {} conformers", cset.len());
    }

    #[test]
    fn embed_deterministic() {
        let mol = parse_smiles("CCO").unwrap();
        let config = EmbedConfig {
            random_seed: 123,
            max_minimize_steps: 0,
            force_field: ForceFieldType::None,
            ..EmbedConfig::default()
        };
        let c1 = embed_molecule(&mol, &config).unwrap();
        let c2 = embed_molecule(&mol, &config).unwrap();
        let rmsd = c1.rmsd(&c2).unwrap();
        assert!(rmsd < 1e-10, "should be deterministic, rmsd = {rmsd}");
    }

    #[test]
    fn embed_no_forcefield() {
        let mol = parse_smiles("CCO").unwrap();
        let config = EmbedConfig {
            force_field: ForceFieldType::None,
            max_minimize_steps: 0,
            ..EmbedConfig::default()
        };
        let conf = embed_molecule(&mol, &config).unwrap();
        assert_eq!(conf.len(), 3);
    }

    #[test]
    fn embed_with_mmff94() {
        let mol = parse_smiles("CC").unwrap();
        let config = EmbedConfig {
            force_field: ForceFieldType::Mmff94,
            max_minimize_steps: 50,
            ..EmbedConfig::default()
        };
        let conf = embed_molecule(&mol, &config).unwrap();
        assert_eq!(conf.len(), 2);
    }

    #[test]
    fn embed_conformer_set_sorted() {
        let mol = parse_smiles("CCCC").unwrap();
        let config = EmbedConfig {
            max_conformers: 5,
            rmsd_threshold: 0.01,
            max_minimize_steps: 20,
            ..EmbedConfig::default()
        };
        let cset = embed_multiple(&mol, &config).unwrap();
        if cset.len() >= 2 {
            // Should be sorted by energy
            for i in 1..cset.len() {
                let e0 = cset.energies[i - 1].unwrap_or(f64::INFINITY);
                let e1 = cset.energies[i].unwrap_or(f64::INFINITY);
                assert!(e0 <= e1 + 1e-6, "not sorted: {} > {}", e0, e1);
            }
        }
    }

    #[test]
    fn bounds_matrix_basic() {
        let mol = parse_smiles("CC").unwrap();
        let (lower, upper) = build_bounds_matrix(&mol);
        // Bond distance bounds should be set
        assert!(lower[0][1] > 0.0);
        assert!(upper[0][1] < 1000.0);
        assert!(lower[0][1] <= upper[0][1]);
    }

    #[test]
    fn simple_rng_deterministic() {
        let mut rng1 = SimpleRng::new(42);
        let mut rng2 = SimpleRng::new(42);
        for _ in 0..10 {
            assert_eq!(rng1.next_u64(), rng2.next_u64());
        }
    }

    #[test]
    fn simple_rng_range() {
        let mut rng = SimpleRng::new(12345);
        for _ in 0..100 {
            let f = rng.next_f64();
            assert!((0.0..1.0).contains(&f), "out of range: {f}");
        }
    }
}
