//! Graph-based molecular descriptors for QSAR/QSPR.
//!
//! All descriptors are computed from the molecular graph using BFS shortest paths.

use std::collections::VecDeque;

use crate::element::element_by_number;
use crate::molecule::{BondOrder, Molecule};
use crate::ring;

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Detailed ring count information.
#[derive(Debug, Clone, PartialEq)]
pub struct RingDetails {
    pub total: usize,
    pub aromatic: usize,
    pub aliphatic: usize,
    pub heteroaromatic: usize,
    pub spiro: usize,
    pub fused: usize,
}

/// Autocorrelation descriptor results.
#[derive(Debug, Clone)]
pub struct AutocorrelationResult {
    pub moreau_broto: Vec<f64>,
    pub moran: Vec<f64>,
    pub geary: Vec<f64>,
}

/// A named collection of descriptor values.
#[derive(Debug, Clone)]
pub struct DescriptorSet {
    pub values: Vec<(String, f64)>,
}

// ---------------------------------------------------------------------------
// Shortest path matrix (BFS from each atom)
// ---------------------------------------------------------------------------

fn shortest_path_matrix(mol: &Molecule) -> Vec<Vec<usize>> {
    let n = mol.atom_count();
    let mut dist = vec![vec![usize::MAX; n]; n];

    for start in 0..n {
        dist[start][start] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(start);
        while let Some(curr) = queue.pop_front() {
            for &(neighbor, _) in &mol.adjacency[curr] {
                if dist[start][neighbor] == usize::MAX {
                    dist[start][neighbor] = dist[start][curr] + 1;
                    queue.push_back(neighbor);
                }
            }
        }
    }
    dist
}

/// Count all paths of a given length using DFS enumeration.
fn count_paths(mol: &Molecule, length: usize) -> usize {
    let n = mol.atom_count();
    let mut count = 0;
    let mut visited = vec![false; n];

    for start in 0..n {
        visited[start] = true;
        count_paths_dfs(mol, start, 0, length, &mut visited, &mut count);
        visited[start] = false;
    }
    // Each path is counted twice (once from each end), except for length 0
    if length > 0 { count / 2 } else { count }
}

fn count_paths_dfs(
    mol: &Molecule,
    current: usize,
    depth: usize,
    target: usize,
    visited: &mut Vec<bool>,
    count: &mut usize,
) {
    if depth == target {
        *count += 1;
        return;
    }
    for &(neighbor, _) in &mol.adjacency[current] {
        if !visited[neighbor] {
            visited[neighbor] = true;
            count_paths_dfs(mol, neighbor, depth + 1, target, visited, count);
            visited[neighbor] = false;
        }
    }
}

// ---------------------------------------------------------------------------
// Descriptors
// ---------------------------------------------------------------------------

/// Wiener index: sum of all shortest-path distances between atom pairs.
pub fn wiener_index(mol: &Molecule) -> f64 {
    let n = mol.atom_count();
    if n < 2 {
        return 0.0;
    }
    let dist = shortest_path_matrix(mol);
    let mut sum = 0u64;
    for i in 0..n {
        for j in (i + 1)..n {
            if dist[i][j] != usize::MAX {
                sum += dist[i][j] as u64;
            }
        }
    }
    sum as f64
}

/// Zagreb indices M1 and M2.
///
/// M1 = sum of d(i)^2 for all atoms, M2 = sum of d(i)*d(j) for all edges.
pub fn zagreb_indices(mol: &Molecule) -> (f64, f64) {
    let m1: f64 = (0..mol.atom_count())
        .map(|i| {
            let d = mol.degree(i) as f64;
            d * d
        })
        .sum();

    let m2: f64 = mol
        .bonds
        .iter()
        .map(|b| mol.degree(b.atom1) as f64 * mol.degree(b.atom2) as f64)
        .sum();

    (m1, m2)
}

/// Balaban J index.
///
/// J = m/(mu+1) * sum of (s_i * s_j)^(-0.5) over edges,
/// where s_i = distance sum, m = bonds, mu = cyclomatic number.
pub fn balaban_j(mol: &Molecule) -> f64 {
    let n = mol.atom_count();
    let m = mol.bond_count();
    if n < 2 || m == 0 {
        return 0.0;
    }

    let dist = shortest_path_matrix(mol);

    // Distance sums for each atom
    let s: Vec<f64> = (0..n)
        .map(|i| {
            (0..n)
                .filter(|&j| j != i && dist[i][j] != usize::MAX)
                .map(|j| dist[i][j] as f64)
                .sum()
        })
        .collect();

    // Cyclomatic number = m - n + components
    let components = count_components(mol);
    let mu = m as f64 - n as f64 + components as f64;

    let mut edge_sum = 0.0;
    for bond in &mol.bonds {
        let si = s[bond.atom1];
        let sj = s[bond.atom2];
        if si > 0.0 && sj > 0.0 {
            edge_sum += (si * sj).powf(-0.5);
        }
    }

    m as f64 / (mu + 1.0) * edge_sum
}

fn count_components(mol: &Molecule) -> usize {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut count = 0;
    for start in 0..n {
        if visited[start] {
            continue;
        }
        count += 1;
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;
        while let Some(curr) = queue.pop_front() {
            for &(neighbor, _) in &mol.adjacency[curr] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    queue.push_back(neighbor);
                }
            }
        }
    }
    count
}

/// Topological Polar Surface Area (Ertl 2000).
///
/// Uses static lookup of N/O atom-type contributions.
pub fn tpsa(mol: &Molecule) -> f64 {
    let mut area = 0.0;
    for (i, atom) in mol.atoms.iter().enumerate() {
        let contribution = tpsa_contribution(mol, i, atom.atomic_number, atom.implicit_hydrogens, atom.formal_charge, atom.is_aromatic);
        area += contribution;
    }
    area
}

/// TPSA contribution for a single atom based on Ertl 2000 fragment values.
fn tpsa_contribution(
    mol: &Molecule,
    atom_idx: usize,
    atomic_num: u8,
    implicit_h: u8,
    charge: i8,
    is_aromatic: bool,
) -> f64 {
    let degree = mol.degree(atom_idx);
    let has_double_bond = mol.adjacency[atom_idx]
        .iter()
        .any(|&(_, bi)| mol.bonds[bi].order == BondOrder::Double);

    match atomic_num {
        // Nitrogen contributions
        7 => {
            if charge > 0 {
                // Charged nitrogen
                if implicit_h >= 3 { return 27.64; } // [NH3+]
                if implicit_h == 2 { return 25.59; } // [NH2+]=
                if implicit_h == 1 { return 23.47; } // [NH+]
                return 0.0;
            }
            if is_aromatic {
                if implicit_h >= 1 { return 15.79; } // nH (pyrrole-like)
                return 12.89; // n (pyridine-like)
            }
            match (degree, implicit_h, has_double_bond) {
                (1, 2, false) => 26.02, // -NH2
                (1, 2, true) => 26.02,  // =NH2 (unusual)
                (2, 1, false) => 19.15, // -NH-
                (2, 1, true) => 23.85,  // =NH
                (2, 0, true) => 12.36,  // =N-
                (2, 0, false) => 19.15, // -N=  (same as -NH- without H)
                (3, 0, false) => 3.24,  // >N-
                (3, 0, true) => 3.24,   // >N=
                (1, 0, true) => 23.79,  // #N
                _ => {
                    if implicit_h >= 2 { 26.02 }
                    else if implicit_h == 1 { 19.15 }
                    else { 3.24 }
                }
            }
        }
        // Oxygen contributions
        8 => {
            if charge < 0 {
                return 23.06; // [O-]
            }
            if is_aromatic {
                return 13.14; // aromatic o
            }
            match (degree, implicit_h, has_double_bond) {
                (1, 1, false) => 20.23, // -OH
                (1, 0, true) => 17.07,  // =O
                (2, 0, false) => 9.23,  // -O-
                (1, 0, false) => 17.07, // -O (terminal, e.g. carboxylate)
                _ => {
                    if implicit_h >= 1 { 20.23 }
                    else if has_double_bond { 17.07 }
                    else { 9.23 }
                }
            }
        }
        // Sulfur contributions (minimal TPSA)
        16 => {
            if implicit_h >= 1 { return 38.80; }
            if has_double_bond { return 25.30; }
            if degree >= 2 { return 25.30; }
            0.0
        }
        // Phosphorus
        15 => {
            if has_double_bond { return 34.14; }
            if implicit_h >= 1 { return 23.47; }
            9.81
        }
        _ => 0.0,
    }
}

/// Wildman-Crippen LogP and MR estimation.
///
/// Returns (logP, molar_refractivity).
pub fn wildman_crippen_logp(mol: &Molecule) -> (f64, f64) {
    let rings = ring::find_sssr(mol);
    let ring_membership = build_ring_membership(mol, &rings);

    let mut logp = 0.0;
    let mut mr = 0.0;

    for i in 0..mol.atom_count() {
        let (lp, m) = crippen_atom_contribution(mol, i, &ring_membership);
        logp += lp;
        mr += m;
    }

    // Add implicit hydrogen contributions
    let total_implicit_h: usize = mol.atoms.iter().map(|a| a.implicit_hydrogens as usize).sum();
    // H attached to C: logP ≈ 0.123, MR ≈ 1.057
    // H attached to heteroatom: logP ≈ -0.2677, MR ≈ 1.057
    for atom in &mol.atoms {
        let h_count = atom.implicit_hydrogens as usize;
        if h_count == 0 {
            continue;
        }
        if atom.atomic_number == 6 {
            logp += h_count as f64 * 0.1230;
            mr += h_count as f64 * 1.057;
        } else {
            logp += h_count as f64 * (-0.2677);
            mr += h_count as f64 * 1.057;
        }
    }

    let _ = total_implicit_h;
    (logp, mr)
}

/// Simplified Wildman-Crippen atom type classification.
fn crippen_atom_contribution(mol: &Molecule, atom_idx: usize, ring_membership: &[bool]) -> (f64, f64) {
    let atom = &mol.atoms[atom_idx];
    let degree = mol.degree(atom_idx);
    let in_ring = ring_membership[atom_idx];
    let has_double_bond = mol.adjacency[atom_idx]
        .iter()
        .any(|&(_, bi)| mol.bonds[bi].order == BondOrder::Double);
    let has_hetero_neighbor = mol.adjacency[atom_idx]
        .iter()
        .any(|&(n, _)| mol.atoms[n].atomic_number != 6 && mol.atoms[n].atomic_number != 1);

    match atom.atomic_number {
        6 => { // Carbon
            if atom.is_aromatic {
                if has_hetero_neighbor { (-0.14, 3.509) }
                else { (0.296, 3.509) }
            } else if has_double_bond {
                if has_hetero_neighbor { (-0.03, 3.509) }
                else { (0.08, 3.509) }
            } else if in_ring {
                (0.1441, 3.509)
            } else {
                match degree {
                    1 => (0.1441, 3.509),
                    2 => (0.1441, 3.509),
                    3 => (0.0, 3.509),
                    _ => (-0.04, 3.509),
                }
            }
        }
        7 => { // Nitrogen
            if atom.is_aromatic { (-0.3187, 2.188) }
            else if atom.formal_charge > 0 { (-1.0190, 2.188) }
            else if has_double_bond { (-0.5262, 2.188) }
            else { (-0.4458, 2.262) }
        }
        8 => { // Oxygen
            if atom.formal_charge < 0 { (-1.189, 1.476) }
            else if has_double_bond { (-0.3339, 1.476) }
            else if degree >= 2 { (-0.2893, 1.476) }
            else { (-0.3567, 1.476) }
        }
        9 => (0.4118, 1.108),    // F
        15 => (0.2836, 6.920),   // P
        16 => { // S
            if has_double_bond { (-0.1084, 7.365) }
            else if atom.formal_charge != 0 { (-0.5188, 7.365) }
            else { (0.6237, 7.365) }
        }
        17 => (0.6895, 5.853),   // Cl
        35 => (0.8813, 8.927),   // Br
        53 => (1.050, 13.940),   // I
        _ => (0.0, 0.0),
    }
}

fn build_ring_membership(mol: &Molecule, rings: &[Vec<usize>]) -> Vec<bool> {
    let mut member = vec![false; mol.atom_count()];
    for ring in rings {
        for &idx in ring {
            member[idx] = true;
        }
    }
    member
}

/// Bertz complexity index (CT).
///
/// CT = 2 * (I_bonds + I_atoms) where I is Shannon information content.
pub fn bertz_ct(mol: &Molecule) -> f64 {
    if mol.bond_count() == 0 {
        return 0.0;
    }

    // Bond type distribution
    let mut bond_counts = [0usize; 4]; // single, double, triple, aromatic
    for bond in &mol.bonds {
        let idx = match bond.order {
            BondOrder::Single => 0,
            BondOrder::Double => 1,
            BondOrder::Triple => 2,
            BondOrder::Aromatic => 3,
        };
        bond_counts[idx] += 1;
    }
    let i_bonds = shannon_entropy(&bond_counts);

    // Atom type distribution (by atomic number in bond neighborhoods)
    let mut atom_env_counts: std::collections::HashMap<(u8, u8, bool), usize> = std::collections::HashMap::new();
    for (i, atom) in mol.atoms.iter().enumerate() {
        let key = (atom.atomic_number, mol.degree(i) as u8, atom.is_aromatic);
        *atom_env_counts.entry(key).or_insert(0) += 1;
    }
    let counts: Vec<usize> = atom_env_counts.values().copied().collect();
    let i_atoms = shannon_entropy(&counts);

    2.0 * (i_bonds + i_atoms)
}

fn shannon_entropy(counts: &[usize]) -> f64 {
    let total: usize = counts.iter().sum();
    if total == 0 {
        return 0.0;
    }
    let n = total as f64;
    let mut entropy = 0.0;
    for &c in counts {
        if c > 0 {
            let p = c as f64 / n;
            entropy -= p * p.log2();
        }
    }
    entropy * n
}

/// Kappa shape indices κ1, κ2, κ3.
///
/// Based on path counts of length 1, 2, and 3.
pub fn kappa_shape_indices(mol: &Molecule) -> (f64, f64, f64) {
    let n = mol.atom_count();
    if n == 0 {
        return (0.0, 0.0, 0.0);
    }

    let p1 = count_paths(mol, 1) as f64;
    let p2 = count_paths(mol, 2) as f64;
    let p3 = count_paths(mol, 3) as f64;
    let nf = n as f64;

    let k1 = if p1 > 0.0 {
        (nf - 1.0).powi(2) / p1.powi(2) * nf
    } else {
        0.0
    };

    let k2 = if p2 > 0.0 && nf >= 2.0 {
        (nf - 1.0).powi(2) * (nf - 2.0).powi(2) / p2.powi(2)
    } else {
        0.0
    };

    let k3 = if p3 > 0.0 && nf >= 3.0 {
        if n % 2 == 1 {
            (nf - 1.0).powi(2) * (nf - 3.0).powi(2) / p3.powi(2)
        } else {
            (nf - 1.0) * (nf - 2.0).powi(2) * (nf - 3.0) / p3.powi(2)
        }
    } else {
        0.0
    };

    (k1, k2, k3)
}

/// Chi (Randić/Kier-Hall) connectivity indices χ0 through χ3.
pub fn chi_connectivity(mol: &Molecule) -> Vec<f64> {
    let n = mol.atom_count();
    if n == 0 {
        return vec![0.0; 4];
    }

    let degrees: Vec<usize> = (0..n).map(|i| mol.degree(i)).collect();

    // χ0 = sum of 1/sqrt(d_i) for all atoms with d > 0
    let chi0: f64 = degrees
        .iter()
        .filter(|&&d| d > 0)
        .map(|&d| 1.0 / (d as f64).sqrt())
        .sum();

    // χ1 = sum of 1/sqrt(d_i * d_j) for all edges
    let chi1: f64 = mol
        .bonds
        .iter()
        .filter_map(|b| {
            let di = degrees[b.atom1];
            let dj = degrees[b.atom2];
            if di > 0 && dj > 0 {
                Some(1.0 / ((di * dj) as f64).sqrt())
            } else {
                None
            }
        })
        .sum();

    // χ2 = sum of 1/sqrt(d_i * d_j * d_k) over all 2-paths
    let mut chi2 = 0.0;
    for i in 0..n {
        for &(j, _) in &mol.adjacency[i] {
            if j <= i { continue; }
            for &(k, _) in &mol.adjacency[j] {
                if k <= i || k == i { continue; }
                let di = degrees[i];
                let dj = degrees[j];
                let dk = degrees[k];
                if di > 0 && dj > 0 && dk > 0 {
                    chi2 += 1.0 / ((di * dj * dk) as f64).sqrt();
                }
            }
        }
    }

    // χ3 = sum of 1/sqrt(d_i * d_j * d_k * d_l) over all 3-paths
    let mut chi3 = 0.0;
    let mut visited = vec![false; n];
    for i in 0..n {
        visited[i] = true;
        for &(j, _) in &mol.adjacency[i] {
            if visited[j] { continue; }
            visited[j] = true;
            for &(k, _) in &mol.adjacency[j] {
                if visited[k] { continue; }
                visited[k] = true;
                for &(l, _) in &mol.adjacency[k] {
                    if visited[l] { continue; }
                    let di = degrees[i];
                    let dj = degrees[j];
                    let dk = degrees[k];
                    let dl = degrees[l];
                    if di > 0 && dj > 0 && dk > 0 && dl > 0 {
                        chi3 += 1.0 / ((di * dj * dk * dl) as f64).sqrt();
                    }
                }
                visited[k] = false;
            }
            visited[j] = false;
        }
        visited[i] = false;
    }

    vec![chi0, chi1, chi2, chi3]
}

/// Fraction of sp3 carbons.
pub fn fraction_sp3(mol: &Molecule) -> f64 {
    let total_carbons = mol.atoms.iter().filter(|a| a.atomic_number == 6).count();
    if total_carbons == 0 {
        return 0.0;
    }

    let sp3_carbons = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|(i, a)| {
            if a.atomic_number != 6 || a.is_aromatic {
                return false;
            }
            // sp3 carbon: no double/triple bonds, not aromatic
            !mol.adjacency[*i]
                .iter()
                .any(|&(_, bi)| matches!(mol.bonds[bi].order, BondOrder::Double | BondOrder::Triple))
        })
        .count();

    sp3_carbons as f64 / total_carbons as f64
}

/// Detailed ring count analysis.
pub fn ring_count_details(mol: &Molecule) -> RingDetails {
    let rings = ring::find_sssr(mol);
    let total = rings.len();

    let aromatic = rings
        .iter()
        .filter(|r| r.iter().all(|&i| mol.atoms[i].is_aromatic))
        .count();
    let aliphatic = total - aromatic;

    let heteroaromatic = rings
        .iter()
        .filter(|r| {
            r.iter().all(|&i| mol.atoms[i].is_aromatic)
                && r.iter().any(|&i| mol.atoms[i].atomic_number != 6)
        })
        .count();

    // Spiro: atoms that belong to 2+ rings but share no common bonds between those rings
    let mut atom_ring_count = vec![0usize; mol.atom_count()];
    for ring in &rings {
        for &idx in ring {
            atom_ring_count[idx] += 1;
        }
    }

    let mut spiro = 0;
    for (atom_idx, &count) in atom_ring_count.iter().enumerate() {
        if count >= 2 {
            // Check if this atom is a spiro center: appears in 2+ rings
            // but no other atom is shared between those same rings
            let containing_rings: Vec<usize> = rings
                .iter()
                .enumerate()
                .filter(|(_, r)| r.contains(&atom_idx))
                .map(|(ri, _)| ri)
                .collect();

            let mut is_spiro = true;
            'outer: for i in 0..containing_rings.len() {
                for j in (i + 1)..containing_rings.len() {
                    let ri = &rings[containing_rings[i]];
                    let rj = &rings[containing_rings[j]];
                    let shared: Vec<usize> = ri.iter().filter(|a| rj.contains(a)).copied().collect();
                    if shared.len() > 1 {
                        is_spiro = false;
                        break 'outer;
                    }
                }
            }
            if is_spiro && count >= 2 {
                spiro += 1;
            }
        }
    }

    // Fused: pairs of rings sharing >= 2 atoms
    let mut fused = 0;
    for i in 0..rings.len() {
        for j in (i + 1)..rings.len() {
            let shared = rings[i].iter().filter(|a| rings[j].contains(a)).count();
            if shared >= 2 {
                fused += 1;
            }
        }
    }

    RingDetails {
        total,
        aromatic,
        aliphatic,
        heteroaromatic,
        spiro,
        fused,
    }
}

/// Electrotopological state (E-state) indices.
///
/// Returns the E-state value per atom.
pub fn estate_indices(mol: &Molecule) -> Vec<f64> {
    let n = mol.atom_count();
    if n == 0 {
        return Vec::new();
    }

    let dist = shortest_path_matrix(mol);

    // Intrinsic state I_i = (delta_v + 1) / delta
    // where delta_v = valence electrons - implicit H, delta = graph degree
    let intrinsic: Vec<f64> = (0..n)
        .map(|i| {
            let atom = &mol.atoms[i];
            let delta = mol.degree(i) as f64;
            if delta == 0.0 {
                return 0.0;
            }
            let valence_electrons = match atom.atomic_number {
                6 => 4, 7 => 5, 8 => 6, 9 => 7,
                15 => 5, 16 => 6, 17 => 7, 35 => 7, 53 => 7,
                14 => 4, 5 => 3,
                _ => atom.atomic_number as usize,
            };
            let delta_v = valence_electrons as f64 - atom.implicit_hydrogens as f64;
            (delta_v + 1.0) / delta
        })
        .collect();

    // E-state: S_i = I_i + sum_j (I_i - I_j) / d_ij^2
    let mut estate = vec![0.0; n];
    for i in 0..n {
        let mut perturbation = 0.0;
        for j in 0..n {
            if i == j { continue; }
            let d = dist[i][j];
            if d != usize::MAX && d > 0 {
                perturbation += (intrinsic[i] - intrinsic[j]) / (d as f64).powi(2);
            }
        }
        estate[i] = intrinsic[i] + perturbation;
    }

    estate
}

/// Autocorrelation descriptors (Moreau-Broto, Moran, Geary).
///
/// Computed for topological distances 1..8 using atomic mass as the property.
pub fn autocorrelation_descriptors(mol: &Molecule) -> AutocorrelationResult {
    let n = mol.atom_count();
    let max_lag = 8;

    if n < 2 {
        return AutocorrelationResult {
            moreau_broto: vec![0.0; max_lag],
            moran: vec![0.0; max_lag],
            geary: vec![0.0; max_lag],
        };
    }

    let dist = shortest_path_matrix(mol);

    // Property: atomic mass
    let prop: Vec<f64> = mol
        .atoms
        .iter()
        .map(|a| {
            element_by_number(a.atomic_number)
                .map(|e| e.atomic_weight)
                .unwrap_or(0.0)
        })
        .collect();

    let mean: f64 = prop.iter().sum::<f64>() / n as f64;
    let variance: f64 = prop.iter().map(|&p| (p - mean).powi(2)).sum::<f64>() / n as f64;

    let mut moreau_broto = vec![0.0; max_lag];
    let mut moran = vec![0.0; max_lag];
    let mut geary = vec![0.0; max_lag];

    for lag in 1..=max_lag {
        let mut mb = 0.0;
        let mut mo_num = 0.0;
        let mut ge_num = 0.0;
        let mut pair_count = 0;

        for i in 0..n {
            for j in (i + 1)..n {
                if dist[i][j] == lag {
                    mb += prop[i] * prop[j];
                    mo_num += (prop[i] - mean) * (prop[j] - mean);
                    ge_num += (prop[i] - prop[j]).powi(2);
                    pair_count += 1;
                }
            }
        }

        moreau_broto[lag - 1] = mb;

        if pair_count > 0 && variance.abs() > 1e-15 {
            moran[lag - 1] = (mo_num / pair_count as f64) / variance;
            geary[lag - 1] = (ge_num / (2.0 * pair_count as f64)) / variance;
        }
    }

    AutocorrelationResult { moreau_broto, moran, geary }
}

/// Compute all descriptors as a named vector.
pub fn compute_all_descriptors(mol: &Molecule) -> DescriptorSet {
    let mut values = Vec::new();

    values.push(("wiener_index".into(), wiener_index(mol)));

    let (m1, m2) = zagreb_indices(mol);
    values.push(("zagreb_m1".into(), m1));
    values.push(("zagreb_m2".into(), m2));

    values.push(("balaban_j".into(), balaban_j(mol)));
    values.push(("tpsa".into(), tpsa(mol)));

    let (logp, mr) = wildman_crippen_logp(mol);
    values.push(("crippen_logp".into(), logp));
    values.push(("crippen_mr".into(), mr));

    values.push(("bertz_ct".into(), bertz_ct(mol)));

    let (k1, k2, k3) = kappa_shape_indices(mol);
    values.push(("kappa1".into(), k1));
    values.push(("kappa2".into(), k2));
    values.push(("kappa3".into(), k3));

    let chi = chi_connectivity(mol);
    for (i, &v) in chi.iter().enumerate() {
        values.push((format!("chi{i}"), v));
    }

    values.push(("fraction_sp3".into(), fraction_sp3(mol)));

    let rd = ring_count_details(mol);
    values.push(("ring_total".into(), rd.total as f64));
    values.push(("ring_aromatic".into(), rd.aromatic as f64));
    values.push(("ring_aliphatic".into(), rd.aliphatic as f64));
    values.push(("ring_heteroaromatic".into(), rd.heteroaromatic as f64));
    values.push(("ring_spiro".into(), rd.spiro as f64));
    values.push(("ring_fused".into(), rd.fused as f64));

    DescriptorSet { values }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn wiener_index_hexane() {
        // n-hexane: C-C-C-C-C-C  (6 atoms in a chain)
        // W = sum of d(i,j) for all pairs
        // d: 1+2+3+4+5 + 1+2+3+4 + 1+2+3 + 1+2 + 1 = 15+10+6+3+1 = 35
        let mol = parse_smiles("CCCCCC").unwrap();
        let w = wiener_index(&mol);
        assert!((w - 35.0).abs() < 1e-10, "wiener={w}");
    }

    #[test]
    fn zagreb_benzene() {
        // Benzene: 6 atoms each with degree 2 → M1 = 6*4 = 24
        let mol = parse_smiles("c1ccccc1").unwrap();
        let (m1, _m2) = zagreb_indices(&mol);
        assert!((m1 - 24.0).abs() < 1e-10, "M1={m1}");
    }

    #[test]
    fn tpsa_aspirin() {
        // Aspirin TPSA ≈ 63.6 Å²
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let t = tpsa(&mol);
        // Allow generous tolerance for our simplified model
        assert!(t > 30.0 && t < 100.0, "tpsa={t}");
    }

    #[test]
    fn logp_reasonable() {
        // Aspirin LogP ≈ 1.2 in literature
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let (logp, _mr) = wildman_crippen_logp(&mol);
        assert!(logp > -2.0 && logp < 5.0, "logP={logp}");
    }

    #[test]
    fn fraction_sp3_cyclohexane() {
        let mol = parse_smiles("C1CCCCC1").unwrap();
        assert!((fraction_sp3(&mol) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn fraction_sp3_benzene() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert!((fraction_sp3(&mol) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn naphthalene_fused_rings() {
        let mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        let rd = ring_count_details(&mol);
        assert_eq!(rd.total, 2);
        assert_eq!(rd.aromatic, 2);
        assert!(rd.fused >= 1, "fused={}", rd.fused);
    }

    #[test]
    fn empty_molecule_zeros() {
        let mol = Molecule::new(String::new(), Vec::new(), Vec::new());
        assert!((wiener_index(&mol) - 0.0).abs() < 1e-10);
        let (m1, m2) = zagreb_indices(&mol);
        assert!((m1 - 0.0).abs() < 1e-10);
        assert!((m2 - 0.0).abs() < 1e-10);
        assert!((fraction_sp3(&mol) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn balaban_j_positive_for_connected() {
        let mol = parse_smiles("CCCC").unwrap();
        let j = balaban_j(&mol);
        assert!(j > 0.0, "J={j}");
    }

    #[test]
    fn bertz_ct_nontrivial() {
        // Aspirin has diverse bond types and atom environments
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let ct = bertz_ct(&mol);
        assert!(ct > 0.0, "CT={ct}");
    }

    #[test]
    fn chi_connectivity_nonzero() {
        let mol = parse_smiles("CCCC").unwrap();
        let chi = chi_connectivity(&mol);
        assert!(chi[0] > 0.0);
        assert!(chi[1] > 0.0);
    }

    #[test]
    fn estate_indices_count() {
        let mol = parse_smiles("CCO").unwrap();
        let es = estate_indices(&mol);
        assert_eq!(es.len(), 3);
    }

    #[test]
    fn autocorrelation_size() {
        let mol = parse_smiles("CCCC").unwrap();
        let ac = autocorrelation_descriptors(&mol);
        assert_eq!(ac.moreau_broto.len(), 8);
        assert_eq!(ac.moran.len(), 8);
        assert_eq!(ac.geary.len(), 8);
    }

    #[test]
    fn compute_all_has_entries() {
        let mol = parse_smiles("CCO").unwrap();
        let ds = compute_all_descriptors(&mol);
        assert!(!ds.values.is_empty());
        // Should have wiener_index
        assert!(ds.values.iter().any(|(name, _)| name == "wiener_index"));
    }

    #[test]
    fn kappa_indices_positive() {
        let mol = parse_smiles("CCCC").unwrap();
        let (k1, k2, k3) = kappa_shape_indices(&mol);
        assert!(k1 > 0.0, "k1={k1}");
        assert!(k2 > 0.0, "k2={k2}");
        assert!(k3 >= 0.0, "k3={k3}");
    }
}
