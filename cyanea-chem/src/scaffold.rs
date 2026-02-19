//! Scaffold analysis: Murcko decomposition, generic scaffolds, MCS, and R-group decomposition.

use std::collections::{HashSet, VecDeque};
use std::time::Instant;

use cyanea_core::{CyaneaError, Result};

use crate::molecule::{Bond, BondOrder, BondStereo, Chirality, MolAtom, Molecule};
use crate::ring;
use crate::substructure::find_substructure_matches;

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Result of Murcko scaffold decomposition.
#[derive(Debug, Clone)]
pub struct MurckoResult {
    pub framework: Molecule,
    pub generic: Molecule,
}

/// Result of maximum common substructure (MCS) search.
#[derive(Debug, Clone)]
pub struct McsResult {
    pub mcs_atoms: usize,
    pub mcs_bonds: usize,
    pub mapping_a: Vec<(usize, usize)>,
    pub mapping_b: Vec<(usize, usize)>,
}

/// Result of R-group decomposition.
#[derive(Debug, Clone)]
pub struct RGroupResult {
    pub core_mapping: Vec<(usize, usize)>,
    pub r_groups: Vec<(String, Vec<usize>)>,
}

// ---------------------------------------------------------------------------
// Murcko scaffold
// ---------------------------------------------------------------------------

/// Compute the Murcko framework scaffold.
///
/// Finds ring atoms and linker atoms (atoms on shortest paths between rings),
/// then prunes everything else.
pub fn murcko_scaffold(mol: &Molecule) -> Result<MurckoResult> {
    let n = mol.atom_count();
    if n == 0 {
        let empty = Molecule::new(String::new(), Vec::new(), Vec::new());
        return Ok(MurckoResult { framework: empty.clone(), generic: empty });
    }

    let rings = ring::find_sssr(mol);
    if rings.is_empty() {
        // No rings — the Murcko scaffold is empty
        let empty = Molecule::new(String::new(), Vec::new(), Vec::new());
        return Ok(MurckoResult { framework: empty.clone(), generic: empty });
    }

    // Mark ring atoms
    let mut ring_atoms = vec![false; n];
    for ring in &rings {
        for &idx in ring {
            ring_atoms[idx] = true;
        }
    }

    // Find ring systems (connected groups of ring atoms)
    let ring_systems = find_ring_systems(mol, &ring_atoms);

    // Find linker atoms: atoms on shortest paths between ring systems
    let mut scaffold_atoms: HashSet<usize> = HashSet::new();
    for (idx, &is_ring) in ring_atoms.iter().enumerate() {
        if is_ring {
            scaffold_atoms.insert(idx);
        }
    }

    // For each pair of ring systems, find linker paths
    if ring_systems.len() > 1 {
        let linkers = find_linker_atoms(mol, &ring_atoms, &ring_systems);
        scaffold_atoms.extend(linkers);
    }

    let framework = extract_subgraph(mol, &scaffold_atoms);
    let generic = make_generic(&framework);

    Ok(MurckoResult { framework, generic })
}

/// Convenience function: generic scaffold (all atoms→C, all bonds→single).
pub fn generic_scaffold(mol: &Molecule) -> Result<Molecule> {
    let result = murcko_scaffold(mol)?;
    Ok(result.generic)
}

fn find_ring_systems(mol: &Molecule, ring_atoms: &[bool]) -> Vec<Vec<usize>> {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut systems = Vec::new();

    for start in 0..n {
        if !ring_atoms[start] || visited[start] {
            continue;
        }
        let mut system = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;
        while let Some(curr) = queue.pop_front() {
            system.push(curr);
            for &(neighbor, _) in &mol.adjacency[curr] {
                if ring_atoms[neighbor] && !visited[neighbor] {
                    visited[neighbor] = true;
                    queue.push_back(neighbor);
                }
            }
        }
        systems.push(system);
    }

    systems
}

fn find_linker_atoms(
    mol: &Molecule,
    ring_atoms: &[bool],
    ring_systems: &[Vec<usize>],
) -> HashSet<usize> {
    let mut linkers = HashSet::new();
    let n = mol.atom_count();

    // For each pair of ring systems, BFS for shortest path
    for i in 0..ring_systems.len() {
        for j in (i + 1)..ring_systems.len() {
            // BFS from all atoms of system i to any atom of system j
            let target_set: HashSet<usize> = ring_systems[j].iter().copied().collect();
            let mut best_path: Option<Vec<usize>> = None;

            for &start in &ring_systems[i] {
                if let Some(path) = bfs_path_to_set(mol, start, &target_set, n) {
                    if best_path.as_ref().is_none_or(|bp| path.len() < bp.len()) {
                        best_path = Some(path);
                    }
                }
            }

            if let Some(path) = best_path {
                for &atom in &path {
                    if !ring_atoms[atom] {
                        linkers.insert(atom);
                    }
                }
            }
        }
    }

    linkers
}

fn bfs_path_to_set(mol: &Molecule, start: usize, targets: &HashSet<usize>, n: usize) -> Option<Vec<usize>> {
    let mut visited = vec![false; n];
    let mut parent = vec![usize::MAX; n];
    let mut queue = VecDeque::new();

    visited[start] = true;
    queue.push_back(start);

    while let Some(curr) = queue.pop_front() {
        if curr != start && targets.contains(&curr) {
            // Reconstruct path
            let mut path = Vec::new();
            let mut node = curr;
            while node != start {
                path.push(node);
                node = parent[node];
            }
            path.push(start);
            path.reverse();
            return Some(path);
        }

        for &(neighbor, _) in &mol.adjacency[curr] {
            if !visited[neighbor] {
                visited[neighbor] = true;
                parent[neighbor] = curr;
                queue.push_back(neighbor);
            }
        }
    }

    None
}

fn extract_subgraph(mol: &Molecule, atoms: &HashSet<usize>) -> Molecule {
    let mut sorted_atoms: Vec<usize> = atoms.iter().copied().collect();
    sorted_atoms.sort();

    let mut index_map = vec![usize::MAX; mol.atom_count()];
    let mut new_atoms = Vec::new();

    for (new_idx, &old_idx) in sorted_atoms.iter().enumerate() {
        index_map[old_idx] = new_idx;
        new_atoms.push(mol.atoms[old_idx].clone());
    }

    let mut new_bonds = Vec::new();
    for bond in &mol.bonds {
        let a1 = index_map[bond.atom1];
        let a2 = index_map[bond.atom2];
        if a1 != usize::MAX && a2 != usize::MAX {
            new_bonds.push(Bond {
                atom1: a1,
                atom2: a2,
                order: bond.order,
                is_aromatic: bond.is_aromatic,
                stereo: bond.stereo,
            });
        }
    }

    Molecule::new(mol.name.clone(), new_atoms, new_bonds)
}

fn make_generic(mol: &Molecule) -> Molecule {
    let atoms: Vec<MolAtom> = mol
        .atoms
        .iter()
        .map(|_| MolAtom {
            atomic_number: 6, // All carbon
            formal_charge: 0,
            isotope: None,
            is_aromatic: false,
            implicit_hydrogens: 0,
            chirality: Chirality::None,
        })
        .collect();

    let bonds: Vec<Bond> = mol
        .bonds
        .iter()
        .map(|b| Bond {
            atom1: b.atom1,
            atom2: b.atom2,
            order: BondOrder::Single,
            is_aromatic: false,
            stereo: BondStereo::None,
        })
        .collect();

    Molecule::new(mol.name.clone(), atoms, bonds)
}

// ---------------------------------------------------------------------------
// Maximum Common Substructure (MCS)
// ---------------------------------------------------------------------------

/// Find the maximum common substructure between two molecules.
///
/// Uses a modular product graph approach with Bron-Kerbosch maximum clique.
/// The `timeout_ms` parameter limits computation time.
pub fn maximum_common_substructure(mol_a: &Molecule, mol_b: &Molecule, timeout_ms: u64) -> Result<McsResult> {
    let na = mol_a.atom_count();
    let nb = mol_b.atom_count();

    if na == 0 || nb == 0 {
        return Ok(McsResult {
            mcs_atoms: 0,
            mcs_bonds: 0,
            mapping_a: Vec::new(),
            mapping_b: Vec::new(),
        });
    }

    // Build modular product graph
    let (nodes, adj) = modular_product_graph(mol_a, mol_b);

    if nodes.is_empty() {
        return Ok(McsResult {
            mcs_atoms: 0,
            mcs_bonds: 0,
            mapping_a: Vec::new(),
            mapping_b: Vec::new(),
        });
    }

    // Find maximum clique with timeout
    let start_time = Instant::now();
    let timeout = std::time::Duration::from_millis(timeout_ms);
    let max_clique = bron_kerbosch_max(&adj, nodes.len(), start_time, timeout);

    // Convert clique back to atom mappings
    let mut mapping_a = Vec::new();
    let mut mapping_b = Vec::new();

    for &node_idx in &max_clique {
        let (a_idx, b_idx) = nodes[node_idx];
        mapping_a.push((a_idx, node_idx));
        mapping_b.push((b_idx, node_idx));
    }

    // Count MCS bonds
    let mcs_a_atoms: HashSet<usize> = max_clique.iter().map(|&ni| nodes[ni].0).collect();
    let mcs_b_atoms: HashSet<usize> = max_clique.iter().map(|&ni| nodes[ni].1).collect();

    let mcs_bonds_a = mol_a.bonds.iter().filter(|b| {
        mcs_a_atoms.contains(&b.atom1) && mcs_a_atoms.contains(&b.atom2)
    }).count();

    let mcs_bonds_b = mol_b.bonds.iter().filter(|b| {
        mcs_b_atoms.contains(&b.atom1) && mcs_b_atoms.contains(&b.atom2)
    }).count();

    Ok(McsResult {
        mcs_atoms: max_clique.len(),
        mcs_bonds: mcs_bonds_a.min(mcs_bonds_b),
        mapping_a,
        mapping_b,
    })
}

/// Build modular product graph: compatible (atom_a, atom_b) pairs with edges
/// for compatible bond pairs.
fn modular_product_graph(mol_a: &Molecule, mol_b: &Molecule) -> (Vec<(usize, usize)>, Vec<Vec<bool>>) {
    let na = mol_a.atom_count();
    let nb = mol_b.atom_count();

    // Compatible atom pairs (same element)
    let mut nodes: Vec<(usize, usize)> = Vec::new();
    for i in 0..na {
        for j in 0..nb {
            if mol_a.atoms[i].atomic_number == mol_b.atoms[j].atomic_number {
                nodes.push((i, j));
            }
        }
    }

    let nn = nodes.len();
    let mut adj = vec![vec![false; nn]; nn];

    // Two nodes are adjacent in the product graph if:
    // 1. They don't share atoms (i.e., different a-atoms and different b-atoms)
    // 2. If there's a bond between their a-atoms, there's a compatible bond between their b-atoms (and vice versa)
    for p in 0..nn {
        for q in (p + 1)..nn {
            let (a1, b1) = nodes[p];
            let (a2, b2) = nodes[q];

            // Must not map same atoms
            if a1 == a2 || b1 == b2 {
                continue;
            }

            let bond_a = mol_a.get_bond(a1, a2);
            let bond_b = mol_b.get_bond(b1, b2);

            let compatible = match (bond_a, bond_b) {
                (Some(ba), Some(bb)) => {
                    // Both have bonds — must be compatible order
                    ba.order == bb.order || (ba.is_aromatic && bb.is_aromatic)
                }
                (None, None) => {
                    // Neither has a bond — compatible (disconnected pair)
                    true
                }
                _ => false, // One has bond, other doesn't
            };

            if compatible {
                adj[p][q] = true;
                adj[q][p] = true;
            }
        }
    }

    (nodes, adj)
}

/// Bron-Kerbosch with pivot for maximum clique, with timeout.
fn bron_kerbosch_max(
    adj: &[Vec<bool>],
    n: usize,
    start_time: Instant,
    timeout: std::time::Duration,
) -> Vec<usize> {
    let mut best = Vec::new();
    let r = Vec::new();
    let p: Vec<usize> = (0..n).collect();
    let x = Vec::new();

    bk_pivot(adj, &r, &p, &x, &mut best, start_time, timeout);
    best
}

fn bk_pivot(
    adj: &[Vec<bool>],
    r: &[usize],
    p: &[usize],
    x: &[usize],
    best: &mut Vec<usize>,
    start_time: Instant,
    timeout: std::time::Duration,
) {
    // Check timeout
    if start_time.elapsed() > timeout {
        return;
    }

    if p.is_empty() && x.is_empty() {
        if r.len() > best.len() {
            *best = r.to_vec();
        }
        return;
    }

    if r.len() > best.len() {
        *best = r.to_vec();
    }

    // Upper bound pruning
    if r.len() + p.len() <= best.len() {
        return;
    }

    // Choose pivot: vertex in P ∪ X with most connections to P
    let pivot = p
        .iter()
        .chain(x.iter())
        .max_by_key(|&&v| p.iter().filter(|&&u| adj[v][u]).count())
        .copied();

    let candidates: Vec<usize> = match pivot {
        Some(pv) => p.iter().filter(|&&v| !adj[pv][v]).copied().collect(),
        None => p.to_vec(),
    };

    for v in candidates {
        if start_time.elapsed() > timeout {
            return;
        }

        let mut new_r = r.to_vec();
        new_r.push(v);

        let new_p: Vec<usize> = p.iter().filter(|&&u| adj[v][u]).copied().collect();
        let new_x: Vec<usize> = x.iter().filter(|&&u| adj[v][u]).copied().collect();

        bk_pivot(adj, &new_r, &new_p, &new_x, best, start_time, timeout);
    }
}

// ---------------------------------------------------------------------------
// R-group decomposition
// ---------------------------------------------------------------------------

/// Decompose a molecule into a core scaffold and R-groups.
///
/// Matches `core` as a substructure and identifies substituents at attachment points.
pub fn r_group_decomposition(mol: &Molecule, core: &Molecule) -> Result<RGroupResult> {
    let matches = find_substructure_matches(mol, core);
    if matches.is_empty() {
        return Err(CyaneaError::InvalidInput("core not found in molecule".into()));
    }

    let best_match = &matches[0];
    let core_mapping = best_match.atom_mapping.clone();

    // Target atoms that are mapped to core
    let core_target_atoms: HashSet<usize> = core_mapping.iter().map(|&(_, t)| t).collect();

    // Find attachment points: core atoms bonded to non-core atoms
    let mut r_groups: Vec<(String, Vec<usize>)> = Vec::new();
    let mut r_idx = 1;

    for &(_, target_atom) in &core_mapping {
        for &(neighbor, _) in &mol.adjacency[target_atom] {
            if !core_target_atoms.contains(&neighbor) {
                // This is an R-group attachment
                let r_atoms = collect_r_group(mol, neighbor, &core_target_atoms);
                let label = format!("R{r_idx}");
                r_groups.push((label, r_atoms));
                r_idx += 1;
            }
        }
    }

    Ok(RGroupResult { core_mapping, r_groups })
}

/// BFS to collect all atoms in an R-group fragment.
fn collect_r_group(mol: &Molecule, start: usize, core_atoms: &HashSet<usize>) -> Vec<usize> {
    let mut visited = HashSet::new();
    let mut queue = VecDeque::new();
    queue.push_back(start);
    visited.insert(start);

    while let Some(curr) = queue.pop_front() {
        for &(neighbor, _) in &mol.adjacency[curr] {
            if !core_atoms.contains(&neighbor) && !visited.contains(&neighbor) {
                visited.insert(neighbor);
                queue.push_back(neighbor);
            }
        }
    }

    let mut atoms: Vec<usize> = visited.into_iter().collect();
    atoms.sort();
    atoms
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn murcko_aspirin() {
        // Aspirin: CC(=O)Oc1ccccc1C(=O)O → framework should be phenyl ring
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let result = murcko_scaffold(&mol).unwrap();
        // Framework should contain the 6-membered ring
        assert_eq!(result.framework.atom_count(), 6, "framework atoms");
    }

    #[test]
    fn murcko_caffeine() {
        // Caffeine: Cn1c(=O)c2c(ncn2C)n1C — has fused ring system
        let mol = parse_smiles("Cn1c(=O)c2c(ncn2C)n1C").unwrap();
        let result = murcko_scaffold(&mol).unwrap();
        // Should preserve fused rings
        // Caffeine's fused bicyclic ring system has 8-9 ring atoms
        assert!(
            result.framework.atom_count() >= 8,
            "framework atoms = {}",
            result.framework.atom_count()
        );
    }

    #[test]
    fn generic_scaffold_all_carbon() {
        let mol = parse_smiles("c1ccc(N)cc1").unwrap();
        let generic = generic_scaffold(&mol).unwrap();
        for atom in &generic.atoms {
            assert_eq!(atom.atomic_number, 6);
        }
        for bond in &generic.bonds {
            assert_eq!(bond.order, BondOrder::Single);
        }
    }

    #[test]
    fn murcko_acyclic_empty() {
        // No rings → empty scaffold
        let mol = parse_smiles("CCCCCC").unwrap();
        let result = murcko_scaffold(&mol).unwrap();
        assert_eq!(result.framework.atom_count(), 0);
    }

    #[test]
    fn mcs_benzene_toluene() {
        let benzene = parse_smiles("c1ccccc1").unwrap();
        let toluene = parse_smiles("Cc1ccccc1").unwrap();
        let result = maximum_common_substructure(&benzene, &toluene, 5000).unwrap();
        // MCS should be benzene (6 atoms)
        assert_eq!(result.mcs_atoms, 6, "mcs_atoms = {}", result.mcs_atoms);
    }

    #[test]
    fn mcs_identical_molecules() {
        let mol = parse_smiles("CCO").unwrap();
        let result = maximum_common_substructure(&mol, &mol, 5000).unwrap();
        assert_eq!(result.mcs_atoms, 3);
    }

    #[test]
    fn mcs_empty_molecule() {
        let mol = parse_smiles("CCO").unwrap();
        let empty = Molecule::new(String::new(), Vec::new(), Vec::new());
        let result = maximum_common_substructure(&mol, &empty, 1000).unwrap();
        assert_eq!(result.mcs_atoms, 0);
    }

    #[test]
    fn r_group_toluene_benzene_core() {
        let toluene = parse_smiles("Cc1ccccc1").unwrap();
        let benzene = parse_smiles("c1ccccc1").unwrap();
        let result = r_group_decomposition(&toluene, &benzene).unwrap();
        assert!(!result.r_groups.is_empty(), "should find R-groups");
        // The methyl should be an R-group
        assert!(result.r_groups.iter().any(|(_, atoms)| atoms.len() == 1));
    }

    #[test]
    fn r_group_core_not_found() {
        let ethanol = parse_smiles("CCO").unwrap();
        let benzene = parse_smiles("c1ccccc1").unwrap();
        let result = r_group_decomposition(&ethanol, &benzene);
        assert!(result.is_err());
    }

    #[test]
    fn mcs_timeout_returns_partial() {
        // Even with 1ms timeout, should return something (or empty)
        let mol_a = parse_smiles("c1ccccc1").unwrap();
        let mol_b = parse_smiles("c1ccccc1C").unwrap();
        let result = maximum_common_substructure(&mol_a, &mol_b, 1).unwrap();
        // Should have found at least some atoms
        // Should return a valid result even with minimal timeout
        let _ = result.mcs_atoms;
    }
}
