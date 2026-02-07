//! Ring detection via smallest set of smallest rings (SSSR).

use std::collections::VecDeque;

use crate::molecule::Molecule;

/// Find the smallest set of smallest rings (SSSR) in a molecule.
///
/// Returns a vector of rings, where each ring is a vector of atom indices.
pub fn find_sssr(mol: &Molecule) -> Vec<Vec<usize>> {
    let n = mol.atom_count();
    if n == 0 || mol.bond_count() == 0 {
        return Vec::new();
    }

    // Expected number of rings = bonds - atoms + connected_components
    let num_components = count_components(mol);
    let expected_rings = mol.bond_count() as isize - mol.atom_count() as isize + num_components as isize;
    if expected_rings <= 0 {
        return Vec::new();
    }

    // Find ring atoms by iteratively pruning terminal atoms
    let ring_atoms = find_ring_atoms(mol);
    if ring_atoms.is_empty() {
        return Vec::new();
    }

    // For each ring bond, find the shortest ring containing it via BFS
    let mut rings: Vec<Vec<usize>> = Vec::new();

    for bond_idx in 0..mol.bonds.len() {
        let bond = &mol.bonds[bond_idx];
        if !ring_atoms[bond.atom1] || !ring_atoms[bond.atom2] {
            continue;
        }

        // BFS from atom1 to atom2 without using this bond
        if let Some(path) = bfs_shortest_path(mol, bond.atom1, bond.atom2, bond_idx, &ring_atoms) {
            // The ring is: path from atom1 to atom2 (which includes both endpoints)
            let mut ring = path;

            // Normalize: start with the smallest index, and choose the direction
            // that gives the lexicographically smallest sequence
            normalize_ring(&mut ring);

            // Check if this ring is already found
            if !rings.iter().any(|r| r == &ring) {
                rings.push(ring);
            }
        }
    }

    // Sort rings by size for deterministic output
    rings.sort_by_key(|r| r.len());

    // Trim to expected SSSR size
    if rings.len() > expected_rings as usize {
        rings.truncate(expected_rings as usize);
    }

    rings
}

/// Check if a specific atom is in any ring.
#[allow(dead_code)]
pub fn is_in_ring(mol: &Molecule, atom_idx: usize) -> bool {
    let ring_atoms = find_ring_atoms(mol);
    atom_idx < ring_atoms.len() && ring_atoms[atom_idx]
}

/// Check if a specific bond is in any ring.
#[allow(dead_code)]
pub fn bond_in_ring(mol: &Molecule, bond_idx: usize) -> bool {
    if bond_idx >= mol.bonds.len() {
        return false;
    }
    let ring_atoms = find_ring_atoms(mol);
    let bond = &mol.bonds[bond_idx];
    ring_atoms[bond.atom1] && ring_atoms[bond.atom2]
}

/// Count connected components.
fn count_components(mol: &Molecule) -> usize {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut components = 0;

    for start in 0..n {
        if visited[start] {
            continue;
        }
        components += 1;
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

    components
}

/// Identify ring atoms by iteratively removing terminal (degree-1) atoms.
fn find_ring_atoms(mol: &Molecule) -> Vec<bool> {
    let n = mol.atom_count();
    let mut degree = vec![0usize; n];
    for i in 0..n {
        degree[i] = mol.adjacency[i].len();
    }

    let mut queue: VecDeque<usize> = VecDeque::new();
    for i in 0..n {
        if degree[i] <= 1 {
            queue.push_back(i);
        }
    }

    let mut removed = vec![false; n];
    while let Some(atom) = queue.pop_front() {
        if removed[atom] {
            continue;
        }
        removed[atom] = true;
        for &(neighbor, _) in &mol.adjacency[atom] {
            if !removed[neighbor] {
                degree[neighbor] -= 1;
                if degree[neighbor] <= 1 {
                    queue.push_back(neighbor);
                }
            }
        }
    }

    removed.iter().map(|&r| !r).collect()
}

/// BFS from `start` to `end` avoiding a specific bond, restricted to ring atoms.
fn bfs_shortest_path(
    mol: &Molecule,
    start: usize,
    end: usize,
    excluded_bond: usize,
    ring_atoms: &[bool],
) -> Option<Vec<usize>> {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut parent = vec![usize::MAX; n];
    let mut queue = VecDeque::new();

    visited[start] = true;
    queue.push_back(start);

    while let Some(curr) = queue.pop_front() {
        if curr == end {
            // Reconstruct path
            let mut path = Vec::new();
            let mut node = end;
            while node != start {
                path.push(node);
                node = parent[node];
            }
            path.push(start);
            path.reverse();
            return Some(path);
        }

        for &(neighbor, bond_idx) in &mol.adjacency[curr] {
            if bond_idx == excluded_bond {
                continue;
            }
            if !visited[neighbor] && ring_atoms[neighbor] {
                visited[neighbor] = true;
                parent[neighbor] = curr;
                queue.push_back(neighbor);
            }
        }
    }

    None
}

/// Normalize a ring so it starts with the smallest index and proceeds
/// in the direction that gives the lexicographically smaller sequence.
fn normalize_ring(ring: &mut Vec<usize>) {
    if ring.is_empty() {
        return;
    }

    // Find the position of the minimum element
    let min_pos = ring
        .iter()
        .enumerate()
        .min_by_key(|&(_, &v)| v)
        .unwrap()
        .0;

    // Rotate so minimum is first
    ring.rotate_left(min_pos);

    // Compare forward vs reverse direction
    let n = ring.len();
    if n > 2 {
        let forward_second = ring[1];
        let reverse_second = ring[n - 1];
        if reverse_second < forward_second {
            // Reverse (keeping first element in place)
            ring[1..].reverse();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn benzene_one_ring() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let rings = find_sssr(&mol);
        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].len(), 6);
    }

    #[test]
    fn naphthalene_two_rings() {
        // Naphthalene: c1ccc2ccccc2c1
        let mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        let rings = find_sssr(&mol);
        assert_eq!(rings.len(), 2);
        for ring in &rings {
            assert_eq!(ring.len(), 6);
        }
    }

    #[test]
    fn acyclic_no_rings() {
        let mol = parse_smiles("CCCC").unwrap();
        let rings = find_sssr(&mol);
        assert_eq!(rings.len(), 0);
        assert!(!is_in_ring(&mol, 0));
    }
}
