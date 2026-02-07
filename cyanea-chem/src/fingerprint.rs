//! Molecular fingerprints and similarity.

use cyanea_core::ContentAddressable;
use sha2::{Digest, Sha256};

use crate::molecule::Molecule;
use crate::ring;

/// A fixed-size bit vector fingerprint.
#[derive(Debug, Clone)]
pub struct Fingerprint {
    bits: Vec<u64>,
    nbits: usize,
}

impl Fingerprint {
    /// Create a new fingerprint of the given size (rounded up to multiple of 64).
    pub fn new(nbits: usize) -> Self {
        let nwords = (nbits + 63) / 64;
        Fingerprint {
            bits: vec![0u64; nwords],
            nbits,
        }
    }

    /// Set a bit at the given position.
    pub fn set_bit(&mut self, pos: usize) {
        let pos = pos % self.nbits;
        let word = pos / 64;
        let bit = pos % 64;
        self.bits[word] |= 1u64 << bit;
    }

    /// Get a bit at the given position.
    pub fn get_bit(&self, pos: usize) -> bool {
        let pos = pos % self.nbits;
        let word = pos / 64;
        let bit = pos % 64;
        (self.bits[word] >> bit) & 1 == 1
    }

    /// Count the number of set bits.
    pub fn count_ones(&self) -> u32 {
        self.bits.iter().map(|w| w.count_ones()).sum()
    }

    /// Number of bits in the fingerprint.
    pub fn nbits(&self) -> usize {
        self.nbits
    }
}

impl ContentAddressable for Fingerprint {
    fn content_hash(&self) -> String {
        let mut hasher = Sha256::new();
        for word in &self.bits {
            hasher.update(word.to_le_bytes());
        }
        hex::encode(hasher.finalize())
    }
}

/// Compute a Morgan (ECFP-like) fingerprint.
///
/// `radius` controls the neighborhood size (2 = ECFP4, 3 = ECFP6).
/// `nbits` is the fingerprint length (commonly 2048).
pub fn morgan_fingerprint(mol: &Molecule, radius: usize, nbits: usize) -> Fingerprint {
    let n = mol.atom_count();
    let mut fp = Fingerprint::new(nbits);

    if n == 0 {
        return fp;
    }

    let ring_atoms = {
        let rings = ring::find_sssr(mol);
        let mut ra = vec![false; n];
        for ring in &rings {
            for &idx in ring {
                ra[idx] = true;
            }
        }
        ra
    };

    // Initial invariants: hash of atom properties
    let mut identifiers: Vec<u64> = Vec::with_capacity(n);
    for (i, atom) in mol.atoms.iter().enumerate() {
        let mut h = fnv1a_init();
        h = fnv1a_update(h, atom.atomic_number as u64);
        h = fnv1a_update(h, mol.degree(i) as u64);
        h = fnv1a_update(h, atom.implicit_hydrogens as u64);
        h = fnv1a_update(h, atom.formal_charge as u64);
        h = fnv1a_update(h, ring_atoms[i] as u64);
        h = fnv1a_update(h, atom.is_aromatic as u64);
        identifiers.push(h);
    }

    // Set bits for radius 0
    for &id in &identifiers {
        fp.set_bit(fold_hash(id, nbits));
    }

    // Iterate for each radius
    for _ in 0..radius {
        let mut new_identifiers = Vec::with_capacity(n);
        for i in 0..n {
            let mut h = fnv1a_init();
            h = fnv1a_update(h, identifiers[i]);

            // Sort neighbor identifiers for determinism
            let mut neighbor_ids: Vec<(u64, u8)> = mol.adjacency[i]
                .iter()
                .map(|&(neighbor, bond_idx)| {
                    (identifiers[neighbor], mol.bonds[bond_idx].order as u8)
                })
                .collect();
            neighbor_ids.sort();

            for (nid, border) in &neighbor_ids {
                h = fnv1a_update(h, *nid);
                h = fnv1a_update(h, *border as u64);
            }

            new_identifiers.push(h);
            fp.set_bit(fold_hash(h, nbits));
        }
        identifiers = new_identifiers;
    }

    fp
}

/// Tanimoto similarity coefficient between two fingerprints.
///
/// Returns 1.0 for identical fingerprints, 0.0 for completely disjoint.
pub fn tanimoto_similarity(fp1: &Fingerprint, fp2: &Fingerprint) -> f64 {
    assert_eq!(fp1.nbits, fp2.nbits, "fingerprints must have the same size");

    let mut and_count = 0u32;
    let mut or_count = 0u32;

    for (w1, w2) in fp1.bits.iter().zip(fp2.bits.iter()) {
        and_count += (w1 & w2).count_ones();
        or_count += (w1 | w2).count_ones();
    }

    if or_count == 0 {
        return 1.0; // Both empty → identical
    }

    and_count as f64 / or_count as f64
}

/// Compute Tanimoto similarity of a query against multiple targets.
pub fn tanimoto_bulk(query: &Fingerprint, targets: &[Fingerprint]) -> Vec<f64> {
    targets.iter().map(|t| tanimoto_similarity(query, t)).collect()
}

// FNV-1a hash functions for deterministic hashing
const FNV_OFFSET: u64 = 0xcbf29ce484222325;
const FNV_PRIME: u64 = 0x100000001b3;

fn fnv1a_init() -> u64 {
    FNV_OFFSET
}

fn fnv1a_update(hash: u64, value: u64) -> u64 {
    let bytes = value.to_le_bytes();
    let mut h = hash;
    for &b in &bytes {
        h ^= b as u64;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h
}

fn fold_hash(hash: u64, nbits: usize) -> usize {
    (hash as usize) % nbits
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn bit_operations() {
        let mut fp = Fingerprint::new(128);
        assert!(!fp.get_bit(42));
        fp.set_bit(42);
        assert!(fp.get_bit(42));
        assert_eq!(fp.count_ones(), 1);
        fp.set_bit(100);
        assert_eq!(fp.count_ones(), 2);
    }

    #[test]
    fn deterministic_fingerprint() {
        let mol = parse_smiles("CCO").unwrap();
        let fp1 = morgan_fingerprint(&mol, 2, 2048);
        let fp2 = morgan_fingerprint(&mol, 2, 2048);
        assert_eq!(fp1.bits, fp2.bits);
        assert_eq!(fp1.content_hash(), fp2.content_hash());
    }

    #[test]
    fn tanimoto_identical_is_one() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let fp = morgan_fingerprint(&mol, 2, 2048);
        let sim = tanimoto_similarity(&fp, &fp);
        assert!((sim - 1.0).abs() < 1e-10);
    }

    #[test]
    fn tanimoto_different_in_range() {
        let mol1 = parse_smiles("CCO").unwrap(); // ethanol
        let mol2 = parse_smiles("CCCO").unwrap(); // propanol
        let fp1 = morgan_fingerprint(&mol1, 2, 2048);
        let fp2 = morgan_fingerprint(&mol2, 2, 2048);
        let sim = tanimoto_similarity(&fp1, &fp2);
        assert!(sim > 0.0 && sim < 1.0, "tanimoto = {sim}");
    }
}
