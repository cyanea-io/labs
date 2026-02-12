# cyanea-chem

Chemistry and small molecule toolkit: molecular graph representation, SMILES/SDF parsing, fingerprints, property calculation, substructure search, canonical SMILES, MACCS keys, and stereochemistry.

## Status: Complete

All planned functionality is implemented including SMILES parsing, SDF parsing, Morgan fingerprints, MACCS fingerprints, Tanimoto similarity, molecular property calculation, substructure matching, canonical SMILES generation, and stereochemistry assignment (R/S, E/Z).

## Public API

### Elements (`element.rs`)

| Type/Function | Description |
|---------------|-------------|
| `Element` | `symbol`, `number`, `mass`, `electronegativity` |
| `element_by_symbol(symbol) -> Option<&Element>` | Lookup by symbol (e.g., "C") |
| `element_by_number(n) -> Option<&Element>` | Lookup by atomic number |

### Molecular graph (`molecule.rs`)

| Type | Description |
|------|-------------|
| `BondOrder` | Enum: `Single`, `Double`, `Triple`, `Aromatic` |
| `MolAtom` | `atomic_number`, `formal_charge`, `aromatic`, `implicit_h_count` |
| `Bond` | `src`, `dst`, `bond_type` |
| `Molecule` | `atoms: Vec<MolAtom>`, `bonds: Vec<Bond>`, `name: Option<String>` |

### SMILES parsing (`smiles.rs`)

| Function | Description |
|----------|-------------|
| `parse_smiles(smiles: &str) -> Result<Molecule>` | Parse SMILES string |
| `parse_smiles_named(smiles, name) -> Result<Molecule>` | Parse with molecule name |

Supports: atoms, bonds, branches, ring closures, aromatic atoms, charges, explicit hydrogen.

### SDF parsing (`sdf.rs`)

| Function | Description |
|----------|-------------|
| `parse_mol_v2000(input) -> Result<Molecule>` | Parse single V2000 mol block |
| `parse_sdf(input) -> Vec<Result<Molecule>>` | Parse multi-molecule SDF |
| `parse_sdf_file(path) -> Result<Vec<Molecule>>` | Parse SDF file (std only) |

### Fingerprints (`fingerprint.rs`)

| Type/Function | Description |
|---------------|-------------|
| `Fingerprint` | Bit-vector molecular fingerprint |
| `morgan_fingerprint(mol, radius, nbits) -> Fingerprint` | Morgan/ECFP circular fingerprint |
| `tanimoto_similarity(fp1, fp2) -> f64` | Tanimoto coefficient [0, 1] |
| `tanimoto_bulk(query, targets) -> Vec<f64>` | Batch similarity search |

### Molecular properties (`properties.rs`)

| Type/Function | Description |
|---------------|-------------|
| `MolecularProperties` | `formula`, `weight`, `hbd`, `hba`, `rotatable_bonds`, `logp_estimate` |
| `compute_properties(mol) -> MolecularProperties` | Compute all properties |
| `molecular_weight(mol) -> f64` | Molecular weight |
| `molecular_formula(mol) -> String` | Hill system formula |
| `hbd_count(mol) -> usize` | Hydrogen bond donors |
| `hba_count(mol) -> usize` | Hydrogen bond acceptors |
| `rotatable_bond_count(mol, rings) -> usize` | Rotatable bonds |

### Substructure search (`substructure.rs`)

| Type/Function | Description |
|---------------|-------------|
| `SubstructureMatch` | `mapped_atoms: Vec<usize>` |
| `has_substructure(target, pattern) -> bool` | Check if pattern exists in target |
| `find_substructure_matches(target, pattern) -> Vec<SubstructureMatch>` | Find all matches |

### Canonical SMILES (`canon.rs`)

| Function | Description |
|----------|-------------|
| `canonical_smiles(mol) -> String` | Generate deterministic canonical SMILES string |

Algorithm: Morgan-like invariant refinement (atomic number, degree, H count, charge, isotope, aromaticity) followed by iterative neighbor summation until convergence, then DFS traversal with canonical atom ordering. Handles disconnected fragments (dot-separated), ring closures (including two-digit `%nn` notation), aromatic atoms, charges, and bracket atoms.

### MACCS fingerprints (`maccs.rs`)

| Function | Description |
|----------|-------------|
| `maccs_fingerprint(mol) -> Fingerprint` | Compute MACCS-like 166-key structural fingerprint |

Each of the 166 bit positions corresponds to a specific structural feature from the MACCS key definitions. Keys are evaluated directly from the molecular graph covering: element presence and count thresholds, ring topology (size, aromaticity, heteroatoms), bond types (single, double, triple, aromatic, specific atom pairs), functional groups (carboxyl, amide, sulfonyl, hydroxyl, amine, aldehyde, ester, nitro, thiol, phosphate, ether, tertiary amine, charged atoms), and count-based thresholds (heavy atoms, ring atoms, bond counts, branching degree).

### Stereochemistry (`stereo.rs`)

| Function | Description |
|----------|-------------|
| `assign_rs(mol, atom_idx) -> Option<char>` | Assign R/S descriptor to a tetrahedral stereocenter |
| `assign_ez(mol, bond_idx) -> Option<char>` | Assign E/Z descriptor to a double bond |

Implements Cahn-Ingold-Prelog (CIP) priority rules with two-level recursive tiebreaking by neighbor atomic numbers. `assign_rs` reads `Chirality` annotations (`@`/`@@` from SMILES) and returns `Some('R')` or `Some('S')` for valid stereocenters with four distinct-priority substituents, or `None` otherwise. `assign_ez` reads `BondStereo` annotations (`/`/`\` from SMILES) on double bonds and returns `Some('E')` (higher-priority groups on opposite sides) or `Some('Z')` (same side), or `None` if stereo markers are absent.

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support (file I/O for SDF) |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |

## Dependencies

- `cyanea-core` -- error types
- `sha2`, `hex` -- fingerprint hashing

## Tests

79 tests across 12 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 44 | Module declarations, re-exports |
| `element.rs` | 113 | Periodic table lookup |
| `molecule.rs` | 212 | Molecular graph representation |
| `smiles.rs` | 576 | SMILES parser |
| `sdf.rs` | 264 | SDF/Mol V2000 parser |
| `fingerprint.rs` | 224 | Morgan fingerprints, Tanimoto similarity |
| `properties.rs` | 255 | Molecular weight, formula, HBD/HBA, rotatable bonds |
| `substructure.rs` | 236 | Substructure matching |
| `ring.rs` | 250 | Ring detection (internal) |
| `canon.rs` | 591 | Canonical SMILES generation (Morgan-like ranking, DFS traversal) |
| `maccs.rs` | 1107 | MACCS 166-key structural fingerprints |
| `stereo.rs` | 456 | Stereochemistry assignment (R/S, E/Z via CIP rules) |
