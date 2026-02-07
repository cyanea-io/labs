# cyanea-chem

Chemistry and small molecule toolkit: molecular graph representation, SMILES/SDF parsing, fingerprints, property calculation, and substructure search.

## Status: Complete

All planned functionality is implemented including SMILES parsing, SDF parsing, Morgan fingerprints, Tanimoto similarity, molecular property calculation, and substructure matching.

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

34 tests across 8 source files.

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
