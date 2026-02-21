# API Reference -- cyanea-struct

Protein and nucleic acid 3D structure analysis: PDB parsing, mmCIF parsing, geometric calculations, secondary structure assignment, structural superposition, contact maps, Ramachandran validation, and B-factor analysis.

## Public API

### Structure types (`types.rs`)

| Type | Description |
|------|-------------|
| `Point3D` | 3D coordinate: `x`, `y`, `z` |
| `Atom` | `id`, `name`, `residue_num`, `alt_loc`, `x`, `y`, `z`, `occupancy`, `temp_factor`, `element` |
| `Residue` | `num`, `code`, `atoms: Vec<Atom>` |
| `Chain` | `id`, `residues: Vec<Residue>` |
| `Structure` | `id`, `chains: Vec<Chain>` |

### PDB parsing (`pdb.rs`)

| Function | Description |
|----------|-------------|
| `parse_pdb(input: &str) -> Result<Structure>` | Parse PDB format string |
| `parse_pdb_file(path) -> Result<Structure>` | Parse PDB file (std only) |

### mmCIF parsing (`mmcif.rs`)

| Function | Description |
|----------|-------------|
| `parse_mmcif(input: &str) -> Result<Structure>` | Parse mmCIF/PDBx format string into a `Structure` |
| `parse_mmcif_loop(lines: &[&str]) -> Vec<BTreeMap<String, String>>` | Parse a `loop_` construct into field-name to value records |

### Geometry (`geometry.rs`)

| Function | Description |
|----------|-------------|
| `distance(a1, a2) -> f64` | Euclidean distance between atoms |
| `distance_points(p1, p2) -> f64` | Distance between points |
| `angle(a1, a2, a3) -> f64` | Bond angle in degrees |
| `angle_points(p1, p2, p3) -> f64` | Angle between points |
| `dihedral(a1, a2, a3, a4) -> f64` | Dihedral angle in degrees |
| `dihedral_points(p1, p2, p3, p4) -> f64` | Dihedral between points |
| `center_of_mass(atoms) -> Point3D` | Center of mass |
| `rmsd(atoms1, atoms2) -> Result<f64>` | Root Mean Square Deviation |

### Secondary structure (`secondary.rs`)

| Type/Function | Description |
|---------------|-------------|
| `SecondaryStructure` | Enum: `Helix`, `Sheet`, `Coil` |
| `SecondaryStructureAssignment` | Per-residue SS assignments |
| `DsspAssignment` | Full DSSP per-residue result (SS state, hydrogen bonds, accessibility) |
| `DsspState` | Enum: `H` (alpha-helix), `B` (beta-bridge), `E` (strand), `G` (3-10 helix), `I` (pi-helix), `T` (turn), `S` (bend), `C` (coil) |
| `assign_secondary_structure(chain) -> Result<SecondaryStructureAssignment>` | Simplified DSSP based on phi/psi angles |
| `dssp(chain) -> Result<Vec<DsspAssignment>>` | Full DSSP secondary structure assignment |
| `backbone_dihedrals(chain) -> Result<Vec<(f64, f64)>>` | Phi/psi Ramachandran angles |

### Superposition (`superposition.rs`)

| Type/Function | Description |
|---------------|-------------|
| `SuperpositionResult` | `rotation_matrix`, `translation_vector`, `rmsd` |
| `kabsch(atoms1, atoms2) -> Result<SuperpositionResult>` | Kabsch optimal alignment |
| `kabsch_points(points1, points2) -> Result<SuperpositionResult>` | Kabsch on raw points |
| `align_structures_by_ca(s1, s2) -> Result<SuperpositionResult>` | Align by CA atoms |

### Contact maps (`contact.rs`)

| Type/Function | Description |
|---------------|-------------|
| `ContactMap` | Residue-residue distance matrix |
| `ContactMap::get(i, j) -> f64` | Distance between residues i and j |
| `compute_contact_map(chain) -> Result<ContactMap>` | CA-only contact map |
| `compute_contact_map_allatom(chain) -> Result<ContactMap>` | All heavy atoms (minimum distance) |

### Ramachandran validation (`ramachandran.rs`)

| Type/Function | Description |
|---------------|-------------|
| `RamachandranRegion` | Enum: `Favored`, `Allowed`, `Outlier`, `Glycine`, `Proline`, `PreProline` |
| `RamachandranRegion::code(&self) -> char` | Single-character code (`F`, `A`, `O`, `G`, `P`, `p`) |
| `validate_ramachandran(phi, psi, residue_type) -> RamachandranRegion` | Classify phi/psi into Ramachandran region |
| `ramachandran_report(structure) -> Result<Vec<(usize, String, f64, f64, RamachandranRegion)>>` | Per-residue Ramachandran report for all chains |

### B-factor analysis (`analysis.rs`)

| Type/Function | Description |
|---------------|-------------|
| `BFactorStats` | Descriptive statistics: `mean`, `std_dev`, `min`, `max`, `median` |
| `residue_bfactors(structure) -> Result<Vec<(usize, String, f64)>>` | Per-residue average B-factors |
| `chain_bfactors(structure) -> Result<Vec<(String, BFactorStats)>>` | Per-chain B-factor statistics |
| `flexibility_score(structure) -> Result<Vec<(usize, f64)>>` | Normalized B-factor Z-scores per residue |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support (file I/O) |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |
| `parallel` | No | Rayon parallelism |

Note: Uses `#![no_std]` with `alloc` -- core algorithms work without std.

## Dependencies

- `cyanea-core` -- error types
- `sha2`, `hex` -- structure hashing

## Tests

76 unit tests + 2 doc tests across 11 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 106 | Module declarations, re-exports |
| `types.rs` | 427 | Point3D, Atom, Residue, Chain, Structure |
| `pdb.rs` | 384 | PDB format parser |
| `mmcif.rs` | 630 | mmCIF/PDBx format parser |
| `geometry.rs` | 204 | Distance, angle, dihedral, RMSD |
| `secondary.rs` | 1123 | Secondary structure assignment (simplified and full DSSP) |
| `superposition.rs` | 247 | Kabsch algorithm |
| `contact.rs` | 324 | Contact map computation |
| `ramachandran.rs` | 402 | Ramachandran plot validation |
| `analysis.rs` | 405 | B-factor analysis and flexibility scoring |
| `linalg.rs` | 295 | Internal 3x3 matrix operations, SVD |
