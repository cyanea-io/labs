# cyanea-struct

Protein and nucleic acid 3D structure analysis: PDB parsing, geometric calculations, secondary structure assignment, structural superposition, and contact maps.

## Status: Complete

All planned functionality is implemented including PDB parsing, distance/angle/dihedral calculation, RMSD, simplified DSSP, Kabsch superposition, and contact map computation.

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
| `assign_secondary_structure(chain) -> Result<SecondaryStructureAssignment>` | Simplified DSSP based on phi/psi angles |
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

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support (file I/O) |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |

Note: Uses `#![no_std]` with `alloc` -- core algorithms work without std.

## Dependencies

- `cyanea-core` -- error types
- `sha2`, `hex` -- structure hashing

## Tests

36 tests across 8 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 96 | Module declarations, re-exports |
| `types.rs` | 426 | Point3D, Atom, Residue, Chain, Structure |
| `pdb.rs` | 383 | PDB format parser |
| `geometry.rs` | 203 | Distance, angle, dihedral, RMSD |
| `secondary.rs` | 455 | Secondary structure assignment (simplified DSSP) |
| `superposition.rs` | 246 | Kabsch algorithm |
| `contact.rs` | 283 | Contact map computation |
| `linalg.rs` | 294 | Internal 3x3 matrix operations, SVD |
