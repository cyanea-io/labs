# Architecture -- cyanea-struct

Internal design documentation for the structural biology crate.

## Structure Hierarchy

The core data model follows a four-level hierarchy:

```
Structure
  +-- chains: Vec<Chain>
        +-- id: String (e.g., "A")
        +-- residues: Vec<Residue>
              +-- num: usize
              +-- code: String (e.g., "ALA")
              +-- atoms: Vec<Atom>
                    +-- id, name, x, y, z, occupancy, temp_factor, element
```

`Point3D` is a standalone coordinate type used by geometry functions. All types support `Clone` and `Debug`. The hierarchy is mutable so parsers can build it incrementally.

## PDB Parser (`pdb.rs`)

Line-based parser that processes ATOM and HETATM records from PDB format files:

- Reads fixed-width columns per the PDB specification (atom name at columns 12-16, coordinates at 30-54, etc.)
- Groups atoms into residues by residue sequence number and insertion code
- Groups residues into chains by chain identifier (column 21)
- Handles MODEL/ENDMDL records for NMR ensembles (only first model is returned)
- Alternate conformations (`alt_loc`) are stored but not filtered; the first conformer is primary
- `parse_pdb_file` is gated behind the `std` feature for filesystem access

## mmCIF Parser (`mmcif.rs`)

Tag-value pair parser for PDBx/mmCIF format:

- Parses `data_` blocks and `loop_` constructs
- `loop_` parsing: reads column headers (`_atom_site.*`), then whitespace-delimited data rows
- Maps mmCIF field names to Structure fields (`_atom_site.Cartn_x` to `x`, etc.)
- Handles quoted strings and multi-line values
- More robust than PDB for large structures and non-standard residues

## DSSP Secondary Structure (`secondary.rs`)

Two levels of secondary structure assignment:

**Simplified DSSP** (`assign_secondary_structure`):
- Computes backbone phi/psi dihedral angles from N, CA, C atoms
- Classifies residues based on angle ranges: alpha-helix (~-60, -45), beta-sheet (~-120, 130), coil (everything else)
- Fast but does not consider hydrogen bonding

**Full DSSP** (`dssp`):
- Implements the Kabsch-Sander algorithm (1983)
- Step 1: Identifies backbone N-H...O=C hydrogen bonds using the electrostatic energy criterion: E = 0.084 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) * 332 kcal/mol. A bond is accepted when E < -0.5 kcal/mol.
- Step 2: Pattern recognition on hydrogen bond networks:
  - Alpha-helix (H): i,i+4 H-bond pattern
  - 3-10 helix (G): i,i+3 H-bond pattern
  - Pi-helix (I): i,i+5 H-bond pattern
  - Beta-strand (E): parallel or antiparallel sheet H-bond ladders
  - Beta-bridge (B): isolated bridge
  - Turn (T): i,i+3/4/5 with H-bond but not helix
  - Bend (S): angle at CA > 70 degrees
  - Coil (C): default
- Returns `DsspAssignment` per residue with state, hydrogen bond partners, and solvent accessibility estimate

## Kabsch Superposition (`superposition.rs`)

SVD-based optimal rigid-body alignment:

1. Compute centroids of both point sets and translate to origin
2. Build the 3x3 cross-covariance matrix H = P1^T * P2
3. Compute SVD of H: H = U * S * V^T using the internal `linalg` module
4. Optimal rotation: R = V * U^T
5. Handle reflection: if det(R) < 0, negate the column of V corresponding to the smallest singular value and recompute R
6. RMSD computed from the aligned coordinates

The `linalg.rs` module provides a self-contained 3x3 SVD implementation (Jacobi iteration) to avoid external linear algebra dependencies.

## Contact Map (`contact.rs`)

Two modes of contact map computation:

**CA-only** (`compute_contact_map`):
- Extracts CA (alpha carbon) atoms from each residue
- Computes O(n^2) pairwise Euclidean distances
- Stores as a symmetric matrix

**All-atom** (`compute_contact_map_allatom`):
- For each residue pair, finds the minimum distance among all heavy atom pairs
- More accurate for detecting true contacts but O(n^2 * m^2) where m is atoms per residue
- Excludes hydrogen atoms

Both return a `ContactMap` with `get(i, j)` for distance lookup. A typical contact cutoff is 8.0 A for CA and 4.5 A for all-atom.

## Ramachandran Validation (`ramachandran.rs`)

Phi/psi dihedral angle classification:

- Computes phi (C_i-1, N_i, CA_i, C_i) and psi (N_i, CA_i, C_i, N_i+1) for each residue
- Classifies into regions based on empirical distributions:
  - **Favored**: core regions of the Ramachandran plot (>98% of validated residues)
  - **Allowed**: outer contour of the plot
  - **Outlier**: outside allowed regions, likely steric clash
- Special handling for glycine (no CB, broader allowed region), proline (constrained phi), and pre-proline residues
- `ramachandran_report` iterates all chains and returns per-residue classifications

## B-Factor Analysis (`analysis.rs`)

Temperature factor analysis for flexibility assessment:

- `residue_bfactors`: averages the `temp_factor` field across all atoms in each residue
- `chain_bfactors`: computes descriptive statistics (mean, std_dev, min, max, median) per chain
- `flexibility_score`: converts B-factors to Z-scores (number of standard deviations from the mean), identifying flexible (high Z) and rigid (low Z) regions
- High B-factor residues often correspond to loops, termini, or crystal packing interfaces

## Module Dependencies

```
types  <--  pdb, mmcif, geometry, secondary, superposition, contact, ramachandran, analysis
geometry  <--  secondary, superposition, ramachandran
linalg  <--  superposition
```

All modules depend on `cyanea-core` for `CyaneaError` and `Result`.
