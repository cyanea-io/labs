# cyanea-chem Architecture

Internal design of the cheminformatics modules.

## Molecule Representation

The `Molecule` struct stores a molecular graph as an atom list, bond list, and adjacency list:

- `atoms: Vec<MolAtom>` -- each atom has `atomic_number`, `formal_charge`, `aromatic`, `implicit_h_count`, `chirality`, `isotope`
- `bonds: Vec<Bond>` -- each bond has `src`, `dst`, `bond_type` (Single/Double/Triple/Aromatic), `stereo` (for E/Z)
- `adjacency: Vec<Vec<(usize, usize)>>` -- per-atom list of `(neighbor_index, bond_index)` for O(1) neighbor lookup

`Molecule` derives `Eq` and `Hash` (for use in sets/maps), which precludes storing `f64` fields like 3D coordinates. Coordinates are stored separately in `Conformer`.

## Ring Perception (SSSR)

Ring detection uses the Hanser algorithm for finding the smallest set of smallest rings (SSSR):

1. Build a path graph from the molecular graph
2. Iteratively remove nodes of degree 2 (merge their two edges)
3. For nodes of degree > 2, enumerate and store all rings formed by merging paths
4. The result is a set of `Vec<usize>` atom index vectors, one per ring

Ring membership is used throughout the codebase: aromaticity perception, rotatable bond counting, scaffold decomposition, MACCS key evaluation, and force field atom typing.

## SMILES Parser

The parser follows the Weininger/OpenSMILES specification:

1. **Tokenization**: Single-pass character-by-character parsing. Bracket atoms `[...]` are parsed with all properties (element, charge, H count, chirality, isotope). Organic subset atoms (B, C, N, O, P, S, F, Cl, Br, I) use implicit valence rules.
2. **Branch handling**: `(` pushes the current atom onto a stack, `)` pops to restore context.
3. **Ring closures**: Digits (and `%nn` two-digit codes) record pending ring-open positions. When a matching closure is found, a bond is created between the two atoms.
4. **Aromaticity**: Lowercase atoms (c, n, o, s) are marked aromatic. After parsing, aromatic bonds are set between adjacent aromatic atoms in rings.
5. **Implicit hydrogens**: Computed from valence rules for organic subset atoms.

## Canonical SMILES

The canonicalization algorithm:

1. **Initial invariants**: Assign each atom a tuple of (atomic_number, degree, H_count, charge, isotope, is_aromatic).
2. **Morgan-like refinement**: Iteratively update each atom's invariant by summing neighbor invariants. Continue until the number of unique invariants stabilizes (convergence).
3. **Canonical ordering**: Sort atoms by their final invariants to break ties deterministically.
4. **DFS traversal**: Starting from the atom with the smallest invariant, perform a DFS traversal of the molecular graph, emitting SMILES tokens. At each branch point, visit children in canonical order. Ring closures are assigned in order of first encounter.

Handles disconnected fragments (dot-separated), bracket atoms, and aromatic notation.

## Fingerprints

### Morgan (ECFP)

1. Assign initial identifiers to each atom (hash of atomic number, degree, H count, charge, ring membership).
2. For each iteration (up to `radius`):
   - For each atom, collect the sorted list of neighbor identifiers
   - Hash the (atom_id, sorted_neighbor_ids) tuple to produce a new identifier
   - Record the identifier and the set of atoms in its environment
3. Fold all collected identifiers into a fixed-length bit vector of size `nbits` using modular hashing.

ECFP4 = radius 2, ECFP6 = radius 3. Hash function uses SHA-256 (from `sha2` crate).

### MACCS 166 Keys

Each of the 166 bit positions is evaluated by a specific structural test:
- Element presence/count thresholds (e.g., "has >= 2 nitrogen atoms")
- Ring tests (e.g., "has 5-membered aromatic ring")
- Bond patterns (e.g., "N-C=O")
- Functional group tests (e.g., "has carboxyl group")

Tests are evaluated directly from the molecular graph and ring perception results.

## CIP Priority Rules (Stereochemistry)

R/S assignment for tetrahedral stereocenters:

1. Identify stereocenters: sp3 atoms with 4 distinct substituents (from SMILES `@`/`@@` annotations).
2. Assign CIP priorities to the 4 neighbors using recursive comparison of atomic numbers along the substituent branches.
3. Orient so the lowest-priority group points away (using the SMILES convention where `@` = counterclockwise and `@@` = clockwise with respect to the first three neighbors).
4. Determine if the remaining three groups go clockwise (R) or counterclockwise (S).

E/Z assignment for double bonds: compare CIP priorities of the two substituents on each end of the double bond using the `/` and `\` bond stereo markers from SMILES.

## SMARTS Pattern Matching

SMARTS extends SMILES with query primitives:

- **Atom expressions**: Logical combinations (AND, OR, NOT) of atom primitives (atomic number, aromatic, degree, H count, charge, ring membership, ring size, connectivity, valence, wildcard).
- **Bond expressions**: Similarly for bonds (single, double, triple, aromatic, ring, any).
- **Recursive SMARTS**: `$(...)` allows embedding a full SMARTS pattern as an atom-level test.

Matching uses a VF2-like subgraph isomorphism algorithm, testing atom and bond expressions at each mapping step.

## Distance Geometry (3D Embedding)

Conformer generation follows the distance geometry approach:

1. **Bounds matrix**: Compute lower and upper bounds for all atom-pair distances from bond lengths, bond angles, and van der Waals radii (using covalent radii table).
2. **Triangle inequality smoothing**: Tighten bounds by iterating: for each triple (i,j,k), enforce l_ij >= l_ik - u_kj and u_ij <= u_ik + u_kj.
3. **Random distance sampling**: For each pair, sample a distance uniformly within [lower, upper].
4. **Metric matrix**: Convert the distance matrix to a metric matrix via double-centering.
5. **SVD/eigendecomposition**: Extract the top 3 eigenvectors, scale by sqrt(eigenvalues) to get 3D coordinates.
6. **ETKDG torsion preferences**: Optionally apply preferred torsion angles from a small knowledge base.
7. **Minimization**: Optionally refine with UFF or MMFF94 force field.

Multiple conformers are generated with different random seeds and deduplicated by RMSD threshold.

## Force Fields

Both UFF and MMFF94 share the same energy functional:

E_total = E_bond + E_angle + E_torsion + E_vdW + E_electrostatic + E_oop

1. **Atom typing**: Assign force field atom types based on element, hybridization, and ring context.
2. **Parameter lookup**: Retrieve equilibrium bond lengths, angles, force constants, and vdW parameters from built-in tables.
3. **Energy evaluation**: Sum contributions from all bonds, angles, proper torsions, improper torsions (out-of-plane), van der Waals (Lennard-Jones 12-6), and electrostatic (Coulombic with Gasteiger charges).
4. **Gradient**: Analytical derivatives of each energy term with respect to atomic coordinates.
5. **Minimization**: Steepest descent (simple, robust) or conjugate gradient (Fletcher-Reeves, faster convergence). Line search by backtracking with Armijo condition.

## Gasteiger Charges

Iterative partial equalization of orbital electronegativity (Gasteiger & Marsili 1980):

1. Initialize all partial charges to formal charges.
2. For each atom, compute electronegativity chi = a + b*q + c*q^2, where (a, b, c) are element/hybridization-dependent parameters.
3. For each bond, transfer charge proportional to the electronegativity difference, damped by a factor that decreases each iteration (damping = 0.5^iteration).
4. Iterate for 6-8 cycles (convergence is typically fast).

Used by the electrostatic term in force field calculations.

## Reactions (SMIRKS)

SMIRKS reactions are parsed as `reactant_SMARTS>>product_SMARTS` with atom-map numbers (`[C:1]`):

1. **Parse**: Split on `>>`, parse each side as a SMARTS pattern, extract atom-map numbers.
2. **Match**: Find all SMARTS matches of the reactant pattern in the input molecule.
3. **Transform**: For each match, apply the product template:
   - Create new atoms/bonds as specified by the product SMARTS
   - Delete atoms/bonds present in the reactant but not in the product
   - Modify bond orders and atom properties according to the mapping
4. **Retrosynthesis**: Apply a library of known disconnection transforms in the reverse direction (product >> reactants).
