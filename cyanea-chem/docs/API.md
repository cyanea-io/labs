# cyanea-chem

Chemistry and small molecule toolkit: molecular graph representation, SMILES/SDF/SMARTS parsing, fingerprints, property calculation, substructure search, canonical SMILES, MACCS keys, stereochemistry, molecular descriptors, drug-likeness filters, scaffold analysis, molecule standardization, 3D conformer generation, force field energy calculations, Gasteiger charges, chemical reactions, and 3D embedding.

## Status: Complete

All planned functionality is implemented across 20 modules covering molecular representation, I/O (SMILES, SDF V2000/V3000, SMARTS), fingerprints (Morgan, MACCS), properties, substructure search, canonical SMILES, stereochemistry (CIP rules), molecular descriptors (topological, physicochemical), drug-likeness (Lipinski, Veber, PAINS, Brenk, QED), scaffold analysis (Murcko, MCS, R-group), standardization pipeline, 3D conformer generation (distance geometry, ETKDG), force fields (UFF, MMFF94), Gasteiger charges, chemical reactions (SMIRKS, retrosynthesis), and 3D coordinate embedding.

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
| `Molecule` | `atoms: Vec<MolAtom>`, `bonds: Vec<Bond>`, `name: Option<String>`, `adjacency` |

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
| `parse_mol_v3000(input) -> Result<Molecule>` | Parse single V3000 mol block |
| `parse_sdf(input) -> Vec<Result<Molecule>>` | Parse multi-molecule SDF |
| `parse_sdf_file(path) -> Result<Vec<Molecule>>` | Parse SDF file (std only) |

### SMARTS pattern matching (`smarts.rs`)

| Type/Function | Description |
|---------------|-------------|
| `AtomPrimitive` | Enum: `AtomicNum`, `Aromatic`, `Aliphatic`, `Degree`, `HCount`, `TotalHCount`, `Charge`, `RingMember`, `RingSize`, `Connectivity`, `Valence`, `Wildcard` |
| `AtomExpr` | Enum: `Prim`, `And(Vec)`, `Or(Vec)`, `Not(Box)`, `Recursive(SmartsPattern)` |
| `BondExpr` | Enum: `Single`, `Double`, `Triple`, `Aromatic`, `Ring`, `Any`, `Not(Box)`, `And(Vec)`, `Or(Vec)` |
| `SmartsPattern` | Parsed SMARTS pattern with atoms, bonds, and topology |
| `parse_smarts(smarts) -> Result<SmartsPattern>` | Parse SMARTS string |
| `smarts_match(mol, pattern) -> bool` | Check if molecule matches SMARTS pattern |
| `smarts_find_all(mol, pattern) -> Vec<SubstructureMatch>` | Find all SMARTS matches with atom mapping |

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

Algorithm: Morgan-like invariant refinement (atomic number, degree, H count, charge, isotope, aromaticity) followed by iterative neighbor summation until convergence, then DFS traversal with canonical atom ordering. Handles disconnected fragments, ring closures, aromatic atoms, charges, and bracket atoms.

### MACCS fingerprints (`maccs.rs`)

| Function | Description |
|----------|-------------|
| `maccs_fingerprint(mol) -> Fingerprint` | Compute MACCS-like 166-key structural fingerprint |

Each of the 166 bit positions corresponds to a specific structural feature from the MACCS key definitions covering element presence, ring topology, bond types, functional groups, and count-based thresholds.

### Stereochemistry (`stereo.rs`)

| Function | Description |
|----------|-------------|
| `assign_rs(mol, atom_idx) -> Option<char>` | Assign R/S descriptor to a tetrahedral stereocenter |
| `assign_ez(mol, bond_idx) -> Option<char>` | Assign E/Z descriptor to a double bond |

Implements Cahn-Ingold-Prelog (CIP) priority rules with recursive tiebreaking.

### Molecular descriptors (`descriptors.rs`)

Graph-based molecular descriptors for QSAR/QSPR.

| Type/Function | Description |
|---------------|-------------|
| `RingDetails` | Ring count breakdown: `total`, `aromatic`, `aliphatic`, `heteroaromatic`, `spiro`, `fused` |
| `AutocorrelationResult` | `moreau_broto`, `moran`, `geary` autocorrelation vectors |
| `DescriptorSet` | Named collection of `(String, f64)` descriptor values |
| `wiener_index(mol) -> usize` | Wiener topological index (sum of shortest paths) |
| `balaban_j(mol) -> f64` | Balaban J index |
| `zagreb_indices(mol) -> (f64, f64)` | First and second Zagreb indices |
| `tpsa(mol) -> f64` | Topological polar surface area |
| `wildman_crippen_logp(mol) -> f64` | Wildman-Crippen LogP estimate |
| `bertz_ct(mol) -> f64` | BertzCT complexity index |
| `kappa_shape(mol) -> (f64, f64, f64)` | Kappa shape indices (1, 2, 3) |
| `chi_connectivity(mol) -> Vec<f64>` | Chi connectivity indices |
| `fsp3(mol) -> f64` | Fraction of sp3-hybridized carbons |
| `ring_counts(mol) -> RingDetails` | Detailed ring count information |
| `estate_indices(mol) -> Vec<f64>` | E-state indices per atom |
| `autocorrelation(mol, property, max_lag) -> AutocorrelationResult` | Moreau-Broto, Moran, Geary autocorrelation descriptors |
| `compute_all_descriptors(mol) -> DescriptorSet` | Batch compute all descriptors |

### Drug-likeness (`druglikeness.rs`)

| Type/Function | Description |
|---------------|-------------|
| `LipinskiResult` | `mw`, `logp`, `hbd`, `hba`, `passes`, `violations` |
| `VeberResult` | `rotatable_bonds`, `tpsa`, `passes` |
| `AlertResult` | Single structural alert match (`name`, `matched`) |
| `AlertFilterResult` | Filter result: `alerts`, `n_hits`, `passes` |
| `QedResult` | `score`, `properties` |
| `DrugLikenessReport` | Comprehensive report (Lipinski, Veber, PAINS, Brenk, lead-likeness, QED) |
| `lipinski(mol) -> LipinskiResult` | Lipinski's Rule of Five |
| `veber(mol) -> VeberResult` | Veber rules (rotatable bonds, TPSA) |
| `pains_filter(mol) -> AlertFilterResult` | PAINS structural alerts |
| `brenk_filter(mol) -> AlertFilterResult` | Brenk structural alerts |
| `lead_likeness(mol) -> bool` | Lead-likeness criteria |
| `qed(mol) -> QedResult` | Quantitative Estimate of Drug-likeness |
| `drug_likeness_report(mol) -> DrugLikenessReport` | Full drug-likeness assessment |

### Scaffold analysis (`scaffold.rs`)

| Type/Function | Description |
|---------------|-------------|
| `MurckoResult` | `framework`, `generic` (Molecule instances) |
| `McsResult` | `mcs_atoms`, `mcs_bonds`, `mapping_a`, `mapping_b` |
| `RGroupResult` | `core_mapping`, `r_groups` |
| `murcko_scaffold(mol) -> Result<MurckoResult>` | Murcko framework scaffold decomposition |
| `generic_scaffold(mol) -> Result<Molecule>` | Generic scaffold (all atoms to carbon, all bonds to single) |
| `mcs(mol_a, mol_b) -> Result<McsResult>` | Maximum Common Substructure search |
| `r_group_decomposition(mol, core) -> Result<RGroupResult>` | R-group decomposition against a core |

### Molecule standardization (`standardize.rs`)

| Type/Function | Description |
|---------------|-------------|
| `StandardizeStep` | Enum: `LargestFragment`, `StripSalts`, `Neutralize`, `CanonicalTautomer` |
| `StandardizeConfig` | Pipeline configuration with ordered steps |
| `standardize(mol, config) -> Result<Molecule>` | Run standardization pipeline |
| `largest_fragment(mol) -> Molecule` | Keep largest connected component by MW |
| `strip_salts(mol) -> Molecule` | Remove small ionic fragments |
| `neutralize(mol) -> Molecule` | Neutralize charges where possible |
| `canonical_tautomer(mol) -> Molecule` | Canonicalize tautomeric form |

### 3D conformer container (`conformer.rs`)

| Type/Function | Description |
|---------------|-------------|
| `Conformer` | 3D coordinate container: `coords: Vec<[f64; 3]>` per atom |
| `ConformerSet` | Collection of conformers with RMSD-based deduplication |
| `Conformer::new(coords) -> Self` | Create conformer from coordinates |
| `Conformer::len() / is_empty()` | Number of atoms |
| `Conformer::distance(i, j) -> f64` | Euclidean distance between atoms |
| `Conformer::angle(i, j, k) -> f64` | Bond angle in radians at atom j |
| `Conformer::dihedral(i, j, k, l) -> f64` | Torsion angle in radians |

### 3D embedding (`embed.rs`)

| Type/Function | Description |
|---------------|-------------|
| `ForceFieldType` | Enum: `None`, `Uff`, `Mmff94` |
| `EmbedConfig` | `max_conformers`, `rmsd_threshold`, `use_torsion_prefs`, `random_seed`, `force_field`, `max_minimize_steps` |
| `embed_molecule(mol, config) -> Result<ConformerSet>` | Generate 3D conformers via distance geometry + ETKDG |

Algorithm: compute distance bounds from connectivity, sample random distances within bounds, convert to metric matrix via eigendecomposition, extract 3D coordinates, optionally apply ETKDG torsion preferences, then minimize with force field.

### Force fields (`forcefield.rs`)

| Type/Function | Description |
|---------------|-------------|
| `EnergyComponents` | Per-term breakdown: `bond_stretch`, `angle_bend`, `torsion`, `van_der_waals`, `electrostatic`, `out_of_plane`, `total` |
| `MinimizeResult` | `conformer`, `initial_energy`, `final_energy`, `n_steps`, `converged`, `energy_components` |
| `MinimizeMethod` | Enum: `SteepestDescent`, `ConjugateGradient` |
| `MinimizeConfig` | `max_steps`, `gradient_threshold`, `method` |
| `uff_energy(mol, conformer) -> Result<EnergyComponents>` | UFF energy calculation |
| `mmff94_energy(mol, conformer) -> Result<EnergyComponents>` | MMFF94 energy calculation |
| `minimize(mol, conformer, config) -> Result<MinimizeResult>` | Energy minimization (steepest descent or conjugate gradient) |

### Gasteiger charges (`gasteiger.rs`)

| Function | Description |
|----------|-------------|
| `gasteiger_charges(mol) -> Result<Vec<f64>>` | Gasteiger-Marsili partial charge calculation |

Iterative partial equalization of orbital electronegativity. Parameters for H, C, N, O, F, Cl, Br, S, P with hybridization-dependent coefficients.

### Chemical reactions (`reaction.rs`)

| Type/Function | Description |
|---------------|-------------|
| `Reaction` | Parsed SMIRKS reaction: `reactant_pattern`, `product_template`, atom maps |
| `ReactionProduct` | Product `molecule` + canonical `smiles` |
| `AtomAtomMapping` | `mapping` pairs + `unmapped_reactant` / `unmapped_product` indices |
| `Disconnection` | Retrosynthetic result: `transform_name`, `smirks`, `precursors` |
| `parse_reaction(smirks) -> Result<Reaction>` | Parse SMIRKS reaction string |
| `apply_reaction(mol, reaction) -> Result<Vec<ReactionProduct>>` | Apply reaction to molecule, return products |
| `enumerate_reactions(mols, reaction) -> Result<Vec<Vec<ReactionProduct>>>` | Enumerate reaction across multiple substrates |
| `atom_atom_mapping(reactant, product, reaction) -> Result<AtomAtomMapping>` | Compute atom-atom mapping |
| `retrosynthetic_disconnection(mol) -> Result<Vec<Disconnection>>` | Apply retrosynthetic transforms to find precursors |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support (file I/O for SDF) |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |
| `parallel` | No | Rayon parallelism |

## Dependencies

- `cyanea-core` -- error types
- `sha2`, `hex` -- fingerprint hashing

## Tests

200 unit tests + 7 doc tests across 22 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | -- | Module declarations, re-exports |
| `element.rs` | -- | Periodic table lookup |
| `molecule.rs` | -- | Molecular graph representation |
| `smiles.rs` | -- | SMILES parser |
| `sdf.rs` | -- | SDF/Mol V2000/V3000 parser |
| `smarts.rs` | -- | SMARTS pattern parser and matcher |
| `fingerprint.rs` | -- | Morgan fingerprints, Tanimoto similarity |
| `maccs.rs` | -- | MACCS 166-key structural fingerprints |
| `properties.rs` | -- | Molecular weight, formula, HBD/HBA, rotatable bonds |
| `substructure.rs` | -- | Substructure matching |
| `ring.rs` | -- | Ring detection SSSR (internal) |
| `canon.rs` | -- | Canonical SMILES generation |
| `stereo.rs` | -- | Stereochemistry (R/S, E/Z via CIP rules) |
| `descriptors.rs` | -- | Topological and physicochemical descriptors |
| `druglikeness.rs` | -- | Lipinski, Veber, PAINS, Brenk, QED |
| `scaffold.rs` | -- | Murcko scaffold, MCS, R-group decomposition |
| `standardize.rs` | -- | Salt stripping, neutralization, tautomer canonicalization |
| `conformer.rs` | -- | 3D conformer container |
| `embed.rs` | -- | Distance geometry + ETKDG 3D embedding |
| `forcefield.rs` | -- | UFF/MMFF94 energy, steepest descent/conjugate gradient minimization |
| `gasteiger.rs` | -- | Gasteiger-Marsili partial charges |
| `reaction.rs` | -- | SMIRKS reactions, enumeration, retrosynthesis |
