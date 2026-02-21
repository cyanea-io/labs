# cyanea-chem

> Chemistry and small molecule toolkit: molecular graphs, SMILES/SDF/SMARTS parsing, fingerprints, properties, descriptors, drug-likeness, scaffolds, 3D conformers, force fields, and chemical reactions.

200 unit tests + 7 doc tests.

## What's Inside

- **Molecular graph** -- `Molecule` with atoms, bonds, ring detection (SSSR)
- **SMILES parsing** -- atoms, bonds, branches, ring closures, aromaticity, charges, explicit hydrogen
- **SDF parsing** -- V2000 and V3000 mol blocks and multi-molecule SDF files
- **SMARTS matching** -- atom/bond query primitives, logical operators (AND/OR/NOT), recursive SMARTS
- **Canonical SMILES** -- deterministic canonical SMILES generation via Morgan-like invariant refinement
- **Morgan fingerprints** -- ECFP-style circular fingerprints with configurable radius and bit width
- **MACCS fingerprints** -- 166-key structural fingerprints covering elements, rings, bonds, functional groups
- **Tanimoto similarity** -- coefficient between fingerprints, including bulk search
- **Molecular properties** -- formula, weight, HBD/HBA, rotatable bonds, LogP estimate
- **Molecular descriptors** -- Wiener, Balaban J, Zagreb, TPSA, Wildman-Crippen LogP, BertzCT, kappa shape, chi connectivity, fsp3, ring counts, E-state, autocorrelation
- **Drug-likeness** -- Lipinski, Veber, PAINS alerts, Brenk alerts, lead-likeness, QED
- **Scaffold analysis** -- Murcko framework, generic scaffold, MCS, R-group decomposition
- **Standardization** -- salt stripping, charge neutralization, tautomer canonicalization, largest fragment
- **Substructure search** -- pattern matching with atom mapping
- **Stereochemistry** -- R/S (tetrahedral) and E/Z (double bond) assignment via CIP priority rules
- **3D conformers** -- distance geometry embedding, ETKDG torsion preferences, multi-conformer enumeration
- **Force fields** -- UFF and MMFF94 energy, steepest descent and conjugate gradient minimization
- **Gasteiger charges** -- iterative partial equalization of electronegativity
- **Chemical reactions** -- SMIRKS parsing, reaction application, atom-atom mapping, retrosynthetic disconnection
- **Element data** -- periodic table lookup by symbol or atomic number

## Quick Start

```toml
[dependencies]
cyanea-chem = "0.1"
```

```rust
use cyanea_chem::{parse_smiles, compute_properties, canonical_smiles};

let mol = parse_smiles("c1ccccc1O").unwrap(); // phenol
let props = compute_properties(&mol);
println!("Formula: {}, MW: {:.1}", props.formula, props.weight);

let canon = canonical_smiles(&mol);
println!("Canonical: {}", canon);
```

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support (SDF file I/O) |
| `wasm` | No | WASM target marker |
| `serde` | No | Serialize/Deserialize derives |
| `parallel` | No | Rayon parallelism |

## Modules

| Module | Description |
|--------|-------------|
| `element` | Periodic table lookup (`Element`, `element_by_symbol`) |
| `molecule` | `Molecule`, `MolAtom`, `Bond`, `BondOrder` |
| `smiles` | SMILES parser |
| `sdf` | SDF/Mol V2000/V3000 parser |
| `smarts` | SMARTS pattern parser and substructure matcher |
| `canon` | Canonical SMILES generation |
| `fingerprint` | Morgan/ECFP fingerprints, Tanimoto similarity |
| `maccs` | MACCS 166-key structural fingerprints |
| `properties` | Molecular weight, formula, HBD/HBA, rotatable bonds |
| `descriptors` | Wiener, Zagreb, TPSA, LogP, BertzCT, kappa, chi, fsp3, E-state, autocorrelation |
| `druglikeness` | Lipinski, Veber, PAINS, Brenk, lead-likeness, QED |
| `scaffold` | Murcko scaffold, MCS, R-group decomposition |
| `standardize` | Salt stripping, neutralization, tautomer canonicalization |
| `substructure` | Substructure matching |
| `stereo` | R/S and E/Z stereochemistry (CIP rules) |
| `conformer` | 3D conformer container |
| `embed` | Distance geometry + ETKDG 3D embedding |
| `forcefield` | UFF/MMFF94 energy, steepest descent/conjugate gradient minimization |
| `gasteiger` | Gasteiger-Marsili partial charge calculation |
| `reaction` | SMIRKS reactions, enumeration, atom-atom mapping, retrosynthesis |

## See Also

- [API Reference](docs/API.md)
- [Usage Guide](docs/GUIDE.md)
- [Architecture](docs/ARCHITECTURE.md)
