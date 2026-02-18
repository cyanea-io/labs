# cyanea-chem

> Chemistry and small molecule toolkit: molecular graphs, SMILES/SDF parsing, fingerprints, properties, and substructure search.

## What's Inside

- **Molecular graph** -- `Molecule` with atoms, bonds, ring detection
- **SMILES parsing** -- atoms, bonds, branches, ring closures, aromaticity, charges, explicit hydrogen
- **SDF parsing** -- V2000 mol block and multi-molecule SDF files
- **Canonical SMILES** -- deterministic canonical SMILES generation via Morgan-like invariant refinement
- **Morgan fingerprints** -- ECFP-style circular fingerprints with configurable radius and bit width
- **MACCS fingerprints** -- 166-key structural fingerprints covering elements, rings, bonds, functional groups
- **Tanimoto similarity** -- coefficient between fingerprints, including bulk search
- **Molecular properties** -- formula, weight, HBD/HBA, rotatable bonds, LogP estimate
- **Substructure search** -- pattern matching with atom mapping
- **Stereochemistry** -- R/S (tetrahedral) and E/Z (double bond) assignment via CIP priority rules
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
| `sdf` | SDF/Mol V2000 parser |
| `canon` | Canonical SMILES generation |
| `fingerprint` | Morgan/ECFP fingerprints, Tanimoto similarity |
| `maccs` | MACCS 166-key structural fingerprints |
| `properties` | Molecular weight, formula, HBD/HBA, rotatable bonds |
| `substructure` | Substructure matching |
| `stereo` | R/S and E/Z stereochemistry (CIP rules) |

## See Also

- [API Reference (STATUS.md)](docs/STATUS.md)
- [Architecture](../ARCHITECTURE.md)
- [Build Guide](../docs/BUILDING.md)
