# cyanea-struct

> Protein and nucleic acid 3D structure analysis: parsing, geometry, secondary structure, superposition, and validation.

## What's Inside

- **PDB parsing** -- ATOM/HETATM records into Structure/Chain/Residue/Atom hierarchy
- **mmCIF parsing** -- PDBx/mmCIF format with `loop_` construct handling
- **Geometry** -- distances, bond angles, dihedral angles, center of mass, RMSD
- **Secondary structure** -- simplified and full DSSP assignment from backbone dihedrals
- **Kabsch superposition** -- optimal rotation/translation alignment with RMSD
- **Contact maps** -- CA-only and all-atom (minimum distance) residue contact matrices
- **Ramachandran validation** -- phi/psi classification (favored/allowed/outlier/glycine/proline)
- **B-factor analysis** -- per-residue/per-chain statistics, flexibility Z-scores

## Quick Start

```toml
[dependencies]
cyanea-struct = "0.1"
```

```rust
use cyanea_struct::{parse_pdb, kabsch, assign_secondary_structure};

let structure = parse_pdb(pdb_text).unwrap();
let chain = &structure.chains[0];
let ss = assign_secondary_structure(chain).unwrap();

println!("{} residues, {} chains", chain.residues.len(), structure.chains.len());
```

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support (file I/O) |
| `wasm` | No | WASM target marker |
| `serde` | No | Serialize/Deserialize derives |
| `parallel` | No | Rayon parallelism |

## Modules

| Module | Description |
|--------|-------------|
| `types` | `Point3D`, `Atom`, `Residue`, `Chain`, `Structure` |
| `pdb` | PDB format parser |
| `mmcif` | mmCIF/PDBx format parser |
| `geometry` | Distance, angle, dihedral, center of mass, RMSD |
| `secondary` | DSSP-based secondary structure assignment |
| `superposition` | Kabsch optimal alignment |
| `contact` | Contact map computation (CA and all-atom) |
| `ramachandran` | Ramachandran region validation |
| `analysis` | B-factor statistics and flexibility scoring |

## See Also

- [API Reference](docs/API.md)
- [Usage Guide](docs/GUIDE.md)
- [Architecture](docs/ARCHITECTURE.md)
- [Build Guide](../docs/BUILDING.md)
