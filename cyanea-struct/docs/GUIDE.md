# Usage Guide -- cyanea-struct

Practical examples for protein and nucleic acid 3D structure analysis.

## PDB File Parsing and Structure Access

```rust
use cyanea_struct::parse_pdb;

let pdb_text = "ATOM      1  N   ALA A   1       1.000  2.000  3.000  1.00 10.00           N\n\
                 ATOM      2  CA  ALA A   1       2.000  3.000  4.000  1.00 12.00           C\n\
                 ATOM      3  C   ALA A   1       3.000  4.000  5.000  1.00 11.00           C\n\
                 END";

let structure = parse_pdb(pdb_text).unwrap();
println!("Structure: {}", structure.id);
println!("Chains: {}", structure.chains.len());

let chain = &structure.chains[0];
println!("Chain {} has {} residues", chain.id, chain.residues.len());

for residue in &chain.residues {
    println!("  Residue {} ({}) - {} atoms", residue.num, residue.code, residue.atoms.len());
}
```

## mmCIF Parsing

```rust
use cyanea_struct::parse_mmcif;

let mmcif_text = "data_example\n\
    loop_\n\
    _atom_site.group_PDB\n\
    _atom_site.id\n\
    _atom_site.type_symbol\n\
    _atom_site.label_atom_id\n\
    _atom_site.label_comp_id\n\
    _atom_site.label_asym_id\n\
    _atom_site.label_seq_id\n\
    _atom_site.Cartn_x\n\
    _atom_site.Cartn_y\n\
    _atom_site.Cartn_z\n\
    _atom_site.occupancy\n\
    _atom_site.B_iso_or_equiv\n\
    ATOM 1 N N ALA A 1 1.000 2.000 3.000 1.00 10.00\n\
    ATOM 2 C CA ALA A 1 2.000 3.000 4.000 1.00 12.00\n";

let structure = parse_mmcif(mmcif_text).unwrap();
println!("Parsed {} chains from mmCIF", structure.chains.len());
```

## Distance, Angle, and Dihedral Calculations

```rust
use cyanea_struct::geometry::{distance, angle, dihedral, center_of_mass};
use cyanea_struct::types::Atom;

// Create atoms with known coordinates
let atom_a = Atom { id: 1, name: "N".into(), residue_num: 1,
    alt_loc: None, x: 0.0, y: 0.0, z: 0.0,
    occupancy: 1.0, temp_factor: 10.0, element: "N".into() };
let atom_b = Atom { id: 2, name: "CA".into(), residue_num: 1,
    alt_loc: None, x: 1.5, y: 0.0, z: 0.0,
    occupancy: 1.0, temp_factor: 12.0, element: "C".into() };
let atom_c = Atom { id: 3, name: "C".into(), residue_num: 1,
    alt_loc: None, x: 1.5, y: 1.5, z: 0.0,
    occupancy: 1.0, temp_factor: 11.0, element: "C".into() };

let d = distance(&atom_a, &atom_b);
println!("N-CA distance: {:.2} A", d);

let a = angle(&atom_a, &atom_b, &atom_c);
println!("N-CA-C angle: {:.1} degrees", a);

let com = center_of_mass(&[atom_a.clone(), atom_b.clone(), atom_c.clone()]);
println!("Center of mass: ({:.2}, {:.2}, {:.2})", com.x, com.y, com.z);
```

## DSSP Secondary Structure Assignment

```rust
use cyanea_struct::{parse_pdb, secondary::assign_secondary_structure};

let structure = parse_pdb(pdb_text).unwrap();
let chain = &structure.chains[0];

let ss = assign_secondary_structure(chain).unwrap();
for (i, assignment) in ss.assignments.iter().enumerate() {
    println!("Residue {}: {:?}", i + 1, assignment);
}
```

For full DSSP with hydrogen bond information:

```rust
use cyanea_struct::{parse_pdb, secondary::dssp};

let structure = parse_pdb(pdb_text).unwrap();
let chain = &structure.chains[0];

let dssp_result = dssp(chain).unwrap();
for assignment in &dssp_result {
    println!("Residue {}: state={:?}", assignment.residue_num, assignment.state);
}
```

## Kabsch Superposition and RMSD

```rust
use cyanea_struct::superposition::{kabsch_points, kabsch};
use cyanea_struct::types::Point3D;

// Superpose two sets of corresponding points
let ref_points = vec![
    Point3D { x: 0.0, y: 0.0, z: 0.0 },
    Point3D { x: 1.0, y: 0.0, z: 0.0 },
    Point3D { x: 0.0, y: 1.0, z: 0.0 },
];
let mobile_points = vec![
    Point3D { x: 0.1, y: 0.1, z: 0.1 },
    Point3D { x: 1.1, y: 0.1, z: 0.1 },
    Point3D { x: 0.1, y: 1.1, z: 0.1 },
];

let result = kabsch_points(&ref_points, &mobile_points).unwrap();
println!("RMSD after superposition: {:.4} A", result.rmsd);
println!("Rotation matrix: {:?}", result.rotation_matrix);
```

For aligning two structures by their CA atoms:

```rust
use cyanea_struct::superposition::align_structures_by_ca;

let result = align_structures_by_ca(&structure1, &structure2).unwrap();
println!("CA RMSD: {:.2} A", result.rmsd);
```

## Contact Map Generation

```rust
use cyanea_struct::{parse_pdb, contact::{compute_contact_map, compute_contact_map_allatom}};

let structure = parse_pdb(pdb_text).unwrap();
let chain = &structure.chains[0];

// CA-only contact map
let cmap = compute_contact_map(chain).unwrap();
let n = chain.residues.len();
for i in 0..n {
    for j in (i + 1)..n {
        let d = cmap.get(i, j);
        if d < 8.0 {
            println!("Contact: residue {} - residue {} ({:.1} A)", i + 1, j + 1, d);
        }
    }
}

// All-atom contact map (minimum heavy atom distance)
let cmap_all = compute_contact_map_allatom(chain).unwrap();
```

## Ramachandran Analysis

```rust
use cyanea_struct::{parse_pdb, ramachandran::{ramachandran_report, validate_ramachandran}};

let structure = parse_pdb(pdb_text).unwrap();

let report = ramachandran_report(&structure).unwrap();
for (resnum, resname, phi, psi, region) in &report {
    println!("Residue {} ({}): phi={:.1}, psi={:.1}, region={:?} [{}]",
        resnum, resname, phi, psi, region, region.code());
}

// Validate a single phi/psi pair
let region = validate_ramachandran(-60.0, -45.0, "ALA");
println!("(-60, -45) for ALA: {:?}", region);
```

## B-Factor Analysis

```rust
use cyanea_struct::{parse_pdb, analysis::{chain_bfactors, flexibility_score, residue_bfactors}};

let structure = parse_pdb(pdb_text).unwrap();

// Per-chain B-factor statistics
let chain_stats = chain_bfactors(&structure).unwrap();
for (chain_id, stats) in &chain_stats {
    println!("Chain {}: mean B={:.1}, std={:.1}, range=[{:.1}, {:.1}]",
        chain_id, stats.mean, stats.std_dev, stats.min, stats.max);
}

// Per-residue B-factors
let res_bfactors = residue_bfactors(&structure).unwrap();
for (resnum, resname, bfactor) in &res_bfactors {
    println!("Residue {} ({}): B = {:.1}", resnum, resname, bfactor);
}

// Flexibility Z-scores (normalized B-factors)
let flex = flexibility_score(&structure).unwrap();
for (resnum, zscore) in &flex {
    if *zscore > 2.0 {
        println!("Residue {} is highly flexible (Z = {:.2})", resnum, zscore);
    }
}
```
