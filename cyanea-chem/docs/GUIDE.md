# cyanea-chem Usage Guide

Practical examples for cheminformatics: SMILES/SDF parsing, fingerprints, properties, substructure search, descriptors, drug-likeness, scaffold analysis, 3D conformers, and chemical reactions.

## SMILES Parsing and Canonical Output

```rust
use cyanea_chem::{parse_smiles, canonical_smiles};

// Parse a SMILES string into a molecular graph
let mol = parse_smiles("c1ccccc1O").unwrap(); // phenol
println!("Atoms: {}, Bonds: {}", mol.atom_count(), mol.bond_count());

// Generate canonical SMILES (deterministic, unique)
let canon = canonical_smiles(&mol);
println!("Canonical: {}", canon);

// Different input SMILES for the same molecule produce the same canonical form
let mol2 = parse_smiles("Oc1ccccc1").unwrap();
assert_eq!(canonical_smiles(&mol2), canon);

// Named molecules
let aspirin = cyanea_chem::parse_smiles_named("CC(=O)Oc1ccccc1C(=O)O", "aspirin").unwrap();
```

## SDF File Reading (V2000 and V3000)

```rust
use cyanea_chem::sdf::*;

// Parse a single V2000 mol block
let mol_block = r#"
  compound1
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
"#;
let mol = parse_mol_v2000(mol_block).unwrap();

// Parse multi-molecule SDF
let molecules = parse_sdf(sdf_contents);
for result in molecules {
    match result {
        Ok(mol) => println!("Parsed: {} atoms", mol.atom_count()),
        Err(e) => eprintln!("Parse error: {}", e),
    }
}

// Read from file (requires std feature)
let mols = parse_sdf_file("compounds.sdf").unwrap();
```

## Molecular Property Calculation

```rust
use cyanea_chem::{parse_smiles, compute_properties, molecular_weight, hbd_count, hba_count};

let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap(); // aspirin
let props = compute_properties(&mol);

println!("Formula: {}", props.formula);
println!("MW: {:.2}", props.weight);
println!("HBD: {}", props.hbd);
println!("HBA: {}", props.hba);
println!("Rotatable bonds: {}", props.rotatable_bonds);
println!("LogP estimate: {:.2}", props.logp_estimate);
```

## Morgan/MACCS Fingerprints and Similarity

```rust
use cyanea_chem::{parse_smiles, morgan_fingerprint, tanimoto_similarity, tanimoto_bulk};
use cyanea_chem::maccs::maccs_fingerprint;

let aspirin = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
let ibuprofen = parse_smiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O").unwrap();
let caffeine = parse_smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C").unwrap();

// Morgan (ECFP4-like) fingerprints
let fp1 = morgan_fingerprint(&aspirin, 2, 2048);
let fp2 = morgan_fingerprint(&ibuprofen, 2, 2048);
let fp3 = morgan_fingerprint(&caffeine, 2, 2048);

let sim = tanimoto_similarity(&fp1, &fp2);
println!("Aspirin-Ibuprofen similarity: {:.3}", sim);

// Bulk similarity search
let targets = vec![fp2.clone(), fp3.clone()];
let sims = tanimoto_bulk(&fp1, &targets);
println!("Similarities: {:?}", sims);

// MACCS 166-key fingerprints
let maccs1 = maccs_fingerprint(&aspirin);
let maccs2 = maccs_fingerprint(&ibuprofen);
let maccs_sim = tanimoto_similarity(&maccs1, &maccs2);
```

## Substructure Search

```rust
use cyanea_chem::{parse_smiles, has_substructure, find_substructure_matches};

let molecule = parse_smiles("c1ccc(cc1)C(=O)O").unwrap(); // benzoic acid
let pattern = parse_smiles("C(=O)O").unwrap(); // carboxyl group

if has_substructure(&molecule, &pattern) {
    println!("Contains carboxyl group");
}

let matches = find_substructure_matches(&molecule, &pattern);
for m in &matches {
    println!("Match at atoms: {:?}", m.mapped_atoms);
}
```

## SMARTS Pattern Search

```rust
use cyanea_chem::smarts::{parse_smarts, smarts_match, smarts_find_all};
use cyanea_chem::parse_smiles;

let mol = parse_smiles("c1ccc(cc1)NC(=O)C").unwrap(); // acetanilide

// SMARTS for amide bond
let pattern = parse_smarts("[NX3][CX3](=O)").unwrap();
let has_amide = smarts_match(&mol, &pattern);
println!("Has amide: {}", has_amide);

// Find all matches
let matches = smarts_find_all(&mol, &pattern);
println!("Found {} amide(s)", matches.len());
```

## Molecular Descriptors

```rust
use cyanea_chem::{parse_smiles, descriptors::*};

let mol = parse_smiles("c1ccccc1O").unwrap(); // phenol

let wi = wiener_index(&mol);
println!("Wiener index: {}", wi);

let (z1, z2) = zagreb_indices(&mol);
println!("Zagreb M1={:.0}, M2={:.0}", z1, z2);

let polar_area = tpsa(&mol);
println!("TPSA: {:.1} A^2", polar_area);

let logp = wildman_crippen_logp(&mol);
println!("LogP: {:.2}", logp);

let sp3 = fsp3(&mol);
println!("Fsp3: {:.2}", sp3);

// All descriptors at once
let all = compute_all_descriptors(&mol);
for (name, value) in &all.values {
    println!("{}: {:.4}", name, value);
}
```

## Drug-Likeness Filters

```rust
use cyanea_chem::{parse_smiles, druglikeness::*};

let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap(); // aspirin

// Lipinski's Rule of Five
let lip = lipinski(&mol);
println!("Passes Ro5: {} (violations: {})", lip.passes, lip.violations);

// Veber rules
let veb = veber(&mol);
println!("Passes Veber: {} (TPSA={:.1}, RotBonds={})",
    veb.passes, veb.tpsa, veb.rotatable_bonds);

// QED score
let q = qed(&mol);
println!("QED: {:.3}", q.score);

// PAINS structural alerts
let pains = pains_filter(&mol);
println!("PAINS clean: {} ({} alerts)", pains.passes, pains.n_hits);

// Full report
let report = drug_likeness_report(&mol);
```

## Scaffold Analysis

```rust
use cyanea_chem::{parse_smiles, scaffold::*};

let mol = parse_smiles("c1ccc2c(c1)cc(cc2)NC(=O)C").unwrap();

// Murcko scaffold
let murcko = murcko_scaffold(&mol).unwrap();
println!("Framework atoms: {}", murcko.framework.atom_count());
println!("Generic scaffold atoms: {}", murcko.generic.atom_count());

// Maximum common substructure
let mol_a = parse_smiles("c1ccccc1CC").unwrap(); // ethylbenzene
let mol_b = parse_smiles("c1ccccc1CO").unwrap(); // benzyl alcohol
let mcs_result = mcs(&mol_a, &mol_b).unwrap();
println!("MCS: {} atoms, {} bonds", mcs_result.mcs_atoms, mcs_result.mcs_bonds);
```

## Molecule Standardization

```rust
use cyanea_chem::{parse_smiles, standardize::*};

// Default pipeline: strip salts -> largest fragment -> neutralize -> canonical tautomer
let mol = parse_smiles("[Na+].OC(=O)c1ccccc1").unwrap(); // sodium benzoate
let clean = standardize(&mol, StandardizeConfig::default()).unwrap();
// Result: benzoic acid (salt stripped, neutralized)

// Custom pipeline
let config = StandardizeConfig {
    steps: vec![
        StandardizeStep::LargestFragment,
        StandardizeStep::Neutralize,
    ],
};
let clean = standardize(&mol, config).unwrap();
```

## 3D Conformer Generation and Optimization

```rust
use cyanea_chem::{parse_smiles, embed::*};

let mol = parse_smiles("CCCC").unwrap(); // butane

// Generate 3D conformer with ETKDG + UFF minimization
let config = EmbedConfig {
    max_conformers: 10,
    rmsd_threshold: 0.5,
    use_torsion_prefs: true,
    force_field: ForceFieldType::Uff,
    max_minimize_steps: 200,
    random_seed: 42,
};
let conformers = embed_molecule(&mol, config).unwrap();
println!("Generated {} unique conformers", conformers.len());

// Access coordinates
for (i, conf) in conformers.iter().enumerate() {
    let d = conf.distance(0, 3); // distance between atoms 0 and 3
    println!("Conformer {}: end-to-end distance = {:.2} A", i, d);
}
```

## Force Field Energy

```rust
use cyanea_chem::{parse_smiles, embed::*, forcefield::*};

let mol = parse_smiles("CCCC").unwrap();
let conformers = embed_molecule(&mol, EmbedConfig::default()).unwrap();
let conf = &conformers[0];

// Calculate UFF energy
let energy = uff_energy(&mol, conf).unwrap();
println!("Total energy: {:.2} kcal/mol", energy.total);
println!("  Bond stretch: {:.2}", energy.bond_stretch);
println!("  Angle bend:   {:.2}", energy.angle_bend);
println!("  Torsion:      {:.2}", energy.torsion);
println!("  vdW:          {:.2}", energy.van_der_waals);

// Minimize with conjugate gradient
let result = minimize(&mol, conf, MinimizeConfig {
    max_steps: 500,
    gradient_threshold: 0.01,
    method: MinimizeMethod::ConjugateGradient,
}).unwrap();
println!("Minimized: {:.2} -> {:.2} in {} steps",
    result.initial_energy, result.final_energy, result.n_steps);
```

## Chemical Reactions (SMIRKS)

```rust
use cyanea_chem::{parse_smiles, reaction::*};

// Parse a SMIRKS reaction (amide bond formation)
let rxn = parse_reaction("[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]").unwrap();

// Apply to a substrate
let acid = parse_smiles("CC(=O)O").unwrap(); // acetic acid
let products = apply_reaction(&acid, &rxn).unwrap();
for p in &products {
    println!("Product: {}", p.smiles);
}

// Retrosynthetic disconnection
let target = parse_smiles("CC(=O)NC").unwrap(); // N-methylacetamide
let disconnections = retrosynthetic_disconnection(&target).unwrap();
for d in &disconnections {
    println!("{}: precursors = {:?}", d.transform_name, d.precursors);
}
```

## Gasteiger Partial Charges

```rust
use cyanea_chem::{parse_smiles, gasteiger::gasteiger_charges};

let mol = parse_smiles("C(=O)O").unwrap(); // formic acid
let charges = gasteiger_charges(&mol).unwrap();
for (i, q) in charges.iter().enumerate() {
    println!("Atom {}: charge = {:.4}", i, q);
}
```
