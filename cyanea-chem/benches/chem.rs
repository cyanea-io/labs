use criterion::{black_box, criterion_group, criterion_main, Criterion};
use cyanea_chem::{morgan_fingerprint, parse_smiles, tanimoto_bulk, tanimoto_similarity, Fingerprint};

/// A set of representative drug-like SMILES strings
const SMILES_SET: &[&str] = &[
    "CCO",                              // ethanol
    "CC(=O)O",                          // acetic acid
    "c1ccccc1",                         // benzene
    "CC(=O)Oc1ccccc1C(=O)O",           // aspirin
    "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C", // testosterone
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",   // caffeine
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  // ibuprofen
    "OC(=O)C1=CC=CC=C1O",              // salicylic acid
    "C1=CC=C(C=C1)O",                  // phenol
    "CC(=O)NC1=CC=C(C=C1)O",           // acetaminophen
    "C(C(=O)O)N",                       // glycine
    "c1ccc2ccccc2c1",                   // naphthalene
    "C1CCCCC1",                         // cyclohexane
    "C(=O)(N)N",                        // urea
    "CC(O)CC",                          // 2-butanol
    "CCCCCCCC",                         // octane
    "c1ccncc1",                         // pyridine
    "C1=CN=CN=C1",                      // pyrimidine
    "c1cc[nH]c1",                       // pyrrole
    "C1=CSC=C1",                        // thiophene
];

fn bench_smiles_parse(c: &mut Criterion) {
    let mut group = c.benchmark_group("smiles_parse");

    // Parse 1k molecules (cycle through the 20 SMILES)
    let smiles_1k: Vec<&str> = SMILES_SET.iter().copied().cycle().take(1000).collect();

    group.bench_function("1k_mols", |b| {
        b.iter(|| {
            for &smi in black_box(&smiles_1k) {
                let _ = parse_smiles(smi);
            }
        })
    });

    group.finish();
}

fn bench_morgan_fp(c: &mut Criterion) {
    let mut group = c.benchmark_group("morgan_fp");

    let mols: Vec<_> = SMILES_SET
        .iter()
        .filter_map(|s| parse_smiles(s).ok())
        .collect();

    // Generate fingerprints for 1k molecules (cycle through parsed mols)
    let mols_1k: Vec<_> = mols.iter().cycle().take(1000).collect();

    group.bench_function("1k_mols_r2_2048", |b| {
        b.iter(|| {
            for mol in black_box(&mols_1k) {
                let _ = morgan_fingerprint(mol, 2, 2048);
            }
        })
    });

    group.finish();
}

fn bench_tanimoto(c: &mut Criterion) {
    let mut group = c.benchmark_group("tanimoto");

    let mols: Vec<_> = SMILES_SET
        .iter()
        .filter_map(|s| parse_smiles(s).ok())
        .collect();

    let fps: Vec<Fingerprint> = mols.iter().map(|m| morgan_fingerprint(m, 2, 2048)).collect();

    // Build 100 query FPs and 1000 target FPs
    let queries: Vec<Fingerprint> = fps.iter().cycle().take(100).cloned().collect();
    let targets: Vec<Fingerprint> = fps.iter().cycle().take(1000).cloned().collect();

    group.bench_function("100x1k_bulk", |b| {
        b.iter(|| {
            for q in black_box(&queries) {
                let _ = tanimoto_bulk(q, black_box(&targets));
            }
        })
    });

    // Single pair
    group.bench_function("single_pair", |b| {
        b.iter(|| tanimoto_similarity(black_box(&fps[0]), black_box(&fps[1])))
    });

    group.finish();
}

criterion_group!(benches, bench_smiles_parse, bench_morgan_fp, bench_tanimoto);
criterion_main!(benches);
