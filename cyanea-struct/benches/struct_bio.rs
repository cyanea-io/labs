use criterion::{black_box, criterion_group, criterion_main, Criterion};
use cyanea_struct::contact::compute_contact_map;
use cyanea_struct::pdb::parse_pdb;
use cyanea_struct::superposition::kabsch_points;
use cyanea_struct::types::Point3D;

/// Generate a synthetic PDB string with `n_residues` residues in chain A.
/// Each residue has 4 backbone atoms (N, CA, C, O) placed in a rough alpha-helix geometry.
fn synthetic_pdb(n_residues: usize) -> String {
    let mut lines = Vec::new();
    lines.push("HEADER                                                        BENCH".to_string());

    let mut serial = 1;
    let residues = ["ALA", "GLY", "VAL", "LEU", "ILE"];
    for i in 0..n_residues {
        let resname = residues[i % residues.len()];
        let resseq = i + 1;
        // Approximate alpha-helix: ~1.5 A rise per residue, 100 degree turn
        let angle = (i as f64) * 100.0_f64.to_radians();
        let rise = i as f64 * 1.5;
        let radius = 2.3;

        let atoms = [
            ("N ", 0.0_f64, 0.0_f64, 0.0_f64),
            ("CA", 1.458, 0.0, 0.0),
            ("C ", 2.009, 1.420, 0.0),
            ("O ", 1.246, 2.390, 0.0),
        ];

        for (name, dx, dy, dz) in &atoms {
            let x = radius * angle.cos() + dx;
            let y = radius * angle.sin() + dy;
            let z = rise + dz;
            lines.push(format!(
                "ATOM  {:>5} {:<4}{} A{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00           {}",
                serial,
                name,
                resname,
                resseq,
                x,
                y,
                z,
                &name[..1]
            ));
            serial += 1;
        }
    }
    lines.push("TER".to_string());
    lines.push("END".to_string());
    lines.join("\n")
}

fn bench_pdb_parse(c: &mut Criterion) {
    let mut group = c.benchmark_group("pdb_parse");

    // ~10k atoms = 2500 residues × 4 atoms
    let pdb_10k = synthetic_pdb(2500);

    group.bench_function("10k_atoms", |b| {
        b.iter(|| parse_pdb(black_box(&pdb_10k)))
    });

    group.finish();
}

fn bench_kabsch(c: &mut Criterion) {
    let mut group = c.benchmark_group("kabsch");

    // Generate 1000 CA-like points
    let mut state: u64 = 42;
    let points_a: Vec<Point3D> = (0..1000)
        .map(|i| {
            let angle = (i as f64) * 100.0_f64.to_radians();
            Point3D {
                x: 2.3 * angle.cos(),
                y: 2.3 * angle.sin(),
                z: i as f64 * 1.5,
            }
        })
        .collect();

    // Slightly perturbed copy
    let points_b: Vec<Point3D> = points_a
        .iter()
        .map(|p| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let noise = (state >> 33) as f64 / (u32::MAX as f64) * 0.5;
            Point3D {
                x: p.x + noise,
                y: p.y - noise * 0.5,
                z: p.z + noise * 0.3,
            }
        })
        .collect();

    group.bench_function("1k_ca_atoms", |b| {
        b.iter(|| kabsch_points(black_box(&points_a), black_box(&points_b)))
    });

    group.finish();
}

fn bench_contact_map(c: &mut Criterion) {
    let mut group = c.benchmark_group("contact_map");

    let pdb = synthetic_pdb(500);
    let structure = parse_pdb(&pdb).unwrap();
    let chain = structure.get_chain('A').unwrap();

    group.bench_function("500_residues", |b| {
        b.iter(|| compute_contact_map(black_box(chain)))
    });

    group.finish();
}

criterion_group!(benches, bench_pdb_parse, bench_kabsch, bench_contact_map);
criterion_main!(benches);
