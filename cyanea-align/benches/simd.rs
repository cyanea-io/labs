//! SIMD backend matrix benchmarks.
//!
//! Measures SIMD vs scalar performance for score-only Smith-Waterman,
//! and banded alignment scaling across different bandwidths.

use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput,
};

use cyanea_align::scoring::{ScoringMatrix, ScoringScheme};
use cyanea_align::simd::banded_sw;
use cyanea_align::{sw_scalar_score, sw_simd_score};

fn dna_scheme() -> ScoringScheme {
    ScoringScheme::Simple(ScoringMatrix::dna_default())
}

fn random_dna(len: usize, seed: u64) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut seq = Vec::with_capacity(len);
    let mut state = seed;
    for _ in 0..len {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        seq.push(bases[((state >> 33) % 4) as usize]);
    }
    seq
}

fn mutate_dna(seq: &[u8], rate: f64, seed: u64) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut out = seq.to_vec();
    let mut state = seed;
    for b in out.iter_mut() {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let r = (state >> 33) as f64 / (u32::MAX as f64);
        if r < rate {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            *b = bases[((state >> 33) % 4) as usize];
        }
    }
    out
}

// =========================================================================
// SIMD vs Scalar — score-only Smith-Waterman
// =========================================================================

fn bench_simd_vs_scalar(c: &mut Criterion) {
    let scoring = dna_scheme();
    let mut group = c.benchmark_group("simd_vs_scalar");

    for &len in &[100, 500, 1_000, 5_000, 10_000] {
        group.throughput(Throughput::Elements(len as u64));

        let query = random_dna(len, 42);
        let target = mutate_dna(&query, 0.10, 137);

        // Scalar baseline
        group.bench_with_input(
            BenchmarkId::new("scalar", len),
            &len,
            |b, _| {
                b.iter(|| {
                    sw_scalar_score(black_box(&query), black_box(&target), &scoring).unwrap()
                })
            },
        );

        // SIMD (NEON on aarch64, AVX2/SSE4.1 on x86_64)
        group.bench_with_input(
            BenchmarkId::new("simd", len),
            &len,
            |b, _| {
                b.iter(|| {
                    sw_simd_score(black_box(&query), black_box(&target), &scoring).unwrap()
                })
            },
        );
    }

    group.finish();
}

// =========================================================================
// Banded alignment scaling — fixed 10,000bp, varying bandwidth
// =========================================================================

fn bench_banded_scaling(c: &mut Criterion) {
    let scoring = dna_scheme();
    let len = 10_000;
    let query = random_dna(len, 42);
    let target = mutate_dna(&query, 0.05, 137);

    let mut group = c.benchmark_group("banded_scaling");
    group.throughput(Throughput::Elements(len as u64));

    for &bw in &[10, 50, 100, 500] {
        group.bench_with_input(
            BenchmarkId::new("banded_sw", bw),
            &bw,
            |b, _| {
                b.iter(|| banded_sw(black_box(&query), black_box(&target), &scoring, bw))
            },
        );
    }

    // Unbanded (full SW, score-only via SIMD)
    group.bench_function("unbanded_simd", |b| {
        b.iter(|| sw_simd_score(black_box(&query), black_box(&target), &scoring).unwrap())
    });

    group.finish();
}

// =========================================================================
// Criterion harness
// =========================================================================

criterion_group!(benches, bench_simd_vs_scalar, bench_banded_scaling);
criterion_main!(benches);
