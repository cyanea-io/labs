use criterion::{black_box, criterion_group, criterion_main, Criterion};
use cyanea_stats::correlation::{pearson, CorrelationMatrix};
use cyanea_stats::descriptive::describe;
use cyanea_stats::reduction::pca;

fn random_f64(n: usize, seed: u64) -> Vec<f64> {
    let mut state = seed;
    (0..n)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 11) as f64 / (1u64 << 53) as f64
        })
        .collect()
}

fn random_matrix(rows: usize, cols: usize, seed: u64) -> Vec<Vec<f64>> {
    let mut state = seed;
    (0..rows)
        .map(|_| {
            (0..cols)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    (state >> 11) as f64 / (1u64 << 53) as f64
                })
                .collect()
        })
        .collect()
}

fn bench_describe(c: &mut Criterion) {
    let mut group = c.benchmark_group("describe");

    let data_100k = random_f64(100_000, 42);
    group.bench_function("100k_values", |b| {
        b.iter(|| describe(black_box(&data_100k)))
    });

    group.finish();
}

fn bench_pca(c: &mut Criterion) {
    let mut group = c.benchmark_group("pca");

    let data = random_matrix(1_000, 100, 42);
    let refs: Vec<&[f64]> = data.iter().map(|r| r.as_slice()).collect();

    group.bench_function("1k_x100_5comp", |b| {
        b.iter(|| pca(black_box(&refs), 5))
    });

    group.finish();
}

fn bench_correlation_matrix(c: &mut Criterion) {
    let mut group = c.benchmark_group("correlation");

    // 500 variables × 500 observations → 500×500 correlation matrix
    let vars: Vec<Vec<f64>> = (0..500)
        .map(|i| random_f64(500, 42 + i))
        .collect();
    let refs: Vec<&[f64]> = vars.iter().map(|v| v.as_slice()).collect();

    group.bench_function("500x500_matrix", |b| {
        b.iter(|| CorrelationMatrix::from_rows(black_box(&refs)))
    });

    group.finish();
}

fn bench_pearson(c: &mut Criterion) {
    let mut group = c.benchmark_group("pearson");

    let x = random_f64(10_000, 42);
    let y = random_f64(10_000, 137);

    group.bench_function("10k_values", |b| {
        b.iter(|| pearson(black_box(&x), black_box(&y)))
    });

    group.finish();
}

criterion_group!(benches, bench_describe, bench_pca, bench_correlation_matrix, bench_pearson);
criterion_main!(benches);
