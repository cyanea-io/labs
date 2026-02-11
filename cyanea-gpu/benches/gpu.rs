use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use cyanea_gpu::{Backend, CpuBackend, DistanceMetricGpu};

/// Deterministic pseudo-random f64 values in [0, 1) using an LCG.
fn random_data(n: usize, seed: u64) -> Vec<f64> {
    let mut state = seed;
    (0..n)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 11) as f64 / (1u64 << 53) as f64
        })
        .collect()
}

fn bench_pairwise_distance_matrix(c: &mut Criterion) {
    let backend = CpuBackend::new();
    let mut group = c.benchmark_group("pairwise_distance_matrix");

    let dim = 10;
    for &n in &[100, 500, 1000] {
        let data = random_data(n * dim, 42);
        let buf = backend.buffer_from_slice(&data).unwrap();

        group.bench_with_input(
            BenchmarkId::new("euclidean", format!("{}x{}", n, dim)),
            &n,
            |b, _| {
                b.iter(|| {
                    backend
                        .pairwise_distance_matrix(
                            black_box(&buf),
                            n,
                            dim,
                            DistanceMetricGpu::Euclidean,
                        )
                        .unwrap()
                })
            },
        );
    }

    group.finish();
}

fn bench_matrix_multiply(c: &mut Criterion) {
    let backend = CpuBackend::new();
    let mut group = c.benchmark_group("matrix_multiply");

    for &size in &[100, 200] {
        let a_data = random_data(size * size, 42);
        let b_data = random_data(size * size, 137);
        let a_buf = backend.buffer_from_slice(&a_data).unwrap();
        let b_buf = backend.buffer_from_slice(&b_data).unwrap();

        group.bench_with_input(
            BenchmarkId::new("square", format!("{}x{}", size, size)),
            &size,
            |b, _| {
                b.iter(|| {
                    backend
                        .matrix_multiply(black_box(&a_buf), black_box(&b_buf), size, size, size)
                        .unwrap()
                })
            },
        );
    }

    group.finish();
}

fn bench_reduce_sum(c: &mut Criterion) {
    let backend = CpuBackend::new();
    let mut group = c.benchmark_group("reduce_sum");

    for &n in &[10_000, 100_000] {
        let data = random_data(n, 42);
        let buf = backend.buffer_from_slice(&data).unwrap();

        group.bench_with_input(BenchmarkId::new("sum", n), &n, |b, _| {
            b.iter(|| backend.reduce_sum(black_box(&buf)).unwrap())
        });
    }

    group.finish();
}

fn bench_elementwise_map(c: &mut Criterion) {
    let backend = CpuBackend::new();
    let mut group = c.benchmark_group("elementwise_map");

    let sqrt_fn: &dyn Fn(f64) -> f64 = &|x: f64| x.sqrt();

    for &n in &[10_000, 100_000] {
        let data = random_data(n, 42);
        let input = backend.buffer_from_slice(&data).unwrap();
        let mut output = backend.buffer_zeros(n).unwrap();

        group.bench_with_input(BenchmarkId::new("sqrt", n), &n, |b, _| {
            b.iter(|| {
                backend
                    .elementwise_map(black_box(&input), black_box(&mut output), sqrt_fn)
                    .unwrap()
            })
        });
    }

    group.finish();
}

fn bench_batch_pairwise(c: &mut Criterion) {
    let backend = CpuBackend::new();
    let mut group = c.benchmark_group("batch_pairwise");

    let n = 50;
    let item_len = 100;
    let data = random_data(n * item_len, 42);
    let buf = backend.buffer_from_slice(&data).unwrap();

    // Use Euclidean distance as the pairwise function
    let euclidean: &dyn Fn(&[f64], &[f64]) -> f64 = &|a: &[f64], b: &[f64]| {
        a.iter()
            .zip(b)
            .map(|(x, y)| (x - y).powi(2))
            .sum::<f64>()
            .sqrt()
    };

    group.bench_function("50_items_100dim", |b| {
        b.iter(|| {
            backend
                .batch_pairwise(black_box(&buf), n, item_len, euclidean)
                .unwrap()
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_pairwise_distance_matrix,
    bench_matrix_multiply,
    bench_reduce_sum,
    bench_elementwise_map,
    bench_batch_pairwise
);
criterion_main!(benches);
