use criterion::{black_box, criterion_group, criterion_main, Criterion};
use cyanea_ml::cluster::{dbscan, kmeans, DbscanConfig, KMeansConfig};
use cyanea_ml::distance::{pairwise_distances, DistanceMetric};
use cyanea_ml::reduction::{pca, tsne, PcaConfig, TsneConfig};

fn random_flat(n: usize, d: usize, seed: u64) -> Vec<f64> {
    let mut state = seed;
    (0..n * d)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 11) as f64 / (1u64 << 53) as f64
        })
        .collect()
}

fn random_points(n: usize, d: usize, seed: u64) -> Vec<Vec<f64>> {
    let mut state = seed;
    (0..n)
        .map(|_| {
            (0..d)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    (state >> 11) as f64 / (1u64 << 53) as f64
                })
                .collect()
        })
        .collect()
}

fn bench_kmeans(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmeans");

    let points = random_points(10_000, 5, 42);
    let refs: Vec<&[f64]> = points.iter().map(|p| p.as_slice()).collect();

    let config = KMeansConfig {
        n_clusters: 5,
        max_iter: 100,
        tolerance: 1e-4,
        seed: 42,
    };

    group.bench_function("10k_pts_k5", |b| {
        b.iter(|| kmeans(black_box(&refs), &config))
    });

    group.finish();
}

fn bench_dbscan(c: &mut Criterion) {
    let mut group = c.benchmark_group("dbscan");

    let points = random_points(5_000, 3, 42);
    let refs: Vec<&[f64]> = points.iter().map(|p| p.as_slice()).collect();

    let config = DbscanConfig {
        eps: 0.1,
        min_samples: 5,
        metric: DistanceMetric::Euclidean,
    };

    group.bench_function("5k_pts", |b| {
        b.iter(|| dbscan(black_box(&refs), &config))
    });

    group.finish();
}

fn bench_tsne(c: &mut Criterion) {
    let mut group = c.benchmark_group("tsne");
    group.sample_size(10); // t-SNE is slow

    let n_features = 10;
    let data = random_flat(500, n_features, 42);

    let config = TsneConfig {
        n_components: 2,
        perplexity: 30.0,
        learning_rate: 200.0,
        n_iter: 250,
        seed: 42,
    };

    group.bench_function("500_pts", |b| {
        b.iter(|| tsne(black_box(&data), n_features, &config))
    });

    group.finish();
}

fn bench_pairwise_distances(c: &mut Criterion) {
    let mut group = c.benchmark_group("pairwise_dist");

    let points = random_points(1_000, 10, 42);
    let refs: Vec<&[f64]> = points.iter().map(|p| p.as_slice()).collect();

    group.bench_function("1k_euclidean", |b| {
        b.iter(|| pairwise_distances(black_box(&refs), DistanceMetric::Euclidean))
    });

    group.finish();
}

fn bench_pca(c: &mut Criterion) {
    let mut group = c.benchmark_group("pca");

    let n_features = 50;
    let data = random_flat(1_000, n_features, 42);

    let config = PcaConfig {
        n_components: 2,
        max_iter: 100,
        tolerance: 1e-6,
    };

    group.bench_function("1k_x50_2comp", |b| {
        b.iter(|| pca(black_box(&data), n_features, &config))
    });

    group.finish();
}

criterion_group!(benches, bench_kmeans, bench_dbscan, bench_tsne, bench_pairwise_distances, bench_pca);
criterion_main!(benches);
