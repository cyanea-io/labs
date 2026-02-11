use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use cyanea_align::msa::progressive_msa;
use cyanea_align::scoring::{ScoringMatrix, ScoringScheme};
use cyanea_align::simd::{banded_nw, banded_sw};
use cyanea_align::{align, align_batch, AlignmentMode};

fn dna_scheme() -> ScoringScheme {
    ScoringScheme::Simple(ScoringMatrix::dna_default())
}

fn random_dna(len: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    // Deterministic pseudo-random for reproducibility
    let mut seq = Vec::with_capacity(len);
    let mut state: u64 = 42;
    for _ in 0..len {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        seq.push(bases[((state >> 33) % 4) as usize]);
    }
    seq
}

fn mutate_dna(seq: &[u8], rate: f64) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut out = seq.to_vec();
    let mut state: u64 = 137;
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

fn bench_nw_sw(c: &mut Criterion) {
    let scoring = dna_scheme();

    let mut group = c.benchmark_group("pairwise");

    for &len in &[100, 1000] {
        let q = random_dna(len);
        let t = mutate_dna(&q, 0.1);

        group.bench_with_input(BenchmarkId::new("nw", len), &len, |b, _| {
            b.iter(|| align(black_box(&q), black_box(&t), AlignmentMode::Global, &scoring))
        });

        group.bench_with_input(BenchmarkId::new("sw", len), &len, |b, _| {
            b.iter(|| align(black_box(&q), black_box(&t), AlignmentMode::Local, &scoring))
        });
    }

    group.finish();
}

fn bench_banded(c: &mut Criterion) {
    let scoring = dna_scheme();
    let mut group = c.benchmark_group("banded");

    for &len in &[1000, 10_000] {
        let q = random_dna(len);
        let t = mutate_dna(&q, 0.05);

        group.bench_with_input(BenchmarkId::new("banded_sw_w50", len), &len, |b, _| {
            b.iter(|| banded_sw(black_box(&q), black_box(&t), &scoring, 50))
        });

        group.bench_with_input(BenchmarkId::new("banded_nw_w50", len), &len, |b, _| {
            b.iter(|| banded_nw(black_box(&q), black_box(&t), &scoring, 50))
        });
    }

    group.finish();
}

fn bench_batch(c: &mut Criterion) {
    let scoring = dna_scheme();
    let mut group = c.benchmark_group("batch");

    let pairs: Vec<(Vec<u8>, Vec<u8>)> = (0..1000)
        .map(|i| {
            let q = random_dna(100);
            let mut t = mutate_dna(&q, 0.1);
            // Vary slightly per pair for realism
            let adjust = (i % 10) as u64;
            if !t.is_empty() {
                t[0] = [b'A', b'C', b'G', b'T'][(adjust % 4) as usize];
            }
            (q, t)
        })
        .collect();
    let pair_refs: Vec<(&[u8], &[u8])> = pairs.iter().map(|(q, t)| (q.as_slice(), t.as_slice())).collect();

    group.bench_function("1000_pairs_100bp", |b| {
        b.iter(|| align_batch(black_box(&pair_refs), AlignmentMode::Local, &scoring))
    });

    group.finish();
}

fn bench_msa(c: &mut Criterion) {
    let scoring = dna_scheme();
    let mut group = c.benchmark_group("msa");

    let base = random_dna(200);
    let seqs: Vec<Vec<u8>> = (0..10)
        .map(|i| {
            let rate = 0.05 + (i as f64 * 0.01);
            mutate_dna(&base, rate)
        })
        .collect();
    let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

    group.bench_function("10_seqs_200bp", |b| {
        b.iter(|| progressive_msa(black_box(&seq_refs), &scoring))
    });

    group.finish();
}

criterion_group!(benches, bench_nw_sw, bench_banded, bench_batch, bench_msa);
criterion_main!(benches);
