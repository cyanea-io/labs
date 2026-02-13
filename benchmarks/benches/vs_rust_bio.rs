//! Head-to-head benchmarks: Cyanea vs rust-bio.
//!
//! Compares Smith-Waterman alignment and FM-index construction/search
//! at multiple sequence lengths with identical scoring parameters.

use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput,
};

// --- Cyanea imports ---
use cyanea_align::scoring::{ScoringMatrix, ScoringScheme};
use cyanea_align::{smith_waterman, sw_simd_score};
use cyanea_seq::FmIndex;

// --- rust-bio imports ---
use bio::alignment::pairwise::Aligner;
use bio::alphabets::dna;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{BackwardSearchResult, FMIndex, FMIndexable};
use bio::data_structures::suffix_array::suffix_array;

// =========================================================================
// Sequence generation — deterministic LCG, ~10% divergence
// =========================================================================

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
// Smith-Waterman benchmarks
// =========================================================================

/// Scoring: match=+2, mismatch=-1, gap_open=-5, gap_extend=-2 (identical for both)
fn bench_smith_waterman(c: &mut Criterion) {
    let cyanea_scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());

    let mut group = c.benchmark_group("smith_waterman");

    for &len in &[100, 500, 1_000, 5_000, 10_000] {
        group.throughput(Throughput::Elements(len as u64));

        let query = random_dna(len, 42);
        let target = mutate_dna(&query, 0.10, 137);

        // Cyanea: full SW with traceback
        group.bench_with_input(
            BenchmarkId::new("cyanea_full", len),
            &len,
            |b, _| {
                b.iter(|| {
                    smith_waterman(black_box(&query), black_box(&target), &cyanea_scoring)
                        .unwrap()
                })
            },
        );

        // Cyanea: SIMD score-only (no traceback)
        group.bench_with_input(
            BenchmarkId::new("cyanea_simd_score", len),
            &len,
            |b, _| {
                b.iter(|| {
                    sw_simd_score(black_box(&query), black_box(&target), &cyanea_scoring)
                        .unwrap()
                })
            },
        );

        // rust-bio: local alignment (pairwise::Aligner)
        // rust-bio uses: gap_open = open + extend for first gap base,
        // then gap_extend for subsequent bases.
        // To match our affine model (open=-5, extend=-2): pass gap_open=-5, gap_extend=-2.
        group.bench_with_input(
            BenchmarkId::new("rust_bio", len),
            &len,
            |b, _| {
                b.iter(|| {
                    let mut aligner = Aligner::with_capacity(
                        query.len(),
                        target.len(),
                        -5i32,  // gap open
                        -2i32,  // gap extend
                        |a: u8, b: u8| {
                            if a == b { 2i32 } else { -1i32 }
                        },
                    );
                    aligner.local(black_box(&query), black_box(&target))
                })
            },
        );
    }

    group.finish();
}

// =========================================================================
// FM-index benchmarks
// =========================================================================

fn bench_fmindex_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("fmindex_build");

    for &len in &[10_000, 100_000, 1_000_000] {
        let text = random_dna(len, 99);

        group.throughput(Throughput::Elements(len as u64));

        // Cyanea: single-step build
        group.bench_with_input(
            BenchmarkId::new("cyanea", len),
            &len,
            |b, _| {
                b.iter(|| FmIndex::build(black_box(&text)))
            },
        );

        // rust-bio: multi-step construction
        let mut text_sentinel = text.clone();
        text_sentinel.push(b'$');
        let alphabet = dna::n_alphabet();

        group.bench_with_input(
            BenchmarkId::new("rust_bio", len),
            &len,
            |b, _| {
                b.iter(|| {
                    let sa = suffix_array(black_box(&text_sentinel));
                    let bwt_result = bwt(&text_sentinel, &sa);
                    let less_arr = less(&bwt_result, &alphabet);
                    let occ = Occ::new(&bwt_result, 3, &alphabet);
                    FMIndex::new(bwt_result, less_arr, occ)
                })
            },
        );
    }

    group.finish();
}

fn bench_fmindex_search(c: &mut Criterion) {
    let text_len = 100_000;
    let text = random_dna(text_len, 99);

    // Pre-build indexes
    let cyanea_fm = FmIndex::build(&text);

    let mut text_sentinel = text.clone();
    text_sentinel.push(b'$');
    let alphabet = dna::n_alphabet();
    let sa = suffix_array(&text_sentinel);
    let bwt_result = bwt(&text_sentinel, &sa);
    let less_arr = less(&bwt_result, &alphabet);
    let occ = Occ::new(&bwt_result, 3, &alphabet);
    let bio_fm = FMIndex::new(bwt_result, less_arr, occ);

    let mut group = c.benchmark_group("fmindex_count");

    for &pat_len in &[20, 50, 100] {
        // Extract a pattern from the text (guaranteed to exist)
        let pattern = text[1000..1000 + pat_len].to_vec();

        // Cyanea: count
        group.bench_with_input(
            BenchmarkId::new("cyanea", pat_len),
            &pat_len,
            |b, _| {
                b.iter(|| cyanea_fm.count(black_box(&pattern)))
            },
        );

        // rust-bio: backward_search interval size
        group.bench_with_input(
            BenchmarkId::new("rust_bio", pat_len),
            &pat_len,
            |b, _| {
                b.iter(|| {
                    let result = bio_fm.backward_search(black_box(pattern.iter()));
                    match result {
                        BackwardSearchResult::Complete(interval) => interval.occ(&sa).len(),
                        BackwardSearchResult::Partial(interval, _) => interval.occ(&sa).len(),
                        BackwardSearchResult::Absent => 0,
                    }
                })
            },
        );
    }

    group.finish();

    // --- Search (return positions) ---
    let mut group = c.benchmark_group("fmindex_search");

    for &pat_len in &[20, 50, 100] {
        let pattern = text[1000..1000 + pat_len].to_vec();

        // Cyanea: search (returns positions)
        group.bench_with_input(
            BenchmarkId::new("cyanea", pat_len),
            &pat_len,
            |b, _| {
                b.iter(|| cyanea_fm.search(black_box(&pattern)))
            },
        );

        // rust-bio: backward_search + occ
        group.bench_with_input(
            BenchmarkId::new("rust_bio", pat_len),
            &pat_len,
            |b, _| {
                b.iter(|| {
                    let result = bio_fm.backward_search(black_box(pattern.iter()));
                    match result {
                        BackwardSearchResult::Complete(interval) => interval.occ(&sa),
                        BackwardSearchResult::Partial(interval, _) => interval.occ(&sa),
                        BackwardSearchResult::Absent => vec![],
                    }
                })
            },
        );
    }

    group.finish();
}

// =========================================================================
// Criterion harness
// =========================================================================

criterion_group!(
    benches,
    bench_smith_waterman,
    bench_fmindex_build,
    bench_fmindex_search,
);
criterion_main!(benches);
