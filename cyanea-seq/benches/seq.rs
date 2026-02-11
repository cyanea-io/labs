use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use cyanea_seq::{DnaSequence, KmerIter};
use std::io::Write;

fn random_dna(len: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut seq = Vec::with_capacity(len);
    let mut state: u64 = 42;
    for _ in 0..len {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        seq.push(bases[((state >> 33) % 4) as usize]);
    }
    seq
}

fn make_fasta(n_seqs: usize, seq_len: usize) -> tempfile::NamedTempFile {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    let bases = [b'A', b'C', b'G', b'T'];
    let mut state: u64 = 42;
    for i in 0..n_seqs {
        writeln!(f, ">seq_{}", i).unwrap();
        for _ in 0..seq_len {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            f.write_all(&[bases[((state >> 33) % 4) as usize]]).unwrap();
        }
        writeln!(f).unwrap();
    }
    f.flush().unwrap();
    f
}

fn make_fastq(n_seqs: usize, seq_len: usize) -> tempfile::NamedTempFile {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    let bases = [b'A', b'C', b'G', b'T'];
    let mut state: u64 = 42;
    for i in 0..n_seqs {
        writeln!(f, "@seq_{}", i).unwrap();
        for _ in 0..seq_len {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            f.write_all(&[bases[((state >> 33) % 4) as usize]]).unwrap();
        }
        writeln!(f).unwrap();
        writeln!(f, "+").unwrap();
        for _ in 0..seq_len {
            f.write_all(b"I").unwrap();
        }
        writeln!(f).unwrap();
    }
    f.flush().unwrap();
    f
}

fn bench_fasta_parse(c: &mut Criterion) {
    let mut group = c.benchmark_group("fasta_parse");

    // ~1MB: 1000 seqs × 1000bp
    let f_1mb = make_fasta(1000, 1000);
    group.bench_function("1MB", |b| {
        b.iter(|| cyanea_seq::parse_fasta_stats(black_box(f_1mb.path().to_str().unwrap())))
    });

    // ~10MB: 10000 seqs × 1000bp
    let f_10mb = make_fasta(10000, 1000);
    group.bench_function("10MB", |b| {
        b.iter(|| cyanea_seq::parse_fasta_stats(black_box(f_10mb.path().to_str().unwrap())))
    });

    group.finish();
}

fn bench_fastq_parse(c: &mut Criterion) {
    let mut group = c.benchmark_group("fastq_parse");

    let f_1mb = make_fastq(3000, 150);
    group.bench_function("1MB", |b| {
        b.iter(|| cyanea_seq::parse_fastq_file(black_box(f_1mb.path().to_str().unwrap())))
    });

    group.finish();
}

fn bench_gc_content(c: &mut Criterion) {
    let mut group = c.benchmark_group("gc_content");

    let seq_10k = DnaSequence::new(&random_dna(10_000)).unwrap();
    group.bench_function("10kb", |b| {
        b.iter(|| black_box(&seq_10k).gc_content())
    });

    group.finish();
}

fn bench_kmer(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer");

    let seq_10k = random_dna(10_000);

    group.bench_with_input(BenchmarkId::new("iter_k21", 10_000), &seq_10k, |b, seq| {
        b.iter(|| {
            let iter = KmerIter::new(black_box(seq), 21).unwrap();
            iter.count()
        })
    });

    group.finish();
}

criterion_group!(benches, bench_fasta_parse, bench_fastq_parse, bench_gc_content, bench_kmer);
criterion_main!(benches);
