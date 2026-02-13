# Cyanea Benchmarks

Competitive benchmarks comparing Cyanea against rust-bio, plus SIMD backend performance matrix.

## Environment

| Parameter | Value |
|-----------|-------|
| Hardware | Apple M4, 10-core (4P + 6E), 16 GB |
| Rust | 1.93.0 (2026-01-19) |
| rust-bio | 3.0.0 |
| Cyanea | 0.1.0 |

## Smith-Waterman Alignment

Scoring: match=+2, mismatch=-1, gap_open=-5, gap_extend=-2. Sequences: deterministic LCG with ~10% divergence.

| Length | Cyanea (full) | Cyanea (SIMD score) | rust-bio | SIMD speedup vs rust-bio |
|-------:|--------------:|--------------------:|---------:|-------------------------:|
| 100bp | 32.7 us | 31.3 us | 49.3 us | 1.6x |
| 500bp | 737.9 us | 278.2 us | 1.67 ms | 6.0x |
| 1,000bp | 2.95 ms | 850.3 us | 7.09 ms | 8.3x |
| 5,000bp | 116.3 ms | 19.5 ms | 178.0 ms | 9.1x |
| 10,000bp | 335.1 ms | 73.7 ms | 900.3 ms | 12.2x |

## FM-Index Build Time

| Text size | Cyanea | rust-bio | Speedup |
|----------:|-------:|---------:|--------:|
| 10kb | 299.1 us | 286.7 us | 0.96x |
| 100kb | 4.04 ms | 4.61 ms | 1.14x |
| 1MB | 60.3 ms | 55.7 ms | 0.92x |

## FM-Index Count

100kb text, pattern extracted from text.

| Pattern length | Cyanea | rust-bio | Speedup |
|---------------:|-------:|---------:|--------:|
| 20bp | 25.6 ns | 133.3 ns | 5.2x |
| 50bp | 87.1 ns | 303.0 ns | 3.5x |
| 100bp | 190.9 ns | 687.4 ns | 3.6x |

## FM-Index Search (positions)

100kb text, pattern extracted from text.

| Pattern length | Cyanea | rust-bio | Speedup |
|---------------:|-------:|---------:|--------:|
| 20bp | 37.2 ns | 129.7 ns | 3.5x |
| 50bp | 98.3 ns | 315.0 ns | 3.2x |
| 100bp | 201.0 ns | 739.5 ns | 3.7x |

## SIMD Speedup (score-only SW)

SIMD backend: NEON (aarch64, Apple M4). 8 x i16 lanes, Farrar striped.

| Length | Scalar | SIMD (NEON) | Speedup |
|-------:|-------:|------------:|--------:|
| 100bp | 14.3 us | 30.8 us | 0.46x |
| 500bp | 370.5 us | 280.4 us | 1.32x |
| 1,000bp | 1.49 ms | 870.3 us | 1.71x |
| 5,000bp | 37.1 ms | 19.6 ms | 1.90x |
| 10,000bp | 148.4 ms | 74.6 ms | 1.99x |

Note: At 100bp the SIMD profile construction overhead dominates. Crossover occurs around ~200bp. At 10,000bp SIMD achieves ~2x scalar throughput.

## Banded Alignment Scaling

10,000bp sequences, ~5% divergence. Banded SW (with traceback) vs unbanded SIMD score-only.

| Bandwidth | Time | Relative to unbanded |
|----------:|-----:|---------------------:|
| 10 | 1.66 ms | 0.023x |
| 50 | 5.60 ms | 0.077x |
| 100 | 11.4 ms | 0.157x |
| 500 | 52.6 ms | 0.726x |
| unbanded (SIMD) | 72.5 ms | 1.0x |

## Reproducing

```bash
# Competitive benchmarks (vs rust-bio)
cargo bench -p cyanea-benchmarks

# SIMD matrix
cargo bench -p cyanea-align -- simd

# Banded scaling
cargo bench -p cyanea-align -- banded_scaling

# HTML reports are generated in target/criterion/
```
