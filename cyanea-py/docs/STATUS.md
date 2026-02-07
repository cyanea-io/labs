# cyanea-py

Python bindings for the Cyanea bioinformatics ecosystem via PyO3. Installable as `pip install cyanea` (via maturin).

## Status: Complete

All bindings are implemented and tested. Wraps cyanea-seq, cyanea-align, cyanea-stats, and cyanea-core.

## Build

Requires maturin and Rust 1.93+.

```bash
# In a virtualenv:
pip install maturin
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release
```

The `PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1` flag is needed for Python 3.14+ (PyO3 0.23 max supported is 3.13).

## Python API

### `cyanea.seq`

**Classes:**

| Class | Description |
|-------|-------------|
| `DnaSequence(data: bytes)` | Validated DNA sequence |
| `RnaSequence(data: bytes)` | Validated RNA sequence |
| `ProteinSequence(data: bytes)` | Validated protein sequence |
| `FastaStats` | FASTA file statistics (frozen, all fields readable) |
| `FastqStats` | FASTQ file statistics (frozen, all fields readable) |
| `FastqRecord` | Single FASTQ record (frozen) |

**DnaSequence methods:** `reverse_complement()`, `transcribe() -> RnaSequence`, `gc_content() -> float`, `kmers(k) -> list[bytes]`, `translate() -> ProteinSequence`, `__len__`, `__bytes__`, `__str__`, `__repr__`, `__eq__`

**RnaSequence methods:** `reverse_complement()`, `translate() -> ProteinSequence`, `reverse_transcribe() -> DnaSequence`, `kmers(k)`, dunders

**ProteinSequence methods:** `molecular_weight() -> float`, `kmers(k)`, dunders

**Functions:**

| Function | Description |
|----------|-------------|
| `fasta_stats(path) -> FastaStats` | Streaming FASTA statistics |
| `parse_fastq(path) -> list[FastqRecord]` | Parse all FASTQ records |
| `fastq_stats(path) -> FastqStats` | Streaming FASTQ statistics |

### `cyanea.align`

| Class/Function | Description |
|----------------|-------------|
| `AlignmentResult` | Frozen: `score`, `aligned_query`, `aligned_target`, `query_start/end`, `target_start/end`, `cigar_string`, `identity`, `matches`, `mismatches`, `gaps`, `length` |
| `align_dna(query, target, *, mode="local", match_score=2, mismatch_score=-1, gap_open=-5, gap_extend=-2)` | DNA alignment |
| `align_protein(query, target, *, mode="global", matrix="blosum62")` | Protein alignment |
| `align_batch(pairs, *, mode="local", ...)` | Batch DNA alignment |

### `cyanea.stats`

| Class/Function | Description |
|----------------|-------------|
| `DescriptiveStats` | Frozen, 15 fields (`count`, `mean`, `median`, `variance`, `std_dev`, ...) |
| `TestResult` | Frozen: `statistic`, `p_value`, `degrees_of_freedom`, `method` |
| `describe(data) -> DescriptiveStats` | Descriptive statistics |
| `pearson(x, y) -> float` | Pearson correlation |
| `spearman(x, y) -> float` | Spearman correlation |
| `t_test(data, *, mu=0.0) -> TestResult` | One-sample t-test |
| `t_test_two_sample(x, y, *, equal_var=False) -> TestResult` | Two-sample t-test |
| `mann_whitney_u(x, y) -> TestResult` | Mann-Whitney U test |
| `bonferroni(p_values) -> list[float]` | Bonferroni correction |
| `benjamini_hochberg(p_values) -> list[float]` | BH FDR correction |

### `cyanea.core`

| Function | Description |
|----------|-------------|
| `sha256(data: bytes) -> str` | SHA-256 hex digest |
| `sha256_file(path: str) -> str` | SHA-256 of a file (streaming) |
| `zstd_compress(data: bytes, *, level=3) -> bytes` | Zstd compression |
| `zstd_decompress(data: bytes) -> bytes` | Zstd decompression |

## Error Mapping

| Rust Error | Python Exception |
|------------|------------------|
| `CyaneaError::Io` | `IOError` |
| `CyaneaError::Parse`, `InvalidInput` | `ValueError` |
| `CyaneaError::Compression`, `Hash`, `Other` | `RuntimeError` |

## Dependencies

- `pyo3` 0.23 (extension-module)
- `cyanea-core`, `cyanea-seq`, `cyanea-io`, `cyanea-align`, `cyanea-stats`

## Tests

No Rust-level tests. Integration tested via Python smoke tests.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 35 | PyO3 module root, submodule registration in `sys.modules` |
| `error.rs` | 28 | `CyaneaError` -> `PyErr` mapping, `IntoPyResult` trait |
| `seq.rs` | 327 | Sequence types and FASTA/FASTQ parsing |
| `align.rs` | 153 | Alignment functions with keyword-only scoring params |
| `stats.rs` | 153 | Statistics and hypothesis testing |
| `core_utils.rs` | 44 | SHA-256 hashing and zstd compression |
