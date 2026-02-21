# cyanea-core Architecture

## Module Map

```
cyanea-core
 +-- error.rs       CyaneaError enum (thiserror 2.x)
 +-- traits.rs      Sequence, ContentAddressable, Compressible, Scored, Annotated, Summarizable
 +-- hash.rs        sha256(), sha256_file() via sha2 crate
 +-- prob.rs        LogProb, PhredProb newtypes
 +-- bitvec.rs      RankSelectBitVec, WaveletMatrix
 +-- fenwick.rs     FenwickTree<T>
 +-- compress.rs    [std] zstd/gzip compress/decompress, auto-detect
 +-- mmap.rs        [std] MappedFile via memmap2
```

No internal module has dependencies on another module within this crate, except that `compress.rs` and `mmap.rs` use `CyaneaError` from `error.rs`. All modules are leaf nodes in the dependency graph.

## Design Decisions

### CyaneaError

Uses `thiserror` 2.x with 6 variants covering all failure modes in the ecosystem:

- `Io(#[from] std::io::Error)` -- file and network operations
- `Parse(String)` -- format parsing (VCF, FASTA, PDB, etc.)
- `InvalidInput(String)` -- invalid arguments or data constraints
- `Compression(String)` -- zstd/gzip failures
- `Hash(String)` -- hashing failures
- `Other(String)` -- catch-all for unclassified errors

The `Result<T>` alias avoids repeating the error type across every function signature in the workspace.

### LogProb and PhredProb

`LogProb(f64)` stores `ln(p)` to avoid floating-point underflow when multiplying many small probabilities (common in HMMs, variant callers, and base quality processing). `PhredProb(f64)` stores `-10 * log10(p)` as used in FASTQ quality encoding.

Key arithmetic: `ln_add` uses the log-sum-exp trick (`ln(a + b) = ln(a) + ln(1 + exp(ln(b) - ln(a)))`) to add probabilities without leaving log-space. `ln_mul` is simply addition of log-values.

Bidirectional `From` conversions between the two types enable seamless use in different contexts without manual conversion.

### RankSelectBitVec

Uses 64-bit blocks (`u64`) with a superblock index for O(1) `rank` queries via hardware popcount. Each superblock covers 512 bits (8 blocks). The `rank1(i)` operation sums the superblock prefix, then counts individual blocks, then uses `count_ones()` with a bit mask for the partial final block.

`select1/select0` use binary search over the superblock index for O(log n) performance, which is sufficient for the bioinformatics workloads (FM-index, wavelet matrix) that drive usage.

### WaveletMatrix

Level-based decomposition for rank/select on alphabets of size sigma. Stores `ceil(log2(sigma))` bitvectors, one per level, each backed by `RankSelectBitVec`. This yields O(log sigma) `access`, `rank`, and `select` operations.

The wavelet matrix (vs. wavelet tree) uses a simpler, cache-friendlier layout: all symbols are partitioned at each level by their bit at that level, with 0-bits going left and 1-bits going right. No explicit tree structure is needed.

### FenwickTree

Generic binary indexed tree parameterized over `T: Copy + Default + Add<Output=T> + Sub<Output=T>`. Supports `i32`, `i64`, `f64`, and other numeric types.

Built from a slice in O(n) using the standard parent-propagation trick (rather than n individual updates). `prefix_sum` and `update` both run in O(log n) by traversing the implicit tree structure encoded in binary representations of indices.

Used internally by `cyanea-align::lcsk` for O(n log n) LCSk++ sparse alignment.

### Compression module

- **zstd**: default compression level 3 (good speed/ratio tradeoff). Levels 1-22 supported.
- **gzip**: via `flate2` with configurable level.
- **Auto-detection**: `detect_algorithm` checks magic bytes (`0x28 0xb5 0x2f 0xfd` for zstd, `0x1f 0x8b` for gzip). `decompress` combines detection with decompression.

### Memory mapping

`MappedFile` is a thin wrapper around `memmap2::Mmap` providing read-only memory-mapped access. Used for large file access (e.g., indexed FASTA reading, BAM random access) where loading the entire file into memory would be impractical.
