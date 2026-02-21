# cyanea-core API Reference

Shared foundation for the Cyanea bioinformatics ecosystem. Defines common traits, error types, content-addressed hashing, compression, memory-mapped file access, log-space probability types, rank/select bitvectors with wavelet matrix, and Fenwick tree.

## Public API

### Traits (`traits.rs`)

| Trait | Description |
|-------|-------------|
| `Sequence` | Biological sequence interface (`len`, `is_empty`, `as_bytes`) |
| `ContentAddressable` | Content-addressed hashing (`content_hash`) |
| `Compressible` | Compress/decompress operations |
| `Scored` | Values with a numeric score |
| `Annotated` | Named and described items (`name`, `description`) |
| `Summarizable` | Summary string representation |

### Error types (`error.rs`)

| Type | Description |
|------|-------------|
| `CyaneaError` | Unified error enum: `Io`, `Parse`, `InvalidInput`, `Compression`, `Hash`, `Other` |
| `Result<T>` | Alias for `std::result::Result<T, CyaneaError>` |

All crates in the workspace use `CyaneaError` via `thiserror` 2.x.

### Hashing (`hash.rs`)

| Function | Description |
|----------|-------------|
| `sha256(data: &[u8]) -> String` | SHA-256 hex digest of in-memory data |
| `sha256_file(path) -> Result<String>` | SHA-256 hex digest of a file (streaming) |

### Compression (`compress.rs`, `std` feature only)

| Function | Description |
|----------|-------------|
| `zstd_compress(data, level) -> Result<Vec<u8>>` | Zstd compression (levels 1-22) |
| `zstd_decompress(data) -> Result<Vec<u8>>` | Zstd decompression |
| `gzip_compress(data, level) -> Result<Vec<u8>>` | Gzip compression |
| `gzip_decompress(data) -> Result<Vec<u8>>` | Gzip decompression |
| `detect_algorithm(data) -> Option<Algorithm>` | Detect compression from magic bytes |
| `decompress(data) -> Result<Vec<u8>>` | Auto-detect and decompress |

### Memory mapping (`mmap.rs`, `std` feature only)

| Type | Description |
|------|-------------|
| `MappedFile` | Read-only memory-mapped file access |

### Probability types (`prob.rs`)

Log-space probability newtypes for numerically stable computation.

| Type | Description |
|------|-------------|
| `LogProb(f64)` | Probability as natural logarithm `ln(p)` |
| `PhredProb(f64)` | Probability as Phred quality score `-10 * log10(p)` |

**LogProb methods:**

| Method | Description |
|--------|-------------|
| `from_prob(p) -> Result<Self>` | Create from raw probability in (0, 1] |
| `to_prob() -> f64` | Convert back to raw probability |
| `ln_add(other) -> Self` | Log-sum-exp (log-space addition) |
| `ln_mul(other) -> Self` | Log-space multiplication |
| `certain() -> Self` | Certain event: ln(1) = 0 |
| `impossible() -> Self` | Impossible event: ln(0) = -inf |

**PhredProb methods:**

| Method | Description |
|--------|-------------|
| `from_phred(q) -> Result<Self>` | Create from Phred score (>= 0) |
| `from_prob(p) -> Result<Self>` | Create from raw probability in (0, 1] |
| `to_phred() -> f64` | The Phred quality score |
| `to_prob() -> f64` | Convert to raw probability |

Bidirectional `From` conversions between `PhredProb` and `LogProb`.

### Rank/select bitvectors and wavelet matrix (`bitvec.rs`)

| Type | Description |
|------|-------------|
| `RankSelectBitVec` | Bitvector with O(1) rank and O(log n) select via u64 blocks + superblock index |
| `WaveletMatrix` | Wavelet matrix over integer alphabet [0, sigma) with O(log sigma) access/rank/select |

**RankSelectBitVec methods:**

| Method | Description |
|--------|-------------|
| `build(bits: &[bool]) -> Self` | Build from boolean slice |
| `get(i) -> bool` | Get bit at position i |
| `rank1(i) -> usize` | Count 1-bits in [0, i) |
| `rank0(i) -> usize` | Count 0-bits in [0, i) |
| `select1(k) -> Option<usize>` | Position of k-th 1-bit (1-indexed) |
| `select0(k) -> Option<usize>` | Position of k-th 0-bit (1-indexed) |
| `len() -> usize` | Total number of bits |
| `count_ones() -> usize` | Total number of 1-bits |
| `count_zeros() -> usize` | Total number of 0-bits |

**WaveletMatrix methods:**

| Method | Description |
|--------|-------------|
| `build(symbols, sigma) -> Result<Self>` | Build from symbol sequence over [0, sigma) |
| `access(i) -> Option<usize>` | Access symbol at position i |
| `rank(c, i) -> usize` | Count occurrences of symbol c in [0, i) |
| `select(c, k) -> Option<usize>` | Position of k-th occurrence of symbol c (1-indexed) |
| `len() -> usize` | Length of indexed sequence |
| `sigma() -> usize` | Alphabet size |

### Fenwick tree (`fenwick.rs`)

Generic Binary Indexed Tree for O(log n) prefix sum and point update queries.

| Type | Description |
|------|-------------|
| `FenwickTree<T>` | Generic Fenwick tree (T: Copy + Default + Add + Sub) |

**FenwickTree methods:**

| Method | Description |
|--------|-------------|
| `new(n) -> Self` | Create tree of size n initialized to zero |
| `from_slice(values) -> Self` | Build from slice in O(n) time |
| `update(i, delta)` | Add delta to element at index i (0-based) |
| `prefix_sum(i) -> T` | Sum of elements in [0, i] (inclusive, 0-based) |
| `range_sum(l, r) -> T` | Sum of elements in [l, r] (inclusive, 0-based) |
| `len() -> usize` | Number of elements |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Enables compression (`zstd`, `flate2`) and memory mapping (`memmap2`) |
| `wasm` | No | WASM target (no-op currently) |
| `serde` | No | Enables `serde::Serialize`/`Deserialize` derives |

## Dependencies

- `thiserror` 2.x -- error derives
- `sha2`, `hex` -- SHA-256 hashing
- `memmap2` -- memory-mapped files (std)
- `zstd` -- Zstd compression (std)
- `flate2` -- Gzip compression (std)

## Tests

58 tests across `hash.rs`, `compress.rs`, `prob.rs`, `bitvec.rs`, and `fenwick.rs`.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 28 | Module declarations, re-exports |
| `error.rs` | 34 | `CyaneaError` enum |
| `traits.rs` | 57 | Core trait definitions |
| `hash.rs` | 82 | SHA-256 hashing |
| `compress.rs` | 134 | Zstd/Gzip compression |
| `mmap.rs` | 83 | Memory-mapped file access |
| `prob.rs` | 247 | LogProb/PhredProb newtypes |
| `bitvec.rs` | 592 | Rank/select bitvectors, wavelet matrix |
| `fenwick.rs` | 221 | Generic Fenwick tree (BIT) |
