# cyanea-core

Shared foundation for the Cyanea bioinformatics ecosystem. Defines common traits, error types, content-addressed hashing, compression, and memory-mapped file access.

## Status: Complete

All planned functionality is implemented and tested.

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

14 tests across `hash.rs` and `compress.rs`.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 22 | Module declarations, re-exports |
| `error.rs` | 34 | `CyaneaError` enum |
| `traits.rs` | 57 | Core trait definitions |
| `hash.rs` | 82 | SHA-256 hashing |
| `compress.rs` | 134 | Zstd/Gzip compression |
| `mmap.rs` | 83 | Memory-mapped file access |
