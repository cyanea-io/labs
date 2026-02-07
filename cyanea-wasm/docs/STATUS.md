# cyanea-wasm

WebAssembly bindings for browser-based bioinformatics. Wraps cyanea-seq, cyanea-align, cyanea-stats, and cyanea-ml into a JSON-based interface suitable for JavaScript/TypeScript consumption.

## Status: Complete

All planned bindings are implemented. Functions accept simple types (strings, numbers) and return JSON strings for easy interop. No wasm-bindgen annotations yet -- the crate provides the logic layer that a thin wasm-bindgen shim would wrap.

## Public API

### Error handling (`error.rs`)

| Function | Description |
|----------|-------------|
| `wasm_ok<T: Serialize>(val) -> String` | Wrap success as `{"ok": value}` |
| `wasm_err(msg) -> String` | Wrap error as `{"error": message}` |
| `wasm_result<T: Serialize>(r) -> String` | Wrap `Result<T>` as JSON |

### Sequence module (`seq.rs`)

| Function | Description |
|----------|-------------|
| `parse_fasta(data: &[u8]) -> String` | Parse FASTA from bytes, return JSON stats |
| `gc_content_json(seq: &str) -> String` | Compute GC content, return JSON |
| `parse_fasta_bytes(data) -> Result<FastaStats>` | Direct Rust-level FASTA parsing |
| `gc_content(seq: &str) -> f64` | Direct Rust-level GC content |

### Alignment module (`align.rs`)

| Function | Description |
|----------|-------------|
| `align_dna(query, target, mode) -> String` | DNA alignment with default scoring, JSON result |
| `align_dna_custom(query, target, mode, match, mismatch, gap_open, gap_extend) -> String` | Custom scoring |
| `align_protein(query, target, mode, matrix) -> String` | Protein alignment, JSON result |

### Statistics module (`stats.rs`)

| Function | Description |
|----------|-------------|
| `describe(data_json: &str) -> String` | Descriptive statistics from JSON array |
| `pearson(x_json, y_json) -> String` | Pearson correlation from JSON arrays |
| `t_test(data_json, mu) -> String` | One-sample t-test from JSON array |

### ML module (`ml.rs`)

| Function | Description |
|----------|-------------|
| `kmer_count(seq, k) -> String` | K-mer counting, JSON result |
| `euclidean_distance(a_json, b_json) -> String` | Euclidean distance from JSON arrays |
| `manhattan_distance(a_json, b_json) -> String` | Manhattan distance |
| `hamming_distance(a, b) -> String` | Byte-level Hamming distance |
| `cosine_similarity(a_json, b_json) -> String` | Cosine similarity |

### Constants

| Constant | Description |
|----------|-------------|
| `VERSION` | Crate version from Cargo.toml |

## Design Decisions

- **JSON-based interface**: All functions accept/return JSON strings for maximum JavaScript interop compatibility. Complex types are serialized via serde.
- **No wasm-bindgen yet**: The crate compiles to a library. A thin wasm-bindgen layer would expose these functions as `#[wasm_bindgen]` exports.
- **Stateless**: All functions are pure -- no global state or initialization needed.

## Feature Flags

None.

## Dependencies

- `cyanea-core`, `cyanea-seq`, `cyanea-io`, `cyanea-align`, `cyanea-stats`, `cyanea-ml`
- `serde`, `serde_json` -- JSON serialization

## Tests

35 tests across 6 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 72 | Module declarations, VERSION constant |
| `error.rs` | 74 | JSON error wrapping |
| `seq.rs` | 162 | FASTA parsing, GC content |
| `align.rs` | 156 | DNA/protein alignment |
| `stats.rs` | 178 | Descriptive statistics, correlation, t-test |
| `ml.rs` | 167 | K-mer counting, distance metrics |
