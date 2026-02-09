# cyanea-wasm

WebAssembly bindings for browser-based bioinformatics. Wraps cyanea-seq, cyanea-align, cyanea-stats, cyanea-ml, and cyanea-core into a JSON-based interface suitable for JavaScript/TypeScript consumption.

## Status: Complete

All bindings are implemented with full `wasm-bindgen` annotations (behind the `wasm` feature flag). Functions accept simple types (strings, numbers) and return JSON strings. Includes sequence manipulation, alignment (all modes + batch), statistics, ML distance metrics, SHA-256 hashing, and zstd compression.

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
| `parse_fasta(data) -> String` | Parse FASTA from string, return JSON stats |
| `gc_content_json(seq) -> String` | Compute GC content, return JSON |
| `reverse_complement(seq) -> String` | DNA reverse complement |
| `transcribe(seq) -> String` | DNA to RNA transcription |
| `translate(seq) -> String` | DNA to protein (standard codon table) |
| `validate(seq, alphabet) -> String` | Validate against "dna"/"rna"/"protein" |
| `parse_fastq(data) -> String` | Parse FASTQ records from string |
| `parse_fasta_bytes(data) -> Result<FastaStats>` | Direct Rust-level FASTA parsing |
| `gc_content(seq) -> f64` | Direct Rust-level GC content |

### Alignment module (`align.rs`)

| Function | Description |
|----------|-------------|
| `align_dna(query, target, mode) -> String` | DNA alignment with default scoring |
| `align_dna_custom(query, target, mode, match, mismatch, gap_open, gap_extend) -> String` | Custom scoring |
| `align_protein(query, target, mode, matrix) -> String` | Protein alignment |
| `align_batch(pairs_json, mode, match, mismatch, gap_open, gap_extend) -> String` | Batch alignment from JSON pairs |

### Statistics module (`stats.rs`)

| Function | Description |
|----------|-------------|
| `describe(data_json) -> String` | Descriptive statistics from JSON array |
| `pearson(x_json, y_json) -> String` | Pearson correlation |
| `spearman(x_json, y_json) -> String` | Spearman rank correlation |
| `t_test(data_json, mu) -> String` | One-sample t-test |
| `t_test_two_sample(x_json, y_json, equal_var) -> String` | Two-sample t-test (Student's or Welch's) |
| `mann_whitney_u(x_json, y_json) -> String` | Mann-Whitney U test |
| `bonferroni(p_json) -> String` | Bonferroni p-value correction |
| `benjamini_hochberg(p_json) -> String` | Benjamini-Hochberg FDR correction |

### ML module (`ml.rs`)

| Function | Description |
|----------|-------------|
| `kmer_count(seq, k) -> String` | K-mer counting |
| `euclidean_distance(a_json, b_json) -> String` | Euclidean distance |
| `manhattan_distance(a_json, b_json) -> String` | Manhattan distance |
| `hamming_distance(a, b) -> String` | Byte-level Hamming distance |
| `cosine_similarity(a_json, b_json) -> String` | Cosine similarity |

### Core utilities module (`core_utils.rs`)

| Function | Description |
|----------|-------------|
| `sha256(data) -> String` | SHA-256 hex digest |
| `zstd_compress(data, level) -> String` | Zstd compression, returns JSON byte array |
| `zstd_decompress(data_json) -> String` | Zstd decompression from JSON byte array |

### Constants

| Constant | Description |
|----------|-------------|
| `VERSION` | Crate version from Cargo.toml |

## Design Decisions

- **JSON-based interface**: All functions accept/return JSON strings for maximum JavaScript interop compatibility. Complex types are serialized via serde.
- **wasm-bindgen ready**: All public functions have `#[cfg_attr(feature = "wasm", wasm_bindgen)]` annotations. Build with `--features wasm` to produce bindgen-compatible exports.
- **Stateless**: All functions are pure -- no global state or initialization needed.

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `wasm` | No | Enables `wasm-bindgen` annotations on all public functions |

## Dependencies

- `cyanea-core` (with `std` for sha256, zstd), `cyanea-seq`, `cyanea-io`, `cyanea-align`, `cyanea-stats`, `cyanea-ml`
- `serde`, `serde_json` -- JSON serialization
- `wasm-bindgen` (optional, behind `wasm` feature)

## Tests

55 unit tests + 1 doc test across 7 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 90 | Module declarations, re-exports, VERSION constant |
| `error.rs` | 74 | JSON error wrapping |
| `seq.rs` | 358 | FASTA/FASTQ parsing, GC content, reverse complement, transcribe, translate, validate |
| `align.rs` | 232 | DNA/protein alignment, batch alignment |
| `stats.rs` | 305 | Descriptive stats, correlation, hypothesis testing, p-value correction |
| `ml.rs` | 175 | K-mer counting, distance metrics |
| `core_utils.rs` | 67 | SHA-256 hashing, zstd compression/decompression |
