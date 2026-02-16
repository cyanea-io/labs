# cyanea-wasm

WebAssembly bindings for browser-based bioinformatics. Wraps cyanea-seq, cyanea-align, cyanea-stats, cyanea-ml, cyanea-chem, cyanea-struct, cyanea-phylo, and cyanea-core into a JSON-based interface suitable for JavaScript/TypeScript consumption.

## Status: Complete

All bindings are implemented with full `wasm-bindgen` annotations (behind the `wasm` feature flag). Functions accept simple types (strings, numbers) and return JSON strings. Includes sequence manipulation, paired-end FASTQ parsing, read trimming, alignment (all modes + batch), statistics, ML distance metrics and UMAP, chemistry (SMILES parsing, fingerprints, similarity, substructure), structural biology (PDB parsing, secondary structure, RMSD), phylogenetics (Newick trees, evolutionary distances, UPGMA/NJ, Robinson-Foulds), SHA-256 hashing, and zstd compression.

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
| `parse_paired_fastq(r1_data, r2_data, validation) -> String` | Parse paired FASTQ from two strings (`validation`: "strict"/"relaxed"/"none") |
| `parse_interleaved_fastq(data, validation) -> String` | Parse interleaved FASTQ (alternating R1/R2) |
| `trim_fastq(data, config_json) -> String` | Trim single-end FASTQ records |
| `trim_paired_fastq(r1_data, r2_data, config_json, orphan_policy) -> String` | Trim paired FASTQ (`orphan_policy`: "drop_both"/"keep_first"/"keep_second") |
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

| Type | Description |
|------|-------------|
| `JsKmerCounts` | Serializable k-mer counts (`k`, `total`, `distinct`, `counts` map) |
| `JsUmapResult` | Serializable UMAP result (`embedding`, `n_samples`, `n_components`, `n_epochs`) |

| Function | Description |
|----------|-------------|
| `kmer_count(seq, k) -> String` | K-mer counting |
| `euclidean_distance(a_json, b_json) -> String` | Euclidean distance |
| `manhattan_distance(a_json, b_json) -> String` | Manhattan distance |
| `hamming_distance(a, b) -> String` | Byte-level Hamming distance |
| `cosine_similarity(a_json, b_json) -> String` | Cosine similarity |
| `umap(data_json, n_features, n_components, n_neighbors, min_dist, n_epochs, metric) -> String` | UMAP dimensionality reduction |

### Chemistry module (`chem.rs`)

| Type | Description |
|------|-------------|
| `JsMolecularProperties` | Serializable molecular properties (`formula`, `molecular_weight`, `atom_count`, `bond_count`, `ring_count`, `rotatable_bonds`, `hbd`, `hba`) |
| `JsFingerprint` | Serializable fingerprint summary (`bits`, `n_bits`, `n_on_bits`) |
| `JsSubstructureResult` | Serializable substructure result (`has_match`, `match_count`) |

| Function | Description |
|----------|-------------|
| `smiles_properties(smiles) -> String` | Parse SMILES and return molecular properties as JSON |
| `canonical(smiles) -> String` | Generate canonical SMILES from input SMILES string |
| `smiles_fingerprint(smiles, radius, n_bits) -> String` | Compute Morgan fingerprint and return set bits as JSON |
| `tanimoto(smiles1, smiles2) -> String` | Tanimoto similarity between two SMILES strings |
| `smiles_substructure(molecule, pattern) -> String` | Check for substructure match between molecule and pattern SMILES |

### Structural biology module (`struct_bio.rs`)

| Type | Description |
|------|-------------|
| `JsStructureInfo` | Serializable structure info (`id`, `chain_count`, `residue_count`, `atom_count`, `chains`) |
| `JsChainInfo` | Serializable chain info (`id`, `residue_count`, `atom_count`) |
| `JsSecondaryStructure` | Serializable secondary structure assignment (`chain_id`, `assignments`) |
| `JsSSAssignment` | Serializable per-residue assignment (`residue_num`, `residue_name`, `structure`) |

| Function | Description |
|----------|-------------|
| `pdb_info(pdb_text) -> String` | Parse PDB text and return structure info as JSON |
| `pdb_secondary_structure(pdb_text) -> String` | Assign secondary structure from PDB text and return JSON |
| `rmsd(coords1_json, coords2_json) -> String` | Compute RMSD between two coordinate sets (JSON arrays of [x,y,z]) |

### Phylogenetics module (`phylo.rs`)

| Type | Description |
|------|-------------|
| `JsTreeInfo` | Serializable tree info (`leaf_count`, `internal_count`, `total_nodes`, `leaf_names`, `newick`) |
| `JsRFDistance` | Serializable Robinson-Foulds distance (`distance`, `normalized`) |

| Function | Description |
|----------|-------------|
| `newick_info(newick) -> String` | Parse a Newick string and return tree info as JSON |
| `evolutionary_distance(seq1, seq2, model) -> String` | Compute evolutionary distance (`model`: "p", "jc", or "k2p") |
| `build_upgma(labels_json, matrix_json) -> String` | Build a UPGMA tree from a JSON distance matrix |
| `build_nj(labels_json, matrix_json) -> String` | Build a Neighbor-Joining tree from a JSON distance matrix |
| `rf_distance(newick1, newick2) -> String` | Robinson-Foulds distance between two Newick trees |

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

- `cyanea-core` (with `std` for sha256, zstd), `cyanea-seq`, `cyanea-io`, `cyanea-align`, `cyanea-stats`, `cyanea-ml`, `cyanea-chem`, `cyanea-struct`, `cyanea-phylo`
- `serde`, `serde_json` -- JSON serialization
- `wasm-bindgen` (optional, behind `wasm` feature)

## Tests

111 tests across 10 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 131 | Module declarations, re-exports, VERSION constant |
| `error.rs` | 74 | JSON error wrapping |
| `seq.rs` | 892 | FASTA/FASTQ parsing, paired-end FASTQ, trimming, GC content, reverse complement, transcribe, translate, validate |
| `align.rs` | 232 | DNA/protein alignment, batch alignment |
| `stats.rs` | 305 | Descriptive stats, correlation, hypothesis testing, p-value correction |
| `ml.rs` | 247 | K-mer counting, distance metrics, UMAP dimensionality reduction |
| `chem.rs` | 214 | SMILES parsing, molecular properties, fingerprints, Tanimoto similarity, substructure search |
| `struct_bio.rs` | 235 | PDB parsing, secondary structure assignment, RMSD computation |
| `phylo.rs` | 339 | Newick parsing, evolutionary distances, UPGMA/NJ tree construction, Robinson-Foulds distance |
| `core_utils.rs` | 67 | SHA-256 hashing, zstd compression/decompression |
