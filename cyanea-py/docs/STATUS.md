# cyanea-py

Python bindings for the Cyanea bioinformatics ecosystem via PyO3. Installable as `pip install cyanea` (via maturin).

## Status: Complete

All bindings are implemented and tested. Wraps cyanea-seq, cyanea-align, cyanea-stats, cyanea-core, cyanea-ml, cyanea-chem, cyanea-struct, cyanea-phylo, and cyanea-io.

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

### `cyanea.ml`

**Classes:**

| Class | Description |
|-------|-------------|
| `UmapResult` | UMAP result: `embedding`, `n_samples`, `n_components`, `n_epochs` |
| `PcaResult` | PCA result (frozen): `components`, `explained_variance`, `explained_variance_ratio`, `transformed`, `mean`, `n_features`, `n_components` |
| `TsneResult` | t-SNE result (frozen): `embedding`, `n_samples`, `n_components`, `kl_divergence` |
| `KMeansResult` | K-Means result (frozen): `centroids`, `labels`, `inertia`, `n_iter`, `n_features` |

**Functions:**

| Function | Description |
|----------|-------------|
| `euclidean_distance(a, b) -> float` | Euclidean distance between two vectors |
| `manhattan_distance(a, b) -> float` | Manhattan distance between two vectors |
| `hamming_distance(a, b) -> int` | Hamming distance between two byte sequences |
| `cosine_similarity(a, b) -> float` | Cosine similarity between two vectors |
| `pairwise_distances(data, n_features, metric="euclidean") -> list[float]` | Pairwise distance matrix (flat n*n) |
| `umap(data, n_features, n_components=2, n_neighbors=15, min_dist=0.1, n_epochs=200, metric="euclidean", seed=42) -> UmapResult` | UMAP dimensionality reduction |
| `pca(data, n_features, n_components=2) -> PcaResult` | PCA dimensionality reduction |
| `tsne(data, n_features, n_components=2, perplexity=30.0, learning_rate=200.0, n_iter=1000, seed=42) -> TsneResult` | t-SNE dimensionality reduction |
| `kmeans(data, n_features, n_clusters=2, max_iter=300, tolerance=1e-4, seed=42) -> KMeansResult` | K-Means clustering |

**NumPy variants (feature `numpy`):**

| Function | Description |
|----------|-------------|
| `pairwise_distances_np(data, n_features, metric="euclidean") -> np.ndarray` | Pairwise distances as NumPy 2D array |
| `umap_np(data, n_features, ...) -> np.ndarray` | UMAP embedding as NumPy 2D array (n_samples x n_components) |
| `pca_np(data, n_features, n_components=2) -> np.ndarray` | PCA transformed data as NumPy 2D array |
| `tsne_np(data, n_features, ...) -> np.ndarray` | t-SNE embedding as NumPy 2D array |

### `cyanea.chem`

**Classes:**

| Class | Description |
|-------|-------------|
| `Molecule(smiles: str)` | Molecular graph parsed from SMILES |
| `MolecularProperties` | Computed properties (frozen): `formula`, `molecular_weight`, `atom_count`, `bond_count`, `ring_count`, `rotatable_bonds`, `hbd`, `hba` |

**Molecule methods:** `atom_count() -> int`, `bond_count() -> int`, `molecular_formula() -> str`, `molecular_weight() -> float`, `canonical_smiles() -> str`, `morgan_fingerprint(radius=2, n_bits=2048) -> list[int]`, `has_substructure(pattern: str) -> bool`, `__len__`, `__str__`, `__repr__`

**Functions:**

| Function | Description |
|----------|-------------|
| `parse_smiles(smiles: str) -> Molecule` | Parse a SMILES string into a Molecule |
| `molecular_properties(smiles: str) -> MolecularProperties` | Compute molecular properties from SMILES |
| `tanimoto(smiles1, smiles2, radius=2, n_bits=2048) -> float` | Tanimoto similarity between two molecules |
| `canonical_smiles(smiles: str) -> str` | Return canonical SMILES for a given SMILES string |

### `cyanea.struct_bio`

**Classes:**

| Class | Description |
|-------|-------------|
| `Structure` | Macromolecular 3D structure parsed from PDB (frozen) |
| `SecondaryStructureAssignment` | Per-residue assignment (frozen): `residue_num`, `residue_name`, `structure` |

**Structure methods:** `chain_count() -> int`, `residue_count() -> int`, `atom_count() -> int`, `id() -> str`, `secondary_structure() -> list[SecondaryStructureAssignment]`, `__str__`, `__repr__`

**Functions:**

| Function | Description |
|----------|-------------|
| `parse_pdb(content: str) -> Structure` | Parse PDB-format text into a Structure |
| `rmsd(coords1, coords2) -> float` | RMSD between two equal-length coordinate sets (list of [x,y,z]) |
| `kabsch_align(moving, target) -> (float, list)` | Kabsch optimal superposition; returns (rmsd, aligned_coords) |

### `cyanea.phylo`

**Classes:**

| Class | Description |
|-------|-------------|
| `PhyloTree` | Rooted phylogenetic tree (frozen) |

**PhyloTree methods:** `leaf_count() -> int`, `leaf_names() -> list[str]`, `to_newick() -> str`, `robinson_foulds(other: PhyloTree) -> int`, `__len__`, `__repr__`

**Functions:**

| Function | Description |
|----------|-------------|
| `parse_newick(newick: str) -> PhyloTree` | Parse a Newick format string into a PhyloTree |
| `evolutionary_distance(seq1, seq2, model="jc") -> float` | Evolutionary distance between two aligned sequences (models: `"p"`, `"jc"`, `"k2p"`) |
| `upgma(labels, matrix) -> PhyloTree` | Build a UPGMA tree from labels and a square distance matrix |
| `neighbor_joining(labels, matrix) -> PhyloTree` | Build a Neighbor-Joining tree from labels and a square distance matrix |

### `cyanea.io`

**Classes:**

| Class | Description |
|-------|-------------|
| `CsvInfo` | CSV metadata (frozen): `row_count`, `column_count`, `columns`, `has_headers` |
| `VcfStats` | VCF summary (frozen): `variant_count`, `snv_count`, `indel_count`, `pass_count`, `chromosomes` |
| `BedStats` | BED summary (frozen): `record_count`, `total_bases`, `chromosomes` |
| `Gff3Stats` | GFF3 summary (frozen): `gene_count`, `transcript_count`, `exon_count`, `protein_coding_count`, `chromosomes` |
| `PySamRecord` | SAM/BAM alignment record (frozen): `qname`, `flag`, `rname`, `pos`, `mapq`, `cigar`, `sequence`, `quality` |
| `PySamStats` | SAM/BAM summary (frozen): `total_reads`, `mapped`, `unmapped`, `avg_mapq`, `avg_length`, `mapq_distribution` |

**Functions:**

| Function | Description |
|----------|-------------|
| `csv_info(path: str) -> CsvInfo` | Parse a CSV file and return metadata |
| `vcf_stats(path: str) -> VcfStats` | Parse a VCF file and return summary statistics |
| `bed_stats(path: str) -> BedStats` | Parse a BED file and return summary statistics |
| `gff3_stats(path: str) -> Gff3Stats` | Parse a GFF3 file and return summary statistics |
| `parse_sam(path: str) -> list[PySamRecord]` | Parse a SAM file and return all alignment records |
| `sam_stats(path: str) -> PySamStats` | Compute summary statistics for a SAM file |
| `parse_bam(path: str) -> list[PySamRecord]` | Parse a BAM file and return all alignment records |
| `bam_stats(path: str) -> PySamStats` | Compute summary statistics for a BAM file |

## Dependencies

- `pyo3` 0.23 (extension-module)
- `numpy` 0.23 (optional, feature-gated)
- `cyanea-core`, `cyanea-seq`, `cyanea-io`, `cyanea-align`, `cyanea-stats`, `cyanea-ml`, `cyanea-chem`, `cyanea-struct`, `cyanea-phylo`

## Tests

No Rust-level tests. Integration tested via Python smoke tests.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 60 | PyO3 module root, submodule registration in `sys.modules` |
| `error.rs` | 28 | `CyaneaError` -> `PyErr` mapping, `IntoPyResult` trait |
| `seq.rs` | 327 | Sequence types and FASTA/FASTQ parsing |
| `align.rs` | 153 | Alignment functions with keyword-only scoring params |
| `stats.rs` | 153 | Statistics and hypothesis testing |
| `core_utils.rs` | 44 | SHA-256 hashing and zstd compression |
| `ml.rs` | 358 | Distance metrics, PCA, t-SNE, UMAP, K-Means, NumPy variants |
| `chem.rs` | 154 | SMILES parsing, molecular properties, fingerprints, Tanimoto similarity |
| `struct_bio.rs` | 159 | PDB parsing, secondary structure, RMSD, Kabsch alignment |
| `phylo.rs` | 179 | Newick parsing, evolutionary distances, UPGMA, Neighbor-Joining |
| `io.rs` | 203 | CSV, VCF, BED, GFF3, SAM, BAM file format parsing |
