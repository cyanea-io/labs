# cyanea-py

> Python bindings for the Cyanea bioinformatics ecosystem via PyO3. Installable as `pip install cyanea`.

## What's Inside

Python classes and functions wrapping 9 domain crates:

- **seq** -- `DnaSequence`, `RnaSequence`, `ProteinSequence` classes with FASTA/FASTQ parsing
- **align** -- `align_dna()`, `align_protein()`, `align_batch()` with keyword-only scoring parameters
- **stats** -- `describe()`, `pearson()`, `spearman()`, t-tests, Mann-Whitney U, p-value correction
- **ml** -- PCA, t-SNE, UMAP, K-means, distance metrics, pairwise distances (with NumPy variants)
- **chem** -- `Molecule` class from SMILES, properties, fingerprints, Tanimoto similarity
- **struct_bio** -- `Structure` class from PDB, secondary structure, RMSD, Kabsch alignment
- **phylo** -- `PhyloTree` from Newick, evolutionary distances, UPGMA/NJ, Robinson-Foulds
- **io** -- CSV, VCF, BED, GFF3, SAM, BAM file parsing with summary statistics
- **core** -- SHA-256, zstd compression/decompression

## Quick Start

### Build

```bash
pip install maturin
cd cyanea-py
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release
```

### Usage

```python
from cyanea import seq, align, stats, ml, chem

# Sequence analysis
dna = seq.DnaSequence(b"ATCGATCGATCG")
print(f"GC: {dna.gc_content():.1%}")
print(f"RevComp: {dna.reverse_complement()}")

# Alignment
result = align.align_dna("ACGTACGT", "ACGTGCGT", mode="local")
print(f"Score: {result.score}, Identity: {result.identity:.1%}")

# Statistics
summary = stats.describe([1.0, 2.0, 3.0, 4.0, 5.0])
print(f"Mean: {summary.mean}, Median: {summary.median}")

# Machine learning
pca = ml.pca([1,2,3,4,5,6,7,8,9,10,11,12], n_features=4, n_components=2)
print(f"Variance explained: {pca.explained_variance_ratio}")

# Chemistry
mol = chem.parse_smiles("c1ccccc1O")
print(f"Formula: {mol.molecular_formula()}, MW: {mol.molecular_weight():.1f}")
```

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `numpy` | No | NumPy array interop (`pca_np`, `umap_np`, `tsne_np`, `pairwise_distances_np`) |
| `h5ad` | No | HDF5 `.h5ad` I/O via cyanea-omics |

## Error Mapping

| Rust `CyaneaError` | Python Exception |
|---------------------|-----------------|
| `Io` | `IOError` |
| `Parse`, `InvalidInput` | `ValueError` |
| `Compression`, `Hash`, `Other` | `RuntimeError` |

## Python Modules

| Module | Description |
|--------|-------------|
| `cyanea.seq` | Sequence types, FASTA/FASTQ parsing |
| `cyanea.align` | DNA/protein alignment |
| `cyanea.stats` | Descriptive stats, hypothesis tests, correlation |
| `cyanea.ml` | PCA, t-SNE, UMAP, K-means, distances |
| `cyanea.chem` | SMILES, molecular properties, fingerprints |
| `cyanea.struct_bio` | PDB parsing, secondary structure, RMSD |
| `cyanea.phylo` | Newick trees, distances, tree construction |
| `cyanea.io` | CSV, VCF, BED, GFF3, SAM, BAM parsing |
| `cyanea.core` | SHA-256, zstd compression |

## See Also

- [API Reference (STATUS.md)](docs/STATUS.md)
- [Architecture](../ARCHITECTURE.md)
- [Build Guide](../docs/BUILDING.md)
