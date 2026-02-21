# cyanea-py

> Python bindings for the Cyanea bioinformatics ecosystem via PyO3. Installable as `pip install cyanea`.

## What's Inside

Python classes and functions wrapping 11 submodules across all domain crates:

- **seq** -- `DnaSequence`, `RnaSequence`, `ProteinSequence` classes, FASTA/FASTQ, RNA folding, protein properties, read simulation, codon usage, assembly stats
- **align** -- `align_dna()`, `align_protein()`, `align_batch()` with keyword-only scoring parameters
- **stats** -- descriptive stats, correlation, hypothesis tests, survival (Kaplan-Meier, Cox PH), population genetics (Fst, Tajima's D), enrichment (GSEA, ORA), diversity, null models
- **ml** -- PCA, t-SNE, UMAP, K-means, random forest, GBDT, distance metrics, HMM, confusion matrix, ROC/PR curves, cross-validation, feature selection (with NumPy variants)
- **chem** -- `Molecule` class from SMILES, properties, fingerprints, Tanimoto similarity
- **struct_bio** -- `Structure` class from PDB, secondary structure, RMSD, Kabsch alignment
- **phylo** -- `PhyloTree` from Newick, evolutionary distances, UPGMA/NJ, Robinson-Foulds, simulation
- **io** -- CSV, VCF, BED, GFF3, SAM, BAM, BLAST XML, bedGraph, GFA, NCBI/UniProt URL builders
- **omics** -- interval operations, coordinate liftover, variant annotation, bisulfite conversion, spatial autocorrelation
- **sc** -- single-cell pipeline: normalize, HVG, neighbors, Leiden/Louvain, diffusion map, DPT, PAGA, markers, Harmony/ComBat/MNN
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
| `cyanea.seq` | Sequence types, FASTA/FASTQ, RNA folding, protein properties, simulation |
| `cyanea.align` | DNA/protein alignment, batch alignment |
| `cyanea.stats` | Statistics, hypothesis tests, survival, popgen, enrichment, diversity |
| `cyanea.ml` | PCA, t-SNE, UMAP, K-means, forests, GBDT, HMM, metrics, CV |
| `cyanea.chem` | SMILES, molecular properties, fingerprints |
| `cyanea.struct_bio` | PDB parsing, secondary structure, RMSD, Kabsch |
| `cyanea.phylo` | Newick trees, distances, construction, simulation |
| `cyanea.io` | CSV, VCF, BED, GFF3, SAM, BAM, BLAST XML, bedGraph, GFA |
| `cyanea.omics` | Interval operations, liftover, variant annotation, methylation, spatial |
| `cyanea.sc` | Single-cell: normalize, HVG, clustering, trajectory, markers, integration |
| `cyanea.core` | SHA-256, zstd compression |

## See Also

- [API Reference](docs/API.md)
- [Usage Guide](docs/GUIDE.md)
- [Internal Architecture](docs/ARCHITECTURE.md)
- [Workspace Architecture](../docs/ARCHITECTURE.md)
- [Build Guide](../docs/BUILDING.md)
- [Bindings Guide](../docs/BINDINGS.md)
