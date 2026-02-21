# cyanea-py Usage Guide

## Installation

```bash
# From source (recommended during development)
pip install maturin
cd cyanea-py
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release

# With NumPy support
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release --features numpy
```

## Quick Start

```python
from cyanea import seq, align, stats

# Sequence analysis
dna = seq.DnaSequence(b"ATCGATCG")
print(f"GC: {dna.gc_content():.1%}")       # 50.0%
print(f"RevComp: {dna.reverse_complement()}")

# Alignment
result = align.align_dna("ACGT", "ACGTT", mode="local")
print(f"Score: {result.score}, CIGAR: {result.cigar_string}")

# Statistics
summary = stats.describe([1.0, 2.0, 3.0, 4.0, 5.0])
print(f"Mean: {summary.mean}, Std: {summary.std_dev:.2f}")
```

## Sequence Analysis

```python
from cyanea import seq

# DNA manipulation
dna = seq.DnaSequence(b"ATCGATCGATCG")
rna = dna.transcribe()          # RnaSequence
protein = dna.translate()       # ProteinSequence
kmers = dna.kmers(3)           # List of 3-mers

# File parsing
fasta = seq.fasta_stats("sequences.fa")
print(f"Sequences: {fasta.count}, Total bases: {fasta.total_length}")

records = seq.parse_fastq("reads.fq")
print(f"Read count: {len(records)}")

# RNA folding
fold = seq.rna_fold_nussinov("GCGCAAUAGCGC")
print(f"Structure: {fold['structure']}, Pairs: {fold['num_pairs']}")

# Protein properties
props = seq.protein_properties("MKWVTFISLLLLFSSAYS")
print(f"MW: {props['molecular_weight']:.1f}, pI: {props['isoelectric_point']:.1f}")
```

## Alignment

```python
from cyanea import align

# Local alignment
result = align.align_dna("ACGTACGT", "CGTAC", mode="local", match_score=2, mismatch_score=-1)
print(f"Score: {result.score}, Identity: {result.identity:.1%}")

# Protein alignment with custom matrix
result = align.align_protein("MKWVT", "MKWIT", mode="global", matrix="blosum62")

# Batch alignment
results = align.align_batch([("ACGT", "ACGTT"), ("GGCC", "GGGCC")], mode="local")
```

## Statistics

```python
from cyanea import stats

# Hypothesis testing
t = stats.t_test_two_sample([1,2,3,4,5], [3,4,5,6,7])
print(f"T={t.statistic:.2f}, p={t.p_value:.4f}")

# Correlation
r = stats.pearson([1,2,3,4,5], [2,4,6,8,10])
print(f"Pearson r: {r}")

# Survival analysis
km = stats.kaplan_meier([1,2,3,5,7,10], [1,1,0,1,0,1])

# Population genetics
fst = stats.fst_hudson([0,1,2,0,1], [2,2,1,2,1])
td = stats.tajimas_d([0,1,0,0,1,2,0,1])

# Enrichment
gsea = stats.gsea_preranked(ranked_genes, gene_set)
```

## Machine Learning

```python
from cyanea import ml

# Dimensionality reduction
pca_result = ml.pca([1,2,3,4,5,6,7,8,9,10,11,12], n_features=4, n_components=2)
print(f"Explained variance: {pca_result.explained_variance_ratio}")

umap_result = ml.umap(data, n_features=100, n_components=2, n_neighbors=15)
tsne_result = ml.tsne(data, n_features=100, n_components=2, perplexity=30.0)

# Clustering
km = ml.kmeans(data, n_features=4, n_clusters=3)
print(f"Labels: {km.labels}, Inertia: {km.inertia:.2f}")

# Classification
predictions = ml.random_forest_classify(train, labels, test, n_features=10, n_trees=100)
cm = ml.confusion_matrix(actual, predicted)
roc = ml.roc_curve(scores, labels)
```

## NumPy Integration

```python
import numpy as np
from cyanea import ml

# NumPy-native variants return np.ndarray directly
coords = ml.pca_np(data, n_features=100, n_components=2)  # np.ndarray (n, 2)
distances = ml.pairwise_distances_np(data, n_features=50)  # np.ndarray (n, n)
embedding = ml.umap_np(data, n_features=100, n_components=2)
```

## Chemistry

```python
from cyanea import chem

mol = chem.parse_smiles("c1ccccc1O")
print(f"Formula: {mol.molecular_formula()}")
print(f"MW: {mol.molecular_weight():.1f}")
print(f"Canonical: {mol.canonical_smiles()}")

fp = mol.morgan_fingerprint(radius=2, n_bits=2048)
print(f"Fingerprint bits: {len(fp)}")

sim = chem.tanimoto("c1ccccc1", "c1ccccc1O")
print(f"Tanimoto: {sim:.3f}")
```

## Structural Biology

```python
from cyanea import struct_bio

structure = struct_bio.parse_pdb(pdb_text)
print(f"Chains: {structure.chain_count()}, Residues: {structure.residue_count()}")

ss = structure.secondary_structure()
for assignment in ss:
    print(f"{assignment.residue_num} {assignment.residue_name}: {assignment.structure}")

rmsd = struct_bio.rmsd(coords1, coords2)
rmsd_aligned, aligned = struct_bio.kabsch_align(moving, target)
```

## Phylogenetics

```python
from cyanea import phylo

tree = phylo.parse_newick("((A:0.1,B:0.2):0.3,C:0.4);")
print(f"Leaves: {tree.leaf_names()}")

dist = phylo.evolutionary_distance("ACGTACGT", "ACGTGCGT", model="jc")
nj_tree = phylo.neighbor_joining(["A","B","C"], [[0,1,2],[1,0,1.5],[2,1.5,0]])

# Simulation
sim = phylo.simulate_coalescent(10, 1000, seed=42)
```

## File I/O

```python
from cyanea import io

# VCF
vcf = io.vcf_stats("variants.vcf")
print(f"Variants: {vcf.variant_count}, SNVs: {vcf.snv_count}")

# SAM/BAM
records = io.parse_sam("alignments.sam")
bam = io.bam_stats("alignments.bam")

# BED/GFF3
bed = io.bed_stats("regions.bed")
gff = io.gff3_stats("annotations.gff3")
```

## Omics

```python
from cyanea import omics

merged = omics.merge_intervals([{"chrom": "chr1", "start": 100, "end": 200},
                                 {"chrom": "chr1", "start": 150, "end": 300}])
vep = omics.annotate_variant(variant, transcripts)
```

## Single-Cell

```python
from cyanea import sc

# Preprocessing
normalized = sc.normalize_total(matrix, target_sum=10000)
log_data = sc.log1p_transform(normalized)
hvg = sc.highly_variable_genes(log_data, n_top=2000)

# Clustering
graph = sc.neighbors(pca_data, n_neighbors=15)
clusters = sc.leiden(graph["adjacency"], resolution=1.0)

# Trajectory
dm = sc.diffusion_map(pca_data, n_components=10)
pseudotime = sc.dpt(dm, root=0)

# Markers
markers = sc.rank_genes(matrix, clusters, method="t-test")

# Integration
corrected = sc.harmony(pca_data, batch_labels, n_components=30)
```

## Error Handling

```python
from cyanea import seq

try:
    dna = seq.DnaSequence(b"ATCGXYZ")
except ValueError as e:
    print(f"Invalid sequence: {e}")

try:
    stats = seq.fasta_stats("/nonexistent.fa")
except IOError as e:
    print(f"File error: {e}")
```
