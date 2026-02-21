# Bindings Guide

Cyanea provides bindings for three non-Rust targets: Python (via PyO3), WebAssembly (via wasm-bindgen), and Elixir NIFs (via Rustler). All bindings wrap the same Rust domain crates, so behavior is identical across targets.

## Python Bindings (`cyanea-py`)

The Python package exposes 11 submodules. Install with maturin:

```bash
cd cyanea-py
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release
```

### `cyanea.seq`

Sequence types and analysis functions.

- **Types**: `DnaSequence`, `RnaSequence`, `ProteinSequence`
- **Functions**: `gc_content`, `reverse_complement`, `transcribe`, `translate`, `parse_fasta`, `parse_fastq`, `suffix_array`, `fm_index`, `kmers`, `minhash`, `trim`, `fold_rna`, `protein_properties`, `codon_usage`, `assembly_stats`, `simulate_reads`

```python
from cyanea import seq

dna = seq.DnaSequence(b"ATCGATCG")
print(dna.gc_content())          # 0.5
print(dna.reverse_complement())  # CGATCGAT

records = seq.parse_fasta(">s1\nACGT\n>s2\nTGCA")
kmers = seq.kmers("ACGTACGT", k=4)
```

### `cyanea.align`

Pairwise and multiple sequence alignment.

- **Functions**: `align_dna`, `align_protein`, `align_local`, `parse_cigar`, `cigar_stats`, `seed_extend`, `progressive_msa`, `poa`

```python
from cyanea import align

result = align.align_dna("ACGTACGT", "ACGACGT", mode="global")
print(result.score, result.cigar_string)

msa = align.progressive_msa(["ACGT", "ACGGT", "ACT"])
```

### `cyanea.stats`

Statistical analysis and hypothesis testing.

- **Descriptive**: `describe`, `pearson`, `spearman`
- **Tests**: `t_test`, `mann_whitney`, `fisher_exact`, `chi_square`, `ks_test`
- **Distributions**: Normal, Poisson, Binomial, Chi-squared, F, Beta, Gamma
- **Dimensionality reduction**: `pca`
- **Population genetics**: `fst`, `tajimas_d`, `hwe_test`, `ld`, `allele_frequencies`
- **Survival analysis**: `kaplan_meier`, `log_rank`, `cox_ph`
- **Enrichment**: `gsea`, `ora`
- **Diversity**: `shannon`, `simpson`, `chao1`
- **Null models**: `permutation_test`, `bootstrap`

```python
from cyanea import stats

summary = stats.describe([1.0, 2.0, 3.0, 4.0])
print(summary.mean, summary.std)

result = stats.t_test([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
print(result.statistic, result.p_value)

fst = stats.fst([[0.1, 0.9], [0.4, 0.6]])
```

### `cyanea.core`

Core utilities: hashing and compression.

- **Functions**: `sha256`, `compress`, `decompress`

```python
from cyanea import core

digest = core.sha256(b"ACGT")
compressed = core.compress(b"ACGTACGTACGT")
original = core.decompress(compressed)
```

### `cyanea.ml`

Machine learning algorithms.

- **Clustering**: `kmeans`, `dbscan`
- **Trees**: `random_forest`, `gbdt`
- **Neighbors**: `knn`
- **Reduction**: `pca`, `tsne`, `umap`
- **Distances**: `pairwise_distances`
- **Metrics**: `confusion_matrix`, `roc_curve`, `auc`
- **Validation**: `cross_validate`
- **Feature selection**: `feature_selection`
- **HMM**: `hmm`

```python
from cyanea import ml

labels = ml.kmeans(data, k=3)
embedding = ml.umap(data, n_components=2)
distances = ml.pairwise_distances(data, metric="euclidean")
```

### `cyanea.chem`

Cheminformatics: molecules, fingerprints, and properties.

- **Types**: `Molecule` (from SMILES)
- **Functions**: `properties`, `fingerprint`, `tanimoto`, `canonical_smiles`, `substructure_search`, `parse_sdf`

```python
from cyanea import chem

mol = chem.Molecule("CCO")
props = chem.properties("CCO")
print(props.molecular_weight, props.logp)

fp1 = chem.fingerprint("CCO", method="morgan", radius=2)
fp2 = chem.fingerprint("CCCO", method="morgan", radius=2)
sim = chem.tanimoto(fp1, fp2)
```

### `cyanea.struct_bio`

Structural biology: PDB/mmCIF analysis.

- **Types**: `Structure` (from PDB string)
- **Functions**: `secondary_structure`, `rmsd`, `kabsch`, `contact_map`, `ramachandran`, `parse_mmcif`

```python
from cyanea import struct_bio

structure = struct_bio.Structure(pdb_string)
ss = struct_bio.secondary_structure(structure)
contacts = struct_bio.contact_map(structure, cutoff=8.0)
```

### `cyanea.phylo`

Phylogenetic tree construction and analysis.

- **Types**: `PhyloTree` (from Newick string)
- **Functions**: `distances`, `upgma`, `nj`, `bootstrap`, `robinson_foulds`, `simulate_evolution`, `simulate_coalescent`

```python
from cyanea import phylo

tree = phylo.PhyloTree("((A:0.1,B:0.2):0.3,C:0.4);")
dist = phylo.distances(["ACGT", "ACGA", "TGCA"])
nj_tree = phylo.nj(dist)
support = phylo.bootstrap(["ACGT", "ACGA", "TGCA"], n=100)
```

### `cyanea.io`

File format parsing.

- **Functions**: `read_csv`, `read_vcf`, `read_bed`, `read_gff3`, `read_sam`, `read_bam`, `blast_xml`, `bedgraph`, `gfa`

```python
from cyanea import io

variants = io.read_vcf("variants.vcf")
regions = io.read_bed("regions.bed")
features = io.read_gff3("annotation.gff3")
records = io.read_sam("alignments.sam")
```

### `cyanea.omics`

Genomic data operations: intervals, annotation, and epigenomics.

- **Functions**: `merge_intervals`, `intersect_intervals`, `liftover`, `annotate_variant`, `bisulfite_convert`, `cpg_islands`, `morans_i`, `gearys_c`

```python
from cyanea import omics

merged = omics.merge_intervals([("chr1", 100, 200), ("chr1", 150, 300)])
hits = omics.intersect_intervals(regions_a, regions_b)
lifted = omics.liftover(intervals, chain)
```

### `cyanea.sc`

Single-cell analysis pipeline.

- **Preprocessing**: `normalize_total`, `log1p_transform`, `highly_variable_genes`
- **Neighbors and clustering**: `neighbors`, `leiden`, `louvain`
- **Trajectory**: `diffusion_map`, `dpt`, `paga`
- **Markers**: `rank_genes`
- **Integration**: `harmony`, `combat`, `mnn`

```python
from cyanea import sc

normalized = sc.normalize_total(counts, target_sum=10000.0)
log_data = sc.log1p_transform(normalized)
hvg = sc.highly_variable_genes(log_data, n_top=2000)
graph = sc.neighbors(pca_data, n_neighbors=15)
clusters = sc.leiden(graph, resolution=1.0)
dm = sc.diffusion_map(graph, n_components=10)
pseudotime = sc.dpt(dm, root_cell=0)
markers = sc.rank_genes(log_data, clusters)
corrected = sc.harmony(pca_data, batch_labels)
```

## WASM/TypeScript Bindings (`cyanea-wasm`)

Published as `@cyanea/bio` on npm. The WASM package wraps 10 domain modules with a JSON-based API.

### Build

```bash
wasm-pack build cyanea-wasm --features wasm
cd cyanea-wasm && npm run build:ts
```

### Modules

| Namespace | Domain |
|-----------|--------|
| `Seq` | Sequences, FASTA/FASTQ, k-mers, FM-index |
| `Align` | Pairwise and multiple alignment |
| `Stats` | Statistical tests, distributions |
| `ML` | Clustering, PCA, t-SNE, UMAP |
| `Chem` | SMILES, fingerprints, properties |
| `StructBio` | PDB/mmCIF, DSSP, Kabsch |
| `Phylo` | Newick trees, NJ, bootstrap |
| `IO` | VCF, BED, GFF3, SAM parsing |
| `Omics` | Intervals, variants, annotation |
| `Core` | SHA-256, compression |

### JSON Envelope

All WASM functions return JSON strings with a consistent envelope:

```json
{"ok": <value>}
{"error": "<message>"}
```

Helper functions in the Rust source:
- `wasm_ok(val)` -- serialize success to `{"ok": val}`
- `wasm_err(msg)` -- serialize error to `{"error": msg}`
- `wasm_result(r)` -- convert `Result<T, E>` to the appropriate envelope

### Usage

```javascript
import { parseNewick, alignDna, describe } from '@cyanea/bio';

const tree = parseNewick('((A:0.1,B:0.2):0.3,C:0.4);');
console.log(JSON.parse(tree));  // { ok: { ... } }

const alignment = alignDna('ACGT', 'ACGTT', 'local');
const stats = describe([1.0, 2.0, 3.0, 4.0]);
```

### TypeScript Types

TypeScript type definitions are in `cyanea-wasm/ts/types.ts`. The wrapper in `cyanea-wasm/ts/index.ts` provides typed functions that parse the JSON envelope and return native objects.

## Elixir NIF Bindings

The NIF crate lives in `../cyanea/native/cyanea_native/` and depends on Cyanea Labs crates via path dependencies.

### Build

```bash
cd ../cyanea
mix compile
```

### Check (Rust only, without BEAM)

```bash
cargo check -p cyanea-native
```

### Elixir Modules

The NIF bindings expose domain functionality through Elixir modules:

- `Cyanea.Seq` -- Sequence types and operations
- `Cyanea.Align` -- Pairwise alignment
- `Cyanea.Stats` -- Statistical functions
- `Cyanea.ML` -- Machine learning
- `Cyanea.Chem` -- Cheminformatics
- `Cyanea.Struct` -- Structural biology
- `Cyanea.Phylo` -- Phylogenetics
- `Cyanea.IO` -- File format parsing

### Usage

```elixir
{:ok, tree} = Cyanea.Phylo.parse_newick("((A:0.1,B:0.2):0.3,C:0.4);")
{:ok, result} = Cyanea.Align.align_dna("ACGT", "ACGTT", :local)
{:ok, stats} = Cyanea.Stats.describe([1.0, 2.0, 3.0, 4.0])
```

## Error Mapping

All bindings convert `CyaneaError` from Rust into target-appropriate error types:

| Rust (`CyaneaError`) | Python Exception | WASM JSON | Elixir |
|-----------------------|-----------------|-----------|--------|
| `Io(io::Error)` | `IOError` | `{"error": "..."}` | `{:error, msg}` |
| `Parse(String)` | `ValueError` | `{"error": "..."}` | `{:error, msg}` |
| `InvalidInput(String)` | `ValueError` | `{"error": "..."}` | `{:error, msg}` |
| `Compression(String)` | `RuntimeError` | `{"error": "..."}` | `{:error, msg}` |
| `Hash(String)` | `RuntimeError` | `{"error": "..."}` | `{:error, msg}` |
| `Other(String)` | `RuntimeError` | `{"error": "..."}` | `{:error, msg}` |

### Python error handling

The `IntoPyResult` trait (defined in `cyanea-py`) converts `Result<T, CyaneaError>` into `PyResult<T>`:

```python
try:
    result = align.align_dna("ACGT", "", mode="global")
except ValueError as e:
    print(f"Invalid input: {e}")
except IOError as e:
    print(f"I/O error: {e}")
```

### WASM error handling

Always check the returned JSON envelope:

```javascript
const result = JSON.parse(alignDna('ACGT', '', 'global'));
if (result.error) {
    console.error(result.error);
} else {
    console.log(result.ok);
}
```

## NumPy Interop

When the `numpy` feature is enabled, ML functions accept and return NumPy arrays directly, avoiding copy overhead for large datasets:

```bash
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release --features numpy
```

```python
import numpy as np
from cyanea import ml

data = np.random.randn(1000, 50).astype(np.float64)

# NumPy-native variants
embedding = ml.pca_numpy(data, n_components=10)    # returns np.ndarray
distances = ml.pairwise_distances_numpy(data)       # returns np.ndarray
labels = ml.kmeans_numpy(data, k=5)                 # returns np.ndarray
tsne_out = ml.tsne_numpy(data, n_components=2)      # returns np.ndarray
umap_out = ml.umap_numpy(data, n_components=2)      # returns np.ndarray
```

Without the `numpy` feature, all functions accept and return Python lists. The NumPy variants are suffixed with `_numpy` and operate on `numpy.ndarray` objects with zero-copy where possible.
