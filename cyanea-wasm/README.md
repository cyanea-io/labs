# cyanea-wasm

> WebAssembly bindings for browser-based bioinformatics. Published as `@cyanea/bio` on npm.

## What's Inside

JSON-based WASM API wrapping 10 domain modules:

- **Sequence** -- FASTA/FASTQ parsing, paired-end FASTQ, read trimming, GC content, reverse complement, transcription, translation, validation
- **Alignment** -- DNA and protein alignment (all modes + batch), custom scoring, BLOSUM/PAM matrices
- **Statistics** -- descriptive stats, correlation, hypothesis tests, survival (Kaplan-Meier, Cox PH), population genetics (Fst, Tajima's D), diversity, null models
- **Machine learning** -- k-mer counting, distance metrics, UMAP, PCA, t-SNE, K-means, random forest, GBDT, HMM, confusion matrix, ROC/PR curves, cross-validation, feature selection
- **Chemistry** -- SMILES parsing, molecular properties, Morgan fingerprints, Tanimoto similarity, canonical SMILES, substructure search
- **Structural biology** -- PDB/mmCIF parsing, secondary structure, RMSD, contact maps, Ramachandran, Kabsch
- **Phylogenetics** -- Newick/NEXUS trees, distances, UPGMA/NJ, Robinson-Foulds, simulation
- **I/O** -- VCF/BED/GFF3/BLAST XML/bedGraph/GFA text parsing, NCBI/UniProt URL builders
- **Omics** -- interval operations, variant annotation, CNV segmentation, methylation, spatial autocorrelation
- **Core utilities** -- SHA-256 hashing, zstd compression/decompression

## Quick Start

```toml
[dependencies]
cyanea-wasm = { version = "0.1", features = ["wasm"] }
```

### JavaScript usage

```javascript
import { alignDna, describe, smilesProperties } from '@cyanea/bio';

const result = JSON.parse(alignDna('ACGTACGT', 'ACGTGCGT', 'local'));
console.log(result.ok.score);

const stats = JSON.parse(describe('[1, 2, 3, 4, 5]'));
console.log(stats.ok.mean);

const props = JSON.parse(smilesProperties('c1ccccc1O'));
console.log(props.ok.formula);
```

### Building

```bash
# Check
cargo check -p cyanea-wasm --features wasm

# Build WASM package
wasm-pack build cyanea-wasm --features wasm

# Build TypeScript wrapper
cd cyanea-wasm && npm run build:ts
```

## Design

- **JSON interface** -- all functions accept simple types (strings, numbers) and return JSON strings
- **Error envelope** -- success: `{"ok": value}`, failure: `{"error": "message"}`
- **Stateless** -- all functions are pure, no global state or initialization
- **wasm-bindgen** -- all public functions have `#[cfg_attr(feature = "wasm", wasm_bindgen)]`

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `wasm` | No | Enables wasm-bindgen annotations |

## Modules

| Module | Description |
|--------|-------------|
| `error` | JSON error wrapping (`wasm_ok`, `wasm_err`, `wasm_result`) |
| `seq` | FASTA/FASTQ, paired-end, trimming, RNA folding, protein properties, simulation, codon usage |
| `align` | DNA/protein alignment, batch, MSA, POA, banded, CIGAR utilities |
| `stats` | Descriptive stats, correlation, hypothesis tests, survival, popgen, diversity, null models |
| `ml` | K-mer counting, distances, PCA, t-SNE, UMAP, K-means, random forest, GBDT, HMM, metrics, CV |
| `chem` | SMILES/SDF, fingerprints (Morgan, MACCS), properties, substructure |
| `struct_bio` | PDB/mmCIF, secondary structure, RMSD, contact maps, Ramachandran, Kabsch |
| `phylo` | Newick/NEXUS, distances, UPGMA/NJ, Robinson-Foulds, simulation |
| `io` | VCF/BED/GFF3/BLAST XML/bedGraph/GFA text parsing, NCBI/UniProt URLs |
| `omics` | Interval ops, variant annotation, CNV, methylation, spatial autocorrelation |
| `core_utils` | SHA-256, zstd compression |

## See Also

- [API Reference](docs/API.md)
- [Usage Guide](docs/GUIDE.md)
- [Internal Architecture](docs/ARCHITECTURE.md)
- [Workspace Architecture](../docs/ARCHITECTURE.md)
- [Build Guide](../docs/BUILDING.md)
- [Bindings Guide](../docs/BINDINGS.md)
