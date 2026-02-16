# Cyanea Labs — Competitive Feature Roadmap

## Tier 1 — Blocking for real-world adoption

- [x] **Paired-end read support** (cyanea-seq)
  - Paired FASTQ reader/writer
  - Read1/Read2 pairing and interleaving/deinterleaving
  - Insert size validation
  - Paired-aware trimming in TrimPipeline

- [ ] **Expose existing Rust features to bindings** (cyanea-wasm, cyanea-py)
  - ~~Python: CIGAR utilities~~ ✓
  - ~~Python: MinHash~~ ✓
  - ~~Python: PCA, t-SNE, K-means~~ ✓
  - ~~WASM: quality trimming~~ ✓
  - ~~WASM: CIGAR utilities~~ ✓
  - Python: genome arithmetic, liftover
  - Python: pileup/depth analysis
  - Python: AnnData/h5ad
  - WASM: PCA, t-SNE, K-means
  - WASM: MinHash
  - WASM: pileup/depth analysis

- [ ] **BAM indexing + indexed region queries** (cyanea-io)
  - .bai/.csi index generation
  - Region-based BAM/CRAM queries (random access)
  - IGV / samtools interop

## Tier 2 — Required for scRNA-seq / omics workflows

- [ ] **10x Genomics format support** (cyanea-omics)
  - 10x HDF5 feature-barcode matrix reader
  - MEX (Matrix Market) format reader
  - Conversion to AnnData

- [ ] **UMI/barcode extraction and deduplication** (cyanea-seq)
  - Barcode extraction from read name / sequence prefix
  - Cell barcode correction / whitelisting
  - UMI deduplication (directional / cluster)

- [ ] **Normalization + differential expression** (cyanea-stats or cyanea-omics)
  - Normalization: TPM, CPM, size factors
  - DE: Wilcoxon rank-sum (single-cell), negative binomial GLM (bulk)
  - Volcano plot data output

## Tier 3 — Differentiation / nice-to-have

- [ ] **CLI tool** (new crate: cyanea-cli)
  - `cyanea trim` — FASTQ quality trimming
  - `cyanea align` — pairwise / batch alignment
  - `cyanea stats` — FASTQ/BAM summary statistics
  - `cyanea convert` — format conversion

- [ ] **Gene set enrichment** (cyanea-omics or cyanea-stats)
  - GMT file parsing
  - Hypergeometric GO enrichment
  - GSEA (rank-based)

- [ ] **Survival analysis** (cyanea-stats)
  - Kaplan-Meier estimator
  - Log-rank test
  - Cox proportional hazards regression

- [ ] **Streaming iterator API** (cyanea-io, cyanea-seq)
  - `Iterator<Item = Record>` for FASTQ, SAM/BAM, VCF
  - Reduce memory footprint for large files
