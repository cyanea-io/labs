# Cyanea Labs — Feature Roadmap

> Capabilities to add across the labs ecosystem. Organized by impact and implementation complexity.
> Items here are *new features* — see `TODO.md` for existing infrastructure work (publishing, CI, docs).

Last updated: 2026-02-17

---

## T1 — High-Impact Sequence & Alignment

### Quality Trimming & Filtering (`cyanea-seq`) — Done
- [x] Sliding-window quality trimming (Trimmomatic-style)
- [x] Leading/trailing low-quality base removal
- [x] Adapter detection and removal (common Illumina adapters)
- [x] Length filtering (min/max)
- [x] Quality score statistics per-base and per-read
- [x] Complexity filtering (entropy-based)
- [x] BWA-style 3' quality trimming
- [x] TrimPipeline builder with batch statistics (TrimReport)

### Low-Complexity & Repeat Masking (`cyanea-seq`) — Done
- [x] DUST algorithm for low-complexity filtering
- [x] SEG algorithm for protein low-complexity regions
- [x] Tandem repeat finder (basic)
- [x] Soft/hard masking output

### Codon Usage & Translation Tables (`cyanea-seq`) — Done
- [x] Codon usage table computation from sequences
- [x] Codon Adaptation Index (CAI)
- [x] Alternative genetic codes (7 NCBI tables: standard, vertebrate/yeast/invertebrate mito, mycoplasma, ciliate, bacterial)
- [x] Synonymous/non-synonymous codon classification
- [x] RSCU, relative adaptiveness, Nei-Gojobori site counting

### X-Drop Alignment (`cyanea-align`) — Done
- [x] X-drop extension for seed-and-extend (terminate early on poor extensions)
- [x] Z-drop variant (MINIMAP2-style, more permissive)
- [x] Configurable drop thresholds and bandwidth
- [x] Seed-extend combining (left + seed + right)

### CIGAR String Utilities (`cyanea-align`) — Done
- [x] CIGAR parsing and validation
- [x] CIGAR ↔ alignment string conversion
- [x] CIGAR arithmetic (split, merge, reverse)
- [x] Alignment statistics from CIGAR (identity, coverage, gaps)
- [x] MD tag generation from CIGAR + sequences

### Spliced Alignment (`cyanea-align`) — Done
- [x] Intron-aware alignment (GT-AG / GC-AG / AT-AC splice sites)
- [x] Configurable splice penalty and canonical/semi-canonical bonuses
- [x] Multi-exon alignment chaining with colinearity enforcement
- [x] Splice site detection, CIGAR with Skip (N) ops

---

## T2 — High-Impact Genomics & Omics

### Genome Arithmetic (`cyanea-omics`) — Done
- [x] Interval intersection (BEDTools-style)
- [x] Interval union, subtraction, complement
- [x] Closest interval queries
- [x] Window/tile generation across genome
- [x] Jaccard similarity between interval sets
- [x] Genomic coordinate liftover (chain file parsing)

### Pileup Generation (`cyanea-io`) — Done
- [x] SAM/BAM → per-position pileup
- [x] Base counts, quality sums at each position
- [x] Variant calling support (Bayesian variant caller with genotype likelihoods)
- [x] Depth-of-coverage calculation
- [x] mpileup-compatible output

### Population Genetics (`cyanea-stats::popgen`) — Done
- [x] Allele frequencies from 012-encoded genotypes
- [x] Hardy-Weinberg equilibrium test
- [x] Fst (fixation index) — Weir-Cockerham and Hudson estimators
- [x] Tajima's D
- [x] Nucleotide diversity (pi, theta)
- [x] Linkage disequilibrium (r², D') with EM haplotype estimation
- [x] Principal components on genotype matrices (eigenanalysis)

### Differential Expression (`cyanea-stats`) — Done
- [x] TPM / FPKM / CPM normalization
- [x] DESeq2-style size factor estimation
- [x] Negative binomial Wald test for differential expression
- [x] Volcano plot data generation (log2FC, p-value, adjusted p-value)
- [x] Gene set enrichment analysis (GSEA) — preranked, score-weighted, permutation-based
- [x] Over-representation analysis (hypergeometric test)

### Survival Analysis (`cyanea-stats::survival`) — Done
- [x] Kaplan-Meier estimator
- [x] Log-rank test
- [x] Cox proportional hazards model (basic)
- [x] Median survival, confidence intervals

---

## T3 — High-Impact RNA & Structure

### RNA Secondary Structure (`cyanea-seq`) — Done
- [x] Nussinov algorithm (maximum base pairs)
- [x] Zuker minimum free energy (nearest-neighbor thermodynamics)
- [x] Dot-bracket notation I/O
- [x] Base pair probability matrix
- [x] Partition function (McCaskill's algorithm)
- [x] Structure comparison (base pair distance, mountain distance)

### Protein Sequence Properties (`cyanea-seq`)
- [x] Hydrophobicity profiles (Kyte-Doolittle, Hopp-Woods)
- [x] Isoelectric point (pI) calculation
- [x] Amino acid composition statistics
- [x] Gravy (grand average of hydropathicity)
- [x] Extinction coefficient estimation
- [x] Intrinsic disorder prediction (IUPred-style simple predictor)
- [x] Secondary structure prediction (GOR, Chou-Fasman)

---

## T4 — Medium-Impact File Formats

### Additional Formats (`cyanea-io`)
- [ ] GenBank/EMBL flat file parser
- [ ] GTF parser (Gene Transfer Format — GFF2 variant)
- [ ] bigWig/bigBed reader (Kent binary formats)
- [ ] BLAST tabular output parser (-outfmt 6/7)
- [ ] MAF (Multiple Alignment Format) parser
- [ ] Pairwise maf (LAST/minimap2 output)
- [ ] FASTQ paired-end interleave/deinterleave
- [ ] BCF (binary VCF) reader

### Parquet Extensions (`cyanea-io`)
- [ ] VCF → Parquet conversion
- [ ] BED → Parquet conversion
- [ ] Expression matrix → Parquet (columnar gene expression)
- [ ] Parquet predicate pushdown for range queries

---

## T5 — Medium-Impact ML & Statistics

### Model Evaluation (`cyanea-ml`) — Done
- [x] ROC curve computation (TPR/FPR at thresholds)
- [x] AUC (area under ROC curve) — trapezoidal
- [x] Precision-recall curve
- [x] F1 score, Matthews correlation coefficient
- [x] Confusion matrix computation
- [x] Cross-validation (k-fold, stratified k-fold, leave-one-out)

### Gradient Boosting (`cyanea-ml`) — Done
- [x] Gradient-boosted decision trees (basic GBDT)
- [x] Feature importance (Gini, permutation)
- [x] Regression and classification modes
- [x] Early stopping

### Feature Selection (`cyanea-ml`) — Done
- [x] Variance threshold filtering
- [x] Mutual information (discrete)
- [x] Recursive feature elimination (with any classifier)
- [x] L1-regularized feature selection (Lasso)

### Enrichment Analysis (`cyanea-stats`) — Done
- [x] Fisher's exact test on gene sets
- [x] Hypergeometric test for overlap significance (ORA)
- [x] Benjamini-Hochberg on enrichment p-values
- [x] GSEA preranked (Subramanian et al. 2005, permutation-based)
- [ ] Gene Ontology term association (if given annotation map)

---

## T6 — Medium-Impact Alignment & Phylo

### Profile HMMs (`cyanea-align`)
- [ ] Profile HMM construction from MSA
- [ ] Viterbi search against profile HMM
- [ ] Forward/backward on profiles
- [ ] E-value estimation (basic)
- [ ] HMMER3-style scoring (optional)

### Phylogenetic Improvements (`cyanea-phylo`)
- [ ] GTR + Gamma model (generalized time-reversible)
- [ ] Ancestral sequence reconstruction (marginal)
- [ ] Tree rerooting and subtree extraction
- [ ] Phylogenetic dating (basic molecular clock)
- [ ] Tree drawing coordinates (rectangular, radial)
- [ ] Consensus tree from bootstrap replicates

---

## T7 — Lower Priority / Advanced

### Assembly & Graph Algorithms (`cyanea-seq` or new crate)
- [ ] De Bruijn graph construction from k-mers
- [ ] Unitig extraction (non-branching paths)
- [ ] Assembly QC metrics (N50, L50, total length, GC)
- [ ] Contig scaffolding statistics

### Metagenomics (`cyanea-seq` or `cyanea-omics`)
- [ ] Alpha diversity (Shannon, Simpson, Chao1, observed species)
- [ ] Beta diversity (Bray-Curtis, UniFrac if trees available)
- [ ] Taxonomic classification framework (LCA-based, k-mer based)
- [ ] OTU/ASV table operations

### Restriction Enzymes (`cyanea-seq`)
- [ ] Restriction enzyme database (common enzymes)
- [ ] Cut site prediction on sequences
- [ ] Fragment size prediction
- [ ] Virtual gel electrophoresis data

### Network Biology (new `cyanea-net` or `cyanea-omics`)
- [ ] Protein-protein interaction graph
- [ ] Network centrality measures (degree, betweenness, closeness)
- [ ] Community detection (Louvain)
- [ ] Regulatory network inference (basic correlation-based)

### Haplotype Analysis (`cyanea-omics` or `cyanea-stats`)
- [ ] Haplotype phasing (basic EM algorithm)
- [ ] Haplotype block detection
- [ ] Haplotype diversity statistics
- [ ] Phase-informative read extraction

### Motif Discovery (`cyanea-seq`)
- [ ] Position weight matrix (PWM) from aligned motifs
- [ ] De novo motif finding (expectation-maximization)
- [ ] Motif scanning with background model
- [ ] Information content and sequence logos data

---

## T8 — GPU Acceleration Extensions

### GPU Kernels (`cyanea-gpu`)
- [ ] GPU-accelerated k-mer counting
- [ ] GPU pairwise distance matrices (larger-than-memory tiling)
- [ ] GPU Smith-Waterman for protein (BLOSUM on device)
- [ ] GPU MinHash sketch computation
- [ ] WebGPU backend (wgpu) for browser GPU compute

---

## Implementation Notes

### Crate Placement Guidelines
- **New capabilities in existing crates** are preferred over new crates
- Only create a new crate when the dependency graph demands it (e.g., `cyanea-rna` if it needs both `cyanea-seq` and `cyanea-struct`)
- Population genetics could be a `popgen` module in `cyanea-stats` or a standalone `cyanea-popgen` crate
- Network biology could go in `cyanea-omics` or a new `cyanea-net`

### Feature Gating
- Large optional deps should be feature-gated (e.g., `hdf5`, `bigwig`)
- Keep `default = ["std"]` convention
- Compute-heavy features behind `parallel` flag where Rayon helps

### Testing Standards
- Unit tests inline, property tests with `proptest` for round-trip invariants
- Fuzz targets for all new parsers
- Criterion benchmarks for performance-critical additions
- Test against known-good outputs from established tools (samtools, bedtools, etc.)

### Binding Propagation
- Each new public API should be exposed through WASM, Python, and NIF bindings
- WASM: add to `cyanea-wasm/src/` with JSON envelope
- Python: add to `cyanea-py/src/` with PyO3 classes
- NIF: add to `cyanea/native/cyanea_native/src/` with rustler exports
- Elixir: add to corresponding domain module (`Cyanea.Seq`, `Cyanea.Stats`, etc.)
