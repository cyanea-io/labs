# Cyanea Labs — Feature Roadmap

> Capabilities to add across the labs ecosystem. Organized by impact and implementation complexity.
> Items here are *new features* — see `TODO.md` for existing infrastructure work (publishing, CI, docs).

Last updated: 2026-02-18

**T1–T8**: Complete (1800+ tests). See below for history.
**T9–T14**: New — competitive gap analysis vs rust-bio, noodles, Biopython, scanpy, RDKit, Bioconductor, etc.

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

## T4 — Medium-Impact File Formats ✅

### Additional Formats (`cyanea-io`) — Done
- [x] GenBank flat file parser
- [x] GTF parser (Gene Transfer Format — GFF2 variant)
- [x] bigWig/bigBed reader (Kent binary formats)
- [x] BLAST tabular output parser (-outfmt 6/7)
- [x] MAF (Multiple Alignment Format) parser — also handles pairwise MAF (LAST/minimap2)
- [x] FASTQ paired-end interleave/deinterleave (in `cyanea-seq`)
- [x] BCF (binary VCF) reader

### Parquet Extensions (`cyanea-io`) — Done
- [x] VCF → Parquet conversion
- [x] BED → Parquet conversion
- [x] Expression matrix → Parquet (columnar gene expression)
- [x] Parquet predicate pushdown for range queries

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
- [x] Gene Ontology term association (if given annotation map)

---

## T6 — Medium-Impact Alignment & Phylo ✅

### Profile HMMs (`cyanea-align`) — Done
- [x] Profile HMM construction from MSA
- [x] Viterbi search against profile HMM
- [x] Forward/backward on profiles
- [x] E-value estimation (basic)
- [x] HMMER3-style scoring (optional)

### Phylogenetic Improvements (`cyanea-phylo`) — Done
- [x] GTR + Gamma model (generalized time-reversible)
- [x] Ancestral sequence reconstruction (marginal)
- [x] Tree rerooting and subtree extraction
- [x] Phylogenetic dating (basic molecular clock)
- [x] Tree drawing coordinates (rectangular, radial)
- [x] Consensus tree from bootstrap replicates

---

## T7 — Lower Priority / Advanced ✅

### Assembly & Graph Algorithms (`cyanea-seq`) — Done
- [x] De Bruijn graph construction from k-mers
- [x] Unitig extraction (non-branching paths)
- [x] Assembly QC metrics (N50, L50, N90, L90, GC, auN)
- [x] Contig scaffolding statistics

### Metagenomics (`cyanea-stats`, `cyanea-omics`, `cyanea-seq`) — Done
- [x] Alpha diversity (Shannon, Simpson, Chao1, observed species)
- [x] Beta diversity (Bray-Curtis dissimilarity, distance matrices)
- [x] Taxonomic classification framework (LCA-based, k-mer based)
- [x] OTU/ASV table operations (filter, rarefy, collapse, merge)

### Restriction Enzymes (`cyanea-seq`) — Done
- [x] Restriction enzyme database (20 common enzymes)
- [x] Cut site prediction on sequences (IUPAC-aware)
- [x] Fragment size prediction
- [x] In-silico digestion with multiple enzymes

### Network Biology (`cyanea-omics`) — Done
- [x] Weighted graph with correlation matrix construction
- [x] Network centrality measures (degree, betweenness, closeness)
- [x] Community detection (Louvain)
- [x] Regulatory network inference (correlation-based)

### Haplotype Analysis (`cyanea-omics`) — Done
- [x] Haplotype phasing (EM algorithm)
- [x] Haplotype block detection (LD-based)
- [x] Haplotype diversity statistics
- [x] Rarefaction curves for diversity estimation

### Motif Discovery (`cyanea-seq`) — Done
- [x] Position weight matrix (PWM) from aligned motifs
- [x] De novo motif finding (expectation-maximization, MEME-style)
- [x] Motif scanning with background model (both strands)
- [x] Information content and consensus sequences

---

## T8 — GPU Acceleration Extensions ✅

### GPU Kernels (`cyanea-gpu`) — Done
- [x] GPU-accelerated k-mer counting
- [x] GPU pairwise distance matrices (larger-than-memory tiling)
- [x] GPU Smith-Waterman for protein (BLOSUM on device)
- [x] GPU MinHash sketch computation
- [x] WebGPU backend (wgpu) for browser GPU compute

---

## T9 — Indexed Access & BAM/VCF Operations

> **Why**: Without indexed random access, any BAM/VCF over a few GB is impractical. VCF write support is needed for any pipeline that produces variants. These are the #1 operational gaps vs pysam, samtools, bcftools.

### Indexed Random Access (`cyanea-io`)
- [ ] BAI index reading (BAM index) — O(log n) region queries on BAM files
- [ ] TBI index reading (Tabix) — region queries on VCF/BED/GFF
- [ ] CSI index reading (coordinate-sorted index for large genomes)
- [ ] `fetch(chrom, start, end)` API for indexed BAM and VCF
- [ ] BGZF virtual offset support for seeking within compressed data

### BAM/CRAM Operations (`cyanea-io`)
- [ ] BAM coordinate sort (in-memory + external merge sort for large files)
- [ ] BAI index creation from sorted BAM
- [ ] BAM merge (combine multiple BAM files with header reconciliation)
- [ ] Mark PCR duplicates (position + CIGAR based)
- [ ] Fixmate (fill in mate information from name-sorted BAM)
- [ ] idxstats (per-chromosome mapped/unmapped counts from index)
- [ ] flagstat (QC pass/fail counts by flag)
- [ ] Comprehensive alignment statistics (insert size distribution, GC bias, error rates by cycle)

### VCF/BCF Writing & Manipulation (`cyanea-io`)
- [ ] VCF writer (header construction, record formatting, INFO/FORMAT field serialization)
- [ ] BCF writer (binary VCF output)
- [ ] Variant normalization (left-alignment, multiallelic splitting/joining)
- [ ] VCF merging (combine samples from multiple VCFs)
- [ ] VCF filtering with expression language
- [ ] VCF annotation (add INFO fields from external sources)
- [ ] VCF set operations (intersection, complement, concordance)
- [ ] VCF statistics (Ti/Tv ratio, per-sample stats, site frequency spectrum)

### Interval Tree Data Structure (`cyanea-omics`)
- [ ] Augmented interval tree (balanced BST with max-endpoint propagation)
- [ ] O(log n + k) overlap queries (vs current linear scan)
- [ ] Nearest/preceding/following interval queries
- [ ] Bulk loading for static interval sets (cache-oblivious layout)
- [ ] Replace or augment `IntervalSet` with tree-backed implementation
- [ ] Coverage as RLE (run-length encoded) vectors for memory-efficient genome-wide coverage

---

## T10 — Single-Cell Pipeline Completeness

> **Why**: Cyanea has AnnData/h5ad/Zarr infrastructure but lacks the algorithms to make it useful. Leiden clustering, HVG selection, and pseudotime are the standard scRNA-seq workflow (scanpy/Seurat). Without these, the container format is an empty shell.

### Preprocessing (`cyanea-omics` or `cyanea-stats`)
- [ ] Highly variable gene (HVG) selection — Seurat v3 method (variance-stabilizing transformation)
- [ ] HVG selection — Cell Ranger method (mean/dispersion)
- [ ] Library size normalization (normalize_total + log1p)
- [ ] Regress out confounding variables (e.g., mitochondrial fraction, cell cycle scores)
- [ ] Scrublet-style doublet detection
- [ ] Gene signature scoring (score_genes)

### Graph Construction & Clustering (`cyanea-ml` or `cyanea-omics`)
- [ ] kNN graph construction from PCA embeddings (with UMAP connectivities)
- [ ] Leiden community detection on kNN graphs
- [ ] Louvain community detection on kNN graphs
- [ ] Resolution parameter for cluster granularity control
- [ ] Modularity optimization with quality metrics

### Trajectory & Pseudotime (`cyanea-ml` or `cyanea-omics`)
- [ ] Diffusion maps (nonlinear dimensionality reduction)
- [ ] Diffusion pseudotime (DPT) — ordering cells along trajectories
- [ ] PAGA (Partition-based Graph Abstraction) — coarse-grained trajectory inference
- [ ] RNA velocity (spliced/unspliced ratio-based future state prediction)

### Marker Genes & Differential Expression
- [ ] Cluster-vs-rest marker gene detection (rank_genes_groups equivalent)
- [ ] Multiple test methods: t-test, Wilcoxon, logistic regression
- [ ] Marker gene filtering (log2FC threshold, pct expressed, adjusted p-value)

### Batch Correction & Integration
- [ ] Harmony integration (iterative PCA correction)
- [ ] ComBat batch correction (parametric empirical Bayes)
- [ ] MNN (Mutual Nearest Neighbors) correction
- [ ] Integration quality metrics (kBET, LISI)

---

## T11 — Microbiome & Ecological Statistics

> **Why**: Cyanea has alpha diversity and Bray-Curtis but lacks the statistical tests and ordination methods that define standard microbiome analysis (scikit-bio, vegan/R, phyloseq). UniFrac + PCoA + PERMANOVA together unlock the entire microbiome analysis workflow.

### Phylogenetic Diversity (`cyanea-stats` or `cyanea-phylo`)
- [ ] Unweighted UniFrac distance (presence/absence on phylogenetic tree)
- [ ] Weighted UniFrac distance (abundance-weighted branch lengths)
- [ ] Generalized UniFrac (alpha parameter for weighting control)
- [ ] Faith's phylogenetic diversity (PD)

### Ordination Methods (`cyanea-stats`)
- [ ] PCoA (Principal Coordinates Analysis) — eigendecomposition of distance matrices
- [ ] CCA (Canonical Correspondence Analysis) — constrained ordination
- [ ] RDA (Redundancy Analysis) — linear constrained ordination
- [ ] NMDS (Non-metric Multidimensional Scaling)
- [ ] Procrustes analysis (compare ordinations)

### Multivariate Tests (`cyanea-stats`)
- [ ] PERMANOVA (Permutational MANOVA) — test group centroids via permutation
- [ ] ANOSIM (Analysis of Similarities) — rank-based multivariate test
- [ ] Mantel test — correlation between two distance matrices
- [ ] BIOENV — environmental variable selection explaining community patterns
- [ ] AMOVA (Analysis of Molecular Variance)

### Diversity Extensions (`cyanea-stats`)
- [ ] Alpha rarefaction curves (species richness vs sampling depth)
- [ ] Beta diversity: Jaccard (presence/absence), weighted Jaccard
- [ ] Diversity profile (Hill numbers / effective number of species)

---

## T12 — Cheminformatics Depth

> **Why**: cyanea-chem can parse molecules and compute fingerprints, but drug discovery workflows require SMARTS, 200+ descriptors, 3D conformers, and reaction chemistry. These are standard in RDKit; adding them moves cyanea-chem from "structural parsing" to "computational chemistry."

### SMARTS & Advanced Matching (`cyanea-chem`)
- [ ] SMARTS pattern language parser
- [ ] SMARTS-based substructure search (atom/bond properties, recursive SMARTS)
- [ ] SMARTS logical operators (AND, OR, NOT on atom primitives)
- [ ] Reaction SMARTS (SMIRKS) for transformation rules

### Molecular Descriptors (`cyanea-chem`)
- [ ] Topological descriptors (Wiener index, Balaban J, Zagreb indices)
- [ ] TPSA (Topological Polar Surface Area)
- [ ] Wildman-Crippen LogP/MR (atom contribution method)
- [ ] BertzCT complexity index
- [ ] Kappa shape indices (1κ, 2κ, 3κ)
- [ ] Chi connectivity indices
- [ ] Fraction sp3 carbons
- [ ] Ring count details (aliphatic, aromatic, heteroaromatic, spiro, fused)
- [ ] EState descriptors
- [ ] Autocorrelation descriptors (Moreau-Broto, Moran, Geary)
- [ ] Descriptor batch computation (all descriptors at once for QSAR matrices)

### 3D Coordinate Generation (`cyanea-chem`)
- [ ] Distance geometry embedding (initial 3D coords from connectivity)
- [ ] ETKDG-style conformer generation (torsion angle preferences)
- [ ] Multiple conformer enumeration (rotatable bond sampling)
- [ ] MMFF94 energy calculation
- [ ] UFF energy calculation
- [ ] Basic energy minimization (steepest descent / conjugate gradient)

### Molecule Standardization (`cyanea-chem`)
- [ ] Salt/fragment stripping (remove counterions, solvents)
- [ ] Charge neutralization
- [ ] Tautomer canonicalization
- [ ] Largest fragment extraction
- [ ] Standardization pipeline (composable steps)

### Drug-Likeness Filters (`cyanea-chem`)
- [ ] Lipinski Rule of 5
- [ ] Veber rules (rotatable bonds + TPSA)
- [ ] PAINS (Pan Assay Interference) structural alerts
- [ ] Brenk structural alerts
- [ ] Lead-likeness criteria
- [ ] Drug-likeness score (QED — Quantitative Estimate of Drug-likeness)

### Scaffold Analysis (`cyanea-chem`)
- [ ] Murcko scaffold decomposition (framework + sidechains)
- [ ] Generic scaffold generation
- [ ] Maximum Common Substructure (MCS) between molecule pairs
- [ ] R-group decomposition (core + variable groups)

### Chemical Reactions (`cyanea-chem`)
- [ ] Reaction SMILES (SMIRKS) application to molecules
- [ ] Reaction enumeration (virtual library generation)
- [ ] Atom-atom mapping in reactions
- [ ] Retrosynthetic disconnection (single-step)

---

## T13 — Phylogenetics & Alignment Maturity

> **Why**: NNI-only tree search is weak for real phylogenetics. Missing protein models (LG, WAG, JTT) means ML phylo is nucleotide-only. SPR search + model selection + protein models would bring cyanea-phylo to IQ-TREE/RAxML competitive parity. Additional alignment formats round out interoperability.

### Tree Search (`cyanea-phylo`)
- [ ] SPR (Subtree Pruning and Regrafting) moves
- [ ] TBR (Tree Bisection and Reconnection) moves
- [ ] Parsimony ratchet for starting tree estimation
- [ ] Stochastic NNI with simulated annealing
- [ ] Lazy SPR evaluation (partial likelihood recomputation)

### Model Selection (`cyanea-phylo`)
- [ ] AIC (Akaike Information Criterion) for model comparison
- [ ] BIC (Bayesian Information Criterion)
- [ ] ModelFinder-style automatic model selection (test all models, rank by IC)
- [ ] LRT (Likelihood Ratio Test) for nested models

### Protein Substitution Models (`cyanea-phylo`)
- [ ] LG (Le & Gascuel 2008)
- [ ] WAG (Whelan & Goldman 2001)
- [ ] JTT (Jones, Taylor & Thornton 1992)
- [ ] Dayhoff (Dayhoff, Schwartz & Orcutt 1978)
- [ ] Rate matrix + frequency vector loading from file
- [ ] Protein + Gamma rate variation

### Bayesian Phylogenetics (`cyanea-phylo`)
- [ ] MCMC sampler over tree topologies and branch lengths
- [ ] Metropolis-Hastings proposals (NNI, SPR on tree; scale on branches)
- [ ] Coalescent priors (constant population, exponential growth)
- [ ] Relaxed molecular clocks (uncorrelated lognormal)
- [ ] Convergence diagnostics (ESS, trace plots data)
- [ ] Posterior summary (MAP tree, credible intervals on node ages)

### Gene Tree / Species Tree (`cyanea-phylo`)
- [ ] ASTRAL-style species tree estimation from gene trees
- [ ] Gene/species tree reconciliation (duplication/loss/ILS)
- [ ] Concordance factors (gene and site)

### Alignment Format I/O (`cyanea-io` or `cyanea-seq`)
- [ ] Stockholm format (Pfam/Rfam alignments, used by HMMER)
- [ ] Clustal format (read/write)
- [ ] PHYLIP format (interleaved and sequential)
- [ ] EMBL sequence format
- [ ] PIR/NBRF format

### Motif Format I/O (`cyanea-seq`)
- [ ] MEME motif format (read/write)
- [ ] TRANSFAC format (read/write)
- [ ] JASPAR format (read/write)
- [ ] Motif-to-motif comparison (Pearson correlation of PWM columns)

---

## T14 — Variant Annotation, Ecosystem & Emerging

> **Why**: Variant effect prediction (VEP), database clients, and specialized -omics analyses. These are the long tail of features that complete the ecosystem — individually less urgent, collectively they determine whether Cyanea is "production-ready" or "research-grade."

### Variant Effect Prediction (`cyanea-omics`)
- [ ] Coding consequence classification (missense, nonsense, frameshift, splice site, synonymous, UTR)
- [ ] Codon change annotation from VCF + GFF/GTF transcript models
- [ ] Amino acid change notation (p.V600E style)
- [ ] SIFT-style functional impact scoring (sequence conservation-based)
- [ ] Splice site disruption scoring

### Database & API Clients (`cyanea-io` or new `cyanea-fetch`)
- [ ] NCBI Entrez client (efetch, esearch, elink for sequences, annotations, publications)
- [ ] UniProt REST client (protein records, batch retrieval)
- [ ] KEGG REST client (pathways, compounds, reactions)
- [ ] htsget client (GA4GH streaming protocol for BAM/VCF)
- [ ] refget client (GA4GH reference sequence retrieval)

### Copy Number & Structural Variants (`cyanea-omics`)
- [ ] Read depth-based CNV detection (circular binary segmentation)
- [ ] B-allele frequency segmentation
- [ ] SV breakpoint detection from split/discordant reads
- [ ] CNV call merging and annotation

### Methylation Analysis (`cyanea-omics` or `cyanea-seq`)
- [ ] Bisulfite sequencing alignment support (C→T conversion-aware)
- [ ] Methylation calling from bisulfite BAM
- [ ] Differential methylation regions (DMRs)
- [ ] CpG island annotation

### Spatial Transcriptomics (`cyanea-omics`)
- [ ] Spatial neighbors graph (Delaunay triangulation, k-nearest spatial)
- [ ] Spatial autocorrelation (Moran's I, Geary's C)
- [ ] Co-occurrence analysis
- [ ] Ligand-receptor interaction scoring

### Simulation (`cyanea-seq`, `cyanea-stats`)
- [ ] Sequence evolution simulation (along a phylogenetic tree with substitution model)
- [ ] Population genetics simulation (Wright-Fisher, coalescent)
- [ ] Read simulator (Illumina error profile, coverage model)
- [ ] Null model generation for statistical testing

### Additional I/O
- [ ] BLAST XML output parser (legacy but common)
- [ ] ABI chromatogram / .ab1 file parser (Sanger sequencing)
- [ ] bedGraph/Wiggle output (for genome browser visualization)
- [ ] GFA (Graphical Fragment Assembly) format for pangenomics

---

## Implementation Notes

### Crate Placement Guidelines
- **New capabilities in existing crates** are preferred over new crates
- Only create a new crate when the dependency graph demands it
- T9 indexed access goes in `cyanea-io` (already depends on noodles which has index support)
- T10 single-cell algorithms split between `cyanea-ml` (graph clustering, trajectory) and `cyanea-omics` (preprocessing, HVG)
- T11 microbiome stats go in `cyanea-stats` (ordination, multivariate tests) and `cyanea-phylo` (UniFrac)
- T12 cheminformatics depth stays in `cyanea-chem`
- T14 database clients could be a new `cyanea-fetch` crate or go in `cyanea-io` with feature gating
- T14 variant annotation goes in `cyanea-omics` (already has variant types + gene annotations)

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
