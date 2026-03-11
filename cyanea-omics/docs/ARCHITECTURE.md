# cyanea-omics Architecture

Internal design of the omics data structures and analysis modules.

## Coordinate System

All genomic coordinates use 0-based, half-open intervals `[start, end)`, matching BED format and most bioinformatics tools. `GenomicInterval` carries a chromosome name (`String`), start/end (`u64`), and strand (`Strand` enum).

## IntervalTree

The `IntervalTree<T>` is a static augmented binary search tree stored in an implicit array layout (children of node `i` at `2i+1` and `2i+2`). Build once, query many times.

- **Construction**: Sort intervals by start coordinate, then recursively build a balanced BST into the array. After placement, propagate `max_end` values bottom-up so each node stores the maximum end coordinate in its subtree.
- **Query**: For a query range `[q_start, q_end)`, recurse into the left subtree only if `left.max_end > q_start`, and into the right subtree only if the node's start < q_end. This achieves O(log n + k) time where k is the number of results.
- **No dynamic updates**: The tree is immutable after construction. For mutable workloads, use `IntervalSet` instead.

## Coverage Vectors

`RleCoverage` stores depth values as `(depth: u32, run_length: u64)` pairs. Built from intervals using a sweep-line approach: collect `+1` events at interval starts and `-1` events at interval ends, sort, sweep to compute depth at each position, then merge consecutive runs with the same depth. Memory usage is proportional to the number of distinct depth transitions, not genome size.

## ExpressionMatrix

Dense row-major layout: `data[gene_idx][sample_idx]` stored as `Vec<Vec<f64>>`. Feature and sample names stored alongside. Suited for bulk RNA-seq where the matrix is typically small enough to be fully dense.

## SparseMatrix

The `SparseMatrix` stores entries in COO format (coordinate list of `(row, col, value)` triplets). This is flexible for incremental construction and random access patterns common in single-cell data.

- **COO format**: Direct storage of non-zero entries. O(nnz) for element lookup, efficient for iteration and modification.
- **CSR conversion**: `to_csr()` produces compressed sparse row format `(data, indices, indptr)` for efficient row-wise access and matrix-vector products.
- **Column operations**: `column_sums()`, `column_means()`, `col_nnz()` iterate over all stored entries (O(nnz)).
- **Row scaling**: `scale_rows()` multiplies each entry by a per-row factor, used for library-size normalization.

## AnnData Container

`AnnData` mirrors the Python anndata container:

- `X: MatrixData` -- primary expression matrix (Dense or Sparse)
- `obs: HashMap<String, ColumnData>` -- per-observation annotations (cell metadata)
- `var: HashMap<String, ColumnData>` -- per-variable annotations (gene metadata)
- `obsm: HashMap<String, Vec<Vec<f64>>>` -- observation embeddings (PCA, UMAP)
- `varm: HashMap<String, Vec<Vec<f64>>>` -- variable embeddings
- `obsp: HashMap<String, SparseMatrix>` -- pairwise observation matrices (kNN graphs)
- `layers: HashMap<String, MatrixData>` -- alternative expression matrices (raw, normalized)
- `uns: HashMap<String, String>` -- unstructured metadata

`ColumnData` is an enum supporting `Strings`, `Numeric`, and `Categorical` (codes + categories) columns, allowing typed metadata.

### Lazy Backends

The `h5ad` and `zarr` modules provide serialization. Both read the full container into memory on load (no lazy/chunked access). For h5ad, CSR sparse matrices are stored using the `X/data`, `X/indices`, `X/indptr` convention. For zarr, the same structure is replicated using zarrs 0.18 directory layout.

## Variant Annotation

The annotation pipeline works in three stages:

1. **TranscriptModel construction**: Build an `IntervalTree` from gene/transcript coordinates for fast overlap queries.
2. **Overlap query**: For each variant position, query the interval tree to find overlapping transcripts.
3. **Consequence prediction**: For each overlapping transcript, determine the variant's location (exon, intron, UTR, splice site) and compute the coding consequence:
   - Extract the reference codon and mutant codon using the transcript's CDS coordinates
   - Translate both codons using the standard genetic code
   - Classify as missense, nonsense, synonymous, frameshift, etc.
   - Generate HGVS notation (c. and p.)
   - Compute splice disruption scores based on distance from splice sites

## Genome Arithmetic

All operations work on `&[GenomicInterval]` slices. Internally, intervals are grouped by chromosome using a `BTreeMap`, then sorted by start position within each group. Operations like intersect and subtract use a merge-sort-like sweep over the two sorted lists per chromosome. The `StrandMode` enum controls whether strand matching is ignored, required to be same, or required to be opposite.

## Single-Cell Pipeline

The single-cell pipeline (feature-gated behind `single-cell`) follows a modular design where each step operates on `&mut AnnData`:

1. **Preprocessing** (`sc_preprocess`): Normalize, log-transform, select HVGs, regress confounders. Results stored in `X` or `layers["raw"]`.
2. **Neighbor graph** (`sc_cluster`): Compute PCA on HVGs, build kNN graph, store in `obsp["connectivities"]`.
3. **Clustering** (`sc_cluster`): Leiden/Louvain on the kNN graph, store assignments in `obs["leiden"]`.
4. **Trajectory** (`sc_trajectory`): Diffusion map from the graph, DPT from a root cell, PAGA for coarse-grained topology.
5. **Markers** (`sc_markers`): One-vs-rest differential expression per cluster using t-test, Wilcoxon, or logistic regression.
6. **Integration** (`sc_integrate`): Harmony iteratively corrects PCA embeddings in `obsm["X_pca"]`. ComBat corrects the expression matrix. MNN corrects via mutual nearest neighbors.

Each module reads and writes standard AnnData slots, so steps can be composed in any order.

## CNV and Methylation

- **CBS** (Circular Binary Segmentation): Recursively finds the split point that maximizes the t-statistic between segment means, using permutation testing for significance. Segments are then optionally merged based on log2-ratio similarity.
- **Methylation**: CpG sites carry read-level methylation counts. DMR detection uses a sliding window approach, testing beta-value differences between groups with Welch's t-test. CpG island detection uses the Gardiner-Garden & Frommer criteria (GC% > 50%, observed/expected CpG ratio > 0.6, length > 200bp).

## Spatial Transcriptomics

- **Neighbor graphs**: Delaunay triangulation uses divide-and-conquer on the 2D point set. k-nearest neighbor graphs compute all pairwise distances and select the k closest per point.
- **Moran's I**: Computed as a weighted covariance of deviations from the mean, normalized by the variance, with analytic p-value from the normal approximation.
- **Ligand-receptor**: For each ligand-receptor pair, compute the mean product of ligand expression in sending cells and receptor expression in receiving neighbors. Significance by permuting cell labels.

## ACMG/AMP Variant Classification

Implements the Richards et al. 2015 combining rules for the 5-tier ACMG/AMP classification system. The `AcmgEvidence` builder collects criteria with their pathogenic/benign direction and evidence strength. On `classify()`, criteria are counted by category (PVS, PS, PM, PP for pathogenic; BA, BS, BP for benign) and the combining rules are applied in priority order: benign rules first (BA1 is stand-alone), then pathogenic, then likely pathogenic, defaulting to VUS.

- **Auto-evidence** (`auto_evidence`): A simplified auto-classifier that inspects variant properties to assign basic criteria: PVS1 for null variants (frameshift indels) in known LOF genes, PM2 for absent/rare population frequency (<0.0001), PP3/BP4 from in silico predictions, BA1 for common variants (>5% AF), and BS1 for moderately common variants (>1% AF). Not sufficient for clinical use without manual review.
- **ClinVar matching**: Exact-match lookup by (chrom, position, ref_allele, alt_allele) against a vector of annotations. `parse_clinvar_tsv` reads a simple 9-column tab-separated format for building annotation databases.

## Pharmacogenomics

Star allele calling follows the CPIC (Clinical Pharmacogenetics Implementation Consortium) model:

1. **Allele definition**: Each `StarAllele` carries a set of defining variants (chrom, position, ref, alt) and an activity score (0.0 = no function through 2.0 = increased function).
2. **Calling**: `call_star_alleles` builds a set of observed variant keys, then scores each allele definition by the fraction of its defining variants that are present. Alleles with >50% match fraction are candidates, sorted by match quality. The top two are selected as the diplotype; if fewer than two match, the reference allele (*1, activity 1.0) fills in.
3. **Phenotype**: The combined activity score (sum of both alleles) maps to a metabolizer phenotype via CPIC thresholds: >2.25 ultrarapid, >2.0 rapid, 1.25-2.0 normal, 0.25-1.0 intermediate, <0.25 poor.
4. **Drug interactions**: `lookup_drug_interactions` filters the `PgxDatabase` interaction list by gene name and phenotype to retrieve clinical recommendations with evidence levels and source guidelines.

The `PgxDatabase` stores allele definitions keyed by gene name (HashMap) and a flat list of drug-gene interactions, designed for extensibility with custom gene panels.

## Clinical Genomics

Three clinical assays commonly used in oncology and transplant medicine:

### HLA Typing

`HlaAllele` stores a gene locus (e.g., "A", "B", "DRB1") and an allele designation at arbitrary resolution (e.g., "02:01:01:01"). Resolution accessors (`two_digit`, `four_digit`, `full_name`) parse the colon-delimited fields. `HlaTypingResult` maps each typed gene to a diplotype. Compatibility scoring (`hla_compatibility`) sums shared alleles across specified loci between donor and recipient, with a maximum of 2 matches per locus. Allele matching is by exact string comparison at the stored resolution. `parse_hla_typing` reads a simple tab-delimited format (one locus per line: GENE, ALLELE1, ALLELE2).

### Tumor Mutational Burden

`compute_tmb` counts somatic mutations in a variant set, divides by the assessed coding region size in megabases, and classifies the result using standard clinical thresholds: <6 mut/Mb (Low), 6-19 (Intermediate), >=20 (High, the FDA-approved threshold for pembrolizumab). The `count_indels` flag controls whether insertions and deletions are included (some targeted panels count only SNVs). Variants are classified by their `VariantType`: SNVs and MNVs always count; insertions, deletions, and complex variants are conditional.

### Microsatellite Instability

MSI analysis compares observed repeat counts at microsatellite loci against reference values. `bethesda_markers()` provides the standard 5 mononucleotide markers (BAT25, BAT26, NR21, NR24, MONO27) with their genomic coordinates and reference repeat counts. `call_msi` determines instability at each locus by checking whether the observed count differs from reference by more than a configurable shift threshold (typically 2-3 repeats). The final status follows standard criteria: MSS (0 unstable), MSI-Low (1 unstable), MSI-High (>=2 unstable or >=30% unstable fraction).
