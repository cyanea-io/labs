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

## Hi-C / 3D Genome Analysis

`hic.rs` provides data structures and algorithms for chromosome conformation capture (Hi-C) analysis: contact matrices, TAD calling, A/B compartment detection, chromatin loop calling, and file format I/O.

### ContactMatrix

The `ContactMatrix` stores Hi-C contacts as a dense `Vec<Vec<f64>>` matrix. For intra-chromosomal (cis) matrices, `set()` and `add_contact()` automatically maintain symmetry (setting both `[i][j]` and `[j][i]`). The matrix supports both square (intra-chromosomal) and rectangular (inter-chromosomal) layouts via separate `n_bins1`/`n_bins2` dimensions.

Two balancing methods normalize systematic biases:

- **ICE** (Imakaev et al. 2012): Iteratively scales rows and columns so that all marginal sums converge to a uniform value. Each iteration computes row sums, divides both the row and corresponding column by `sqrt(row_sum)`, and accumulates the bias factor. Converges when the maximum relative change in bias drops below the tolerance.

- **KR** (Knight-Ruiz): Computes a target marginal (`sqrt(mean_marginal)`) and iteratively adjusts per-bin weights so that the weighted row sums match the target. After convergence, applies the weight products `w_i * w_j` to the matrix. Better suited for sparse matrices than ICE.

**Observed/expected** (`observed_expected`): For each diagonal offset `d`, computes the mean contact count across all entries at that distance, then divides each cell by its expected value. This removes the distance-decay signal, highlighting structural features like TADs and compartments.

### TAD Calling

TADs are called using the insulation score method (Crane et al. 2015):

1. **Insulation score**: For each bin `i`, compute the mean contact count in the square window `[i-w, i) x [i, i+w)` on the contact matrix, then take `ln()` (clamped to -10 to handle zeros). Low insulation indicates a boundary where contacts do not cross.

2. **Z-score normalization**: Insulation scores are standardized to z-scores across all bins.

3. **Boundary detection**: Local minima in the z-scored insulation profile that fall below `boundary_threshold` are identified as TAD boundaries. A point is a local minimum if its z-score is <= both neighbors.

4. **TAD construction**: Consecutive boundary pairs define TADs. The matrix start (bin 0) and end (bin n) are implicitly included as boundaries. TADs smaller than `min_tad_size` are discarded.

### A/B Compartments

Compartment calling follows the Lieberman-Aiden method:

1. Compute the O/E matrix (distance normalization).
2. Compute the Pearson correlation matrix of O/E rows -- bins with similar interaction patterns have high correlation.
3. Extract the first eigenvector (PC1) via power iteration (100 iterations) on the correlation matrix.
4. Assign compartments: positive PC1 = A (active/open), negative PC1 = B (inactive/closed).

If GC content (or gene density) is provided, the eigenvector is oriented so that A compartments correlate positively with GC content, resolving the sign ambiguity inherent in eigenvector extraction.

### Loop Calling

Loop calling uses a HiCCUPS-style local enrichment approach:

1. For each pixel `(i, j)` within the allowed distance range `[min_distance, max_distance]`, compute the donut background: the mean contact count in a ring of radius `background_window` around the pixel, excluding the 3x3 center.

2. Compute enrichment as `observed / expected`. Filter by `min_enrichment`.

3. Compute a p-value using a normal approximation to the Poisson distribution: `z = (observed - expected) / sqrt(expected)`, then `p = erfc(z / sqrt(2)) / 2`. The complementary error function uses the Abramowitz-Stegun polynomial approximation.

4. Loops passing the p-value threshold are collected and sorted by enrichment (descending).

### File Format I/O

- **parse_cool_text**: Parses tab-delimited contact exports from cooler tools (7 columns: chrom1, start1, end1, chrom2, start2, end2, count). Genomic positions are converted to bin indices by dividing by the resolution. Lines starting with `#` are skipped.

- **parse_pairs**: Parses the 4DN consortium pairs format (readID, chr1, pos1, chr2, pos2, strand1, strand2). Returns `(chrom1, pos1, chrom2, pos2)` tuples. Header lines (starting with `#`) are skipped.

- **contacts_to_matrix**: Builds a dense `ContactMatrix` from a slice of `SparseContact` records, calling `add_contact` for each entry (which handles symmetry).

- **write_pairs**: Serializes sparse contacts back to 4DN pairs text format, converting bin indices to genomic positions by multiplying by the resolution.

## CNV and Methylation

- **CBS** (Circular Binary Segmentation): Recursively finds the split point that maximizes the t-statistic between segment means, using permutation testing for significance. Segments are then optionally merged based on log2-ratio similarity.
- **Methylation**: CpG sites carry read-level methylation counts. DMR detection uses a sliding window approach, testing beta-value differences between groups with Welch's t-test. CpG island detection uses the Gardiner-Garden & Frommer criteria (GC% > 50%, observed/expected CpG ratio > 0.6, length > 200bp).

## Microarray Analysis

`microarray.rs` provides normalization, summarization, and differential analysis for expression and methylation microarray data.

### Expression Microarray Pipeline

- **Quantile normalization**: Sorts each sample's values, replaces each rank with the mean across samples at that rank, then restores original ordering. This forces all samples to have identical value distributions, removing systematic technical variation while preserving relative differences within samples.

- **RMA** (Robust Multi-array Average): A two-step pipeline combining background correction (log2 transform with offset to stabilize low-intensity values) followed by quantile normalization. This is the standard preprocessing for Affymetrix GeneChip data.

- **Median polish**: Tukey's iterative median polish for probe set summarization. Alternately subtracts row medians and column medians until convergence (or `max_iter`), decomposing the probe-by-sample matrix into an overall effect, row effects (probe affinities), and column effects (sample expression levels). The column effects become the summarized expression values per sample.

- **limma differential expression**: Implements moderated t-statistics with empirical Bayes variance moderation (Smyth 2004). For each gene, fits a linear model comparing two groups, then shrinks per-gene variance estimates toward a common prior estimated across all genes. This borrowing of information across genes improves power for small sample sizes. P-values are computed from the moderated t-statistic and corrected for multiple testing using Benjamini-Hochberg.

### Methylation Microarray Pipeline

- **Beta values**: Computed as `M / (M + U + offset)` where M and U are methylated and unmethylated signal intensities. The offset (typically 100) prevents division instability when both signals are near zero. Beta values range from 0 (unmethylated) to 1 (fully methylated) and are biologically interpretable but heteroscedastic.

- **M-values**: The logit transform `log2(beta / (1 - beta))` of beta values. M-values have better statistical properties (approximately homoscedastic) and are preferred for differential analysis, though less interpretable. Conversion functions are provided in both directions.

- **SWAN normalization**: Subset-quantile Within Array Normalization corrects the technical bias between Infinium I and Infinium II probe designs on Illumina BeadChip arrays. TypeI probes use two beads (methylated and unmethylated channels), while TypeII probes use a single bead, leading to systematic intensity differences. SWAN performs separate quantile normalization within each design type to equalize their distributions.

- **Differential methylation**: Converts beta values to M-values, then applies the same empirical Bayes moderated t-test framework as limma. Reports delta-beta (difference in mean beta between groups) for biological interpretability alongside the statistical results from M-value analysis. BH correction controls the false discovery rate.

## Spatial Transcriptomics

- **Neighbor graphs**: Delaunay triangulation uses divide-and-conquer on the 2D point set. k-nearest neighbor graphs compute all pairwise distances and select the k closest per point.
- **Moran's I**: Computed as a weighted covariance of deviations from the mean, normalized by the variance, with analytic p-value from the normal approximation.
- **Ligand-receptor**: For each ligand-receptor pair, compute the mean product of ligand expression in sending cells and receptor expression in receiving neighbors. Significance by permuting cell labels.

### Spatial Platforms

`VisiumData`, `MerfishData`, and `SlideseqData` are typed containers for the three major spatial transcriptomics technologies. Each constructor validates dimension consistency (counts rows vs. barcodes/cell_ids, counts columns vs. genes, coordinate array lengths). Platform-specific fields capture technology details: Visium carries `array_positions` (hex grid row/col), `in_tissue` flags, and `VisiumScaleFactors` for coordinate system conversion; MERFISH supports optional `cell_volumes` (3D segmentation), `fov_ids`, and `blank_counts` for error rate estimation; Slide-seq stores bead coordinates in microns.

The `*_to_spatial_points` converter functions bridge platform containers to the core `SpatialPoint` type used by the spatial analysis functions (`morans_i`, `knn_graph`, etc.), mapping each platform's native coordinates to `(x, y, index)` tuples.

### Spatial Segmentation

Three complementary segmentation strategies for assigning transcripts or pixels to cells:

1. **Voronoi segmentation**: Assigns each transcript to the nearest seed (nucleus centroid) via brute-force nearest-neighbor search. An optional `max_radius` parameter clips assignments beyond a threshold, leaving distant transcripts unassigned. Cell boundaries are approximated by the convex hull of assigned transcripts (Andrew's monotone chain algorithm), with area computed by the shoelace formula.

2. **Nucleus expansion**: A simplified Baysor/Cellpose-style approach. For each seed, the effective radius is `min(max_radius, (nearest_neighbor_distance - min_gap) / 2)`, preventing cell overlap. Transcripts are assigned to the nearest seed only if within that seed's effective radius. This naturally handles variable cell density -- seeds in sparse regions expand to `max_radius`, while closely packed seeds shrink their territories.

3. **Watershed grid**: Operates on a 2D intensity image (e.g., DAPI). Seeds are placed at detected nucleus positions. A BFS flood-fill processes pixels in order of decreasing intensity, propagating labels from seeds outward via 4-connected neighbors. Returns a label grid where each pixel carries a cell ID (1-based; no explicit boundary/background detection in this simplified variant).

### Spatial Domains

Domain detection uses a spatially-aware k-means variant:

1. **Normalization**: Expression is z-scored per gene; spatial coordinates are min-max normalized to [0, 1].
2. **Combined distance**: `d = (1 - alpha) * expr_distance + alpha * spatial_distance`, where `alpha = spatial_weight` (0 = expression only, 1 = space only). Expression distance is Euclidean on z-scored values, normalized by `sqrt(n_genes)`.
3. **K selection**: If `n_domains` is not specified, an automatic heuristic uses `sqrt(n_spots / 2)` clamped to [2, 20].
4. **Iteration**: Standard k-means assignment + update loop with early convergence termination.

**HMRF smoothing** refines domain labels post-clustering using Iterated Conditional Modes (ICM). The energy function for assigning spot `i` to domain `k` is: `E = ||x_i - mu_k||^2 + beta * (number of neighbors with label != k)`. The `beta` parameter controls the spatial smoothness penalty -- higher values produce more spatially coherent domains at the cost of expression fidelity. Cluster means are recomputed each iteration. Convergence is detected when no labels change.

**SVG detection** (`find_spatially_variable_genes`) computes Moran's I per gene on a pre-built neighbor graph, using the standard formula `I = (N/W) * sum(w_ij * z_i * z_j) / sum(z_i^2)`. Significance is assessed by permutation testing (Fisher-Yates shuffle of expression values, xorshift64 PRNG). Multiple testing correction uses Benjamini-Hochberg with monotonicity enforcement (step-up procedure). Results are sorted by adjusted p-value.

### Spatial Cell-Cell Communication

Implements a CellChat-style analysis pipeline with spatial awareness:

1. **L-R database**: `demo_lr_database()` provides 12 curated ligand-receptor pairs across 8 signaling pathways (CXCL, CCL, NOTCH, WNT, VEGF, TGFb, COLLAGEN, FN1, PD-L1, EGF, HGF). Multi-subunit complexes are supported (e.g., WNT3A/FZD1+LRP6, TGFB1/TGFBR1+TGFBR2). Three interaction types: Secreted Signaling, Cell-Cell Contact, ECM-Receptor.

2. **Communication probability**: For each L-R pair and each (source_type, target_type) pair, the probability is `P = mean(L_i * R_j * w_ij)` across all source cell `i` and target cell `j` pairs (excluding self-interactions). `L_i` and `R_j` are geometric means of subunit expression (handles multi-subunit complexes -- if any subunit is zero, the geometric mean is zero). `w_ij = exp(-d^2 / (2 * sigma^2))` is a Gaussian spatial weight that downweights distant cell pairs.

3. **Significance**: Permutation testing shuffles cell type labels (not expression values), recomputing the probability under the null. The p-value is `(count_permuted >= observed + 1) / (n_permutations + 1)`.

4. **Pathway aggregation**: `aggregate_pathways` groups L-R pair results by pathway, summing probabilities and counting significant pairs (below `p_threshold`).

### Spatial Deconvolution

Two approaches for estimating cell type composition at each spatial spot:

1. **NNLS deconvolution** (`nnls_deconvolve`): Given cell type reference signatures (gene sets with expected expression weights), estimates per-spot proportions via non-negative least squares using coordinate descent. For each cell type `t`, the optimal proportion is `dot(residual, reference_t) / dot(reference_t, reference_t)`, clamped to >= 0. Iterates until convergence (max 200 iterations). Proportions are normalized to sum to 1. Residual error (RMSE) quantifies reconstruction quality. Genes present in at least one signature are used; genes absent from a signature get weight 0.

2. **Enrichment scoring** (`score_enrichment`): Per-spot, per-signature scoring based on mean z-scored expression of signature genes. Z-scores are computed per gene across all spots. Significance is assessed by comparing the observed score against random gene sets of the same size (permutation test, xorshift64 PRNG). Unlike deconvolution, enrichment does not estimate proportions -- it provides a relative enrichment score and p-value, useful for identifying hotspots of cell type activity.

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
