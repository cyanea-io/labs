# cyanea-py API Reference

Python bindings for the Cyanea bioinformatics ecosystem via PyO3. Installable as `pip install cyanea` (via maturin).

## Status: Complete

Wraps all domain crates. 11 Python submodules: seq, align, stats, core, ml, chem, struct_bio, phylo, io, omics, sc.

## `cyanea.seq`

**Classes:**

| Class | Description |
|-------|-------------|
| `DnaSequence(data: bytes)` | Validated DNA sequence |
| `RnaSequence(data: bytes)` | Validated RNA sequence |
| `ProteinSequence(data: bytes)` | Validated protein sequence |
| `FastaStats` | FASTA file statistics (frozen) |
| `FastqStats` | FASTQ file statistics (frozen) |
| `FastqRecord` | Single FASTQ record (frozen) |

**DnaSequence methods:** `reverse_complement()`, `transcribe()`, `gc_content()`, `kmers(k)`, `translate()`, `__len__`, `__bytes__`, `__str__`, `__repr__`, `__eq__`

**RnaSequence methods:** `reverse_complement()`, `translate()`, `reverse_transcribe()`, `kmers(k)`

**ProteinSequence methods:** `molecular_weight()`, `kmers(k)`

**Functions:**

| Function | Description |
|----------|-------------|
| `fasta_stats(path) -> FastaStats` | Streaming FASTA statistics |
| `parse_fastq(path) -> list[FastqRecord]` | Parse all FASTQ records |
| `fastq_stats(path) -> FastqStats` | Streaming FASTQ statistics |
| `rna_fold_nussinov(seq) -> dict` | RNA folding (Nussinov algorithm) |
| `rna_fold_zuker(seq) -> dict` | RNA folding (Zuker MFE) |
| `protein_properties(seq) -> dict` | Protein physicochemical properties |
| `simulate_reads(ref_seq, num_reads, read_length, error_rate, seed) -> list` | Simulate Illumina reads |
| `codon_usage(seq) -> dict` | Codon usage table and CAI |
| `assembly_stats(lengths) -> dict` | Assembly QC (N50, L50, etc.) |

## `cyanea.align`

| Class/Function | Description |
|----------------|-------------|
| `AlignmentResult` | Frozen: `score`, `aligned_query`, `aligned_target`, `query_start/end`, `target_start/end`, `cigar_string`, `identity`, `matches`, `mismatches`, `gaps`, `length` |
| `align_dna(query, target, *, mode, match_score, mismatch_score, gap_open, gap_extend)` | DNA alignment |
| `align_protein(query, target, *, mode, matrix)` | Protein alignment |
| `align_batch(pairs, *, mode, ...)` | Batch DNA alignment |

## `cyanea.stats`

| Class/Function | Description |
|----------------|-------------|
| `DescriptiveStats` | Frozen, 15 fields |
| `TestResult` | Frozen: `statistic`, `p_value`, `degrees_of_freedom`, `method` |
| `describe(data) -> DescriptiveStats` | Descriptive statistics |
| `pearson(x, y) -> float` | Pearson correlation |
| `spearman(x, y) -> float` | Spearman correlation |
| `t_test(data, *, mu) -> TestResult` | One-sample t-test |
| `t_test_two_sample(x, y, *, equal_var) -> TestResult` | Two-sample t-test |
| `mann_whitney_u(x, y) -> TestResult` | Mann-Whitney U test |
| `bonferroni(p_values) -> list[float]` | Bonferroni correction |
| `benjamini_hochberg(p_values) -> list[float]` | BH FDR correction |
| `kaplan_meier(times, events) -> dict` | Kaplan-Meier estimator |
| `log_rank(times1, events1, times2, events2) -> dict` | Log-rank test |
| `cox_ph(times, events, covariates) -> dict` | Cox proportional hazards |
| `permutation_null(data, n_permutations) -> list` | Permutation null model |
| `bootstrap_null(data, n_bootstraps) -> list` | Bootstrap null model |
| `alpha_diversity(counts, method) -> float` | Alpha diversity |
| `fst_hudson(pop1, pop2) -> float` | Hudson Fst |
| `tajimas_d(genotypes) -> float` | Tajima's D |
| `gsea_preranked(ranked_genes, gene_set) -> dict` | GSEA preranked |
| `ora(gene_list, gene_set, background_size) -> dict` | Over-representation analysis |

## `cyanea.core`

| Function | Description |
|----------|-------------|
| `sha256(data: bytes) -> str` | SHA-256 hex digest |
| `sha256_file(path: str) -> str` | SHA-256 of a file |
| `zstd_compress(data: bytes, *, level) -> bytes` | Zstd compression |
| `zstd_decompress(data: bytes) -> bytes` | Zstd decompression |

## `cyanea.ml`

**Classes:**

| Class | Description |
|-------|-------------|
| `UmapResult` | `embedding`, `n_samples`, `n_components`, `n_epochs` |
| `PcaResult` | `components`, `explained_variance`, `explained_variance_ratio`, `transformed`, `mean` |
| `TsneResult` | `embedding`, `n_samples`, `n_components`, `kl_divergence` |
| `KMeansResult` | `centroids`, `labels`, `inertia`, `n_iter` |
| `GradientBoostedTrees` | GBDT model with `predict()` |

**Functions:**

| Function | Description |
|----------|-------------|
| `euclidean_distance(a, b) -> float` | Euclidean distance |
| `manhattan_distance(a, b) -> float` | Manhattan distance |
| `hamming_distance(a, b) -> int` | Hamming distance |
| `cosine_similarity(a, b) -> float` | Cosine similarity |
| `pairwise_distances(data, n_features, metric) -> list` | Pairwise distances |
| `umap(data, n_features, ...) -> UmapResult` | UMAP |
| `pca(data, n_features, n_components) -> PcaResult` | PCA |
| `tsne(data, n_features, ...) -> TsneResult` | t-SNE |
| `kmeans(data, n_features, n_clusters, ...) -> KMeansResult` | K-means |
| `random_forest_classify(train, labels, test, n_features, n_trees) -> list` | Random forest |
| `gbdt_regression(train, targets, test, n_features, n_trees, learning_rate) -> list` | GBDT |
| `confusion_matrix(actual, predicted) -> dict` | Confusion matrix |
| `roc_curve(scores, labels) -> dict` | ROC curve |
| `pr_curve(scores, labels) -> dict` | PR curve |
| `cross_validate(data, labels, n_features, k_folds) -> dict` | K-fold CV |
| `variance_threshold(data, n_features, threshold) -> list` | Variance threshold |
| `mutual_info(data, labels, n_features, k) -> list` | Mutual information |
| `hmm_viterbi(obs, states, trans, emit, init) -> list` | HMM Viterbi |

**NumPy variants** (feature `numpy`): `pairwise_distances_np`, `umap_np`, `pca_np`, `tsne_np`

## `cyanea.chem`

**Classes:**

| Class | Description |
|-------|-------------|
| `Molecule(smiles: str)` | Molecular graph from SMILES |
| `MolecularProperties` | Computed properties (frozen) |

**Molecule methods:** `atom_count()`, `bond_count()`, `molecular_formula()`, `molecular_weight()`, `canonical_smiles()`, `morgan_fingerprint(radius, n_bits)`, `has_substructure(pattern)`

**Functions:**

| Function | Description |
|----------|-------------|
| `parse_smiles(smiles) -> Molecule` | Parse SMILES |
| `molecular_properties(smiles) -> MolecularProperties` | Molecular properties |
| `tanimoto(smiles1, smiles2, radius, n_bits) -> float` | Tanimoto similarity |
| `canonical_smiles(smiles) -> str` | Canonical SMILES |

## `cyanea.struct_bio`

**Classes:**

| Class | Description |
|-------|-------------|
| `Structure` | 3D structure from PDB (frozen) |
| `SecondaryStructureAssignment` | Per-residue DSSP assignment |

**Structure methods:** `chain_count()`, `residue_count()`, `atom_count()`, `id()`, `secondary_structure()`

**Functions:**

| Function | Description |
|----------|-------------|
| `parse_pdb(content) -> Structure` | Parse PDB text |
| `rmsd(coords1, coords2) -> float` | RMSD |
| `kabsch_align(moving, target) -> (float, list)` | Kabsch superposition |

## `cyanea.phylo`

**Classes:**

| Class | Description |
|-------|-------------|
| `PhyloTree` | Rooted phylogenetic tree |

**PhyloTree methods:** `leaf_count()`, `leaf_names()`, `to_newick()`, `robinson_foulds(other)`

**Functions:**

| Function | Description |
|----------|-------------|
| `parse_newick(newick) -> PhyloTree` | Parse Newick string |
| `evolutionary_distance(seq1, seq2, model) -> float` | Evolutionary distance |
| `upgma(labels, matrix) -> PhyloTree` | UPGMA tree |
| `neighbor_joining(labels, matrix) -> PhyloTree` | NJ tree |
| `simulate_evolution(newick, root_seq, model, params) -> dict` | Sequence evolution |
| `simulate_coalescent(n_samples, pop_size, seed) -> PhyloTree` | Coalescent simulation |
| `simulate_coalescent_growth(n_samples, pop_size, growth_rate, seed) -> PhyloTree` | Coalescent with growth |

## `cyanea.io`

**Classes:**

| Class | Description |
|-------|-------------|
| `CsvInfo` | CSV metadata (frozen) |
| `VcfStats` | VCF summary (frozen) |
| `BedStats` | BED summary (frozen) |
| `Gff3Stats` | GFF3 summary (frozen) |
| `PySamRecord` | SAM/BAM record (frozen) |
| `PySamStats` | SAM/BAM summary (frozen) |

**Functions:**

| Function | Description |
|----------|-------------|
| `csv_info(path) -> CsvInfo` | CSV metadata |
| `vcf_stats(path) -> VcfStats` | VCF statistics |
| `bed_stats(path) -> BedStats` | BED statistics |
| `gff3_stats(path) -> Gff3Stats` | GFF3 statistics |
| `parse_sam(path) -> list[PySamRecord]` | Parse SAM |
| `sam_stats(path) -> PySamStats` | SAM statistics |
| `parse_bam(path) -> list[PySamRecord]` | Parse BAM |
| `bam_stats(path) -> PySamStats` | BAM statistics |
| `parse_blast_xml(path) -> list` | Parse BLAST XML |
| `write_bedgraph(regions, path)` | Write bedGraph |
| `parse_gfa(path) -> dict` | Parse GFA |
| `ncbi_efetch_url(db, id, rettype, retmode) -> str` | NCBI URL builder |
| `uniprot_url(accession) -> str` | UniProt URL builder |

## `cyanea.omics`

| Function | Description |
|----------|-------------|
| `merge_intervals(intervals) -> list` | Merge overlapping intervals |
| `intersect_intervals(a, b) -> list` | Interval intersection |
| `subtract_intervals(a, b) -> list` | Interval subtraction |
| `liftover(chain, chrom, start, end) -> dict` | Coordinate liftover |
| `annotate_variant(variant, transcripts) -> dict` | Variant effect prediction |
| `bisulfite_convert(seq) -> str` | Bisulfite conversion |
| `cpg_islands(seq) -> list` | CpG island detection |
| `morans_i(values, weights) -> float` | Moran's I |
| `gearys_c(values, weights) -> float` | Geary's C |

## `cyanea.sc`

| Function | Description |
|----------|-------------|
| `normalize_total(matrix, target_sum) -> list` | Library size normalization |
| `log1p_transform(matrix) -> list` | Log1p transformation |
| `highly_variable_genes(matrix, n_top) -> list` | HVG selection |
| `neighbors(data, n_neighbors) -> dict` | kNN graph construction |
| `leiden(adjacency, resolution) -> list` | Leiden clustering |
| `louvain(adjacency, resolution) -> list` | Louvain clustering |
| `diffusion_map(data, n_components) -> dict` | Diffusion map |
| `dpt(diffusion_map, root) -> list` | Diffusion pseudotime |
| `paga(adjacency, clusters) -> dict` | PAGA trajectory |
| `rank_genes(matrix, clusters, method) -> dict` | Marker gene detection |
| `harmony(data, batch_labels, n_components) -> list` | Harmony integration |
| `combat(matrix, batch_labels) -> list` | ComBat batch correction |
| `mnn(datasets) -> list` | MNN correction |

## Error Mapping

| Rust Error | Python Exception |
|------------|------------------|
| `CyaneaError::Io` | `IOError` |
| `CyaneaError::Parse`, `InvalidInput` | `ValueError` |
| `CyaneaError::Compression`, `Hash`, `Other` | `RuntimeError` |

## Dependencies

- `pyo3` 0.23 (extension-module)
- `numpy` 0.23 (optional, feature-gated)
- All domain crates: core, seq, io, align, stats, ml, chem, struct, phylo, omics
