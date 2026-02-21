# cyanea-wasm API Reference

WebAssembly bindings for browser-based bioinformatics. Wraps all domain crates into a JSON-based interface for JavaScript/TypeScript consumption.

## Status: Complete

223 tests across 11 source files. All bindings implemented with `wasm-bindgen` annotations (behind `wasm` feature flag). Functions accept simple types (strings, numbers) and return JSON strings.

## Error Handling (`error.rs`)

| Function | Description |
|----------|-------------|
| `wasm_ok<T: Serialize>(val) -> String` | Wrap success as `{"ok": value}` |
| `wasm_err(msg) -> String` | Wrap error as `{"error": message}` |
| `wasm_result<T: Serialize>(r) -> String` | Wrap `Result<T>` as JSON |

## Sequence Module (`seq.rs`)

| Function | Description |
|----------|-------------|
| `parse_fasta(data) -> String` | Parse FASTA from string, return JSON stats |
| `gc_content_json(seq) -> String` | Compute GC content |
| `reverse_complement(seq) -> String` | DNA reverse complement |
| `transcribe(seq) -> String` | DNA to RNA transcription |
| `translate(seq) -> String` | DNA to protein (standard codon table) |
| `validate(seq, alphabet) -> String` | Validate against "dna"/"rna"/"protein" |
| `parse_fastq(data) -> String` | Parse FASTQ records from string |
| `parse_paired_fastq(r1_data, r2_data, validation) -> String` | Parse paired FASTQ |
| `parse_interleaved_fastq(data, validation) -> String` | Parse interleaved FASTQ |
| `trim_fastq(data, config_json) -> String` | Trim single-end FASTQ records |
| `trim_paired_fastq(r1_data, r2_data, config_json, orphan_policy) -> String` | Trim paired FASTQ |
| `rna_fold_nussinov(seq) -> String` | RNA folding via Nussinov algorithm |
| `rna_fold_zuker(seq) -> String` | RNA folding via Zuker MFE |
| `protein_properties(seq) -> String` | Protein physicochemical properties |
| `simulate_reads(ref_seq, num_reads, read_length, error_rate, seed) -> String` | Simulate Illumina reads |
| `codon_usage(seq) -> String` | Codon usage table and CAI |
| `assembly_stats_json(lengths_json) -> String` | Assembly QC (N50, L50, etc.) |

## Alignment Module (`align.rs`)

| Function | Description |
|----------|-------------|
| `align_dna(query, target, mode) -> String` | DNA alignment with default scoring |
| `align_dna_custom(query, target, mode, match, mismatch, gap_open, gap_extend) -> String` | Custom scoring |
| `align_protein(query, target, mode, matrix) -> String` | Protein alignment |
| `align_batch(pairs_json, mode, match, mismatch, gap_open, gap_extend) -> String` | Batch alignment |
| `progressive_msa(sequences_json, match, mismatch, gap_open, gap_extend) -> String` | Progressive MSA |
| `poa_consensus(sequences_json, match, mismatch, gap_open, gap_extend) -> String` | POA consensus |
| `banded_align(query, target, mode, bandwidth, match, mismatch, gap_open, gap_extend) -> String` | Banded alignment |
| `parse_cigar(cigar_str) -> String` | Parse CIGAR string to JSON ops |
| `cigar_stats(cigar_str) -> String` | CIGAR statistics (identity, gaps, etc.) |

## Statistics Module (`stats.rs`)

| Function | Description |
|----------|-------------|
| `describe(data_json) -> String` | Descriptive statistics |
| `pearson(x_json, y_json) -> String` | Pearson correlation |
| `spearman(x_json, y_json) -> String` | Spearman rank correlation |
| `t_test(data_json, mu) -> String` | One-sample t-test |
| `t_test_two_sample(x_json, y_json, equal_var) -> String` | Two-sample t-test |
| `mann_whitney_u(x_json, y_json) -> String` | Mann-Whitney U test |
| `bonferroni(p_json) -> String` | Bonferroni correction |
| `benjamini_hochberg(p_json) -> String` | Benjamini-Hochberg FDR |
| `shannon_index(counts_json) -> String` | Shannon diversity index |
| `kaplan_meier(times_json, events_json) -> String` | Kaplan-Meier survival estimator |
| `log_rank_test(times1_json, events1_json, times2_json, events2_json) -> String` | Log-rank test |
| `cox_ph(times_json, events_json, covariates_json) -> String` | Cox proportional hazards |
| `permutation_null(data_json, n_permutations) -> String` | Permutation null model |
| `bootstrap_null(data_json, n_bootstraps) -> String` | Bootstrap null model |
| `wright_fisher(pop_size, n_generations, initial_freq) -> String` | Wright-Fisher simulation |
| `alpha_diversity(counts_json, method) -> String` | Alpha diversity (Shannon, Simpson, Chao1) |
| `beta_diversity_bray_curtis(a_json, b_json) -> String` | Bray-Curtis dissimilarity |
| `fst_hudson(pop1_json, pop2_json) -> String` | Hudson Fst |
| `tajimas_d(genotypes_json) -> String` | Tajima's D |
| `nucleotide_diversity(genotypes_json) -> String` | Nucleotide diversity (pi) |

## ML Module (`ml.rs`)

| Function | Description |
|----------|-------------|
| `kmer_count(seq, k) -> String` | K-mer counting |
| `euclidean_distance(a_json, b_json) -> String` | Euclidean distance |
| `manhattan_distance(a_json, b_json) -> String` | Manhattan distance |
| `hamming_distance(a, b) -> String` | Hamming distance |
| `cosine_similarity(a_json, b_json) -> String` | Cosine similarity |
| `umap(data_json, n_features, n_components, n_neighbors, min_dist, n_epochs, metric) -> String` | UMAP |
| `pca(data_json, n_features, n_components) -> String` | PCA |
| `tsne(data_json, n_features, n_components, perplexity, learning_rate, n_iter) -> String` | t-SNE |
| `kmeans(data_json, n_features, n_clusters, max_iter) -> String` | K-means clustering |
| `random_forest_classify(train_json, labels_json, test_json, n_features, n_trees) -> String` | Random forest |
| `random_forest_importance(train_json, labels_json, n_features, n_trees) -> String` | Feature importance |
| `gbdt_regression(train_json, targets_json, test_json, n_features, n_trees, learning_rate) -> String` | GBDT |
| `hmm_viterbi(obs_json, states, trans_json, emit_json, init_json) -> String` | HMM Viterbi |
| `hmm_forward(obs_json, states, trans_json, emit_json, init_json) -> String` | HMM forward |
| `confusion_matrix(actual_json, predicted_json) -> String` | Confusion matrix |
| `roc_curve_json(scores_json, labels_json) -> String` | ROC curve |
| `pr_curve_json(scores_json, labels_json) -> String` | PR curve |
| `cross_validate(data_json, labels_json, n_features, k_folds) -> String` | K-fold CV |
| `variance_threshold(data_json, n_features, threshold) -> String` | Variance threshold selection |
| `mutual_info_selection(data_json, labels_json, n_features, k) -> String` | Mutual information selection |

## Chemistry Module (`chem.rs`)

| Function | Description |
|----------|-------------|
| `smiles_properties(smiles) -> String` | Molecular properties from SMILES |
| `canonical(smiles) -> String` | Canonical SMILES |
| `smiles_fingerprint(smiles, radius, n_bits) -> String` | Morgan fingerprint |
| `tanimoto(smiles1, smiles2) -> String` | Tanimoto similarity |
| `smiles_substructure(molecule, pattern) -> String` | Substructure match |
| `parse_sdf_text(sdf_text) -> String` | Parse SDF V2000/V3000 text |
| `maccs_fingerprint(smiles) -> String` | MACCS 166-key fingerprint |
| `maccs_tanimoto(smiles1, smiles2) -> String` | MACCS Tanimoto similarity |

## Structural Biology Module (`struct_bio.rs`)

| Function | Description |
|----------|-------------|
| `pdb_info(pdb_text) -> String` | Parse PDB, return structure info |
| `pdb_secondary_structure(pdb_text) -> String` | DSSP secondary structure |
| `rmsd(coords1_json, coords2_json) -> String` | RMSD between coordinate sets |
| `contact_map(pdb_text, cutoff, mode) -> String` | Contact map (CA or all-atom) |
| `ramachandran_analysis(pdb_text) -> String` | Ramachandran phi/psi analysis |
| `parse_mmcif_text(mmcif_text) -> String` | Parse mmCIF text |
| `kabsch_align(moving_json, target_json) -> String` | Kabsch superposition |

## Phylogenetics Module (`phylo.rs`)

| Function | Description |
|----------|-------------|
| `newick_info(newick) -> String` | Parse Newick, return tree info |
| `evolutionary_distance(seq1, seq2, model) -> String` | Evolutionary distance (p/jc/k2p) |
| `build_upgma(labels_json, matrix_json) -> String` | UPGMA tree construction |
| `build_nj(labels_json, matrix_json) -> String` | Neighbor-Joining tree |
| `rf_distance(newick1, newick2) -> String` | Robinson-Foulds distance |
| `parse_nexus(nexus_text) -> String` | Parse NEXUS format |
| `write_nexus(newick) -> String` | Write NEXUS from Newick |
| `simulate_evolution(newick, root_seq, model, params_json) -> String` | Sequence evolution simulation |
| `simulate_coalescent(n_samples, pop_size, seed) -> String` | Coalescent simulation |

## I/O Module (`io.rs`)

| Function | Description |
|----------|-------------|
| `pileup_from_sam(sam_text) -> String` | Generate pileup from SAM text |
| `depth_from_sam(sam_text) -> String` | Per-position depth from SAM |
| `parse_vcf_text(vcf_text) -> String` | Parse VCF text |
| `parse_bed_text(bed_text) -> String` | Parse BED text |
| `parse_gff3_text(gff3_text) -> String` | Parse GFF3 text |
| `parse_blast_xml(xml_text) -> String` | Parse BLAST XML output |
| `write_bedgraph(regions_json) -> String` | Write bedGraph from JSON |
| `parse_gfa_text(gfa_text) -> String` | Parse GFA text |
| `ncbi_efetch_url(db, id, rettype, retmode) -> String` | NCBI Entrez fetch URL builder |
| `uniprot_url(accession) -> String` | UniProt REST URL builder |

## Omics Module (`omics.rs`)

| Function | Description |
|----------|-------------|
| `merge_intervals(intervals_json) -> String` | Merge overlapping intervals |
| `intersect_intervals(a_json, b_json) -> String` | Interval intersection |
| `subtract_intervals(a_json, b_json) -> String` | Interval subtraction |
| `complement_intervals(intervals_json, chrom_size) -> String` | Interval complement |
| `closest_intervals(a_json, b_json) -> String` | Closest interval queries |
| `window_intervals(chrom, start, end, size, step) -> String` | Window/tile generation |
| `jaccard_intervals(a_json, b_json) -> String` | Jaccard similarity |
| `liftover_interval(chain_json, chrom, start, end) -> String` | Coordinate liftover |
| `annotate_variant(variant_json, transcripts_json) -> String` | Variant effect prediction |
| `cbs_segment(values_json) -> String` | Circular binary segmentation |
| `bisulfite_convert(seq) -> String` | Bisulfite conversion |
| `cpg_islands(seq) -> String` | CpG island detection |
| `morans_i(values_json, weights_json) -> String` | Moran's I spatial autocorrelation |
| `gearys_c(values_json, weights_json) -> String` | Geary's C spatial autocorrelation |

## Core Utilities Module (`core_utils.rs`)

| Function | Description |
|----------|-------------|
| `sha256(data) -> String` | SHA-256 hex digest |
| `zstd_compress(data, level) -> String` | Zstd compression |
| `zstd_decompress(data_json) -> String` | Zstd decompression |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `wasm` | No | Enables `wasm-bindgen` annotations on all public functions |

## Dependencies

- `cyanea-core`, `cyanea-seq`, `cyanea-io`, `cyanea-align`, `cyanea-stats`, `cyanea-ml`, `cyanea-chem`, `cyanea-struct`, `cyanea-phylo`, `cyanea-omics`
- `serde`, `serde_json` -- JSON serialization
- `wasm-bindgen` (optional, behind `wasm` feature)
