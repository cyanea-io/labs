# cyanea-meta API Reference

Metagenomics analysis: taxonomy, profiling, diversity, composition, functional annotation, binning, and assembly QC.

## Public API

### Error types (`error.rs`)

| Type | Description |
|------|-------------|
| `MetaError` | Metagenomics-specific errors (InvalidInput, Parse, Taxonomy, EmptyInput, DimensionMismatch) |
| `Result<T>` | Alias for `std::result::Result<T, MetaError>` |

### Taxonomy (`taxonomy.rs`)

| Type | Description |
|------|-------------|
| `TaxonRank` | Enum: Domain, Phylum, Class, Order, Family, Genus, Species, Strain, Unranked |
| `TaxonNode` | Node in taxonomy tree: taxid, parent_id, rank, name |
| `TaxonomyDB` | In-memory taxonomy tree with k-mer index for classification |

**TaxonomyDB methods:**

| Method | Description |
|--------|-------------|
| `new(root_id: u32) -> Self` | Create taxonomy with given root taxon ID |
| `add_node(taxid, parent_id, rank, name) -> Result<()>` | Insert a node into the tree |
| `add_reference(sequence, taxid) -> Result<()>` | Index reference sequence k-mers for classification |
| `classify_sequence(read) -> Result<Option<u32>>` | Classify a read by k-mer LCA (Kraken-style) |
| `get_lineage(taxid) -> Result<Vec<u32>>` | Full lineage from taxid to root |
| `get_rank(taxid) -> Result<TaxonRank>` | Get rank for a taxon |
| `get_name(taxid) -> Result<String>` | Get name for a taxon |
| `lca(taxid1, taxid2) -> Result<u32>` | Lowest common ancestor of two taxa |

### Profiling (`profile.rs`)

| Type | Description |
|------|-------------|
| `TaxonomicProfile` | Taxon ID → (read count, relative abundance) mapping |

| Function | Description |
|----------|-------------|
| `profile_from_classifications(classifications: &[u32]) -> Result<TaxonomicProfile>` | Aggregate classified reads into a profile |
| `reestimate_abundance(profile, target_rank) -> Result<TaxonomicProfile>` | Bracken-style abundance re-estimation (placeholder) |
| `filter_profile(profile, min_abundance) -> Result<TaxonomicProfile>` | Remove taxa below abundance threshold |
| `normalize_profile(profile) -> Result<TaxonomicProfile>` | Normalize to relative abundances (sum to 1.0) |
| `merge_profiles(profiles: &[&TaxonomicProfile]) -> Result<TaxonomicProfile>` | Merge multiple profiles by summing counts |

**TaxonomicProfile methods:**

| Method | Description |
|--------|-------------|
| `total_count() -> u64` | Sum of all read counts |
| `richness() -> usize` | Number of distinct taxa |
| `as_vec() -> Vec<(u32, u64, f64)>` | Sorted (taxid, count, abundance) tuples |

### Diversity (`diversity.rs`)

| Type | Description |
|------|-------------|
| `AlphaDiversity` | Shannon, Simpson, Chao1, ACE, Fisher's alpha, observed species |
| `BetaDiversityMatrix` | Pairwise dissimilarity matrix with sample labels |

| Function | Description |
|----------|-------------|
| `alpha_diversity(counts: &[u64]) -> Result<AlphaDiversity>` | Compute all alpha diversity metrics |
| `beta_diversity_matrix(samples: &[&[u64]]) -> Result<BetaDiversityMatrix>` | Pairwise Bray-Curtis dissimilarity |
| `rarefaction_curve(counts, steps) -> Result<Vec<(usize, f64)>>` | Expected species richness at subsampled depths |
| `rarefy(counts, depth) -> Result<Vec<u64>>` | Subsample to even depth |

### Composition (`composition.rs`)

| Type | Description |
|------|-------------|
| `CompositionTransform` | CLR-transformed sample with original and transformed values |
| `DifferentialAbundanceResult` | Per-taxon: mean_diff, t_statistic, p_value, q_value |
| `AncemResult` | ANCOM result: taxon_index, w_statistic, is_significant |

| Function | Description |
|----------|-------------|
| `clr_transform(abundances: &[Vec<f64>]) -> Result<Vec<CompositionTransform>>` | Centered log-ratio transform |
| `ilr_transform(abundances: &[Vec<f64>]) -> Result<Vec<Vec<f64>>>` | Isometric log-ratio transform |
| `differential_abundance(group1, group2, taxa) -> Result<Vec<DifferentialAbundanceResult>>` | ALDEx2-style CLR + Welch's t-test |
| `ancom(group1, group2) -> Result<Vec<AncemResult>>` | ANCOM W-statistic differential abundance |

### Functional annotation (`functional.rs`)

| Type | Description |
|------|-------------|
| `FunctionalProfile` | Function ID → abundance mapping |

| Function | Description |
|----------|-------------|
| `map_to_functions(profile, taxon_functions) -> Result<FunctionalProfile>` | Map taxa to functional categories |
| `pathway_abundance(func_profile, pathway_map) -> Result<FunctionalProfile>` | Aggregate gene families to pathways |
| `functional_diversity(profile) -> Result<f64>` | Shannon entropy of functional categories |

### Binning (`binning.rs`)

| Type | Description |
|------|-------------|
| `Contig` | Contig with id, sequence, coverage, gc_content, length, tnf |
| `Bin` | Collection of contigs with completeness and contamination estimates |

| Function | Description |
|----------|-------------|
| `tetranucleotide_frequency(sequence) -> Result<Vec<f64>>` | 256-dimensional TNF vector |
| `bin_contigs(contigs, n_clusters) -> Result<Vec<Bin>>` | K-means binning by coverage + TNF |
| `assess_bin_quality(bin) -> Result<(f64, f64)>` | Estimate (completeness, contamination) |
| `filter_bins(bins, min_completeness, max_contamination) -> Result<Vec<Bin>>` | Quality-based filtering |

### Assembly QC (`assembly.rs`)

| Type | Description |
|------|-------------|
| `AssemblyStats` | N50, L50, N90, L90, total_length, contig_count, longest, shortest, gc_content, mean_length, median_length, aun |

| Function | Description |
|----------|-------------|
| `assembly_stats(contigs: &[&[u8]]) -> Result<AssemblyStats>` | Comprehensive assembly metrics |
| `nx_values(contigs, x) -> Result<(usize, usize)>` | Compute Nx and Lx for custom percentile |
| `coverage_depth(contigs, reads_per_contig) -> Result<Vec<f64>>` | Per-contig coverage estimate |
| `filter_contigs(contigs, min_length, min_coverage) -> Result<Vec<&[u8]>>` | Filter by length/coverage |
