# cyanea-epi API Reference

Epigenomics analysis: peak calling, signal tracks, motif discovery, chromatin states, differential binding, nucleosome positioning, and ATAC-seq QC.

## Public API

### Error types (`error.rs`)

| Type | Description |
|------|-------------|
| `EpiError` | Epigenomics-specific errors (InvalidInput, Parse, EmptyInput, DimensionMismatch, InsufficientData) |
| `Result<T>` | Alias for `std::result::Result<T, EpiError>` |

### Peak calling (`peaks.rs`)

| Type | Description |
|------|-------------|
| `Peak` | Called peak: chrom, start, end, summit, score, p_value, q_value, fold_enrichment, name |
| `PeakCallParams` | Parameters: bandwidth, q_value_cutoff, min_length, max_gap, fragment_size |
| `PeakSet` | Collection of peaks with set operations |
| `PeakStats` | Summary: count, total_bp, median_width, mean_fold_enrichment |

| Function | Description |
|----------|-------------|
| `call_peaks(treatment, control, params) -> Result<Vec<Peak>>` | MACS2-style narrow peak calling with Poisson test and BH correction |
| `call_broad_peaks(treatment, control, params) -> Result<Vec<Peak>>` | Broad peak calling for histone marks (two-pass seed + link) |

**Peak methods:**

| Method | Description |
|--------|-------------|
| `len() -> u64` | Peak width in bp |
| `contains(position) -> bool` | Position within peak |
| `overlaps(other) -> bool` | Two peaks overlap |

**PeakSet methods:**

| Method | Description |
|--------|-------------|
| `new(peaks) -> Self` | Create from peak vector |
| `merge(max_gap) -> PeakSet` | Merge overlapping/nearby peaks |
| `intersect(other) -> PeakSet` | Peaks overlapping another set |
| `subtract(other) -> PeakSet` | Peaks not overlapping another set |
| `filter_by_score(min_score) -> PeakSet` | Filter by score threshold |
| `stats() -> PeakStats` | Summary statistics |

### Signal pileup (`pileup.rs`)

| Type | Description |
|------|-------------|
| `TagPileup` | Per-chromosome coverage vectors with total read count |

| Function | Description |
|----------|-------------|
| `build_pileup(reads, fragment_size) -> TagPileup` | Build from aligned reads (chrom, start, frag_size) |
| `normalize_pileup(pileup, method) -> Result<TagPileup>` | Normalize: "cpm" or "rpkm" |
| `smooth_pileup(pileup, bandwidth) -> TagPileup` | Gaussian kernel smoothing |
| `pileup_correlation(p1, p2) -> Option<f64>` | Pearson correlation between two pileups |
| `fingerprint(pileup) -> Vec<(f64, f64)>` | Cumulative enrichment plot data (deepTools-style) |

**TagPileup methods:**

| Method | Description |
|--------|-------------|
| `get(chrom, position) -> u32` | Coverage at a position |
| `mean_coverage(chrom) -> f64` | Mean coverage for a chromosome |

### Motif discovery (`motifs.rs`)

| Type | Description |
|------|-------------|
| `Motif` | Name, PWM (Vec<[f64; 4]>), consensus sequence |
| `MotifMatch` | Position, strand (+/-), score, p_value, matched_seq |
| `DiscoveryParams` | motif_width, n_motifs, n_sequences, background_freq |

| Function | Description |
|----------|-------------|
| `scan_sequence(seq, motif, threshold) -> Vec<MotifMatch>` | Scan sequence for PWM matches (both strands) |
| `discover_motifs(sequences, params) -> Result<Vec<Motif>>` | K-mer enrichment-based motif discovery |
| `parse_meme(content) -> Result<Vec<Motif>>` | Parse MEME format motif file |
| `write_meme(motifs) -> String` | Write motifs in MEME format |
| `compare_motifs(m1, m2) -> f64` | Pairwise similarity (PCC of aligned columns) |
| `motif_enrichment(motif, target_seqs, bg_seqs) -> Result<(f64, f64)>` | Fisher's exact test for enrichment |

### Chromatin states (`chromatin.rs`)

| Type | Description |
|------|-------------|
| `ChromatinState` | State: id, name, emission_probs, color |
| `ChromHMMModel` | Learned model: states, transition_matrix, marks |
| `ChromHMMParams` | n_states, bin_size, max_iter, convergence |
| `ChromatinSegmentation` | Regions with assigned states |

| Function | Description |
|----------|-------------|
| `learn_chromatin_states(marks_matrix, params) -> Result<ChromHMMModel>` | ChromHMM-like EM learning (multivariate Bernoulli HMM) |
| `segment_genome(model, marks_matrix) -> Result<ChromatinSegmentation>` | Viterbi decoding of state assignments |
| `state_enrichment(segmentation, annotations) -> Vec<Vec<f64>>` | Enrichment of states at genomic annotations |

### Differential binding (`differential.rs`)

| Type | Description |
|------|-------------|
| `DiffResult` | region, log2_fc, p_value, q_value, mean_count_1, mean_count_2 |

| Function | Description |
|----------|-------------|
| `differential_peaks(counts, conditions, method) -> Result<Vec<DiffResult>>` | DESeq2-style differential analysis: size factors + Welch's t-test + BH correction |
| `count_reads_in_peaks(peaks, reads) -> Vec<u64>` | Count reads overlapping each peak |
| `ma_plot_data(results) -> Vec<(f64, f64)>` | (mean_expression, log2FC) pairs for MA plot |

### Nucleosome positioning (`nucleosome.rs`)

| Type | Description |
|------|-------------|
| `NucleosomePosition` | center, occupancy, fuzziness, score |
| `NucleosomeParams` | fragment_lower, fragment_upper, smooth_bandwidth, min_occupancy |

| Function | Description |
|----------|-------------|
| `call_nucleosomes(pileup, params) -> Vec<NucleosomePosition>` | Detect nucleosome positions from MNase/ATAC signal |
| `nfr_score(positions, tss_list) -> f64` | Nucleosome-free region score around TSSs |
| `periodicity(coverage) -> f64` | Detect ~147bp / ~200bp periodicity via autocorrelation |

### ATAC-seq QC (`accessibility.rs`)

| Type | Description |
|------|-------------|
| `InsertSizeMetrics` | nfr_ratio, mono_nuc_ratio, median_size, periodicity |
| `AtacQcResult` | tss_enrichment, frip, insert_metrics, total_reads |

| Function | Description |
|----------|-------------|
| `tss_enrichment(pileup, tss_positions, flank) -> f64` | TSS enrichment score (ENCODE QC metric) |
| `fragment_size_distribution(sizes, bin_size) -> Vec<(u64, u64)>` | Histogram of fragment sizes |
| `insert_size_metrics(sizes) -> InsertSizeMetrics` | NFR ratio, mono-nucleosome ratio, periodicity |
| `frip(reads_in_peaks, total_reads) -> f64` | Fraction of reads in peaks |
| `atacqc(pileup, peaks, tss, fragment_sizes) -> Result<AtacQcResult>` | Comprehensive ATAC-seq QC |
