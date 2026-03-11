# cyanea-omics Usage Guide

Practical examples for working with genomic data structures, single-cell analysis, variant annotation, spatial transcriptomics (including platform-specific containers, segmentation, domain detection, cell-cell communication, and deconvolution), Hi-C / 3D genome analysis, ACMG variant classification, pharmacogenomics, and clinical genomics.

## Genomic Intervals and Interval Sets

```rust
use cyanea_omics::{GenomicInterval, IntervalSet, Strand};

// Create intervals (0-based, half-open)
let a = GenomicInterval::new("chr1", 100, 200, Strand::Forward);
let b = GenomicInterval::new("chr1", 150, 300, Strand::Forward);
let c = GenomicInterval::new("chr1", 500, 600, Strand::Forward);

// Build an interval set and merge overlapping intervals
let set = IntervalSet::from(vec![a, b, c]);
let merged = set.merge();
// merged contains: [chr1:100-300, chr1:500-600]

// Query overlaps
let query = GenomicInterval::new("chr1", 250, 550, Strand::Forward);
let overlapping = set.overlapping(&query);
```

## Interval Tree for Fast Overlap Queries

```rust
use cyanea_omics::interval_tree::{Interval, IntervalTree};

// Build a tree of gene annotations
let intervals = vec![
    Interval::new(1000, 2000, "gene_A"),
    Interval::new(1500, 3000, "gene_B"),
    Interval::new(5000, 6000, "gene_C"),
];
let tree = IntervalTree::from_unsorted(intervals);

// O(log n + k) overlap query
let hits = tree.query(1800, 2500);
for hit in &hits {
    println!("{}: {}-{}", hit.data, hit.start, hit.end);
}

// Count without allocating
let n = tree.count(1000, 3000);
```

## Genome Arithmetic

BEDTools-style operations on interval sets.

```rust
use cyanea_omics::genome_arithmetic::*;
use cyanea_omics::{GenomicInterval, Strand};

let peaks = vec![
    GenomicInterval::new("chr1", 100, 500, Strand::Unstranded),
    GenomicInterval::new("chr1", 800, 1200, Strand::Unstranded),
];
let genes = vec![
    GenomicInterval::new("chr1", 300, 900, Strand::Forward),
];

// Intersect: regions present in both sets
let overlap = intersect(&peaks, &genes, StrandMode::Ignore);

// Subtract: peaks minus gene regions
let non_genic = subtract(&peaks, &genes, StrandMode::Ignore);

// Union: merge all overlapping intervals
let all = union(&[peaks.clone(), genes.clone()].concat());

// Complement against genome
let genome = genome_info(&[("chr1", 10000)]);
let gaps = complement(&peaks, &genome);

// Jaccard similarity
let stats = jaccard(&peaks, &genes).unwrap();
println!("Jaccard: {:.3}", stats.jaccard);

// Find closest gene for each peak
let closest = closest(&peaks, &genes, StrandMode::Ignore);
```

## Expression Matrices (Dense and Sparse)

```rust
use cyanea_omics::{ExpressionMatrix, SparseMatrix};

// Dense matrix: 3 genes x 4 samples
let data = vec![
    vec![10.0, 20.0, 15.0, 25.0],  // gene1
    vec![5.0,  8.0,  3.0,  12.0],  // gene2
    vec![0.0,  1.0,  0.0,  2.0],   // gene3
];
let genes = vec!["TP53".into(), "BRCA1".into(), "MYC".into()];
let samples = vec!["S1".into(), "S2".into(), "S3".into(), "S4".into()];
let matrix = ExpressionMatrix::new(data, genes, samples).unwrap();

println!("Shape: {:?}", matrix.shape()); // (3, 4)

// Sparse matrix for single-cell data
let mut sparse = SparseMatrix::new(1000, 20000);
sparse.insert(0, 42, 3.5);
sparse.insert(0, 100, 1.2);
println!("Density: {:.4}", sparse.density());

// Convert between formats
let csr = sparse.to_csr();
let dense = sparse.to_dense();
```

## AnnData Container

```rust
use cyanea_omics::single_cell::{AnnData, MatrixData, ColumnData};

// Create from dense data: 100 cells x 500 genes
let x = MatrixData::Dense(vec![vec![0.0; 500]; 100]);
let obs_names: Vec<String> = (0..100).map(|i| format!("cell_{}", i)).collect();
let var_names: Vec<String> = (0..500).map(|i| format!("gene_{}", i)).collect();

let mut adata = AnnData::new(x, obs_names, var_names).unwrap();

// Add cell-level annotations
adata.add_obs("cell_type".into(), vec!["T-cell".into(); 100]);
adata.add_obs_numeric("n_counts".into(), vec![5000.0; 100]);

// Add gene-level annotations
adata.add_var_numeric("mean_expression".into(), vec![1.5; 500]);

// Add embeddings
adata.add_obsm("X_pca".into(), vec![vec![0.0; 50]; 100]);

// Subset
let subset = adata.subset_obs(&[0, 1, 2, 3, 4]).unwrap();
println!("Subset: {} cells", subset.n_obs());

// QC metrics
let qc = adata.qc_metrics();
```

## Variant Annotation Pipeline

```rust
use cyanea_omics::annotation::{Gene, Transcript, Exon, GeneType};
use cyanea_omics::variant::Variant;
use cyanea_omics::variant_annotation::{TranscriptModel, annotate_variant};
use cyanea_omics::genomic::Strand;

// Build transcript model from gene annotations
let genes = vec![/* genes loaded from GFF/GTF */];
let model = TranscriptModel::from_genes(&genes);

// Annotate a variant
let variant = Variant {
    chrom: "chr17".into(),
    position: 7674220,
    ref_allele: "C".into(),
    alt_alleles: vec!["T".into()],
    quality: Some(99.0),
    filter: Default::default(),
    info: Default::default(),
};

let annotations = annotate_variant(&variant, &model);
for ann in &annotations {
    println!(
        "{} {} {:?} {}",
        ann.gene, ann.transcript, ann.consequence, ann.amino_acid_change
    );
}
```

## Single-Cell Analysis Pipeline

Full Scanpy-style workflow: normalize, select HVGs, cluster, infer trajectory, find markers, and integrate batches.

```rust
use cyanea_omics::single_cell::{AnnData, MatrixData};
use cyanea_omics::sc_preprocess::*;
use cyanea_omics::sc_cluster::*;
use cyanea_omics::sc_trajectory::*;
use cyanea_omics::sc_markers::*;
use cyanea_omics::sc_integrate::*;

// Start with raw count data
let mut adata = /* load or create AnnData */;

// 1. Normalize: library-size normalization + log1p
normalize_total(&mut adata, NormalizeConfig {
    target_sum: Some(10_000.0),
    log_transform: true,
    save_raw: true,
}).unwrap();

// 2. Select highly variable genes
highly_variable_genes(&mut adata, HvgConfig {
    n_top_genes: 2000,
    method: HvgMethod::SeuratV3,
    ..Default::default()
}).unwrap();

// 3. Build neighbor graph from PCA
neighbors(&mut adata, NeighborsConfig {
    n_neighbors: 15,
    n_pcs: Some(50),
    ..Default::default()
}).unwrap();

// 4. Cluster with Leiden
leiden(&mut adata, ClusterConfig {
    resolution: 1.0,
    key_added: "leiden".into(),
    ..Default::default()
}).unwrap();

// 5. Diffusion pseudotime
diffusion_map(&mut adata, DiffusionConfig {
    n_components: 10,
    ..Default::default()
}).unwrap();

dpt(&mut adata, DptConfig {
    root_cell: 0,
    ..Default::default()
}).unwrap();

// 6. Find marker genes per cluster
let markers = rank_genes_groups(&adata, MarkerConfig {
    method: MarkerMethod::Wilcoxon,
    cluster_key: "leiden".into(),
    n_genes: 100,
    ..Default::default()
}).unwrap();

for (cluster, genes) in &markers.markers {
    println!("Cluster {}: top gene = {} (log2FC = {:.2})",
        cluster, genes[0].gene_name, genes[0].log2_fold_change);
}

// 7. Batch integration with Harmony
harmony(&mut adata, HarmonyConfig {
    batch_key: "batch".into(),
    ..Default::default()
}).unwrap();

// Evaluate integration quality
let metrics = integration_metrics(&adata, MetricsConfig {
    batch_key: "batch".into(),
    label_key: "leiden".into(),
    ..Default::default()
}).unwrap();
println!("kBET accept: {:.2}, iLISI: {:.2}", metrics.kbet_accept_rate, metrics.mean_ilisi);
```

## Spatial Transcriptomics

```rust
use cyanea_omics::spatial::*;

// Build spatial neighbor graph from spot coordinates
let points: Vec<SpatialPoint> = coordinates.iter().enumerate()
    .map(|(i, &(x, y))| SpatialPoint { x, y, index: i })
    .collect();
let graph = knn_graph(&points, 6).unwrap();

// Moran's I: test spatial autocorrelation of a gene
let gene_expr = vec![/* expression values per spot */];
let result = morans_i(&gene_expr, &graph).unwrap();
println!("Moran's I = {:.3}, p = {:.4}", result.morans_i, result.p_value);

// Geary's C
let gc = gearys_c(&gene_expr, &graph).unwrap();

// Ligand-receptor interaction analysis
let pairs = vec![("LIG1".into(), "REC1".into())];
let lr = ligand_receptor(&expression_matrix, &graph, &pairs, 1000).unwrap();
```

## Spatial Platform Data: Visium, MERFISH, Slide-seq

```rust
use cyanea_omics::spatial_platforms::*;

// Load Visium data with validation
let visium = VisiumData::new(
    gene_names,
    barcodes,
    count_matrix,      // counts[spot][gene]
    spot_coordinates,  // (x, y) pixel coords
    array_positions,   // (row, col) on Visium grid
    in_tissue_flags,   // which spots are under tissue
).unwrap();

// QC: filter to tissue spots and inspect
let tissue_only = visium.filter_tissue();
println!("{} tissue spots, {} genes", tissue_only.n_spots(), tissue_only.n_genes());
let totals = tissue_only.total_counts_per_spot();
let detected = tissue_only.genes_detected_per_spot();

// Convert to SpatialPoint for downstream spatial analysis
let points = visium_to_spatial_points(&tissue_only);
let graph = cyanea_omics::spatial::knn_graph(&points, 6).unwrap();

// MERFISH: subcellular-resolution imaging data
let merfish = MerfishData::new(
    gene_panel,       // typically 100-500 genes
    cell_ids,
    cell_centroids,   // (x, y) in µm
    cell_by_gene,     // counts[cell][gene]
).unwrap();
// Estimate false positive rate from blank barcodes
merfish.blank_counts = Some(blank_barcode_counts);
if let Some(fpr) = merfish.estimated_fpr() {
    println!("Mean FPR: {:.4}", fpr.iter().sum::<f64>() / fpr.len() as f64);
}

// Slide-seq: bead-based near-cellular resolution
let slideseq = SlideseqData::new(genes, barcodes, bead_coords, counts).unwrap();
let filtered = slideseq.filter_min_counts(100.0); // remove low-quality beads
```

## Spatial Cell Segmentation

```rust
use cyanea_omics::spatial_segmentation::*;

// Voronoi segmentation: assign transcripts to nearest nucleus
let nuclei = vec![(100.0, 200.0), (300.0, 150.0), (500.0, 400.0)];
let transcript_positions = vec![/* (x, y) per detected transcript */];
let result = voronoi_segmentation(&nuclei, &transcript_positions, Some(50.0)).unwrap();
println!("{} cells, {} assigned, {} unassigned",
    result.cells.len(), result.assigned_transcripts, result.unassigned_transcripts);

// Nucleus expansion: grow circles from seeds, stop at neighbors
let params = ExpansionParams {
    max_radius: 15.0,  // max expansion in µm
    min_gap: 1.0,      // min gap between cell boundaries
};
let result = expansion_segmentation(&nuclei, &transcript_positions, &params).unwrap();
for cell in &result.cells {
    println!("Cell {}: area={:.1}, transcripts={}", cell.cell_id, cell.area, cell.n_transcripts);
}

// Watershed on DAPI intensity grid
let dapi_grid = vec![vec![0.0; 512]; 512]; // row-major intensity image
let detected_nuclei = vec![(100, 200), (300, 150)]; // (row, col) seeds
let labels = watershed_grid(&dapi_grid, 512, 512, &detected_nuclei).unwrap();
// labels[row][col] gives cell ID (1-based), 0 = background
```

## Spatial Domain Detection

```rust
use cyanea_omics::spatial_domains::*;

// Detect tissue domains via spatially-aware k-means
let params = DomainParams {
    spatial_weight: 0.5,    // balance expression vs. space (0=expr only, 1=space only)
    n_domains: Some(5),     // or None for auto (sqrt(n/2))
    n_neighbors: 6,
    max_iter: 50,
    tolerance: 1e-4,
};
let result = detect_domains(&expression, &coords, &params).unwrap();
println!("{} domains detected", result.n_domains);
for domain in &result.domains {
    println!("Domain {}: {} spots", domain.domain_id, domain.members.len());
}

// HMRF smoothing: refine noisy labels with spatial consistency
let neighbors: Vec<Vec<usize>> = /* build kNN neighbor lists */;
let smoothed_labels = hmrf_smooth(
    &result.labels,    // initial labels from k-means
    &expression,       // spot x gene matrix
    &neighbors,        // per-spot neighbor indices
    1.0,               // beta: higher = smoother
    20,                // max ICM iterations
).unwrap();

// Find spatially variable genes (SVGs) via Moran's I
let svgs = find_spatially_variable_genes(
    &expression,
    &gene_names,
    &neighbors,
    999,   // permutations
    42,    // seed
).unwrap();
// Results sorted by BH-adjusted p-value
for svg in svgs.iter().take(10) {
    println!("{}: Moran's I = {:.3}, padj = {:.4}", svg.gene_name, svg.morans_i, svg.adjusted_p_value);
}
```

## Spatial Cell-Cell Communication

```rust
use cyanea_omics::spatial_cellchat::*;

// Use the built-in curated L-R database (12 pairs, 8 pathways)
let lr_db = demo_lr_database();

// Analyze communication with spatial distance weighting
let params = CommParams {
    distance_sigma: 100.0,   // Gaussian kernel sigma
    n_permutations: 100,
    p_threshold: 0.05,
    min_pct: 0.1,            // min fraction expressing
    seed: 42,
};
let results = analyze_communication(
    &expression,     // cell x gene
    &gene_names,
    &cell_types,     // per-cell type label
    &coords,         // (x, y) per cell
    &lr_db,
    &params,
).unwrap();

// Filter significant interactions
for r in results.iter().filter(|r| r.p_value < 0.05) {
    println!("{}: {} -> {} (prob={:.3}, p={:.3}, pathway={})",
        r.lr_pair, r.source, r.target, r.probability, r.p_value, r.pathway);
}

// Aggregate to pathway level
let pathways = aggregate_pathways(&results, 0.05);
for pw in &pathways {
    println!("{} {} -> {}: strength={:.2}, {} significant pairs",
        pw.pathway, pw.source, pw.target, pw.strength, pw.n_significant);
}
```

## Spatial Deconvolution

```rust
use cyanea_omics::spatial_deconvolution::*;

// Define cell type signatures (reference profiles)
let signatures = vec![
    CellTypeSignature {
        cell_type: "T_cell".into(),
        genes: vec!["CD3D".into(), "CD3E".into(), "CD8A".into()],
        weights: vec![10.0, 8.0, 6.0],
    },
    CellTypeSignature {
        cell_type: "Macrophage".into(),
        genes: vec!["CD68".into(), "CD163".into(), "CSF1R".into()],
        weights: vec![9.0, 7.0, 5.0],
    },
];

// NNLS deconvolution: estimate cell type proportions per spot
let result = nnls_deconvolve(&expression, &gene_names, &signatures).unwrap();
for spot in &result.spots {
    println!("Spot {}: T_cell={:.1}%, Macrophage={:.1}%, residual={:.3}",
        spot.spot_idx,
        spot.proportions[0] * 100.0,
        spot.proportions[1] * 100.0,
        spot.residual);
}
println!("Mean residual: {:.3}", result.mean_residual);

// Enrichment scoring with permutation p-values
let enrichment = score_enrichment(&expression, &gene_names, &signatures, 999, 42).unwrap();
for e in enrichment.iter().filter(|e| e.p_value < 0.05) {
    println!("Spot {} enriched for {} (score={:.2}, p={:.3})",
        e.spot_idx, e.cell_type, e.score, e.p_value);
}
```

## Hi-C / 3D Genome Analysis

```rust
use cyanea_omics::hic::*;

// Create a contact matrix from scratch
let mut matrix = ContactMatrix::new("chr1", 10000, 100); // 100 bins at 10kb resolution

// Add contacts (automatically symmetric for intra-chromosomal)
matrix.add_contact(5, 20, 150.0);
matrix.add_contact(10, 15, 200.0);
assert!((matrix.get(20, 5) - 150.0).abs() < 1e-10); // symmetric

// Or parse from .cool text export
let cool_data = "chr1\t0\t10000\tchr1\t10000\t20000\t150\n\
                  chr1\t0\t10000\tchr1\t20000\t30000\t50\n\
                  chr1\t10000\t20000\tchr1\t20000\t30000\t120\n";
let contacts = parse_cool_text(cool_data, 10000).unwrap();
let matrix = contacts_to_matrix(&contacts, "chr1", 10000, 100);

// Parse 4DN pairs format
let pairs_data = "## pairs format v1.0\n\
                   #columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n\
                   read_1\tchr1\t1000\tchr1\t50000\t+\t-\n\
                   read_2\tchr1\t2000\tchr1\t30000\t+\t+\n";
let pairs = parse_pairs(pairs_data).unwrap();
println!("{} contacts parsed", pairs.len());

// ICE balancing (iterative correction)
let mut matrix = ContactMatrix::from_dense("chr1", 10000, raw_data).unwrap();
let bias = matrix.ice_balance(50, 1e-6);

// KR (Knight-Ruiz) balancing
let weights = matrix.kr_balance(100, 1e-6);

// Observed/expected (distance normalization)
let oe = matrix.observed_expected();
// Within-TAD entries will have O/E > 1, between-TAD < 1
```

### TAD Calling

```rust
use cyanea_omics::hic::*;

// Call TADs using the insulation score method (Crane et al. 2015)
let params = TadParams {
    window_size: 10,             // insulation window in bins
    boundary_threshold: -0.5,    // z-score threshold for boundaries
    min_tad_size: 3,             // minimum TAD size in bins
};
let tads = call_tads(&matrix, &params).unwrap();
for tad in &tads {
    println!("TAD: {}:{}-{} ({} bins, boundary_score={:.2})",
        tad.chrom, tad.start_bp, tad.end_bp,
        tad.end_bin - tad.start_bin, tad.boundary_score);
}

// Get raw insulation scores for custom analysis
let scores = insulation_scores(&matrix, 10);
```

### A/B Compartment Calling

```rust
use cyanea_omics::hic::*;

// Call compartments via PC1 of the O/E correlation matrix
let result = call_compartments(&matrix, None).unwrap();
for (i, comp) in result.compartments.iter().enumerate() {
    println!("Bin {}: {:?} (eigenvector={:.3})", i, comp, result.eigenvector[i]);
}

// Use GC content to orient the eigenvector (A = high GC)
let gc_content: Vec<f64> = (0..matrix.n_bins1)
    .map(|i| if i < 50 { 0.55 } else { 0.40 })
    .collect();
let result = call_compartments(&matrix, Some(&gc_content)).unwrap();
```

### Loop Calling

```rust
use cyanea_omics::hic::*;

// Call chromatin loops (HiCCUPS-style local enrichment)
let params = LoopParams {
    background_window: 5,    // donut background window size in bins
    min_enrichment: 1.5,     // minimum enrichment over background
    min_distance: 5,         // min genomic distance in bins
    max_distance: 500,       // max genomic distance in bins
    p_threshold: 0.01,       // p-value threshold
};
let loops = call_loops(&matrix, &params).unwrap();
// Results sorted by enrichment (descending)
for lp in loops.iter().take(10) {
    println!("Loop: {}:{}-{} (enrichment={:.1}x, p={:.2e})",
        lp.chrom, lp.anchor1_bp, lp.anchor2_bp, lp.enrichment, lp.p_value);
}

// Write contacts back to 4DN pairs format
let output = write_pairs(&contacts, "chr1", 10000);
```

## Microarray Expression Analysis

```rust
use cyanea_omics::microarray::*;

// Quantile normalization across samples
let mut data = vec![
    vec![5.0, 2.0, 3.0],   // sample 1
    vec![4.0, 1.0, 4.0],   // sample 2
    vec![3.0, 4.0, 6.0],   // sample 3
];
quantile_normalize(&mut data).unwrap();
// Each sample now has identical rank distributions

// Full RMA pipeline (background correction + quantile normalization)
let probe_intensities = vec![
    vec![100.0, 200.0, 150.0],  // probe set 1
    vec![50.0, 80.0, 60.0],     // probe set 2
];
let normalized = rma_normalize(&probe_intensities).unwrap();

// Median polish for probe set summarization
let probes = vec![
    vec![5.0, 6.0, 7.0],   // probe 1 across samples
    vec![4.0, 5.0, 8.0],   // probe 2 across samples
    vec![6.0, 7.0, 6.0],   // probe 3 across samples
];
let summary = median_polish(&probes, 100).unwrap();
println!("Summarized expression: {:?}", summary);

// Differential expression with limma (moderated t-test)
let expression = vec![
    vec![2.5, 3.1, 2.8, 5.2, 6.0, 5.5],   // gene 1
    vec![1.0, 1.2, 0.9, 1.1, 1.0, 1.3],   // gene 2
];
let gene_names = vec!["BRCA1".into(), "GAPDH".into()];
let groups = vec![0, 0, 0, 1, 1, 1]; // control vs. treatment
let results = limma_diff_expr(&expression, &gene_names, &groups).unwrap();
for r in &results {
    println!("{}: log2FC={:.2}, p_adj={:.4}", r.gene_name, r.log2_fold_change, r.adjusted_p_value);
}
```

## Methylation Microarray Analysis

```rust
use cyanea_omics::microarray::*;

// Compute beta values from methylated (M) and unmethylated (U) signals
let methylated = vec![1000.0, 5000.0, 200.0];
let unmethylated = vec![4000.0, 500.0, 9800.0];
let betas = compute_beta(&methylated, &unmethylated, 100.0).unwrap();
// beta = M / (M + U + offset), range [0, 1]
println!("Beta values: {:?}", betas);

// Convert between beta values and M-values
let m_values = beta_to_m_value(&betas);  // log2(beta / (1 - beta))
let back_to_beta = m_value_to_beta(&m_values);

// SWAN normalization to correct Infinium I/II probe design bias
let beta_values = vec![
    vec![0.2, 0.8, 0.5],   // sample 1
    vec![0.3, 0.7, 0.6],   // sample 2
];
let design_types = vec![InfiniumType::TypeI, InfiniumType::TypeII, InfiniumType::TypeI];
let swan_normalized = swan_normalize(&beta_values, &design_types).unwrap();

// Differential methylation analysis
let probe_ids = vec!["cg00001".into(), "cg00002".into(), "cg00003".into()];
let groups = vec![0, 0, 0, 1, 1, 1];
let beta_matrix = vec![
    vec![0.1, 0.12, 0.11, 0.8, 0.85, 0.78],  // probe 1: hypermethylated in group 1
    vec![0.5, 0.52, 0.48, 0.51, 0.49, 0.53],  // probe 2: no change
    vec![0.9, 0.88, 0.91, 0.2, 0.22, 0.18],   // probe 3: hypomethylated in group 1
];
let dm_results = diff_methylation(&beta_matrix, &probe_ids, &groups).unwrap();
for r in &dm_results {
    println!("{}: delta_beta={:.2}, p_adj={:.4}", r.probe_id, r.delta_beta, r.adjusted_p_value);
}
```

## CNV Detection

```rust
use cyanea_omics::cnv::*;

// Circular binary segmentation on log2 ratio data
let log2_ratios = vec![0.1, 0.15, 0.12, 0.8, 0.9, 0.85, 0.1, 0.05];
let positions = vec![1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000];

let segments = cbs_segment(&log2_ratios, &positions, "chr1", Default::default()).unwrap();
for seg in &segments {
    println!("{}:{}-{} CN={} log2R={:.2}",
        seg.chrom, seg.start, seg.end, seg.copy_number, seg.log2_ratio);
}

// Merge adjacent similar segments
let merged = merge_segments(&segments, 10000, 0.3);
```

## Methylation Analysis

```rust
use cyanea_omics::methylation::*;

// Bisulfite conversion simulation
let reference = b"ACGTCGATCG";
let converted = bisulfite_convert(reference);

// Find CpG islands in a sequence
let islands = find_cpg_islands(reference, "chr1", Default::default()).unwrap();

// Detect differentially methylated regions
let group1_sites = vec![/* CpgSite vec for condition 1 */];
let group2_sites = vec![/* CpgSite vec for condition 2 */];
let dmrs = find_dmrs(&group1_sites, &group2_sites, Default::default()).unwrap();
```

## ACMG/AMP Variant Classification

```rust
use cyanea_omics::acmg::*;
use cyanea_omics::variant::Variant;

// Build evidence manually using the builder pattern
let evidence = AcmgEvidence::new()
    .pvs1("Frameshift in BRCA1 — known LOF mechanism")
    .strong_pathogenic("PS1", "Same amino acid change as established pathogenic variant")
    .moderate_pathogenic("PM2", "Absent from gnomAD");

let result = evidence.classify();
println!("{}", result.classification.as_str()); // "Pathogenic"
println!("PVS={}, PS={}, PM={}", result.pvs_count, result.ps_count, result.pm_count);

// Auto-assign evidence from variant properties
let variant = Variant::new("chr17", 43092919, b"ACGTG".to_vec(), vec![b"A".to_vec()]).unwrap();
let auto_ev = auto_evidence(&variant, None, true, Some(true));
let auto_result = auto_ev.classify();
// PVS1 (frameshift in LOF gene) + PM2 (absent) + PP3 (in silico) → Likely Pathogenic or Pathogenic

// Match against ClinVar annotations
let clinvar_tsv = "chrom\tpos\tref\talt\tsignificance\treview_status\tconditions\tsubmitter_count\tstar_rating\n\
                   chr17\t43092919\tA\tG\tPathogenic\texpert panel\tBreast cancer\t5\t3\n";
let db = parse_clinvar_tsv(clinvar_tsv).unwrap();

let query = Variant::new("chr17", 43092919, b"A".to_vec(), vec![b"G".to_vec()]).unwrap();
if let Some(ann) = match_clinvar(&query, &db) {
    println!("{} ({}, {} stars)", ann.significance, ann.review_status, ann.star_rating);
}
```

## Pharmacogenomics

```rust
use cyanea_omics::pharmacogenomics::*;
use cyanea_omics::variant::Variant;

// Use the built-in demo CYP2D6 database (*4, *10, *17 alleles)
let db = demo_cyp2d6_database();

// Call star alleles from observed variants
let variants = vec![
    Variant::new("chr22", 42128945, b"C".to_vec(), vec![b"T".to_vec()]).unwrap(), // CYP2D6*4
];
let call = call_star_alleles("CYP2D6", &variants, &db).unwrap();
println!("Diplotype: {}", call.diplotype);  // e.g., "*4/*1"
println!("Activity score: {}", call.activity_score);
println!("Phenotype: {}", call.phenotype.as_str());

// No variants → reference *1/*1 (Normal Metabolizer)
let ref_call = call_star_alleles("CYP2D6", &[], &db).unwrap();
assert_eq!(ref_call.diplotype, "*1/*1");

// Look up drug interaction recommendations
let recs = lookup_drug_interactions(&db, "CYP2D6", MetabolizerPhenotype::PoorMetabolizer);
for rec in &recs {
    println!("{}: {} [{}] ({})", rec.drug, rec.recommendation, rec.evidence_level, rec.source);
}

// Build a custom PGx database
let mut custom_db = PgxDatabase::new();
custom_db.add_allele(StarAllele {
    gene: "CYP2C19".into(),
    allele: "*2".into(),
    defining_variants: vec![("chr10".into(), 94781859, b"G".to_vec(), b"A".to_vec())],
    activity_score: 0.0,
    function: AlleleFunction::NoFunction,
});
```

## Clinical Genomics: HLA Typing

```rust
use cyanea_omics::clinical::*;

// Parse HLA typing from tab-separated input
let donor_typing = parse_hla_typing("A\t02:01\t03:01\nB\t07:02\t44:02\nDRB1\t04:01\t15:01\n").unwrap();
let recipient_typing = parse_hla_typing("A\t02:01\t68:01\nB\t07:02\t44:02\nDRB1\t04:01\t07:01\n").unwrap();

// Check compatibility across HLA-A, -B, -DRB1
let match_count = hla_compatibility(&donor_typing, &recipient_typing, &["A", "B", "DRB1"]);
println!("Matched alleles: {} / 6", match_count);

// Inspect allele resolution
let (a1, a2) = donor_typing.diplotype("A").unwrap();
println!("{} / {}", a1.four_digit(), a2.four_digit()); // "A*02:01 / A*03:01"
```

## Clinical Genomics: Tumor Mutational Burden

```rust
use cyanea_omics::clinical::*;
use cyanea_omics::variant::Variant;

// Compute TMB from somatic variants
let somatic_variants: Vec<Variant> = (0..600)
    .map(|i| Variant::new("chr1", i * 100, b"A".to_vec(), vec![b"G".to_vec()]).unwrap())
    .collect();

let tmb = compute_tmb(&somatic_variants, 30.0, false).unwrap();
println!("TMB: {:.1} mut/Mb ({})", tmb.tmb, tmb.category.as_str());
// TMB: 20.0 mut/Mb (TMB-High) — above FDA pembrolizumab threshold
```

## Clinical Genomics: Microsatellite Instability

```rust
use cyanea_omics::clinical::*;
use std::collections::HashMap;

// Use the standard 5 Bethesda markers
let markers = bethesda_markers();

// Simulate tumor with instability at BAT25 and BAT26
let mut observed: HashMap<String, usize> = markers.iter()
    .map(|m| (m.name.clone(), m.reference_count))
    .collect();
*observed.get_mut("BAT25").unwrap() = 20; // shifted by 5 from reference 25
*observed.get_mut("BAT26").unwrap() = 20; // shifted by 6 from reference 26

let result = call_msi(&observed, &markers, 2);
println!("{}: {}/{} unstable loci ({:.0}%)",
    result.status.as_str(), result.unstable_loci, result.total_loci,
    result.instability_fraction * 100.0);
// MSI-H: 2/5 unstable loci (40%)
```

## Liftover

```rust
use cyanea_omics::liftover::*;
use cyanea_omics::{GenomicInterval, Strand};

// Parse chain file
let chain = parse_chain_file(&chain_file_contents).unwrap();

// Lift over a single interval
let interval = GenomicInterval::new("chr1", 100000, 100500, Strand::Forward);
match liftover(&interval, &chain) {
    LiftoverResult::Mapped(new_iv) => {
        println!("Mapped to {}:{}-{}", new_iv.chrom, new_iv.start, new_iv.end);
    }
    LiftoverResult::Partial { mapped, fraction_mapped } => {
        println!("Partially mapped ({:.0}%)", fraction_mapped * 100.0);
    }
    LiftoverResult::Unmapped => {
        println!("Could not map interval");
    }
}
```
