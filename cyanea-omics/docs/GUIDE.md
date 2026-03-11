# cyanea-omics Usage Guide

Practical examples for working with genomic data structures, single-cell analysis, variant annotation, spatial transcriptomics, ACMG variant classification, pharmacogenomics, and clinical genomics.

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
