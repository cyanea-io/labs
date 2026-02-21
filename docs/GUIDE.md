# Cross-Crate Usage Guide

This guide shows how to combine multiple Cyanea crates to build common bioinformatics pipelines. Each example includes the Cargo.toml dependencies, required feature flags, and a complete Rust code snippet.

## 1. FASTQ to Alignment to Variants

Parse FASTQ reads, trim adapters, align to a reference, and extract variant positions.

**Crates**: cyanea-seq, cyanea-align, cyanea-io, cyanea-omics

```toml
[dependencies]
cyanea-seq = "0.1"
cyanea-align = "0.1"
cyanea-io = { version = "0.1", features = ["sam", "vcf"] }
cyanea-omics = "0.1"
```

```rust
use cyanea_seq::fastq::parse_fastq_records;
use cyanea_seq::trim::{TrimConfig, trim_quality};
use cyanea_align::seed_extend::seed_and_extend;
use cyanea_align::scoring::Scoring;
use cyanea_io::sam::format_sam_record;
use cyanea_omics::variant::Variant;
use cyanea_omics::genomic::GenomicInterval;

fn reads_to_variants(fastq_data: &str, reference: &[u8]) -> Vec<Variant> {
    // 1. Parse and trim reads
    let records = parse_fastq_records(fastq_data).unwrap();
    let config = TrimConfig { min_quality: 20, min_length: 30, ..Default::default() };
    let trimmed: Vec<_> = records.iter()
        .filter_map(|r| trim_quality(r, &config).ok())
        .collect();

    // 2. Align each read to the reference
    let scoring = Scoring::default();
    let mut variants = Vec::new();
    for read in &trimmed {
        let result = seed_and_extend(read.sequence(), reference, 11, &scoring).unwrap();

        // 3. Extract mismatches as variant candidates
        for (qpos, rpos) in result.mismatches() {
            variants.push(Variant {
                chrom: "chr1".to_string(),
                pos: rpos as u64,
                ref_allele: String::from_utf8(vec![reference[rpos]]).unwrap(),
                alt_allele: String::from_utf8(vec![read.sequence()[qpos]]).unwrap(),
                ..Default::default()
            });
        }
    }
    variants
}
```

## 2. Single-Cell Analysis Pipeline

Load a count matrix, normalize, find highly variable genes, reduce dimensions, cluster, and compute trajectory.

**Crates**: cyanea-omics (with `single-cell` feature), cyanea-ml

```toml
[dependencies]
cyanea-omics = { version = "0.1", features = ["single-cell"] }
cyanea-ml = "0.1"
```

```rust
use cyanea_omics::sc_preprocess::{normalize_total, log1p_transform, highly_variable_genes};
use cyanea_omics::sc_cluster::{build_knn_graph, leiden};
use cyanea_omics::sc_trajectory::{diffusion_map, dpt};
use cyanea_omics::sc_markers::rank_genes_groups;
use cyanea_omics::sc_integrate::harmony;
use cyanea_ml::reduction::pca;

fn single_cell_pipeline(counts: &[Vec<f64>], batch_labels: &[usize]) {
    // 1. Preprocessing
    let normalized = normalize_total(counts, 10_000.0).unwrap();
    let log_data = log1p_transform(&normalized);
    let hvg_mask = highly_variable_genes(&log_data, 2000).unwrap();
    let filtered = select_features(&log_data, &hvg_mask);

    // 2. Dimensionality reduction
    let pca_result = pca(&filtered, 50).unwrap();

    // 3. Batch correction (if multiple batches)
    let corrected = harmony(&pca_result.embeddings, batch_labels, 20).unwrap();

    // 4. Neighborhood graph and clustering
    let graph = build_knn_graph(&corrected, 15).unwrap();
    let clusters = leiden(&graph, 1.0).unwrap();

    // 5. Trajectory inference
    let dm = diffusion_map(&graph, 10).unwrap();
    let pseudotime = dpt(&dm, 0).unwrap();  // root cell = 0

    // 6. Marker gene detection
    let markers = rank_genes_groups(&log_data, &clusters).unwrap();
}
```

## 3. Cheminformatics Pipeline

Parse a molecule from SMILES, compute properties, generate fingerprints, assess drug-likeness, decompose scaffolds, generate 3D conformers, and optimize geometry.

**Crates**: cyanea-chem

```toml
[dependencies]
cyanea-chem = "0.1"
```

```rust
use cyanea_chem::smiles::parse_smiles;
use cyanea_chem::properties::molecular_properties;
use cyanea_chem::fingerprint::{morgan_fingerprint, tanimoto_similarity};
use cyanea_chem::descriptors::compute_descriptors;
use cyanea_chem::druglikeness::lipinski_violations;
use cyanea_chem::scaffold::murcko_scaffold;
use cyanea_chem::conformer::generate_conformer;
use cyanea_chem::forcefield::{UffForceField, minimize_energy};

fn chem_pipeline(smiles: &str) {
    // 1. Parse and canonicalize
    let mol = parse_smiles(smiles).unwrap();

    // 2. Compute molecular properties
    let props = molecular_properties(&mol).unwrap();
    println!("MW: {:.1}, LogP: {:.2}, HBA: {}, HBD: {}",
        props.molecular_weight, props.logp,
        props.h_bond_acceptors, props.h_bond_donors);

    // 3. Fingerprints and similarity
    let fp1 = morgan_fingerprint(&mol, 2, 2048).unwrap();
    let mol2 = parse_smiles("CCO").unwrap();
    let fp2 = morgan_fingerprint(&mol2, 2, 2048).unwrap();
    let sim = tanimoto_similarity(&fp1, &fp2);

    // 4. Drug-likeness assessment
    let violations = lipinski_violations(&props);
    println!("Lipinski violations: {}", violations);

    // 5. Scaffold decomposition
    let scaffold = murcko_scaffold(&mol).unwrap();

    // 6. 3D conformer generation and optimization
    let conformer = generate_conformer(&mol, 42).unwrap();
    let ff = UffForceField::new(&mol, &conformer).unwrap();
    let optimized = minimize_energy(&ff, 500).unwrap();
    println!("Final energy: {:.2} kcal/mol", optimized.energy);
}
```

## 4. Phylogenetics Pipeline

Build a multiple sequence alignment, construct a distance matrix, build a tree with neighbor joining, assess bootstrap support, optimize with maximum likelihood, date the tree, and render it.

**Crates**: cyanea-align, cyanea-phylo (with `ml` feature)

```toml
[dependencies]
cyanea-align = "0.1"
cyanea-phylo = { version = "0.1", features = ["ml"] }
```

```rust
use cyanea_align::msa::progressive_msa;
use cyanea_align::scoring::Scoring;
use cyanea_phylo::distance::jukes_cantor_distances;
use cyanea_phylo::reconstruct::neighbor_joining;
use cyanea_phylo::bootstrap::bootstrap_support;
use cyanea_phylo::likelihood::optimize_branch_lengths;
use cyanea_phylo::models::GtrModel;
use cyanea_phylo::dating::root_to_tip_dating;
use cyanea_phylo::drawing::newick_string;
use cyanea_phylo::tree_search::nni_search;
use cyanea_phylo::model_selection::select_model_aic;

fn phylo_pipeline(sequences: &[&[u8]], names: &[&str]) {
    // 1. Multiple sequence alignment
    let scoring = Scoring::default();
    let msa = progressive_msa(sequences, &scoring).unwrap();

    // 2. Distance matrix
    let distances = jukes_cantor_distances(&msa).unwrap();

    // 3. Neighbor-joining tree
    let mut tree = neighbor_joining(&distances, names).unwrap();

    // 4. Model selection
    let best_model = select_model_aic(&msa, &tree).unwrap();

    // 5. ML optimization with NNI tree search
    let model = GtrModel::default();
    optimize_branch_lengths(&mut tree, &msa, &model).unwrap();
    nni_search(&mut tree, &msa, &model).unwrap();

    // 6. Bootstrap support (100 replicates)
    let support = bootstrap_support(&msa, names, 100, &model).unwrap();

    // 7. Molecular dating
    let dated = root_to_tip_dating(&tree, &tip_dates).unwrap();

    // 8. Export
    let nwk = newick_string(&tree).unwrap();
    println!("{}", nwk);
}
```

## 5. Population Genetics

Parse a VCF file, compute allele frequencies, test for Hardy-Weinberg equilibrium, calculate Fst between populations, run Tajima's D, compute linkage disequilibrium, and perform PCA on genotypes.

**Crates**: cyanea-io (with `vcf` feature), cyanea-stats

```toml
[dependencies]
cyanea-io = { version = "0.1", features = ["vcf"] }
cyanea-stats = "0.1"
```

```rust
use cyanea_io::vcf::parse_vcf;
use cyanea_stats::popgen::{
    allele_frequencies, hwe_exact_test, weir_cockerham_fst,
    tajimas_d, linkage_disequilibrium, genotype_pca,
};

fn popgen_pipeline(vcf_data: &str, populations: &[Vec<usize>]) {
    // 1. Parse VCF
    let variants = parse_vcf(vcf_data).unwrap();

    // 2. Allele frequencies per population
    for (i, pop) in populations.iter().enumerate() {
        let freqs = allele_frequencies(&variants, pop).unwrap();
        println!("Pop {}: {} variants", i, freqs.len());
    }

    // 3. Hardy-Weinberg equilibrium test
    for variant in &variants {
        let hwe = hwe_exact_test(variant).unwrap();
        if hwe.p_value < 0.001 {
            println!("HWE violation at {}:{}", variant.chrom, variant.pos);
        }
    }

    // 4. Pairwise Fst between populations
    let fst = weir_cockerham_fst(&variants, &populations[0], &populations[1]).unwrap();
    println!("Fst = {:.4}", fst);

    // 5. Tajima's D in sliding windows
    let td = tajimas_d(&variants, &populations[0]).unwrap();
    println!("Tajima's D = {:.4}", td);

    // 6. Linkage disequilibrium
    let ld = linkage_disequilibrium(&variants, 100_000).unwrap();

    // 7. Population structure via PCA
    let pca = genotype_pca(&variants, 10).unwrap();
    println!("PC1 variance explained: {:.1}%", pca.variance_explained[0] * 100.0);
}
```

## 6. Microbiome Analysis

Compute alpha diversity, build a distance matrix, run PCoA ordination, test for group differences with PERMANOVA and ANOSIM, and calculate phylogenetic diversity with UniFrac.

**Crates**: cyanea-stats, cyanea-ml, cyanea-phylo

```toml
[dependencies]
cyanea-stats = "0.1"
cyanea-ml = "0.1"
cyanea-phylo = { version = "0.1", features = ["ml"] }
```

```rust
use cyanea_stats::diversity::{shannon, simpson, chao1};
use cyanea_ml::distance::pairwise_distances;
use cyanea_stats::ordination::pcoa;
use cyanea_stats::multivariate::{permanova, anosim};
use cyanea_phylo::unifrac::{unweighted_unifrac, weighted_unifrac};
use cyanea_phylo::newick::parse_newick;

fn microbiome_pipeline(
    otu_table: &[Vec<f64>],
    groups: &[usize],
    tree_nwk: &str,
) {
    // 1. Alpha diversity per sample
    for (i, sample) in otu_table.iter().enumerate() {
        let h = shannon(sample).unwrap();
        let s = simpson(sample).unwrap();
        let c = chao1(sample).unwrap();
        println!("Sample {}: Shannon={:.2}, Simpson={:.2}, Chao1={:.0}", i, h, s, c);
    }

    // 2. Beta diversity (Bray-Curtis distance matrix)
    let distances = pairwise_distances(otu_table, "bray_curtis").unwrap();

    // 3. Ordination (PCoA)
    let coords = pcoa(&distances, 3).unwrap();
    println!("PCoA axis 1 explains {:.1}% variance",
        coords.eigenvalues[0] / coords.eigenvalues.iter().sum::<f64>() * 100.0);

    // 4. Statistical testing
    let perm = permanova(&distances, groups, 999).unwrap();
    println!("PERMANOVA: F={:.2}, p={:.4}", perm.f_statistic, perm.p_value);

    let anos = anosim(&distances, groups, 999).unwrap();
    println!("ANOSIM: R={:.3}, p={:.4}", anos.r_statistic, anos.p_value);

    // 5. Phylogenetic diversity (UniFrac)
    let tree = parse_newick(tree_nwk).unwrap();
    let uf = unweighted_unifrac(otu_table, &tree).unwrap();
    let wuf = weighted_unifrac(otu_table, &tree).unwrap();
    println!("UniFrac distances computed: {}x{}", uf.len(), uf[0].len());
}
```
