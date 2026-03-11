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

## 7. Shotgun Metagenomics Pipeline

Build a taxonomy, classify reads, profile community composition, compute diversity, and perform differential abundance testing.

**Crates**: cyanea-meta, cyanea-seq

```toml
[dependencies]
cyanea-meta = "0.1"
cyanea-seq = "0.1"
```

```rust
use cyanea_meta::taxonomy::{TaxonomyDB, TaxonRank};
use cyanea_meta::profile::{profile_from_classifications, normalize_profile, filter_profile};
use cyanea_meta::diversity::{alpha_diversity, beta_diversity_matrix, rarefaction_curve};
use cyanea_meta::composition::{clr_transform, differential_abundance};
use cyanea_meta::assembly::assembly_stats;

fn metagenomics_pipeline(
    reads: &[&[u8]],
    contigs: &[&[u8]],
    sample1_counts: &[u64],
    sample2_counts: &[u64],
) {
    // 1. Build taxonomy and classify reads
    let mut db = TaxonomyDB::new(1);
    db.add_node(2, 1, TaxonRank::Domain, "Bacteria").unwrap();
    // ... add more nodes, index references ...

    let mut classifications = Vec::new();
    for read in reads {
        if let Ok(Some(taxid)) = db.classify_sequence(read) {
            classifications.push(taxid);
        }
    }

    // 2. Build taxonomic profile
    let profile = profile_from_classifications(&classifications).unwrap();
    let normalized = normalize_profile(&profile).unwrap();
    let filtered = filter_profile(&normalized, 0.001).unwrap(); // >0.1% abundance

    // 3. Alpha diversity
    let ad = alpha_diversity(sample1_counts).unwrap();
    println!("Shannon: {:.3}, Chao1: {:.1}", ad.shannon, ad.chao1);

    // 4. Beta diversity across samples
    let beta = beta_diversity_matrix(&[sample1_counts, sample2_counts]).unwrap();
    println!("Bray-Curtis: {:.3}", beta.get(0, 1));

    // 5. Rarefaction curve
    let steps = vec![100, 500, 1000, 5000];
    let curve = rarefaction_curve(sample1_counts, &steps).unwrap();

    // 6. Assembly QC
    let stats = assembly_stats(contigs).unwrap();
    println!("N50: {} bp, {} contigs, GC: {:.1}%",
        stats.n50, stats.contig_count, stats.gc_content * 100.0);
}
```

## 8. ChIP-seq / ATAC-seq Epigenomics Pipeline

Build signal tracks, call peaks, discover motifs, learn chromatin states, and run differential binding analysis.

**Crates**: cyanea-epi

```toml
[dependencies]
cyanea-epi = "0.1"
```

```rust
use cyanea_epi::pileup::{build_pileup, normalize_pileup, smooth_pileup, pileup_correlation};
use cyanea_epi::peaks::{call_peaks, PeakCallParams, PeakSet};
use cyanea_epi::motifs::{discover_motifs, scan_sequence, DiscoveryParams, write_meme};
use cyanea_epi::chromatin::{learn_chromatin_states, segment_genome, ChromHMMParams};
use cyanea_epi::differential::differential_peaks;
use cyanea_epi::accessibility::{tss_enrichment, atacqc};

fn epigenomics_pipeline(
    chip_reads: &[(String, u64, u64)],
    input_reads: &[(String, u64, u64)],
    peak_sequences: &[&[u8]],
) {
    // 1. Build and normalize pileups
    let treatment = build_pileup(chip_reads, 200);
    let control = build_pileup(input_reads, 200);
    let normalized = normalize_pileup(&treatment, "cpm").unwrap();
    let smoothed = smooth_pileup(&normalized, 150.0);

    // 2. Call peaks
    let params = PeakCallParams::default();
    let peaks = call_peaks(&treatment, &control, &params).unwrap();
    let peak_set = PeakSet::new(peaks);
    let stats = peak_set.stats();
    println!("{} peaks, median width {} bp", stats.count, stats.median_width);

    // 3. Motif discovery at peak summits
    let disc_params = DiscoveryParams {
        motif_width: 8,
        n_motifs: 10,
        ..Default::default()
    };
    let motifs = discover_motifs(peak_sequences, &disc_params).unwrap();
    let meme_output = write_meme(&motifs);
    println!("Discovered {} motifs", motifs.len());

    // 4. Chromatin state learning (from multiple histone marks)
    let marks = vec![
        vec![true, true, false, false],   // active promoter bin
        vec![false, false, true, false],   // repressed bin
        // ...
    ];
    let hmm_params = ChromHMMParams { n_states: 5, ..Default::default() };
    let model = learn_chromatin_states(&marks, &hmm_params).unwrap();
    let segmentation = segment_genome(&model, &marks).unwrap();

    // 5. Differential binding between conditions
    let counts = vec![
        vec![100.0, 110.0, 20.0, 25.0],
        vec![50.0, 45.0, 55.0, 48.0],
    ];
    let conditions = vec!["treated", "treated", "control", "control"];
    let diff = differential_peaks(&counts, &conditions, "wald").unwrap();
    let significant: Vec<_> = diff.iter().filter(|d| d.q_value < 0.05).collect();
    println!("{} differentially bound regions (FDR < 0.05)", significant.len());
}
```

## 9. Proteomics: MS/MS Database Search Pipeline

Parse mass spectra, digest a protein database, search for peptide-spectrum matches, control FDR, infer proteins, and quantify.

**Crates**: cyanea-proteomics

**Cargo.toml**:

```toml
[dependencies]
cyanea-proteomics = { path = "../cyanea-proteomics" }
```

```rust
use cyanea_proteomics::mgf::parse_mgf;
use cyanea_proteomics::peptide::{digest, DigestConfig, Protease};
use cyanea_proteomics::search::{search_all, generate_decoys, SearchConfig};
use cyanea_proteomics::fdr::{filter_fdr, fdr_summary};
use cyanea_proteomics::protein::{infer_proteins, ProteinEntry};
use cyanea_proteomics::quantification::{spectral_counting, quantify_tmt, TmtPlex, median_normalize};
use cyanea_proteomics::mztab::write_mztab_proteins;

fn main() {
    // 1. Parse MGF spectra
    let mgf_text = std::fs::read_to_string("experiment.mgf").unwrap();
    let spectra = parse_mgf(&mgf_text).unwrap();
    println!("{} spectra loaded", spectra.len());

    // 2. In-silico digest of protein database
    let protein_sequences = vec![
        ("P12345", b"MAAAKPEPTIDEKSEQENCER".as_slice()),
        ("Q67890", b"MDDDEKFFFFERGGGGK".as_slice()),
    ];

    let digest_config = DigestConfig {
        protease: Protease::Trypsin,
        max_missed_cleavages: 2,
        min_length: 6,
        max_length: 50,
        min_mass: 400.0,
        max_mass: 5000.0,
    };

    let mut all_peptides = Vec::new();
    for (_, seq) in &protein_sequences {
        all_peptides.extend(digest(seq, &digest_config).unwrap());
    }
    println!("{} target peptides", all_peptides.len());

    // 3. Generate decoy database (reversed sequences)
    let decoy_peptides = generate_decoys(&all_peptides);

    // 4. Search spectra against target and decoy databases
    let search_config = SearchConfig {
        fragment_tolerance: 0.02,   // 20 mDa for high-res MS2
        precursor_tolerance: 10.0,  // Da
        max_fragment_charge: 2,
        min_matched_ions: 4,
    };

    let target_psms = search_all(&spectra, &all_peptides, &search_config).unwrap();
    let decoy_psms = search_all(&spectra, &decoy_peptides, &search_config).unwrap();
    println!("{} target PSMs, {} decoy PSMs", target_psms.len(), decoy_psms.len());

    // 5. FDR control at 1%
    let summary = fdr_summary(&target_psms, &decoy_psms).unwrap();
    println!("PSMs at 1% FDR: {}, at 5% FDR: {}", summary.passing_1pct, summary.passing_5pct);

    let filtered_psms = filter_fdr(&target_psms, &decoy_psms, 0.01).unwrap();

    // 6. Protein inference (parsimony)
    let protein_db: Vec<ProteinEntry> = protein_sequences.iter()
        .map(|(acc, seq)| ProteinEntry::new(*acc, *seq))
        .collect();

    let protein_groups = infer_proteins(&filtered_psms, &protein_db).unwrap();
    for g in &protein_groups {
        println!("Protein {:?}: {} unique peptides, {:.0}% coverage",
            g.accessions, g.unique_peptides.len(), g.coverage * 100.0);
    }

    // 7. Label-free quantification (spectral counting)
    let mut quants = spectral_counting(&protein_groups, &filtered_psms);
    median_normalize(&mut quants);
    for q in &quants {
        println!("{}: raw={:.0}, normalized={:.3}",
            q.accessions.join(","), q.raw_value, q.normalized_value);
    }

    // 8. TMT quantification (if using TMT labeling)
    let tmt_quants = quantify_tmt(&spectra, TmtPlex::Tmt10, 0.01);
    println!("{} spectra with TMT reporter ions", tmt_quants.len());

    // 9. Export results in mzTab format
    let mztab = write_mztab_proteins(&protein_groups, Some(&quants));
    println!("{}", mztab);
}
```

## 10. Network and Pathway Biology

Build biological networks, compute centrality and community structure, infer gene regulatory networks, score pathway topology, and analyze pathway crosstalk.

**Crates**: cyanea-network

```toml
[dependencies]
cyanea-network = "0.1"
```

```rust
use cyanea_network::graph::{Graph, EdgeWeight};
use cyanea_network::centrality::{degree_centrality, betweenness_centrality, pagerank};
use cyanea_network::community::{louvain, label_propagation};
use cyanea_network::ppi::{build_ppi_network, hub_proteins, ppi_enrichment};
use cyanea_network::grn::{correlation_network, mutual_information_network, clr_network};
use cyanea_network::pathway::{score_topology, pathway_crosstalk};
use cyanea_network::io::{parse_graphml, write_gexf, parse_gmt, write_sif};

fn network_pipeline(
    interactions: &[(String, String, f64)],
    expression_matrix: &[Vec<f64>],
    gene_names: &[&str],
    gmt_data: &str,
) {
    // 1. Build a PPI network from interaction data
    let ppi = build_ppi_network(interactions).unwrap();
    println!("{} nodes, {} edges", ppi.node_count(), ppi.edge_count());

    // 2. Centrality analysis — find hub proteins
    let degree = degree_centrality(&ppi).unwrap();
    let betweenness = betweenness_centrality(&ppi).unwrap();
    let pr = pagerank(&ppi, 0.85, 100).unwrap();
    let hubs = hub_proteins(&ppi, &degree, 0.95).unwrap();
    println!("{} hub proteins (top 5% degree)", hubs.len());

    // 3. Community detection (Louvain modularity optimization)
    let communities = louvain(&ppi, 1.0).unwrap();
    println!("{} communities, modularity={:.3}",
        communities.n_communities, communities.modularity);

    // 4. Label propagation (faster alternative)
    let lp = label_propagation(&ppi, 100).unwrap();
    println!("{} communities via LP", lp.n_communities);

    // 5. Gene regulatory network inference from expression data
    let corr_net = correlation_network(expression_matrix, gene_names, 0.7).unwrap();
    let mi_net = mutual_information_network(expression_matrix, gene_names, 0.5).unwrap();
    let clr_net = clr_network(expression_matrix, gene_names).unwrap();
    println!("CLR network: {} edges", clr_net.edge_count());

    // 6. Pathway topology scoring
    let gene_sets = parse_gmt(gmt_data).unwrap();
    for gs in &gene_sets {
        let score = score_topology(&ppi, &gs.genes).unwrap();
        println!("Pathway '{}': topology score={:.3}", gs.name, score);
    }

    // 7. Pathway crosstalk analysis
    let crosstalk = pathway_crosstalk(&ppi, &gene_sets).unwrap();
    for ct in crosstalk.iter().take(10) {
        println!("{} <-> {}: {} shared edges, score={:.3}",
            ct.pathway_a, ct.pathway_b, ct.shared_edges, ct.score);
    }

    // 8. PPI enrichment test (are input genes more connected than expected?)
    let test_genes: Vec<&str> = gene_names[..20].to_vec();
    let enrich = ppi_enrichment(&ppi, &test_genes, 1000).unwrap();
    println!("PPI enrichment p={:.4}", enrich.p_value);

    // 9. I/O — read/write network formats
    let graphml_data = std::fs::read_to_string("network.graphml").unwrap();
    let imported = parse_graphml(&graphml_data).unwrap();
    let gexf_output = write_gexf(&ppi).unwrap();
    let sif_output = write_sif(&ppi).unwrap();
}
```

## 11. Long-Read Sequencing Analysis

Analyze PacBio HiFi and Oxford Nanopore long reads: compute read statistics, self-correct reads, call structural variants from split alignments, run nanopore-specific QC, and detect CpG methylation.

**Crates**: cyanea-seq

```toml
[dependencies]
cyanea-seq = "0.1"
```

```rust
use cyanea_seq::longread::{
    LongRead, LongReadPlatform, longread_stats, self_correct,
    simple_consensus, LongReadSimConfig, simulate_long_reads,
    trim_adapters, common_adapters,
};
use cyanea_seq::sv::{
    SplitAlignment, SvCallConfig, call_svs, cluster_svs,
    svs_from_cigar, sv_summary,
};
use cyanea_seq::nanopore::{
    parse_signal_metadata, parse_methylation_calls,
    aggregate_methylation, nanopore_qc,
};

fn longread_pipeline(
    reads: &[LongRead],
    split_alignments: &[SplitAlignment],
    nanopore_metadata_lines: &[&str],
    methylation_text: &str,
) {
    // 1. Long-read statistics (N50, mean/median length, quality)
    let stats = longread_stats(reads).unwrap();
    println!("N50: {} bp, mean length: {} bp, mean quality: {:.1}",
        stats.n50, stats.mean_length as u64, stats.mean_quality);

    // 2. Self-correct a read using k-mer consensus
    let corrected = self_correct(&reads[0], 9).unwrap();
    println!("Corrected {} -> {} bases, {:.1}% identity",
        reads[0].len(), corrected.sequence.len(),
        corrected.identity * 100.0);

    // 3. Consensus from overlapping reads
    let seqs: Vec<&[u8]> = reads.iter().map(|r| r.sequence.as_slice()).collect();
    let consensus = simple_consensus(&seqs).unwrap();
    println!("Consensus: {} bp", consensus.len());

    // 4. Trim platform-specific adapters
    let adapters = common_adapters();
    for read in reads {
        let trimmed = trim_adapters(read, &adapters, 0.8).unwrap();
        println!("{}: {} -> {} bp after adapter trimming",
            read.id, read.len(), trimmed.len());
    }

    // 5. Simulate long reads from a reference (for benchmarking)
    let reference = b"ATCGATCGATCG";
    let sim_config = LongReadSimConfig {
        platform: LongReadPlatform::PacBioHiFi,
        num_reads: 100,
        mean_length: 15000,
        ..Default::default()
    };
    let simulated = simulate_long_reads(reference, &sim_config).unwrap();
    println!("Simulated {} reads", simulated.len());

    // 6. Call structural variants from split alignments
    let sv_config = SvCallConfig::default();
    let mut svs = call_svs(split_alignments, &sv_config).unwrap();
    println!("Called {} raw SVs", svs.len());

    // 7. Cluster nearby SVs and summarize
    let merged = cluster_svs(&mut svs, &sv_config);
    let summary = sv_summary(&merged);
    println!("{} INS, {} DEL, {} INV, {} DUP, {} BND",
        summary.insertions, summary.deletions, summary.inversions,
        summary.duplications, summary.breakends);

    // 8. Extract SVs from CIGAR strings (large indels)
    let cigar_svs = svs_from_cigar("chr1", 1000, "50M500D30M2000I20M", 50).unwrap();
    println!("{} SVs from CIGAR", cigar_svs.len());

    // 9. Nanopore-specific QC
    let qc = nanopore_qc(reads).unwrap();
    println!("Nanopore QC: {} reads, N50={}, median Q={:.1}",
        qc.total_reads, qc.n50, qc.median_quality);

    // 10. Parse signal metadata from POD5/FAST5
    for line in nanopore_metadata_lines {
        let meta = parse_signal_metadata(line).unwrap();
        println!("Read {} on channel {}, duration {:.1}s",
            meta.read_id, meta.channel, meta.duration);
    }

    // 11. CpG methylation from modified base calls
    let meth_calls = parse_methylation_calls(methylation_text).unwrap();
    let sites = aggregate_methylation(&meth_calls).unwrap();
    for site in &sites {
        println!("{}:{} methylation: {:.1}% ({} reads)",
            site.chrom, site.position,
            site.methylation_frequency * 100.0, site.coverage);
    }
}
```

## 12. Using Sample Datasets

Use `cyanea-datasets` to get pre-built sample data for prototyping and testing. All data is generated in-memory -- no external files needed.

**Crates**: cyanea-datasets

```toml
[dependencies]
cyanea-datasets = "0.1"
```

```rust
use cyanea_datasets::genomics;
use cyanea_datasets::alignment;
use cyanea_datasets::epigenomics;
use cyanea_datasets::single_cell;
use cyanea_datasets::chemistry;
use cyanea_datasets::phylogenetics;
use cyanea_datasets::metagenomics;
use cyanea_datasets::structural_biology;

fn explore_sample_data() {
    // Genomics: sample sequences, FASTA/FASTQ data, k-mer sets
    let seqs = genomics::sample_sequences();
    let fastq = genomics::sample_fastq();

    // Alignment: pre-aligned sequence pairs, scoring matrices
    let pair = alignment::sample_alignment_pair();
    let msa = alignment::sample_msa();

    // Epigenomics: ChIP-seq peaks, signal tracks
    let peaks = epigenomics::sample_peaks();

    // Single-cell: count matrices, cell metadata
    let counts = single_cell::sample_count_matrix();
    let metadata = single_cell::sample_cell_metadata();

    // Chemistry: SMILES strings, molecular structures
    let molecules = chemistry::sample_molecules();

    // Phylogenetics: Newick trees, distance matrices
    let tree = phylogenetics::sample_newick_tree();

    // Metagenomics: OTU tables, taxonomy
    let otu = metagenomics::sample_otu_table();

    // Structural biology: PDB coordinates, contact maps
    let structure = structural_biology::sample_pdb();
}
```

## 13. Protocol Templates

Access 16 structured protocol templates (10 wet lab, 6 dry lab) for common bioinformatics and laboratory workflows. Each protocol includes categorized steps, difficulty level, estimated time, required materials, tips, cautions, and renders to markdown.

**Crates**: cyanea-datasets

```toml
[dependencies]
cyanea-datasets = "0.1"
```

```rust
use cyanea_datasets::protocols::{
    all_protocols, wet_lab_protocols, dry_lab_protocols,
    elisa, gwas_pipeline, ProtocolCategory,
};

fn explore_protocols() {
    // 1. List all 16 protocol templates
    let all = all_protocols();
    println!("{} protocols available", all.len());  // 16

    // 2. Filter by category: wet lab (10) or dry lab (6)
    let wet = wet_lab_protocols();
    let dry = dry_lab_protocols();
    assert_eq!(wet.len(), 10);
    assert_eq!(dry.len(), 6);

    // 3. Inspect a protocol's metadata
    let protocol = elisa();
    println!("{} (slug: {})", protocol.title, protocol.slug);
    println!("Difficulty: {:?}, Time: {}", protocol.difficulty, protocol.estimated_time);
    println!("Requirements: {}", protocol.requirements.len());

    // 4. Access structured steps
    for step in &protocol.steps {
        println!("  Step {}: {}", step.number, step.title);
        if let Some(dur) = step.duration {
            println!("    Duration: {}", dur);
        }
        for tip in &step.tips {
            println!("    Tip: {}", tip);
        }
        if let Some(caution) = step.caution {
            println!("    Caution: {}", caution);
        }
    }

    // 5. Render a protocol to markdown
    let md = protocol.to_markdown();
    assert!(md.contains("# ELISA"));
    assert!(md.contains("## Steps"));

    // 6. Unique slugs for URL routing
    let gwas = gwas_pipeline();
    assert_eq!(gwas.slug, "gwas-pipeline");
    assert_eq!(gwas.category, ProtocolCategory::DryLab);
}
```

## 14. Clinical Genomics Pipeline

Classify variants using ACMG/AMP guidelines, match against ClinVar, call pharmacogenomics star alleles, determine metabolizer phenotypes, check drug-gene interactions, type HLA alleles, and compute tumor biomarkers (TMB, MSI).

**Crates**: cyanea-omics

```toml
[dependencies]
cyanea-omics = "0.1"
```

```rust
use cyanea_omics::acmg::{auto_evidence, AcmgClassification, AcmgEvidence, parse_clinvar_tsv, match_clinvar};
use cyanea_omics::pharmacogenomics::{
    call_star_alleles, activity_to_phenotype, lookup_drug_interactions,
    demo_cyp2d6_database,
};
use cyanea_omics::clinical::{compute_tmb, call_msi, bethesda_markers, parse_hla_typing, hla_compatibility};

fn clinical_pipeline() {
    // 1. ACMG/AMP variant classification
    //    Provide evidence criteria, get 5-tier classification
    let evidence = AcmgEvidence::new();
    // evidence.add(AcmgCriterion::PS1, EvidenceStrength::Strong);
    let classification = AcmgClassification::from_evidence(&evidence);
    println!("Classification: {:?}", classification.class);  // Pathogenic, Likely Pathogenic, VUS, etc.

    // 2. Auto-evidence from variant properties
    let auto = auto_evidence("chr17", 7674220, "C", "T", Some(0.0001), None);
    println!("Auto-gathered {} criteria", auto.criteria.len());

    // 3. ClinVar matching
    let clinvar_tsv = "...";  // ClinVar TSV data
    let db = parse_clinvar_tsv(clinvar_tsv).unwrap();
    let hit = match_clinvar("chr17", 7674220, "C", "T", &db);

    // 4. Pharmacogenomics: star allele calling with demo CYP2D6 database
    let pgx_db = demo_cyp2d6_database();
    let variants_for_gene = vec![];  // variant positions
    let calls = call_star_alleles(&variants_for_gene, &pgx_db);

    // 5. Metabolizer phenotype from activity scores
    let phenotype = activity_to_phenotype(1.5);  // normal metabolizer
    println!("Phenotype: {:?}", phenotype);

    // 6. Drug-gene interaction lookup
    let interactions = lookup_drug_interactions("CYP2D6", &pgx_db);
    for dgi in &interactions {
        println!("{}: {}", dgi.drug, dgi.recommendation);
    }

    // 7. HLA typing and compatibility
    let typing = parse_hla_typing("A*02:01/A*03:01");
    let compat = hla_compatibility(&typing, &typing);

    // 8. Tumor mutational burden
    let tmb = compute_tmb(150, 30.0);  // 150 somatic mutations, 30 Mb panel
    println!("TMB: {:.1} mut/Mb, category: {:?}", tmb.mutations_per_mb, tmb.category);

    // 9. Microsatellite instability
    let loci = bethesda_markers();
    let msi = call_msi(&loci, 3);  // 3 of 5 unstable
    println!("MSI status: {:?}", msi.status);
}
```

## 15. Spatial Biology Analysis

Analyze spatial transcriptomics data: load platform-specific data (Visium, MERFISH, Slide-seq), segment cells, detect spatial domains, find spatially variable genes, infer cell-cell communication, and deconvolve spot-level expression.

**Crates**: cyanea-omics

```toml
[dependencies]
cyanea-omics = "0.1"
```

```rust
use cyanea_omics::spatial_platforms::{VisiumData, visium_to_spatial_points};
use cyanea_omics::spatial_segmentation::{voronoi_segmentation, expansion_segmentation, ExpansionParams};
use cyanea_omics::spatial_domains::{detect_domains, hmrf_smooth, find_spatially_variable_genes, DomainParams};
use cyanea_omics::spatial_cellchat::{analyze_communication, demo_lr_database, CommParams};
use cyanea_omics::spatial_deconvolution::{nnls_deconvolve, score_enrichment, CellTypeSignature};
use cyanea_omics::spatial::{knn_spatial_neighbors, morans_i, SpatialPoint};

fn spatial_pipeline() {
    // Example data: 5 spots with 3 genes
    let coords = vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (0.0, 1.0), (1.0, 1.0)];
    let expression = vec![
        vec![5.0, 2.0, 8.0], vec![3.0, 7.0, 1.0], vec![6.0, 1.0, 9.0],
        vec![4.0, 6.0, 2.0], vec![5.0, 5.0, 5.0],
    ];
    let gene_names = vec!["GENE1", "GENE2", "GENE3"];

    // 1. Voronoi cell segmentation
    let seeds: Vec<(f64, f64)> = coords.clone();
    let transcripts = coords.clone(); // in practice, individual molecule locations
    let cells = voronoi_segmentation(&seeds, &transcripts, Some(50.0));
    println!("{} cells segmented", cells.len());

    // 2. Spatial domain detection (spatially-aware k-means + HMRF smoothing)
    let params = DomainParams { k: 2, spatial_weight: 0.5, max_iter: 50, seed: 42 };
    let domains = detect_domains(&expression, &coords, &params);
    let neighbors: Vec<Vec<usize>> = vec![vec![1,3], vec![0,2,4], vec![1], vec![0,4], vec![1,3]];
    let smoothed = hmrf_smooth(&domains.labels, &expression, &neighbors, 1.0, 10);

    // 3. Spatially variable genes (Moran's I with permutation test)
    let svgs = find_spatially_variable_genes(&expression, &gene_names, &neighbors, 99, 42);
    for svg in &svgs {
        println!("{}: Moran's I = {:.3}, p = {:.3}", svg.gene, svg.morans_i, svg.p_value);
    }

    // 4. Cell-cell communication (CellChat-style L-R analysis)
    let lr_pairs = demo_lr_database();
    let cell_types = vec!["TypeA", "TypeA", "TypeB", "TypeB", "TypeA"];
    let comm_params = CommParams { max_distance: 200.0, n_permutations: 100, seed: 42 };
    let results = analyze_communication(
        &expression, &gene_names, &cell_types, &coords, &lr_pairs, &comm_params,
    );

    // 5. NNLS deconvolution
    let signatures = vec![
        CellTypeSignature { name: "Tcell".into(), gene_names: vec!["GENE1".into()], expression: vec![8.0] },
        CellTypeSignature { name: "Bcell".into(), gene_names: vec!["GENE2".into()], expression: vec![7.0] },
    ];
    let deconv = nnls_deconvolve(&expression, &gene_names, &signatures);
    for spot in &deconv.proportions {
        println!("Proportions: {:.2}/{:.2}", spot[0], spot[1]);
    }
}
```

## 16. Microarray Analysis Pipeline

Parse microarray data files (Affymetrix CEL, GenePix GPR, Illumina IDAT), normalize expression values with RMA, run limma-style differential expression analysis, and analyze Illumina methylation arrays.

**Crates**: cyanea-io, cyanea-omics

```toml
[dependencies]
cyanea-io = "0.1"
cyanea-omics = "0.1"
```

```rust
use cyanea_io::microarray::{parse_cel_v3, parse_gpr, parse_idat, CelFile, GprFile, IdatFile};
use cyanea_omics::microarray::{
    rma_normalize, quantile_normalize, median_polish,
    limma_diff_expr, DiffExprResult,
    compute_beta, beta_to_m_value, swan_normalize, diff_methylation,
    InfiniumType,
};

fn microarray_pipeline(cel_text: &str, gpr_text: &str, idat_bytes: &[u8]) {
    // 1. Parse Affymetrix CEL v3 file (probe-level intensities)
    let cel = parse_cel_v3(cel_text).unwrap();
    println!("{} probes, algorithm: {}", cel.intensities.len(), cel.header.algorithm);

    // 2. Parse GenePix GPR file (two-color microarray)
    let gpr = parse_gpr(gpr_text).unwrap();
    println!("{} spots in {} blocks", gpr.spots.len(), gpr.header.len());

    // 3. Parse Illumina IDAT file (binary BeadArray)
    let idat = parse_idat(idat_bytes).unwrap();
    println!("{} probes, {} fields", idat.n_probes, idat.fields.len());

    // 4. RMA normalization (log2 background correction + quantile normalization)
    let probe_data = vec![vec![100.0, 200.0, 50.0], vec![150.0, 180.0, 60.0]];
    let rma = rma_normalize(&probe_data);
    println!("RMA: {} samples normalized", rma.len());

    // 5. Quantile normalization (rank-based cross-sample normalization)
    let mut data = vec![vec![5.0, 2.0, 3.0], vec![4.0, 1.0, 6.0]];
    quantile_normalize(&mut data);

    // 6. Limma-style differential expression (empirical Bayes moderated t-test)
    let expression = vec![
        vec![5.2, 5.1, 3.0, 2.9],  // gene1: 2 control, 2 treated
        vec![2.0, 2.1, 2.0, 2.1],  // gene2: no change
    ];
    let gene_names = vec!["GeneA".into(), "GeneB".into()];
    let groups = vec![0, 0, 1, 1];  // 0=control, 1=treated
    let results = limma_diff_expr(&expression, &gene_names, &groups);
    for r in &results {
        println!("{}: logFC={:.2}, adj_p={:.4}", r.gene, r.log_fc, r.adj_p_value);
    }

    // 7. Methylation analysis: beta values, M-values, SWAN normalization
    let methylated = 800.0;
    let unmethylated = 200.0;
    let beta = compute_beta(methylated, unmethylated);  // 0.8
    let m_val = beta_to_m_value(beta);                  // log2(0.8/0.2)

    let beta_values = vec![vec![0.1, 0.9, 0.5], vec![0.2, 0.8, 0.4]];
    let design_types = vec![InfiniumType::TypeI, InfiniumType::TypeII, InfiniumType::TypeI];
    let normalized_betas = swan_normalize(&beta_values, &design_types);

    // 8. Differential methylation
    let probe_ids = vec!["cg001".into(), "cg002".into(), "cg003".into()];
    let dm_results = diff_methylation(&beta_values, &probe_ids, &groups);
}
```

## 17. Flow Cytometry (FCS) Analysis

Parse FCS files (versions 2.0, 3.0, 3.1), inspect parameter metadata, extract channel data, and write FCS output.

**Crates**: cyanea-io

```toml
[dependencies]
cyanea-io = "0.1"
```

```rust
use cyanea_io::fcs::{parse_fcs, write_fcs, fcs_stats, FcsFile};

fn fcs_pipeline(fcs_data: &[u8]) {
    // 1. Parse an FCS file (supports FCS 2.0, 3.0, 3.1)
    let fcs = parse_fcs(fcs_data).unwrap();
    println!("FCS version: {}", fcs.version);
    println!("{} events, {} parameters", fcs.n_events(), fcs.n_parameters());

    // 2. Inspect parameter metadata
    for param in &fcs.parameters {
        println!("  {}: bits={}, range={}, amplification={:?}",
            param.name, param.bits, param.range, param.amplification);
        if let Some(ref stain) = param.stain {
            println!("    stain: {}", stain);
        }
    }

    // 3. Access TEXT segment keywords
    if let Some(date) = fcs.keyword("$DATE") {
        println!("Acquisition date: {}", date);
    }

    // 4. Extract channel data by name or index
    let fsc_data = fcs.parameter_data_by_name("FSC-A");
    let ssc_data = fcs.parameter_data(1);  // by index

    // 5. Summary statistics
    let stats = fcs.stats();
    println!("{} parameters: {:?}", stats.n_parameters, stats.parameter_names);

    // 6. Write an FCS 3.1 file (round-trip)
    let output = write_fcs(&fcs).unwrap();
    let fcs2 = parse_fcs(&output).unwrap();
    assert_eq!(fcs2.n_events(), fcs.n_events());
}
```

## 18. Metabolomics

Identify metabolites by exact mass matching, generate theoretical isotope patterns, predict retention times, and run KEGG pathway enrichment on detected metabolites.

**Crates**: cyanea-chem

```toml
[dependencies]
cyanea-chem = "0.1"
```

```rust
use cyanea_chem::metabolomics::{
    match_by_mass, positive_adducts, negative_adducts, calc_mz,
    isotope_pattern, isotope_cosine_score,
    predict_rt,
    pathway_enrichment, demo_metabolite_database, demo_metabolic_pathways,
};

fn metabolomics_pipeline() {
    // 1. Match observed m/z against a metabolite database
    let db = demo_metabolite_database();  // 20 common metabolites
    let adducts = positive_adducts();     // [M+H]+, [M+Na]+, [M+K]+, etc.

    let observed_mz = 181.0707;  // glucose [M+H]+
    let matches = match_by_mass(observed_mz, &db, &adducts, 10.0);  // 10 ppm tolerance
    for m in &matches {
        println!("{} ({}): {:.1} ppm, adduct {}, calc m/z {:.4}",
            m.metabolite.name, m.metabolite.formula,
            m.ppm_error, m.adduct, m.calc_mz);
    }

    // 2. Generate theoretical isotope pattern from formula
    let pattern = isotope_pattern("C6H12O6", 4).unwrap();  // glucose, M+0 through M+3
    for peak in &pattern {
        println!("  M+{:.3}: abundance {:.3}", peak.mass_offset, peak.abundance);
    }

    // Compare observed vs theoretical isotope patterns (cosine score)
    let observed = vec![1.0, 0.068, 0.003];
    let theoretical: Vec<f64> = pattern.iter().take(3).map(|p| p.abundance).collect();
    let score = isotope_cosine_score(&observed, &theoretical);
    println!("Isotope match score: {:.3}", score);  // > 0.99 = good match

    // 3. Predict reversed-phase HPLC retention time
    let rt = predict_rt(0.5, 180.0, 110.0);  // logP, MW, PSA
    println!("Predicted RT: {:.1} +/- {:.1} min", rt.rt_minutes, rt.error_margin);

    // 4. KEGG pathway enrichment (hypergeometric test)
    let pathways = demo_metabolic_pathways();  // 6 core KEGG pathways
    let matched_ids = vec!["C00031", "C00022", "C00186", "C00074", "C00024"];  // glycolysis hits
    let enriched = pathway_enrichment(&matched_ids, &pathways, 20);
    for pw in &enriched {
        println!("{}: {}/{} hits, p={:.2e}, impact={:.2}",
            pw.pathway_name, pw.hits, pw.total, pw.p_value, pw.impact);
    }
}
```

## 19. CRISPR Guide Design and Screen Analysis

Score guide RNAs (Rule Set 2 and CFD), search for off-targets with mismatches, run MAGeCK-style CRISPR screen analysis, and predict base editing outcomes.

**Crates**: cyanea-omics

```toml
[dependencies]
cyanea-omics = "0.1"
```

```rust
use cyanea_omics::crispr::{
    score_guide_rs2, cfd_score, find_off_targets, count_mismatches,
    analyze_screen, predict_editing, BaseEditor,
};

fn crispr_pipeline() {
    // 1. Score guide RNA on-target activity (Rule Set 2)
    //    Requires 30-nt context: 4nt upstream + 20nt spacer + 3nt PAM + 3nt downstream
    let context = b"ACGTGCATGCTAGCTAGCGATGGNNGGATT";
    let rs2_score = score_guide_rs2(context).unwrap();
    println!("Rule Set 2 score: {:.3}", rs2_score);  // 0-1, higher = better

    // 2. CFD off-target scoring (cutting frequency determination)
    let guide     = b"ACGTACGTACGTACGTACGT";
    let off_target = b"ACGTACGTACGTACGTACGA";  // 1 mismatch at position 20
    let cfd = cfd_score(guide, off_target).unwrap();
    println!("CFD score: {:.3}", cfd);  // 1.0 = perfect match, lower = less likely to cut
    println!("Mismatches: {}", count_mismatches(guide, off_target));

    // 3. Genome-wide off-target search (finds NGG PAM sites with up to N mismatches)
    let genome = b"NNNNNACGTACGTACATACGTACGTNGGNNNNNN";  // 1 mismatch near guide
    let off_targets = find_off_targets(guide, genome, "chr1", 3).unwrap();
    for ot in &off_targets {
        println!("{}:{} strand={}, mm={}, CFD={:.3}",
            ot.chrom, ot.position, ot.strand, ot.mismatches, ot.cfd_score);
    }

    // 4. MAGeCK-style CRISPR screen analysis (robust rank aggregation)
    let screen_data = vec![
        ("g1".into(), "TP53".into(), -3.5),   // guide, gene, log2FC
        ("g2".into(), "TP53".into(), -2.8),
        ("g3".into(), "TP53".into(), -4.1),
        ("g4".into(), "GFP".into(),   0.1),   // control gene
        ("g5".into(), "GFP".into(),  -0.2),
        ("g6".into(), "GFP".into(),   0.3),
    ];
    let results = analyze_screen(&screen_data, true);  // true = negative selection
    for r in &results {
        println!("{}: median_lfc={:.2}, RRA={:.4}, FDR={:.4}, class={}",
            r.gene, r.median_lfc, r.rra_score, r.fdr, r.classification);
    }

    // 5. Base editing outcome prediction (CBE: C→T, ABE: A→G)
    let spacer = b"AAGCACGTACGTACGTACGT";
    let cbe_outcomes = predict_editing(BaseEditor::Cbe, spacer).unwrap();
    for o in &cbe_outcomes {
        println!("pos {}: {}→{}, efficiency={:.2}, in_window={}",
            o.position, o.ref_base, o.alt_base, o.efficiency, o.in_window);
    }

    let abe_outcomes = predict_editing(BaseEditor::Abe, spacer).unwrap();
    for o in &abe_outcomes {
        println!("pos {}: {}→{}, efficiency={:.2}", o.position, o.ref_base, o.alt_base, o.efficiency);
    }
}
```

## 20. Hi-C / 3D Genome Analysis

Parse Hi-C contact data, build contact matrices with ICE/KR balancing, call TADs via insulation scores, identify A/B compartments, and detect chromatin loops.

**Crates**: cyanea-omics

```toml
[dependencies]
cyanea-omics = "0.1"
```

```rust
use cyanea_omics::hic::{
    ContactMatrix, contacts_to_matrix, parse_pairs, write_pairs, parse_cool_text,
    call_tads, TadParams, insulation_scores,
    call_compartments, Compartment,
    call_loops, LoopParams,
};

fn hic_pipeline() {
    // 1. Build a contact matrix from sparse contacts
    let n = 20;
    let mut matrix = ContactMatrix::new("chr1", n, 10000);  // 20 bins, 10 kb resolution

    // Add contacts (symmetric for cis)
    matrix.add_contact(0, 1, 5.0);
    matrix.add_contact(1, 2, 8.0);
    matrix.add_contact(5, 6, 12.0);
    println!("{}x{} matrix, resolution {} bp", matrix.size, matrix.size, matrix.resolution);

    // 2. ICE or KR balancing (iterative correction)
    let mut balanced = matrix.clone();
    balanced.ice_balance(50, 1e-5);   // 50 iterations, tolerance 1e-5
    // Or: balanced.kr_balance(100, 1e-6);

    // 3. Observed/Expected normalization
    let oe = balanced.observed_expected();

    // 4. Call TADs via insulation scores
    let tad_params = TadParams { window_size: 5, min_size: 3, zscore_threshold: -0.5 };
    let tads = call_tads(&balanced, &tad_params);
    for tad in &tads {
        println!("TAD: bins {}-{}, mean insulation {:.3}", tad.start, tad.end, tad.mean_insulation);
    }

    // Raw insulation scores
    let scores = insulation_scores(&balanced, 5);

    // 5. A/B compartment calling (PCA on correlation of O/E matrix)
    let gc_content: Option<&[f64]> = None;  // optional GC for orientation
    let compartments = call_compartments(&balanced, gc_content);
    let n_a = compartments.compartments.iter().filter(|c| **c == Compartment::A).count();
    let n_b = compartments.compartments.iter().filter(|c| **c == Compartment::B).count();
    println!("{} A bins, {} B bins", n_a, n_b);

    // 6. Loop detection (donut background enrichment)
    let loop_params = LoopParams {
        donut_inner: 2, donut_outer: 5,
        min_distance: 5, p_threshold: 0.01, min_enrichment: 1.5,
    };
    let loops = call_loops(&balanced, &loop_params);
    for lp in &loops {
        println!("Loop ({},{}) enrichment={:.2}, p={:.2e}", lp.bin1, lp.bin2, lp.enrichment, lp.p_value);
    }

    // 7. Parse/write .pairs format
    let pairs_text = "## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\nr1\tchr1\t1000\tchr1\t5000\t+\t-\n";
    let contacts = parse_pairs(pairs_text).unwrap();
    let output = write_pairs(&contacts, "test");

    // 8. Parse cooler text dump
    let cool_text = "chr1\t0\t10000\tchr1\t10000\t20000\t5.0\n";
    let sparse = parse_cool_text(cool_text, 10000).unwrap();
    let mat = contacts_to_matrix(&sparse.contacts, "chr1", 20, 10000);
}
```
