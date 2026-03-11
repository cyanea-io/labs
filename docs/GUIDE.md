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

## 13. Clinical Genomics Pipeline

Classify variants using ACMG/AMP guidelines, match against ClinVar, call pharmacogenomics star alleles, determine metabolizer phenotypes, check drug-gene interactions, type HLA alleles, and compute tumor biomarkers (TMB, MSI).

**Crates**: cyanea-omics, cyanea-io (with `vcf` feature)

```toml
[dependencies]
cyanea-omics = "0.1"
cyanea-io = { version = "0.1", features = ["vcf"] }
```

```rust
use cyanea_io::vcf::parse_vcf;
use cyanea_omics::acmg::{classify_variant, AcmgEvidence, ClinVarDB, clinvar_match};
use cyanea_omics::pharmacogenomics::{
    call_star_alleles, metabolizer_phenotype, drug_gene_interactions,
    StarAlleleDB,
};
use cyanea_omics::clinical::{hla_type, tumor_mutational_burden, microsatellite_instability};

fn clinical_pipeline(
    vcf_data: &str,
    clinvar_db: &ClinVarDB,
    pgx_db: &StarAlleleDB,
    tumor_variants: &[cyanea_omics::variant::Variant],
    normal_variants: &[cyanea_omics::variant::Variant],
    msi_loci: &[(String, u64, &[u8])],
) {
    // 1. Parse variants
    let variants = parse_vcf(vcf_data).unwrap();

    // 2. ACMG/AMP variant classification
    for variant in &variants {
        let evidence = AcmgEvidence::gather(variant, clinvar_db).unwrap();
        let classification = classify_variant(&evidence).unwrap();
        println!("{}:{} {} -> {} ({} criteria met)",
            variant.chrom, variant.pos, variant.alt_allele,
            classification.significance, classification.criteria.len());
    }

    // 3. ClinVar matching
    for variant in &variants {
        if let Some(entry) = clinvar_match(variant, clinvar_db).unwrap() {
            println!("{}:{} ClinVar {}: {} (review: {})",
                variant.chrom, variant.pos, entry.accession,
                entry.significance, entry.review_status);
        }
    }

    // 4. Pharmacogenomics: star allele calling
    let star_alleles = call_star_alleles(&variants, pgx_db).unwrap();
    for sa in &star_alleles {
        println!("Gene {}: diplotype {} / {}",
            sa.gene, sa.allele1, sa.allele2);
    }

    // 5. Metabolizer phenotype prediction
    for sa in &star_alleles {
        let phenotype = metabolizer_phenotype(sa).unwrap();
        println!("{}: {} (activity score: {:.1})",
            sa.gene, phenotype.phenotype, phenotype.activity_score);
    }

    // 6. Drug-gene interaction lookup
    for sa in &star_alleles {
        let interactions = drug_gene_interactions(&sa.gene, pgx_db).unwrap();
        for dgi in &interactions {
            println!("{} + {}: {} — {}",
                dgi.drug, dgi.gene, dgi.recommendation, dgi.evidence_level);
        }
    }

    // 7. HLA typing
    let hla_results = hla_type(&variants).unwrap();
    for hla in &hla_results {
        println!("HLA-{}: {} / {}", hla.gene, hla.allele1, hla.allele2);
    }

    // 8. Tumor mutational burden (TMB)
    let tmb = tumor_mutational_burden(tumor_variants, normal_variants, 30.0).unwrap();
    println!("TMB: {:.1} mutations/Mb ({})",
        tmb.mutations_per_mb, tmb.classification);

    // 9. Microsatellite instability (MSI)
    let msi = microsatellite_instability(msi_loci, tumor_variants).unwrap();
    println!("MSI score: {:.2}, status: {} ({}/{} unstable loci)",
        msi.score, msi.status, msi.unstable_loci, msi.total_loci);
}
```

## 14. Spatial Biology Analysis

Analyze spatial transcriptomics data: load platform-specific data (Visium, MERFISH, Slide-seq), segment cells, detect spatial domains, find spatially variable genes, infer cell-cell communication, and deconvolve spot-level expression.

**Crates**: cyanea-omics (with `single-cell` feature), cyanea-stats

```toml
[dependencies]
cyanea-omics = { version = "0.1", features = ["single-cell"] }
cyanea-stats = "0.1"
```

```rust
use cyanea_omics::spatial::{
    VisiumData, MerfishData, SlideSeqData,
    voronoi_segmentation, nucleus_expansion, watershed_segmentation,
    spatial_kmeans, hmrf_smooth,
    morans_i, spatially_variable_genes,
    LigandReceptorDB, spatial_communication, CellChatConfig,
    nnls_deconvolution, enrichment_scoring,
};

fn spatial_pipeline(
    coordinates: &[(f64, f64)],
    counts: &[Vec<f64>],
    gene_names: &[&str],
    cell_types: &[&str],
    reference_profiles: &[Vec<f64>],
) {
    // 1. Load platform-specific spatial data
    let visium = VisiumData::from_coordinates_and_counts(coordinates, counts).unwrap();

    // 2. Cell segmentation (choose one strategy)
    let voronoi_cells = voronoi_segmentation(coordinates).unwrap();
    let expanded = nucleus_expansion(coordinates, 15.0).unwrap();  // 15 um radius
    let watershed_cells = watershed_segmentation(coordinates, counts).unwrap();

    // 3. Spatial domain detection
    //    Spatially-aware k-means incorporates coordinate proximity
    let domains = spatial_kmeans(coordinates, counts, 8, 0.5).unwrap();  // 8 domains, alpha=0.5
    //    HMRF smoothing refines cluster boundaries using neighbor context
    let smoothed = hmrf_smooth(&domains, coordinates, 5, 0.3).unwrap();  // 5 iters, beta=0.3
    println!("{} spatial domains detected", smoothed.n_domains);

    // 4. Spatially variable gene (SVG) detection via Moran's I
    let svg_results = spatially_variable_genes(coordinates, counts, gene_names).unwrap();
    for svg in svg_results.iter().filter(|s| s.fdr < 0.05).take(20) {
        println!("SVG: {} (Moran's I={:.3}, FDR={:.2e})", svg.gene, svg.morans_i, svg.fdr);
    }

    // Single-gene Moran's I
    let mi = morans_i(coordinates, &counts[0]).unwrap();
    println!("Moran's I = {:.4}, p = {:.4}", mi.statistic, mi.p_value);

    // 5. Cell-cell communication (CellChat-style)
    //    Uses a ligand-receptor database with multi-subunit complex support
    let lr_db = LigandReceptorDB::builtin();
    let config = CellChatConfig {
        max_distance: 200.0,   // um
        min_expression: 0.1,
        ..Default::default()
    };
    let interactions = spatial_communication(
        coordinates, counts, gene_names, cell_types, &lr_db, &config,
    ).unwrap();
    for inter in interactions.iter().take(10) {
        println!("{} -> {} via {}-{}: score={:.3}",
            inter.sender, inter.receiver,
            inter.ligand, inter.receptor, inter.score);
    }

    // 6. Deconvolution of spot-level expression
    //    NNLS: non-negative least squares for cell-type proportions
    let nnls_props = nnls_deconvolution(counts, reference_profiles).unwrap();
    println!("NNLS deconvolution: {} spots x {} cell types",
        nnls_props.len(), nnls_props[0].len());

    //    Enrichment scoring: rank-based cell-type enrichment per spot
    let enrichment = enrichment_scoring(counts, reference_profiles, gene_names).unwrap();
    for (i, scores) in enrichment.iter().enumerate().take(5) {
        println!("Spot {}: top cell type = {} (score={:.3})",
            i, cell_types[scores.top_type_index], scores.top_score);
    }
}
```

## 15. Microarray Analysis Pipeline

Parse microarray data files (Affymetrix CEL, GenePix GPR, Illumina IDAT), normalize expression values with RMA, run limma-style differential expression analysis, and analyze Illumina methylation arrays.

**Crates**: cyanea-io, cyanea-omics

```toml
[dependencies]
cyanea-io = { version = "0.1", features = ["cel", "gpr", "idat"] }
cyanea-omics = "0.1"
```

```rust
use cyanea_io::cel::parse_cel;
use cyanea_io::gpr::parse_gpr;
use cyanea_io::idat::parse_idat;
use cyanea_omics::microarray::{
    rma_normalize, limma_differential_expression, LimmaConfig,
    methylation_array_analysis, MethylationArrayConfig,
};

fn microarray_pipeline(
    cel_data: &[u8],
    gpr_data: &str,
    idat_data: &[u8],
) {
    // 1. Parse Affymetrix CEL file (probe-level intensities)
    let cel = parse_cel(cel_data).unwrap();
    println!("{} probes, {} cells", cel.num_probes(), cel.num_cells());

    // 2. Parse GenePix GPR file (two-color microarray)
    let gpr = parse_gpr(gpr_data).unwrap();
    println!("{} spots, {} blocks", gpr.num_spots(), gpr.num_blocks());

    // 3. Parse Illumina IDAT file (BeadArray intensities)
    let idat = parse_idat(idat_data).unwrap();
    println!("{} probes in IDAT", idat.num_probes());

    // 4. RMA normalization (background correction + quantile normalization + summarization)
    let expression_matrix = vec![cel.intensities()];  // multiple CEL files in practice
    let normalized = rma_normalize(&expression_matrix).unwrap();
    println!("RMA normalized: {} genes x {} samples",
        normalized.num_genes(), normalized.num_samples());

    // 5. Limma-style differential expression
    let config = LimmaConfig {
        group_labels: vec!["control", "treated", "control", "treated"],
        contrast: ("treated", "control"),
        ..Default::default()
    };
    let de_results = limma_differential_expression(&normalized, &config).unwrap();
    let significant: Vec<_> = de_results.iter()
        .filter(|r| r.adj_p_value < 0.05 && r.log_fold_change.abs() > 1.0)
        .collect();
    println!("{} DE genes (|logFC| > 1, FDR < 0.05)", significant.len());

    // 6. Illumina methylation array analysis (450K/EPIC)
    let meth_config = MethylationArrayConfig::default();
    let meth_results = methylation_array_analysis(&idat, &meth_config).unwrap();
    println!("{} CpG sites, {} DMPs (FDR < 0.05)",
        meth_results.num_sites(),
        meth_results.differentially_methylated(0.05).len());
}
```

## 16. Hi-C / 3D Genome Analysis

Parse Hi-C contact matrices, call topologically associating domains (TADs), identify A/B compartments, and detect chromatin loops.

**Crates**: cyanea-omics

```toml
[dependencies]
cyanea-omics = "0.1"
```

```rust
use cyanea_omics::hic::{
    ContactMatrix, load_contact_matrix, normalize_contacts, NormMethod,
    call_tads, TadConfig, ab_compartments, CompartmentConfig,
    detect_loops, LoopConfig,
};

fn hic_pipeline(
    matrix_data: &[Vec<f64>],
    resolution: u64,
    chrom: &str,
) {
    // 1. Load and normalize Hi-C contact matrix
    let matrix = ContactMatrix::from_dense(matrix_data, chrom, resolution).unwrap();
    let normalized = normalize_contacts(&matrix, NormMethod::KR).unwrap();
    println!("Contact matrix: {}x{} at {} bp resolution",
        normalized.dim(), normalized.dim(), resolution);

    // 2. Call topologically associating domains (TADs)
    let tad_config = TadConfig {
        min_size: 3,           // minimum 3 bins
        max_size: 200,         // maximum 200 bins
        ..Default::default()
    };
    let tads = call_tads(&normalized, &tad_config).unwrap();
    println!("{} TADs called, median size: {} bins",
        tads.len(), tads.median_size());

    for tad in tads.iter().take(5) {
        println!("TAD {}:{}-{} ({} bins, insulation score={:.3})",
            chrom, tad.start, tad.end, tad.n_bins, tad.insulation_score);
    }

    // 3. A/B compartment identification via eigenvector decomposition
    let comp_config = CompartmentConfig::default();
    let compartments = ab_compartments(&normalized, &comp_config).unwrap();
    let n_a = compartments.iter().filter(|c| c.is_a()).count();
    let n_b = compartments.iter().filter(|c| c.is_b()).count();
    println!("{} A compartments, {} B compartments", n_a, n_b);

    for comp in compartments.iter().take(5) {
        println!("Bin {}: compartment {} (eigenvector={:.4})",
            comp.bin_index, comp.label(), comp.eigenvector_value);
    }

    // 4. Chromatin loop detection
    let loop_config = LoopConfig {
        min_distance: 5,       // minimum 5 bins apart
        fdr_threshold: 0.05,
        ..Default::default()
    };
    let loops = detect_loops(&normalized, &loop_config).unwrap();
    println!("{} loops detected (FDR < {})", loops.len(), loop_config.fdr_threshold);

    for lp in loops.iter().take(5) {
        println!("Loop {}:{}-{} (observed/expected={:.2}, FDR={:.2e})",
            chrom, lp.anchor1, lp.anchor2,
            lp.observed_over_expected, lp.fdr);
    }
}
```
