# cyanea-meta Usage Guide

## Installation

```toml
[dependencies]
cyanea-meta = { version = "0.1", features = ["std"] }
```

## Building a Taxonomy and Classifying Reads

```rust
use cyanea_meta::taxonomy::{TaxonomyDB, TaxonRank};

// Build a small taxonomy tree
let mut db = TaxonomyDB::new(1); // root taxid = 1
db.add_node(2, 1, TaxonRank::Domain, "Bacteria").unwrap();
db.add_node(1224, 2, TaxonRank::Phylum, "Proteobacteria").unwrap();
db.add_node(1236, 1224, TaxonRank::Class, "Gammaproteobacteria").unwrap();
db.add_node(91347, 1236, TaxonRank::Order, "Enterobacterales").unwrap();
db.add_node(543, 91347, TaxonRank::Family, "Enterobacteriaceae").unwrap();
db.add_node(561, 543, TaxonRank::Genus, "Escherichia").unwrap();
db.add_node(562, 561, TaxonRank::Species, "Escherichia coli").unwrap();

// Index reference sequences (associates k-mers with taxon IDs)
db.add_reference(b"ACGTACGTACGTACGTACGTACGTACGTACGTACG", 562).unwrap();

// Classify a metagenomic read
if let Ok(Some(taxid)) = db.classify_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACG") {
    let lineage = db.get_lineage(taxid).unwrap();
    println!("Classified to taxid {}, lineage: {:?}", taxid, lineage);
}

// Lowest common ancestor
let ancestor = db.lca(562, 561).unwrap();
assert_eq!(ancestor, 561); // Escherichia genus
```

## Taxonomic Profiling

```rust
use cyanea_meta::profile::*;

// Build a profile from classified reads
let classifications = vec![562, 562, 562, 561, 561, 1224];
let profile = profile_from_classifications(&classifications).unwrap();
println!("Taxa: {}, Total reads: {}", profile.richness(), profile.total_count());

// Normalize to relative abundances
let normalized = normalize_profile(&profile).unwrap();

// Filter low-abundance taxa
let filtered = filter_profile(&normalized, 0.1).unwrap();

// Merge profiles from multiple samples
let merged = merge_profiles(&[&profile, &profile]).unwrap();
```

## Alpha and Beta Diversity

```rust
use cyanea_meta::diversity::*;

// Alpha diversity from species abundance counts
let counts = vec![100, 50, 25, 10, 5, 3, 2, 1, 1, 1];
let ad = alpha_diversity(&counts).unwrap();
println!("Shannon: {:.3}", ad.shannon);
println!("Simpson: {:.3}", ad.simpson);
println!("Chao1: {:.1}", ad.chao1);
println!("Observed species: {}", ad.observed_species);

// Rarefaction curve
let steps = vec![10, 50, 100, 150, 198];
let curve = rarefaction_curve(&counts, &steps).unwrap();
for (depth, expected_species) in &curve {
    println!("depth={}: {:.1} species", depth, expected_species);
}

// Rarefy to even depth
let rarefied = rarefy(&counts, 50).unwrap();

// Beta diversity (Bray-Curtis) between samples
let sample1 = vec![100, 50, 25];
let sample2 = vec![80, 60, 30];
let sample3 = vec![10, 5, 200];
let beta = beta_diversity_matrix(&[&sample1, &sample2, &sample3]).unwrap();
println!("BC(1,2): {:.3}", beta.get(0, 1));
```

## Compositional Data Analysis

```rust
use cyanea_meta::composition::*;

// CLR transform (handles compositional data properly)
let samples = vec![
    vec![100.0, 50.0, 25.0, 10.0],
    vec![80.0, 60.0, 30.0, 15.0],
];
let clr = clr_transform(&samples).unwrap();

// Differential abundance (ALDEx2-style)
let group1 = vec![
    vec![100.0, 50.0, 25.0],
    vec![110.0, 45.0, 30.0],
];
let group2 = vec![
    vec![20.0, 150.0, 25.0],
    vec![25.0, 140.0, 28.0],
];
let taxa = vec!["Species_A".into(), "Species_B".into(), "Species_C".into()];
let results = differential_abundance(&group1, &group2, &taxa).unwrap();
for r in &results {
    println!("{}: log2FC={:.2}, q={:.4}", r.taxon, r.mean_diff, r.q_value);
}

// ANCOM
let ancom_results = ancom(&group1, &group2).unwrap();
for r in &ancom_results {
    println!("Taxon {}: W={}, significant={}", r.taxon_index, r.w_statistic, r.is_significant);
}
```

## Metagenomic Binning

```rust
use cyanea_meta::binning::*;

// Compute tetranucleotide frequency
let tnf = tetranucleotide_frequency(b"ACGTACGTACGTACGTACGTACGTACGT").unwrap();
assert_eq!(tnf.len(), 256);

// Create contigs and bin by TNF + coverage
let contig1 = Contig::new("contig_1", b"ACGTACGTACGTACGT...", 50.0).unwrap();
let contig2 = Contig::new("contig_2", b"TGCATGCATGCATGCA...", 48.0).unwrap();
let bins = bin_contigs(&[contig1, contig2], 2).unwrap();

// Filter bins by quality
let good_bins = filter_bins(&bins, 50.0, 10.0).unwrap();
```

## Assembly Quality Control

```rust
use cyanea_meta::assembly::*;

let contigs: Vec<&[u8]> = vec![
    b"ACGTACGTACGTACGTACGTACGTACGT",  // 28 bp
    b"ACGTACGT",                        // 8 bp
    b"ACGTACGTACGTACGT",                // 16 bp
];

let stats = assembly_stats(&contigs).unwrap();
println!("Total length: {} bp", stats.total_length);
println!("Contigs: {}", stats.contig_count);
println!("N50: {} bp", stats.n50);
println!("L50: {}", stats.l50);
println!("GC: {:.1}%", stats.gc_content * 100.0);
println!("Longest: {} bp", stats.longest);

// Custom Nx value
let (n75, l75) = nx_values(&contigs, 0.75).unwrap();
println!("N75: {} bp, L75: {}", n75, l75);
```
