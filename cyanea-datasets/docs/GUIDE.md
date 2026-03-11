# cyanea-datasets Usage Guide

## Overview

`cyanea-datasets` provides small, self-contained demo datasets for every domain in the Cyanea bioinformatics ecosystem. All data is generated in-memory with no external files. Designed for notebooks ("open and run"), integration tests, and tutorials.

## 1. Genomics

### Sequence data

```rust
use cyanea_datasets::genomics;

// E. coli 16S rRNA gene (1,542 bp)
let (name, seq) = genomics::ecoli_16s_rrna();
println!("{}: {} bp", name, seq.len());

// SARS-CoV-2 Spike RBD — nucleotide (~670 bp)
let (name, seq) = genomics::sars_cov2_spike_rbd();
println!("{}: {} bp", name, seq.len());

// SARS-CoV-2 Spike RBD — protein (~222 aa)
let (name, seq) = genomics::sars_cov2_spike_rbd_protein();
println!("{}: {} aa", name, seq.len());
```

### Variants and VCF

```rust
use cyanea_datasets::genomics;

// Structured variant records
let variants = genomics::demo_variants_chr22();
for v in &variants {
    println!("{} {} {} {}/{} AF={:.4}", v.chrom, v.pos, v.id, v.ref_a, v.alt, v.af);
}

// Format individual variants as VCF lines
let line = variants[0].to_vcf_line();
println!("{}", line);

// Full VCF string (header + 10 data lines)
let vcf = genomics::demo_vcf_chr22();
assert!(vcf.starts_with("##fileformat=VCFv4.2"));
```

### Gene annotations

```rust
use cyanea_datasets::genomics;

let tp53 = genomics::tp53_gene();
println!("{} ({}:{}-{}, strand {}, {} exons)",
    tp53.name, tp53.chrom, tp53.start, tp53.end, tp53.strand, tp53.exon_count);

let brca1 = genomics::brca1_gene();
println!("{}: {}", brca1.name, brca1.description);
```

## 2. Alignment

### Multiple alignment input

```rust
use cyanea_datasets::alignment;

// SARS-CoV-2 spike variants for MSA demo
let seqs = alignment::spike_alignment_seqs();
for (name, seq) in &seqs {
    println!("{}: {} aa", name, seq.len());
}

// Cross-species hemoglobin alpha (~141 aa each)
let hba = alignment::hemoglobin_alpha();
assert_eq!(hba.len(), 5);
```

### Pre-computed alignment

```rust
use cyanea_datasets::alignment;

let aln = alignment::demo_pairwise_alignment();
println!("{} vs {}", aln.seq_a_name, aln.seq_b_name);
println!("Score: {}, Identity: {:.1}%, Gaps: {}",
    aln.score, aln.identity * 100.0, aln.gaps);
// Both aligned sequences have equal length (including gaps)
assert_eq!(aln.aligned_a.len(), aln.aligned_b.len());
```

### CIGAR strings

```rust
use cyanea_datasets::alignment;

let cigars = alignment::demo_cigar_strings();
for (label, cigar) in &cigars {
    println!("{}: {}", label, cigar);
}
// perfect_match: 100M
// with_insertion: 50M2I48M
// spliced_rna: 75M5000N25M
```

## 3. Epigenomics

### ChIP-seq peaks

```rust
use cyanea_datasets::epigenomics;

let peaks = epigenomics::demo_chipseq_peaks();
for p in &peaks {
    println!("{} {}:{}-{} signal={:.1} pvalue={:.2e}",
        p.name, p.chrom, p.start, p.end, p.signal, p.pvalue);
}

// Format as narrowPeak or BED
let np_line = peaks[0].to_narrowpeak();
let bed_line = peaks[0].to_bed();

// Full narrowPeak output (10 lines)
let narrowpeak = epigenomics::demo_narrowpeak();
```

### CpG methylation

```rust
use cyanea_datasets::epigenomics;

let sites = epigenomics::demo_cpg_methylation();
for s in &sites {
    println!("{} pos={} beta={:.3} ({}/{})",
        s.chrom, s.pos, s.beta, s.methylated, s.total);
}
// Shows a gradient from unmethylated to methylated across a promoter region
```

## 4. Single-Cell

### PBMC expression matrix

```rust
use cyanea_datasets::single_cell;

let data = single_cell::demo_pbmc_50();
println!("{} genes x {} cells", data.num_genes(), data.num_cells());

// Cell types: 20 T cells, 15 B cells, 15 Monocytes
for ct in &["T_cell", "B_cell", "Monocyte"] {
    let count = data.cell_types.iter().filter(|t| t.as_str() == *ct).count();
    println!("{}: {} cells", ct, count);
}
```

### Gene expression queries

```rust
use cyanea_datasets::single_cell;

let data = single_cell::demo_pbmc_50();

// Look up expression for a specific gene
if let Some(expr) = data.gene_expression("CD3D") {
    let t_mean: f64 = expr[0..20].iter().sum::<f64>() / 20.0;
    let b_mean: f64 = expr[20..35].iter().sum::<f64>() / 15.0;
    println!("CD3D mean: T cells={:.1}, B cells={:.1}", t_mean, b_mean);
}

// QC metrics
let total_counts = data.total_counts_per_cell();
let genes_detected = data.genes_per_cell();
println!("Median UMI/cell: {:.0}", total_counts[25]);
println!("Median genes/cell: {}", genes_detected[25]);
```

## 5. Chemistry

### Drug library

```rust
use cyanea_datasets::chemistry;

let drugs = chemistry::fda_approved_drugs();
for mol in &drugs {
    let status = if mol.lipinski_compliant() { "PASS" } else { "FAIL" };
    println!("{}: MW={:.1}, LogP={:.1}, Lipinski={} ({})",
        mol.name, mol.molecular_weight, mol.logp, status, mol.therapeutic_area);
}

// Count violations for a specific drug
let doxorubicin = &drugs[9];
println!("{}: {} Lipinski violations", doxorubicin.name, doxorubicin.lipinski_violations());
```

### SDF output

```rust
use cyanea_datasets::chemistry;

let sdf = chemistry::demo_sdf_summary();
let record_count = sdf.matches("$$$$").count();
println!("SDF records: {}", record_count); // 12
```

## 6. Phylogenetics

### Tree data

```rust
use cyanea_datasets::phylogenetics;

// Newick format (12 species)
let tree = phylogenetics::primate_newick();
println!("Newick: {}...{}", &tree[..40], &tree[tree.len()-10..]);

// NEXUS format
let nexus = phylogenetics::primate_nexus();
assert!(nexus.contains("#NEXUS"));

// Species list
let species = phylogenetics::primate_species();
println!("Species: {:?}", species);
```

### Sequence data for tree building

```rust
use cyanea_datasets::phylogenetics;

// Cytochrome b protein sequences (6 primates, aligned)
let seqs = phylogenetics::cytochrome_b_primates();
for (name, seq) in &seqs {
    println!("{}: {} aa", name, seq.len());
}

// FASTA format output
let fasta = phylogenetics::cytochrome_b_fasta();
println!("{}", fasta);
```

### Distance matrix

```rust
use cyanea_datasets::phylogenetics;

let (names, matrix) = phylogenetics::primate_distance_matrix();
println!("Distance matrix ({}x{}):", names.len(), names.len());
for (i, name) in names.iter().enumerate() {
    print!("{:15}", name);
    for j in 0..names.len() {
        print!(" {:.3}", matrix[i][j]);
    }
    println!();
}
```

## 7. Metagenomics

### OTU table analysis

```rust
use cyanea_datasets::metagenomics;

let otu = metagenomics::demo_otu_table();
println!("{} taxa x {} samples", otu.taxa.len(), otu.samples.len());

// Total counts per sample
let totals = otu.sample_totals();
for (sample, total) in otu.samples.iter().zip(totals.iter()) {
    println!("{}: {} reads", sample, total);
}

// Relative abundance
let rel = otu.relative_abundance();
for (i, taxon) in otu.taxa.iter().enumerate() {
    println!("{}: healthy_mean={:.3}, ibd_mean={:.3}",
        taxon,
        rel[i][0..3].iter().sum::<f64>() / 3.0,
        rel[i][3..6].iter().sum::<f64>() / 3.0);
}

// Shannon diversity — healthy > IBD
let div = otu.shannon_diversity();
let healthy_mean: f64 = div[0..3].iter().sum::<f64>() / 3.0;
let ibd_mean: f64 = div[3..6].iter().sum::<f64>() / 3.0;
println!("Shannon diversity: healthy={:.3}, IBD={:.3}", healthy_mean, ibd_mean);
```

### Taxonomy

```rust
use cyanea_datasets::metagenomics;

let taxa = metagenomics::demo_taxonomy();
for t in &taxa {
    println!("{} -> {} -> {} -> {}",
        t.phylum, t.class, t.family, t.species);
}
```

## 8. Structural Biology

### Protein structure

```rust
use cyanea_datasets::structural;

let structure = structural::insulin_chain_b();
println!("{} ({}, chain {}, {} residues)",
    structure.name, structure.pdb_id, structure.chain, structure.len());

// Cα distance between adjacent residues (~3.8 A)
let d = structure.ca_distance(0, 1);
println!("Residues 1-2 distance: {:.1} A", d);
```

### Contact map

```rust
use cyanea_datasets::structural;

let structure = structural::insulin_chain_b();
let contacts = structure.contact_map(8.0);

let mut num_contacts = 0;
for i in 0..structure.len() {
    for j in i+1..structure.len() {
        if contacts[i][j] {
            num_contacts += 1;
        }
    }
}
println!("Contacts (8 A threshold): {}", num_contacts);
```

### PDB output

```rust
use cyanea_datasets::structural;

let structure = structural::insulin_chain_b();
let pdb = structure.to_pdb();
println!("{}", &pdb[..200]); // First few ATOM lines
```

### Ramachandran angles

```rust
use cyanea_datasets::structural;

let angles = structural::demo_ramachandran();
for (label, phi, psi) in &angles {
    println!("{}: phi={:.0}, psi={:.0}", label, phi, psi);
}
// alpha_helix_1: phi=-57, psi=-47
// beta_sheet_1: phi=-120, psi=130
```

## 9. Combining Datasets

Datasets are designed to work together across cyanea crates.

```rust
use cyanea_datasets::{genomics, alignment, phylogenetics};

// Use genomics sequences as input for alignment
let (_, spike_nt) = genomics::sars_cov2_spike_rbd();
let (_, spike_aa) = genomics::sars_cov2_spike_rbd_protein();
println!("Spike RBD: {} nt -> {} aa", spike_nt.len(), spike_aa.len());

// Use alignment sequences with cyanea-align
let seqs = alignment::hemoglobin_alpha();
// Pass to cyanea_align::msa::align_multiple(...)

// Use phylogenetics sequences and distance matrix with cyanea-phylo
let (names, dist) = phylogenetics::primate_distance_matrix();
// Pass to cyanea_phylo::nj::neighbor_joining(...)

// Use VCF data with cyanea-io
let vcf = genomics::demo_vcf_chr22();
// Parse with cyanea_io::vcf::parse_vcf(...)
```
