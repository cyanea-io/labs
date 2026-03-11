# cyanea-epi Usage Guide

## Installation

```toml
[dependencies]
cyanea-epi = { version = "0.1", features = ["std"] }
```

## Building a Pileup and Calling Peaks

```rust
use cyanea_epi::pileup::*;
use cyanea_epi::peaks::*;

// Aligned reads as (chrom, start, fragment_size)
let treatment_reads = vec![
    ("chr1".into(), 1000u64, 200u64),
    ("chr1".into(), 1010, 200),
    ("chr1".into(), 1020, 200),
    ("chr1".into(), 1030, 200),
    ("chr1".into(), 1040, 200),
    // ... enriched region
];
let control_reads = vec![
    ("chr1".into(), 500u64, 200u64),
    ("chr1".into(), 1500, 200),
];

// Build pileups
let treatment = build_pileup(&treatment_reads, 200);
let control = build_pileup(&control_reads, 200);

// Call narrow peaks (for TFs, H3K4me3, H3K27ac)
let params = PeakCallParams::default();
let peaks = call_peaks(&treatment, &control, &params).unwrap();
for peak in &peaks {
    println!("{}:{}-{} q={:.2e} fold={:.1}",
        peak.chrom, peak.start, peak.end, peak.q_value, peak.fold_enrichment);
}

// Call broad peaks (for H3K27me3, H3K36me3)
let broad = call_broad_peaks(&treatment, &control, &params).unwrap();
```

## Peak Set Operations

```rust
use cyanea_epi::peaks::{PeakSet, Peak};

let peaks1 = PeakSet::new(vec![/* ... */]);
let peaks2 = PeakSet::new(vec![/* ... */]);

// Merge nearby peaks
let merged = peaks1.merge(100); // merge peaks within 100 bp

// Set operations
let common = peaks1.intersect(&peaks2);
let unique = peaks1.subtract(&peaks2);

// Filter and summarize
let significant = merged.filter_by_score(50.0);
let stats = significant.stats();
println!("{} peaks, median width {} bp", stats.count, stats.median_width);
```

## Signal Normalization and Correlation

```rust
use cyanea_epi::pileup::*;

let pileup = build_pileup(&reads, 200);

// Normalize to CPM (counts per million)
let normalized = normalize_pileup(&pileup, "cpm").unwrap();

// Gaussian smoothing
let smoothed = smooth_pileup(&normalized, 150.0);

// Replicate correlation
let rep1 = build_pileup(&rep1_reads, 200);
let rep2 = build_pileup(&rep2_reads, 200);
if let Some(r) = pileup_correlation(&rep1, &rep2) {
    println!("Replicate Pearson r: {:.3}", r);
}

// Fingerprint plot data (enrichment QC)
let fp = fingerprint(&pileup);
// fp is Vec<(fraction_of_genome, fraction_of_signal)>
```

## Motif Discovery and Scanning

```rust
use cyanea_epi::motifs::*;

// Build a PWM manually (A, C, G, T frequencies per position)
let motif = Motif::new(
    "AP-1".into(),
    vec![
        [0.1, 0.1, 0.1, 0.7],  // T
        [0.1, 0.1, 0.7, 0.1],  // G
        [0.7, 0.1, 0.1, 0.1],  // A
        [0.1, 0.7, 0.1, 0.1],  // C
        [0.1, 0.1, 0.1, 0.7],  // T
        [0.1, 0.7, 0.1, 0.1],  // C
        [0.7, 0.1, 0.1, 0.1],  // A
    ],
);

// Scan a sequence
let matches = scan_sequence(b"ACGTGACTCAGTACGT", &motif, 5.0);
for m in &matches {
    println!("pos={} strand={} score={:.2}", m.position, m.strand, m.score);
}

// Discover motifs from peak sequences
let peak_seqs: Vec<&[u8]> = vec![b"ATGACTCAGT", b"CCTGACTCAG", b"AATGACTCAA"];
let params = DiscoveryParams { motif_width: 7, n_motifs: 5, ..Default::default() };
let discovered = discover_motifs(&peak_seqs, &params).unwrap();

// Parse and write MEME format
let meme_str = write_meme(&discovered);
let parsed = parse_meme(&meme_str).unwrap();

// Motif enrichment test
let (odds_ratio, p_value) = motif_enrichment(&motif, &peak_seqs, &background_seqs).unwrap();
```

## Chromatin State Learning

```rust
use cyanea_epi::chromatin::*;

// Binary mark presence matrix: rows = genomic bins, cols = histone marks
// marks: [H3K4me3, H3K27ac, H3K27me3, H3K36me3]
let marks = vec![
    vec![true, true, false, false],   // active promoter
    vec![true, true, false, false],
    vec![false, false, true, false],  // repressed
    vec![false, false, true, false],
    vec![false, false, false, true],  // transcribed
];

let params = ChromHMMParams { n_states: 3, ..Default::default() };
let model = learn_chromatin_states(&marks, &params).unwrap();

// Segment the genome
let segmentation = segment_genome(&model, &marks).unwrap();
for region in &segmentation.regions {
    println!("bins {}-{}: state {}", region.start, region.end, region.state_id);
}

// Enrichment at annotations (e.g., TSS positions as bin indices)
let annotations = vec![vec![0usize, 1], vec![4]]; // two annotation types
let enrichments = state_enrichment(&segmentation, &annotations);
```

## Differential Binding Analysis

```rust
use cyanea_epi::differential::*;

// Read counts at peaks: rows = peaks, cols = samples
// Conditions: ["treated", "treated", "control", "control"]
let counts = vec![
    vec![100.0, 120.0, 20.0, 25.0],  // peak 1: enriched in treatment
    vec![50.0, 45.0, 55.0, 48.0],    // peak 2: unchanged
    vec![10.0, 15.0, 80.0, 90.0],    // peak 3: enriched in control
];
let conditions = vec!["treated", "treated", "control", "control"];

let results = differential_peaks(&counts, &conditions, "wald").unwrap();
for r in &results {
    println!("peak {}: log2FC={:.2}, q={:.4}", r.region, r.log2_fc, r.q_value);
}

// MA plot data
let ma = ma_plot_data(&results);
// ma is Vec<(mean_expression, log2_fold_change)>
```

## Nucleosome Analysis

```rust
use cyanea_epi::nucleosome::*;
use cyanea_epi::pileup::build_pileup;

// MNase-seq pileup
let pileup = build_pileup(&mnase_reads, 147);

let params = NucleosomeParams::default();
let positions = call_nucleosomes(&pileup, &params);
for nuc in &positions {
    println!("center={}, occupancy={:.1}, fuzziness={:.1}",
        nuc.center, nuc.occupancy, nuc.fuzziness);
}

// NFR score around TSSs
let tss_list = vec![10000, 20000, 30000];
let score = nfr_score(&positions, &tss_list);
println!("NFR score: {:.3}", score);
```

## ATAC-seq Quality Control

```rust
use cyanea_epi::accessibility::*;
use cyanea_epi::pileup::build_pileup;

// Fragment size analysis
let fragment_sizes: Vec<u64> = vec![120, 135, 150, 180, 200, 300, 350];
let dist = fragment_size_distribution(&fragment_sizes, 50);
let metrics = insert_size_metrics(&fragment_sizes);
println!("NFR ratio: {:.2}, mono-nuc ratio: {:.2}", metrics.nfr_ratio, metrics.mono_nuc_ratio);

// TSS enrichment score
let pileup = build_pileup(&atac_reads, 200);
let tss_positions = vec![10000, 20000, 30000];
let tss_score = tss_enrichment(&pileup, &tss_positions, 2000);
println!("TSS enrichment: {:.1}", tss_score);  // ENCODE: >7 is good

// FRiP
let frip_score = frip(500000, 2000000);
println!("FRiP: {:.1}%", frip_score * 100.0);

// Comprehensive QC
let qc = atacqc(&pileup, &peaks, &tss_positions, &fragment_sizes).unwrap();
println!("TSS: {:.1}, FRiP: {:.1}%, NFR: {:.2}",
    qc.tss_enrichment, qc.frip * 100.0, qc.insert_metrics.nfr_ratio);
```
