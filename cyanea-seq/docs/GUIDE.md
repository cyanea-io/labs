# cyanea-seq Usage Guide

## Installation

```toml
[dependencies]
cyanea-seq = { version = "0.1", features = ["std", "minhash"] }
```

## DNA/RNA/Protein Sequence Creation

```rust
use cyanea_seq::{DnaSequence, RnaSequence, ProteinSequence};

// Create validated sequences (returns error on invalid characters)
let dna = DnaSequence::new(b"ACGTACGT").unwrap();
let rna = RnaSequence::new(b"ACGUACGU").unwrap();
let protein = ProteinSequence::new(b"MKWVTFISLLFLFSSAYS").unwrap();

// Basic properties
println!("DNA length: {}", dna.len());
println!("GC content: {:.1}%", dna.gc_content() * 100.0);

// Reverse complement
let rc = dna.reverse_complement();
assert_eq!(rc.as_bytes(), b"ACGTACGT");

// Transcription and translation
let rna_from_dna = dna.transcribe();
let protein_from_dna = dna.translate().unwrap();

// RNA reverse transcription
let dna_from_rna = rna.reverse_transcribe();

// Protein molecular weight
let mw = protein.molecular_weight();
println!("Molecular weight: {:.1} Da", mw);
```

## FASTA/FASTQ Parsing

```rust
use cyanea_seq::fasta::parse_fasta_stats;
use cyanea_seq::fastq::{parse_fastq_file, parse_fastq_stats, FastqRecord};

// Streaming FASTA statistics (memory efficient)
let stats = parse_fasta_stats("genome.fa").unwrap();
println!("{} sequences, {} total bases", stats.sequence_count, stats.total_bases);
println!("GC: {:.1}%, avg length: {:.0}", stats.gc_content * 100.0, stats.avg_length);

// Parse all FASTQ records into memory
let records = parse_fastq_file("reads.fq").unwrap();
for rec in &records {
    println!("{}: {} bp", rec.id, rec.seq.len());
}

// Streaming FASTQ statistics (for large files)
let fq_stats = parse_fastq_stats("reads.fq").unwrap();
println!("Mean Q: {:.1}, Q20: {:.1}%, Q30: {:.1}%",
    fq_stats.mean_quality,
    fq_stats.q20_fraction * 100.0,
    fq_stats.q30_fraction * 100.0);
```

## K-mer Iteration and 2-bit Encoding

```rust
use cyanea_seq::DnaSequence;
use cyanea_seq::twobit::TwoBitSequence;

// K-mer iteration
let seq = DnaSequence::new(b"ACGTACGT").unwrap();
let kmers: Vec<&[u8]> = seq.kmers(4).unwrap().collect();
// ["ACGT", "CGTA", "GTAC", "TACG", "ACGT"]

// 2-bit encoding (4x compression, ACGT only)
let encoded = TwoBitSequence::encode(b"ACGTACGT").unwrap();
assert_eq!(encoded.len(), 8);
assert_eq!(encoded.get(0), Some(b'A'));

// Extract k-mer as integer encoding (for hash tables)
let kmer_int = encoded.kmer(0, 4).unwrap(); // ACGT as u64

// Complement
let comp = encoded.complement(); // TGCATGCA
```

## FM-Index Construction and Search

```rust
use cyanea_seq::fm_index::FmIndex;

// Build FM-Index from reference sequence
let reference = b"ACGTACGTACGTACGT";
let fm = FmIndex::build(reference);

// Exact pattern search
let positions = fm.search(b"ACGT");
println!("Found {} occurrences: {:?}", positions.len(), positions);

// Count-only (more efficient when positions not needed)
let count = fm.count(b"ACGT");
println!("{} occurrences of ACGT", count);
```

## Pattern Matching

```rust
use cyanea_seq::pattern::{horspool, kmp, shift_and, myers_bitparallel, ukkonen};

let text = b"ACGTACGTACGTACGT";
let pattern = b"ACGT";

// Exact matching (choose based on pattern length and use case)
let hits_horspool = horspool(text, pattern);   // O(n/m) avg, best for long patterns
let hits_kmp = kmp(text, pattern);             // O(n+m), stable worst case
let hits_bitpar = shift_and(text, pattern);    // O(n), fast for patterns <= 64

// Approximate matching (find pattern with up to k errors)
let approx = myers_bitparallel(text, b"ACGX", 1);
// Returns Vec<(end_position, edit_distance)>
for (end_pos, dist) in &approx {
    println!("Match ending at {} with {} edits", end_pos, dist);
}

// Ukkonen for longer patterns
let approx2 = ukkonen(text, b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGX", 2);
```

## Quality Trimming Pipeline

```rust
use cyanea_seq::trim::TrimPipeline;
use cyanea_seq::fastq::parse_fastq_file;

let records = parse_fastq_file("reads.fq").unwrap();

// Build a trimming pipeline (Trimmomatic-style)
let pipeline = TrimPipeline::new()
    .illumina_adapters()           // Remove common Illumina adapters
    .leading(3)                     // Trim leading bases below Q3
    .trailing(3)                    // Trim trailing bases below Q3
    .sliding_window(4, 15.0)       // Cut when 4-base window mean < Q15
    .min_length(36)                 // Discard reads shorter than 36 bp
    .min_mean_quality(20.0);        // Discard reads with mean Q < 20

// Process batch with statistics
let report = pipeline.process_batch_with_stats(&records);
println!("Kept {}/{} reads ({:.1}%)",
    report.kept, records.len(),
    report.kept as f64 / records.len() as f64 * 100.0);
```

## MinHash Sketching for Similarity

```rust
use cyanea_seq::minhash::{MinHash, FracMinHash};

// Build bottom-k MinHash sketches
let sketch_a = MinHash::from_sequence(b"ACGTACGTACGT", 4, 100).unwrap();
let sketch_b = MinHash::from_sequence(b"ACGTACGTTTTT", 4, 100).unwrap();

// Estimate similarity
let jaccard = sketch_a.jaccard(&sketch_b).unwrap();
let containment = sketch_a.containment(&sketch_b).unwrap();
let ani = sketch_a.ani(&sketch_b).unwrap();
println!("Jaccard: {:.3}, Containment: {:.3}, ANI: {:.4}", jaccard, containment, ani);

// FracMinHash for variable-size genomes
let frac = FracMinHash::from_sequence(b"ACGTACGTACGT", 4, 1000).unwrap();
```

## RNA Secondary Structure Prediction

```rust
use cyanea_seq::rna_structure::{nussinov, zuker_mfe, mccaskill, RnaSecondaryStructure};

let seq = b"GGGAAACCC";

// Nussinov: maximize base pair count (simple model)
let nuss = nussinov(seq, 3).unwrap();
println!("{} base pairs", nuss.max_pairs);
println!("{}", nuss.structure.to_dot_bracket());

// Zuker MFE: minimum free energy with Turner parameters
let mfe = zuker_mfe(seq).unwrap();
println!("MFE: {:.1} kcal/mol", mfe.energy);
println!("{}", mfe.structure.to_dot_bracket());

// McCaskill: base pair probabilities
let part = mccaskill(seq, 310.15).unwrap();  // 37 C
println!("P(0,8 paired): {:.3}", part.pair_probability(0, 8));
```

## ORF Finding

```rust
use cyanea_seq::orf::{find_orfs_both_strands, OrfResult};

let seq = b"ATGATGATGATGATGATGATGATGATGAAATTTCCCGGGTAA";

// Find ORFs in all 6 reading frames (min 30 bp)
let orfs = find_orfs_both_strands(seq, 30);
for orf in &orfs {
    println!("ORF: {}..{} frame={} strand={:?} len={}",
        orf.start, orf.end, orf.frame, orf.strand, orf.sequence.len());
}
```

## Codon Usage Analysis

```rust
use cyanea_seq::codon::{translate_codon, translate_sequence};

// Single codon translation
let aa = translate_codon(b"ATG"); // Some(b'M')

// Full coding sequence translation
let protein = translate_sequence(b"ATGAAATTTCCCGGGTAA");
// Returns amino acid bytes, stops at first stop codon
```

## Long-Read Analysis

```rust
use cyanea_seq::longread::{
    LongRead, LongReadPlatform, longread_stats,
    self_correct, simple_consensus, simulate_long_reads,
    common_adapters, trim_adapters, LongReadSimConfig,
};

// Create a long read
let read = LongRead::new(
    "read_001",
    b"ACGTACGTACGTACGT".to_vec(),
    vec![30u8; 16],
    LongReadPlatform::PacBioHiFi,
).unwrap();

println!("Length: {} bp", read.len());
println!("Mean Q: {:.1}", read.mean_quality());
println!("Error rate: {:.4}", read.error_rate());
println!("GC: {:.1}%", read.gc_content() * 100.0);

// Compute stats for a collection of reads
let reads = vec![read];
let stats = longread_stats(&reads).unwrap();
println!("N50: {} bp, {} total reads", stats.n50, stats.num_reads);

// Self-correction using k-mer consensus
let corrected = self_correct(&reads[0], 11).unwrap();
println!("{} bases corrected", corrected.corrections);

// Majority-vote consensus from overlapping reads
let overlapping: Vec<&[u8]> = vec![b"ACGTACGT", b"ACGTACGT", b"ACGTTCGT"];
let consensus = simple_consensus(&overlapping).unwrap();

// Simulate PacBio HiFi reads from a reference
let config = LongReadSimConfig {
    platform: LongReadPlatform::PacBioHiFi,
    mean_length: 15000,
    length_std: 3000,
    coverage: 30.0,
    num_passes: 10,
    seed: 42,
};
let sim_reads = simulate_long_reads(b"ACGT...", &config).unwrap();

// Trim platform-specific adapters
let adapters = common_adapters();
let trimmed = trim_adapters(&reads[0], &adapters, 1);
```

## Structural Variant Detection

```rust
use cyanea_seq::sv::{
    call_svs, svs_from_cigar, sv_summary,
    SplitAlignment, SvCallConfig, SvType,
};

// Detect SVs from split alignments
let alignments = vec![
    SplitAlignment {
        read_name: "read_1".into(),
        chrom: "chr1".into(),
        ref_start: 1000, ref_end: 2000,
        query_start: 0, query_end: 1000,
        is_reverse: false, mapq: 60,
    },
    SplitAlignment {
        read_name: "read_1".into(),
        chrom: "chr1".into(),
        ref_start: 5000, ref_end: 6000,
        query_start: 1000, query_end: 2000,
        is_reverse: false, mapq: 60,
    },
];

let config = SvCallConfig { min_support: 1, ..Default::default() };
let svs = call_svs(&alignments, &config).unwrap();
for sv in &svs {
    println!("{} {}:{}-{} len={}", sv.sv_type.as_str(), sv.chrom, sv.start, sv.end, sv.length);
}

// Extract large indels from CIGAR strings
let cigar = vec![('M', 1000), ('D', 500), ('M', 500), ('I', 200), ('M', 1000)];
let cigar_svs = svs_from_cigar("chr1", 10000, &cigar, 50);

// Summarize SV calls
let summary = sv_summary(&svs);
println!("{} total SVs: {} DEL, {} INS", summary.total, summary.deletions, summary.insertions);
```

## Nanopore Analysis

```rust
use cyanea_seq::nanopore::{
    parse_signal_metadata, parse_methylation_calls,
    aggregate_methylation, nanopore_qc,
};

// Parse signal metadata from pod5/slow5 output
let line = "read_id=abc123 run_id=run1 flow_cell_id=FC001 channel=42 start_time=100.5 duration=5.2 sampling_rate=4000 num_samples=20800 median_signal=85.3";
let meta = parse_signal_metadata(line).unwrap();
println!("Read {} on channel {}, {:.1}s duration", meta.read_id, meta.channel, meta.duration);

// Parse methylation calls (modBAM-style tab-separated)
let calls_text = "read1\tchr1\t1000\t+\t5mC\t0.95\nread2\tchr1\t1000\t+\t5mC\t0.85\n";
let calls = parse_methylation_calls(calls_text).unwrap();

// Aggregate per-read calls into per-site summaries
let sites = aggregate_methylation(&calls, 0.5); // threshold = 0.5
for site in &sites {
    println!("{}:{} {} coverage={} freq={:.2}",
        site.chrom, site.position, site.mod_type.as_str(),
        site.coverage, site.frequency);
}

// Nanopore run QC
let metadata = vec![meta];
let lengths = vec![5000usize];
let qualities = vec![15.0f64];
let qc = nanopore_qc(&metadata, &lengths, &qualities, 512, 10.0).unwrap();
println!("N50: {}, pore occupancy: {:.1}%, speed: {:.0} bp/s",
    qc.n50, qc.pore_occupancy * 100.0, qc.translocation_speed);
```
