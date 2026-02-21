# cyanea-io Usage Guide

## Installation

```toml
[dependencies]
cyanea-io = { version = "0.1", features = ["vcf", "bed", "sam", "bam", "genbank"] }
```

Enable only the format-specific features you need. Each feature pulls in only its required dependencies.

## CSV Parsing

```rust
use cyanea_io::csv::{parse_csv_info, csv_preview};

// Extract metadata
let info = parse_csv_info("data.csv").unwrap();
println!("{} rows x {} columns, delimiter: {:?}",
    info.rows, info.columns, info.delimiter);

// Preview first 5 rows
let preview = csv_preview("data.csv", 5).unwrap();
println!("{}", preview);
```

## VCF Reading and Writing

```rust
use cyanea_io::vcf::{parse_vcf, vcf_stats};

// Parse all variants
let variants = parse_vcf("variants.vcf").unwrap();
for v in &variants {
    println!("{}:{} {} -> {:?} (QUAL: {:?})",
        v.chrom, v.position, v.ref_allele, v.alt_alleles, v.quality);
}

// Streaming statistics (memory-efficient for large files)
let stats = vcf_stats("variants.vcf").unwrap();
println!("{} variants: {} SNVs, {} indels, {} PASS",
    stats.variant_count, stats.snv_count, stats.indel_count, stats.pass_count);
```

## BED File Operations

```rust
use cyanea_io::bed::{parse_bed, parse_bed_intervals, bed_stats};

// Parse BED records
let records = parse_bed("regions.bed").unwrap();
for rec in &records {
    println!("{}:{}-{} {:?}", rec.chrom, rec.start, rec.end, rec.name);
}

// Parse directly into GenomicInterval for interval arithmetic
let intervals = parse_bed_intervals("regions.bed").unwrap();

// Statistics
let stats = bed_stats("regions.bed").unwrap();
println!("{} regions covering {} bases across {} chromosomes",
    stats.record_count, stats.total_bases, stats.chromosomes.len());
```

## GFF3/GTF Annotation Parsing

```rust
use cyanea_io::gff::{parse_gff3, gff3_stats};
use cyanea_io::gtf::{parse_gtf, gtf_stats};

// GFF3: hierarchical Gene -> Transcript -> Exon
let genes = parse_gff3("annotations.gff3").unwrap();
for gene in &genes {
    println!("Gene: {} ({} transcripts)", gene.name, gene.transcripts.len());
}

// GFF3 statistics
let stats = gff3_stats("annotations.gff3").unwrap();
println!("{} genes, {} transcripts, {} exons",
    stats.gene_count, stats.transcript_count, stats.exon_count);

// GTF parsing (same Gene output type)
let gtf_genes = parse_gtf("annotations.gtf").unwrap();
```

## SAM/BAM Parsing and Statistics

```rust
use cyanea_io::sam::{parse_sam, sam_stats, pair_sam_records, paired_sam_stats};
use cyanea_io::bam::{parse_bam, bam_stats};

// Parse SAM records
let records = parse_sam("reads.sam").unwrap();
println!("{} records", records.len());

// Compute alignment statistics
let stats = sam_stats(&records);
println!("{} mapped / {} total, mean MAPQ: {:.1}",
    stats.mapped, stats.total_reads, stats.avg_mapq);

// Paired-end analysis
let pairs = pair_sam_records(&records);
let paired_stats = paired_sam_stats(&records);
println!("{} pairs, {:.1}% proper, mean insert: {:.0}",
    paired_stats.paired_count,
    paired_stats.proper_pair_count as f64 / paired_stats.paired_count as f64 * 100.0,
    paired_stats.avg_insert_size);

// BAM parsing (same SamRecord output)
let bam_records = parse_bam("reads.bam").unwrap();
let bam_stats = bam_stats("reads.bam").unwrap();
```

## BAM Operations (Sort, Index, Merge, Mark Duplicates)

```rust
use cyanea_io::bam::parse_bam;
use cyanea_io::bam_ops::{
    coordinate_sort, queryname_sort, merge_sorted,
    mark_duplicates, flagstat, depth_of_coverage,
};

let mut records = parse_bam("reads.bam").unwrap();

// Sort by coordinate
let ref_order = vec!["chr1".to_string(), "chr2".to_string()];
coordinate_sort(&mut records, &ref_order);

// Mark duplicates
let dup_report = mark_duplicates(&mut records);
println!("{} duplicates out of {} reads ({:.1}% rate)",
    dup_report.duplicates, dup_report.total, dup_report.rate * 100.0);

// flagstat-equivalent
let flags = flagstat(&records);

// Per-position depth of coverage
let depths = depth_of_coverage(&records, "chr1");
```

## Indexed BAM/VCF Random Access

```rust
use cyanea_io::indexed_bam::IndexedBamReader;
use cyanea_io::indexed_vcf::IndexedVcfReader;

// Open indexed BAM (auto-detects .bai)
let mut bam = IndexedBamReader::open("reads.bam").unwrap();
let records = bam.query("chr1:1000-2000").unwrap();
println!("{} reads in region", records.len());

// Open indexed VCF (auto-detects .tbi)
let mut vcf = IndexedVcfReader::open("variants.vcf.gz").unwrap();
let variants = vcf.query("chr1:1000-2000").unwrap();
println!("{} variants in region", variants.len());
```

## Parquet Columnar Analytics

```rust
use cyanea_io::parquet::{
    parquet_info, write_variants_parquet, read_variants_parquet,
    read_variants_parquet_region, write_expression_parquet, read_expression_parquet,
};

// Write variants to Parquet (columnar, compressed)
let variants = vec![/* ... */];
write_variants_parquet(&variants, "variants.parquet").unwrap();

// Read back (optionally with region filter using predicate pushdown)
let all = read_variants_parquet("variants.parquet").unwrap();
let region = read_variants_parquet_region("variants.parquet", "chr1", 1000, 2000).unwrap();

// File metadata
let info = parquet_info("data.parquet").unwrap();
println!("{} rows, {} columns, {} row groups",
    info.num_rows, info.num_columns, info.num_row_groups);
```

## GenBank/EMBL Sequence Formats

```rust
use cyanea_io::genbank::{parse_genbank, genbank_stats};
use cyanea_io::embl::parse_embl;

// GenBank flat file
let records = parse_genbank("sequence.gb").unwrap();
for rec in &records {
    println!("{}: {} ({} bp, {} features)",
        rec.accession, rec.definition, rec.sequence.len(), rec.features.len());
}

// EMBL format
let embl_input = std::fs::read_to_string("sequence.embl").unwrap();
let embl_records = parse_embl(&embl_input).unwrap();
```

## Alignment Format I/O (Stockholm, Clustal, Phylip)

```rust
use cyanea_io::stockholm::parse_stockholm;
use cyanea_io::clustal::parse_clustal;
use cyanea_io::phylip::parse_phylip;

// Stockholm (Pfam/Rfam/HMMER)
let stk_input = std::fs::read_to_string("alignment.sto").unwrap();
let alignments = parse_stockholm(&stk_input).unwrap();
for aln in &alignments {
    println!("{} sequences", aln.sequences.len());
    // Access GC annotations (e.g., consensus structure)
    if let Some(ss) = aln.gc_annotations.get("SS_cons") {
        println!("Secondary structure: {}", ss);
    }
}

// Clustal
let clu_input = std::fs::read_to_string("alignment.aln").unwrap();
let clustal = parse_clustal(&clu_input).unwrap();
println!("{} sequences aligned", clustal.sequences.len());

// PHYLIP (interleaved)
let phy_input = std::fs::read_to_string("alignment.phy").unwrap();
let phylip = parse_phylip(&phy_input).unwrap();
println!("{} taxa, {} sites", phylip.n_taxa, phylip.n_sites);
```
