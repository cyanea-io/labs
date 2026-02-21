# cyanea-io Architecture

## Module Map

```
cyanea-io
 +-- csv.rs              [csv] CSV metadata, preview
 +-- vcf.rs              [vcf] VCF variant parsing, stats
 +-- vcf_header.rs       [vcf] Structured VCF 4.3 header construction/parsing
 +-- vcf_ops.rs          [vcf] Normalization, multi-allelic split/join, filtering, set ops
 +-- bed.rs              [bed] BED3-BED6 parsing, interval conversion
 +-- bedpe.rs            [bed] BEDPE paired-end intervals
 +-- bedgraph.rs         [bed] bedGraph/Wiggle signal tracks
 +-- gff.rs              [gff] GFF3 hierarchical gene parsing
 +-- gtf.rs              [gtf] GTF (GFF2) gene parsing
 +-- sam.rs              [sam] SAM record parsing, flag helpers, paired-end stats
 +-- pileup.rs           [sam] Pileup generation from SAM, mpileup output
 +-- bgzf.rs             [bam] BGZF block decompression, virtual offsets
 +-- bam.rs              [bam] BAM binary alignment parsing
 +-- bam_ops.rs          [bam] Sort, merge, mark duplicates, flagstat, depth
 +-- indexed_bam.rs      [bam] Random-access BAM via BAI/CSI index (noodles)
 +-- indexed_vcf.rs      [vcf] Random-access VCF via tabix (noodles)
 +-- cram.rs             [cram] CRAM format via noodles
 +-- bcf.rs              [bcf] BCF2 binary VCF reader
 +-- bcf_write.rs        [bcf] BCF2 binary VCF writer
 +-- variant_call.rs     [variant-calling] Bayesian genotype caller
 +-- blast.rs            [blast] BLAST tabular output (-outfmt 6/7)
 +-- blast_xml.rs        [blast] BLAST XML output (-outfmt 5)
 +-- maf.rs              [maf] Multiple Alignment Format
 +-- genbank.rs          [genbank] GenBank flat file
 +-- embl.rs             [genbank] EMBL/ENA format
 +-- stockholm.rs        [genbank] Stockholm MSA format
 +-- clustal.rs          [genbank] ClustalW/Omega format
 +-- phylip.rs           [genbank] PHYLIP interleaved/sequential
 +-- pir.rs              [genbank] PIR/NBRF protein format
 +-- abi.rs              [genbank] ABI Sanger chromatogram binary
 +-- gfa.rs              [genbank] GFA v1 sequence graph
 +-- bigwig.rs           [bigwig] bigWig/bigBed Kent binary format
 +-- parquet.rs          [parquet] Apache Parquet via arrow/parquet crates
 +-- fetch.rs            [fetch] URL builders for NCBI/UniProt/KEGG/htsget/refget
```

## Design Decisions

### Feature Flag Architecture

Every I/O module is behind a Cargo feature flag. This is essential because:

1. **Dependency minimization**: Each format parser may pull in large dependencies (noodles for BAM/CRAM, arrow/parquet for Parquet). Users only pay for formats they use.
2. **Compile time**: A full build with all features is significantly slower than the default (CSV-only) build.
3. **WASM compatibility**: Some formats depend on filesystem or system libraries that are unavailable in WASM.

Feature implications enforce logical groupings:
- `bam` implies `sam` (BAM records are returned as `SamRecord`)
- `bcf` implies `vcf` (BCF records are returned as `Variant`)
- `cram` implies `sam`
- `parquet` implies `vcf` + `bed` (uses `Variant` and `GenomicInterval` types)
- `variant-calling` implies `sam` + `vcf` (reads pileup data, outputs variants)

### BGZF: Block Gzip with Virtual Offsets

BGZF (Blocked GNU Zip Format) is the compression layer for BAM, BCF, and tabix-indexed VCF. Key properties:

- Series of concatenated gzip blocks, each at most 64 KiB uncompressed.
- **Virtual offset** = `(block_offset << 16) | within_block_offset`, enabling random access by seeking to the compressed block and then to the byte within the decompressed block.
- BAI/TBI/CSI indices store virtual offsets for genomic intervals.

The `bgzf.rs` module provides shared decompression used by both `bam.rs` and `bcf.rs`.

### BAI/TBI/CSI Index Reading for Random Access

Indexed readers (`indexed_bam.rs`, `indexed_vcf.rs`) use noodles for index parsing:

- **BAI** (BAM index): bins of virtual offset ranges, linear index of minimum virtual offsets per 16kb window.
- **TBI** (tabix): generalized index for tab-delimited, coordinate-sorted, bgzipped files (VCF, BED, GFF).
- **CSI** (Coordinate-Sorted Index): supports arbitrary bin sizes for very large genomes.

The workflow: parse the index file, look up bins overlapping the query region, seek to each bin's virtual offset in the BGZF file, decompress blocks, and filter records by coordinate.

### Pileup Generation from SAM Records

`pileup.rs` walks SAM records in coordinate order, expanding CIGAR strings to track which reads cover each reference position. For each position, it maintains:

- Base counts (A, C, G, T, N, deletions, insertions)
- Quality score sums per base
- Depth of coverage
- Forward/reverse strand counts

The pileup can be formatted as samtools mpileup output or fed to the variant caller.

### Variant Calling: Bayesian Genotype Likelihoods

`variant_call.rs` implements a diploid Bayesian genotype caller:

1. For each pileup position with sufficient depth, enumerate candidate alleles.
2. Compute genotype likelihoods P(data | genotype) using base quality error probabilities.
3. Apply diploid genotype priors (homozygous ref, heterozygous, homozygous alt).
4. Compute posterior probabilities and report the MAP genotype.
5. Apply quality filters: minimum depth, minimum allele count, strand bias (Fisher exact test).

### Parquet: Arrow Record Batches with Predicate Pushdown

The Parquet integration uses the Apache Arrow ecosystem:

- Variants and intervals are converted to Arrow record batches with typed columns.
- Region queries use row-group statistics (min/max chromosome, min/max position) for predicate pushdown, skipping entire row groups that cannot contain matches.
- Expression matrices use a columnar layout with one column per sample, enabling efficient per-sample access.

### noodles Integration

noodles provides Rust-native BAM, CRAM, BCF, and tabix parsers. cyanea-io uses noodles for:

- **BAM indexed reading**: `noodles-bam::io::IndexedReader` for BAI/CSI-backed random access.
- **CRAM**: `noodles-cram` for reference-based compression format.
- **BCF**: Header parsing and typed-value decoding (though the core BCF reader is custom).
- **Indexed VCF**: `noodles-vcf::io::IndexedReader` for tabix-backed random access.

All noodles records are converted to cyanea-io's own types (`SamRecord`, `Variant`) at the boundary, keeping the public API independent of noodles version.
