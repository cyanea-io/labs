# cyanea-io

Unified file format parsing for bioinformatics and tabular data. Each parser is behind a feature flag to keep the dependency tree minimal.

## Status: Complete

CSV, VCF, BED, BEDPE, GFF3, GTF, SAM, BAM, CRAM, BCF, BLAST, MAF, GenBank, bigWig/bigBed, and Parquet parsers are fully implemented. Each format includes both raw record parsing and summary statistics.

## Public API

### CSV (`csv.rs`, `csv` feature)

| Type/Function | Description |
|---------------|-------------|
| `CsvInfo` | `delimiter`, `rows`, `columns`, `sample_preview` |
| `parse_csv_info(path) -> Result<CsvInfo>` | Extract CSV metadata |
| `csv_preview(path, limit) -> Result<String>` | First N rows as a string |

### VCF (`vcf.rs`, `vcf` feature)

| Type/Function | Description |
|---------------|-------------|
| `parse_vcf(path) -> Result<Vec<Variant>>` | Parse VCF data lines into `cyanea_omics::Variant` records |
| `VcfStats` | `variant_count`, `snv_count`, `indel_count`, `pass_count`, `chromosomes` |
| `vcf_stats(path) -> Result<VcfStats>` | Streaming VCF summary statistics |

### BED (`bed.rs`, `bed` feature)

| Type/Function | Description |
|---------------|-------------|
| `BedRecord` | BED3-BED6 fields: `chrom`, `start`, `end`, `name`, `score`, `strand` |
| `parse_bed(path) -> Result<Vec<BedRecord>>` | Parse BED file |
| `parse_bed_intervals(path) -> Result<Vec<GenomicInterval>>` | Parse into `cyanea_omics::GenomicInterval` |
| `BedStats` | `record_count`, `total_bases`, `chromosomes` |
| `bed_stats(path) -> Result<BedStats>` | Summary statistics |

Coordinate convention: BED uses 0-based half-open internally (matching the BED spec).

### GFF3 (`gff.rs`, `gff` feature)

| Type/Function | Description |
|---------------|-------------|
| `parse_gff3(path) -> Result<Vec<Gene>>` | Hierarchical assembly: Gene -> Transcript -> Exon |
| `GffStats` | `gene_count`, `transcript_count`, `exon_count`, `protein_coding_count` |
| `gff3_stats(path) -> Result<GffStats>` | Summary statistics |

Coordinate conversion: GFF3 1-based closed coordinates are converted to 0-based half-open internally. Handles `##FASTA` sentinel, percent-encoded attributes, and `Parent` attribute hierarchy.

### GTF (`gtf.rs`, `gtf` feature)

| Type/Function | Description |
|---------------|-------------|
| `parse_gtf(path) -> Result<Vec<Gene>>` | Hierarchical assembly from GTF (GFF2) format |
| `GtfStats` | `gene_count`, `transcript_count`, `exon_count`, `cds_count` |
| `gtf_stats(path) -> Result<GtfStats>` | Summary statistics |

GTF attribute syntax: `gene_id "ENSG001"; gene_name "TP53";` (space-separated key-value with quoted values). Hierarchy built from `gene_id`/`transcript_id` attributes. 1-based closed → 0-based half-open coordinate conversion.

### SAM (`sam.rs`, `sam` feature)

| Type/Function | Description |
|---------------|-------------|
| `SamRecord` | `qname`, `flag`, `rname`, `pos`, `mapq`, `cigar`, `rnext`, `pnext`, `tlen`, `sequence`, `quality` |
| `parse_sam(path) -> Result<Vec<SamRecord>>` | Parse SAM text format |
| `SamStats` | `total_reads`, `mapped`, `unmapped`, `avg_mapq`, `mapq_distribution` |
| `sam_stats(records) -> SamStats` | Compute stats from records |
| `sam_stats_from_path(path) -> Result<SamStats>` | Streaming SAM statistics |
| `SamPair` | Paired reads: `r1`, `r2` with `insert_size() -> i64` |
| `PairedSamStats` | `base` (SamStats), `paired_count`, `proper_pair_count`, `singletons`, `avg_insert_size` |
| `pair_sam_records(records) -> Vec<SamPair>` | Group records into mate pairs by qname |
| `filter_proper_pairs(records) -> Vec<&SamRecord>` | Filter records with FLAG 0x2 (proper pair) |
| `paired_sam_stats(records) -> PairedSamStats` | Compute paired-end statistics |

**Flag helpers on `SamRecord`**: `is_paired()`, `is_proper_pair()`, `is_mate_unmapped()`, `is_reverse()`, `is_mate_reverse()`, `is_first_in_pair()`, `is_second_in_pair()`, `is_secondary()`, `is_supplementary()`

### BAM (`bam.rs`, `bam` feature)

| Type/Function | Description |
|---------------|-------------|
| `parse_bam(path) -> Result<Vec<SamRecord>>` | Parse BAM (BGZF-compressed) format |
| `bam_stats(path) -> Result<SamStats>` | Streaming BAM statistics |
| `BamReference` | Reference sequence metadata from BAM header |

### CRAM (`cram.rs`, `cram` feature)

| Type/Function | Description |
|---------------|-------------|
| `CramConfig` | Configuration: `reference_path: Option<PathBuf>` for reference FASTA |
| `parse_cram(path, config) -> Result<Vec<SamRecord>>` | Parse CRAM format with reference |
| `parse_cram_default(path) -> Result<Vec<SamRecord>>` | Parse CRAM without external reference |
| `cram_stats(path, config) -> Result<SamStats>` | CRAM summary statistics |
| `cram_stats_default(path) -> Result<SamStats>` | CRAM stats without external reference |

Reuses `SamRecord` and `SamStats` from the SAM module. The `cram` feature implies `sam`.

### BCF (`bcf.rs`, `bcf` feature)

| Type/Function | Description |
|---------------|-------------|
| `parse_bcf(path) -> Result<Vec<Variant>>` | Parse BCF2 binary VCF into `Variant` records |
| `bcf_stats(path) -> Result<VcfStats>` | BCF summary statistics (reuses `VcfStats`) |

BCF2.1 format with BGZF compression. Header contains VCF text for contig/filter definitions. Records use integer-encoded contig IDs and typed values. The `bcf` feature implies `vcf`.

### BLAST (`blast.rs`, `blast` feature)

| Type/Function | Description |
|---------------|-------------|
| `BlastRecord` | 12-column tabular: `query_id`, `subject_id`, `pct_identity`, `alignment_length`, `mismatches`, `gap_opens`, `query_start`, `query_end`, `subject_start`, `subject_end`, `evalue`, `bit_score` |
| `parse_blast(path) -> Result<Vec<BlastRecord>>` | Parse BLAST `-outfmt 6`/`7` tabular output |
| `BlastStats` | `hit_count`, `unique_queries`, `unique_subjects`, `avg_identity`, `avg_evalue` |
| `blast_stats(path) -> Result<BlastStats>` | Summary statistics |

Supports both `-outfmt 6` (no comments) and `-outfmt 7` (with `#` comment lines).

### MAF (`maf.rs`, `maf` feature)

| Type/Function | Description |
|---------------|-------------|
| `MafBlock` | `score: Option<f64>`, `sequences: Vec<MafSequence>` |
| `MafSequence` | `src`, `start`, `size`, `strand`, `src_size`, `text` |
| `parse_maf(path) -> Result<Vec<MafBlock>>` | Parse MAF alignment blocks |
| `MafStats` | `block_count`, `total_aligned_bases`, `species` |
| `maf_stats(path) -> Result<MafStats>` | Summary statistics |

Handles both UCSC Multiple Alignment Format and pairwise MAF (LAST/minimap2 output).

### GenBank (`genbank.rs`, `genbank` feature)

| Type/Function | Description |
|---------------|-------------|
| `GenbankRecord` | `locus`, `definition`, `accession`, `version`, `organism`, `sequence`, `features` |
| `GenbankFeature` | `feature_type`, `location`, `qualifiers: Vec<(String, String)>` |
| `parse_genbank(path) -> Result<Vec<GenbankRecord>>` | Parse multi-record GenBank flat file |
| `GenbankStats` | `record_count`, `total_bases`, `feature_counts: HashMap<String, u64>` |
| `genbank_stats(path) -> Result<GenbankStats>` | Summary statistics |

Parses LOCUS, DEFINITION, ACCESSION, VERSION, ORGANISM, FEATURES table (with qualifiers), and ORIGIN sequence sections. Multi-record files separated by `//`.

### bigWig/bigBed (`bigwig.rs`, `bigwig` feature)

| Type/Function | Description |
|---------------|-------------|
| `BigWigHeader` | `version`, `zoom_levels`, `chrom_count`, `total_summary` |
| `BigWigSummary` | `bases_covered`, `min_val`, `max_val`, `sum`, `sum_squares` |
| `BigWigInterval` | `chrom`, `start`, `end`, `value` |
| `BigBedRecord` | `chrom`, `start`, `end`, `rest` |
| `read_bigwig_header(path) -> Result<BigWigHeader>` | Read bigWig file header + summary |
| `read_bigwig_intervals(path, chrom, start, end) -> Result<Vec<BigWigInterval>>` | Query intervals in a region |
| `read_bigbed_header(path) -> Result<BigWigHeader>` | Read bigBed file header |
| `read_bigbed_records(path, chrom, start, end) -> Result<Vec<BigBedRecord>>` | Query BED records in a region |

Kent binary formats with B+ tree chromosome index and R-tree spatial index. Supports both little-endian and big-endian files. Data blocks may be zlib-compressed.

### Parquet (`parquet.rs`, `parquet` feature)

| Type/Function | Description |
|---------------|-------------|
| `ParquetInfo` | `num_rows`, `num_columns`, `column_names`, `num_row_groups`, `created_by` |
| `parquet_info(path) -> Result<ParquetInfo>` | Extract Parquet file metadata |
| `write_variants_parquet(variants, path) -> Result<()>` | Write `Vec<Variant>` to Parquet |
| `read_variants_parquet(path) -> Result<Vec<Variant>>` | Read variants from Parquet |
| `read_variants_parquet_region(path, chrom, start, end) -> Result<Vec<Variant>>` | Read variants with range filter |
| `parquet_variant_stats(path) -> Result<VcfStats>` | Variant statistics from Parquet |
| `write_intervals_parquet(intervals, path) -> Result<()>` | Write `Vec<GenomicInterval>` to Parquet |
| `read_intervals_parquet(path) -> Result<Vec<GenomicInterval>>` | Read intervals from Parquet |
| `read_intervals_parquet_region(path, chrom, start, end) -> Result<Vec<GenomicInterval>>` | Read intervals with range filter |
| `parquet_interval_stats(path) -> Result<BedStats>` | Interval statistics from Parquet |
| `ExpressionMatrix` | `gene_ids`, `sample_ids`, `values` (genes × samples) |
| `write_expression_parquet(matrix, path) -> Result<()>` | Write expression matrix to Parquet |
| `read_expression_parquet(path) -> Result<ExpressionMatrix>` | Read expression matrix from Parquet |

The `parquet` feature implies `vcf` and `bed` (reuses `VcfStats`, `BedStats`, `Variant`, `GenomicInterval`). Expression matrix uses columnar layout (one column per sample). Region queries use predicate pushdown with row-group statistics filtering.

### BEDPE (`bedpe.rs`, `bed` feature)

| Type/Function | Description |
|---------------|-------------|
| `BedpeRecord` | Paired-end record: `interval1`, `interval2` (as `GenomicInterval`), `name`, `score` |
| `parse_bedpe(path) -> Result<Vec<BedpeRecord>>` | Parse BEDPE file (6 required + 4 optional columns) |
| `BedpeStats` | `record_count`, `total_span`, `inter_chromosomal`, `intra_chromosomal`, `chromosomes` |
| `bedpe_stats(path) -> Result<BedpeStats>` | Summary statistics |

BEDPE format: 6 required columns (chrom1/start1/end1/chrom2/start2/end2), 4 optional (name/score/strand1/strand2). Skips `#`, `track`, and `browser` header lines. Feature-gated behind `bed`.

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `csv` | Yes | CSV parsing (`csv`, `serde`, `serde_json` deps) |
| `vcf` | No | VCF variant parsing (requires `cyanea-omics`) |
| `bed` | No | BED interval parsing (requires `cyanea-omics`) |
| `gff` | No | GFF3 gene structure parsing (requires `cyanea-omics`) |
| `gtf` | No | GTF (GFF2) gene structure parsing (requires `cyanea-omics`) |
| `sam` | No | SAM text alignment parsing |
| `bam` | No | BAM binary alignment parsing (implies `sam`, requires `flate2`) |
| `cram` | No | CRAM reference-based alignment (implies `sam`, requires noodles) |
| `bcf` | No | BCF binary VCF parsing (implies `vcf`, requires `flate2`) |
| `blast` | No | BLAST tabular output parsing |
| `maf` | No | MAF (Multiple Alignment Format) parsing |
| `genbank` | No | GenBank flat file parsing |
| `bigwig` | No | bigWig/bigBed binary format reading (requires `flate2`) |
| `parquet` | No | Apache Parquet columnar format (implies `vcf`, `bed`) |
| `parallel` | No | Rayon parallelism |
| `wasm` | No | WASM target |

## Dependencies

- `cyanea-core` -- error types
- `cyanea-omics` -- genomic types (`Variant`, `GenomicInterval`, `Gene`), optional
- `csv`, `serde`, `serde_json` -- CSV parsing, optional
- `flate2` -- BAM/BCF/bigWig BGZF/zlib decompression, optional
- `noodles-cram`, `noodles-sam`, `noodles-fasta`, `noodles-core` -- CRAM parsing, optional
- `parquet`, `arrow` -- Apache Parquet/Arrow, optional

## Tests

3 tests with default features (CSV only), 158+ tests with all features enabled: CSV (3), VCF (5), BED (9), BEDPE (10), GFF3 (9), GTF (6), SAM (32), pileup (32), BAM (12), CRAM (7), Parquet (14), BLAST (5), MAF (5), GenBank (6), bigWig (6), BCF (5).
