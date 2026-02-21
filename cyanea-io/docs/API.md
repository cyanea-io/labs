# cyanea-io API Reference

Unified file format parsing for bioinformatics and tabular data. Each parser is behind a feature flag to keep the dependency tree minimal.

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
| `parse_vcf(path) -> Result<Vec<Variant>>` | Parse VCF data lines into `Variant` records |
| `VcfStats` | `variant_count`, `snv_count`, `indel_count`, `pass_count`, `chromosomes` |
| `vcf_stats(path) -> Result<VcfStats>` | Streaming VCF summary statistics |

### VCF header (`vcf_header.rs`, `vcf` feature)

| Type | Description |
|------|-------------|
| `VcfHeader` | Structured VCF 4.3 header with contigs, INFO, FORMAT, FILTER definitions |
| `ContigLine` | Contig definition: `id`, `length` |
| `FieldDef` | INFO/FORMAT field: `id`, `number`, `field_type`, `description` |
| `FilterDef` | FILTER definition: `id`, `description` |

| Method | Description |
|--------|-------------|
| `VcfHeader::new() -> Self` | Create empty header |
| `VcfHeader::add_contig(id, length) -> &mut Self` | Add contig definition |
| `VcfHeader::add_info(field) -> &mut Self` | Add INFO field |
| `VcfHeader::add_format(field) -> &mut Self` | Add FORMAT field |
| `VcfHeader::add_filter(filter) -> &mut Self` | Add FILTER definition |
| `VcfHeader::to_string() -> String` | Serialize to VCF header text |
| `VcfHeader::parse(text) -> Result<Self>` | Parse VCF header text |

### VCF operations (`vcf_ops.rs`, `vcf` feature)

| Function | Description |
|----------|-------------|
| `split_multiallelic(variant) -> Vec<Variant>` | Split multi-allelic into biallelic records |
| `join_biallelic(variants) -> Option<Variant>` | Join biallelic variants at same position |
| `left_align(variant, ref_seq) -> Variant` | Left-align indels against reference |
| `normalize_variant(variant, ref_seq) -> Variant` | Full normalization (left-align + trim) |
| `filter_variants(variants, expression) -> Result<Vec<Variant>>` | Filter with expression parsing |
| `intersect_variants(a, b) -> Vec<Variant>` | Set intersection (matching chrom/pos/ref/alt) |
| `subtract_variants(a, b) -> Vec<Variant>` | Set subtraction |
| `concordance(a, b) -> VariantConcordance` | Concordance statistics (TP, FP, FN) |
| `detailed_vcf_stats(variants) -> DetailedVcfStats` | Ti/Tv, het/hom ratio, per-chrom counts |

### BED (`bed.rs`, `bed` feature)

| Type/Function | Description |
|---------------|-------------|
| `BedRecord` | BED3-BED6 fields: `chrom`, `start`, `end`, `name`, `score`, `strand` |
| `parse_bed(path) -> Result<Vec<BedRecord>>` | Parse BED file |
| `parse_bed_intervals(path) -> Result<Vec<GenomicInterval>>` | Parse into `GenomicInterval` |
| `BedStats` | `record_count`, `total_bases`, `chromosomes` |
| `bed_stats(path) -> Result<BedStats>` | Summary statistics |

### BEDPE (`bedpe.rs`, `bed` feature)

| Type/Function | Description |
|---------------|-------------|
| `BedpeRecord` | Paired-end record: `interval1`, `interval2`, `name`, `score` |
| `parse_bedpe(path) -> Result<Vec<BedpeRecord>>` | Parse BEDPE file |
| `BedpeStats` | `record_count`, `total_span`, `inter_chromosomal`, `intra_chromosomal` |
| `bedpe_stats(path) -> Result<BedpeStats>` | Summary statistics |

### GFF3 (`gff.rs`, `gff` feature)

| Type/Function | Description |
|---------------|-------------|
| `parse_gff3(path) -> Result<Vec<Gene>>` | Hierarchical assembly: Gene -> Transcript -> Exon |
| `GffStats` | `gene_count`, `transcript_count`, `exon_count`, `protein_coding_count` |
| `gff3_stats(path) -> Result<GffStats>` | Summary statistics |

### GTF (`gtf.rs`, `gtf` feature)

| Type/Function | Description |
|---------------|-------------|
| `parse_gtf(path) -> Result<Vec<Gene>>` | Hierarchical assembly from GTF (GFF2) format |
| `GtfStats` | `gene_count`, `transcript_count`, `exon_count`, `cds_count` |
| `gtf_stats(path) -> Result<GtfStats>` | Summary statistics |

### SAM (`sam.rs`, `sam` feature)

| Type/Function | Description |
|---------------|-------------|
| `SamRecord` | Full SAM record: `qname`, `flag`, `rname`, `pos`, `mapq`, `cigar`, etc. |
| `parse_sam(path) -> Result<Vec<SamRecord>>` | Parse SAM text format |
| `SamStats` | `total_reads`, `mapped`, `unmapped`, `avg_mapq`, `mapq_distribution` |
| `sam_stats(records) -> SamStats` | Compute stats from records |
| `sam_stats_from_path(path) -> Result<SamStats>` | Streaming SAM statistics |
| `SamPair` | Paired reads: `r1`, `r2` with `insert_size()` |
| `PairedSamStats` | `base`, `paired_count`, `proper_pair_count`, `singletons`, `avg_insert_size` |
| `pair_sam_records(records) -> Vec<SamPair>` | Group records into mate pairs |
| `filter_proper_pairs(records) -> Vec<&SamRecord>` | Filter proper pairs (FLAG 0x2) |

### Pileup (`pileup.rs`, `sam` feature)

| Type/Function | Description |
|---------------|-------------|
| `PileupColumn` | Per-position column: `chrom`, `pos`, `ref_base`, `depth`, base counts, qualities |
| `Pileup` | Collection of pileup columns across a region |
| `generate_pileup(records, ref_name) -> Result<Pileup>` | Generate pileup from SAM records |
| `pileup_to_mpileup(pileup) -> String` | Format as samtools mpileup output |

### BGZF (`bgzf.rs`, `bam` feature)

| Type/Function | Description |
|---------------|-------------|
| `VirtualOffset` | BGZF virtual offset: block offset (48 bits) + within-block offset (16 bits) |
| `read_bgzf_block(reader) -> Result<Option<Vec<u8>>>` | Read and decompress next block |
| `decompress_bgzf(data) -> Result<Vec<u8>>` | Decompress full BGZF-compressed data |

### BAM (`bam.rs`, `bam` feature)

| Type/Function | Description |
|---------------|-------------|
| `parse_bam(path) -> Result<Vec<SamRecord>>` | Parse BAM (BGZF-compressed) format |
| `bam_stats(path) -> Result<SamStats>` | Streaming BAM statistics |
| `BamReference` | Reference sequence metadata from BAM header |

### BAM operations (`bam_ops.rs`, `bam` feature)

| Type/Function | Description |
|---------------|-------------|
| `SortOrder` | Enum: `Coordinate`, `Queryname`, `Unsorted` |
| `coordinate_sort(records, ref_order)` | Sort by reference + position |
| `queryname_sort(records)` | Sort by query name |
| `merge_sorted(a, b, ref_order) -> Vec<SamRecord>` | Merge two sorted record lists |
| `mark_duplicates(records) -> DuplicateReport` | Mark PCR/optical duplicates |
| `DuplicateReport` | Duplicate counts: `total`, `duplicates`, `optical`, `rate` |
| `depth_of_coverage(records, ref_name) -> Vec<u32>` | Per-position depth |
| `flagstat(records) -> FlagstatReport` | samtools flagstat-equivalent |

### Indexed BAM (`indexed_bam.rs`, `bam` feature)

| Type | Description |
|------|-------------|
| `IndexedBamReader` | Random-access BAM reader using BAI/CSI index |

| Method | Description |
|--------|-------------|
| `open(bam_path) -> Result<Self>` | Open BAM with auto-detected `.bai` index |
| `open_with_index(bam_path, bai_path) -> Result<Self>` | Open with explicit index path |
| `query(region) -> Result<Vec<SamRecord>>` | Fetch records overlapping a genomic region |
| `references() -> &[BamReference]` | Reference sequence metadata |

### Indexed VCF (`indexed_vcf.rs`, `vcf` feature)

| Type | Description |
|------|-------------|
| `IndexedVcfReader` | Random-access VCF reader using tabix (.tbi) index |

| Method | Description |
|--------|-------------|
| `open(vcf_gz_path) -> Result<Self>` | Open bgzipped VCF with auto-detected `.tbi` index |
| `open_with_index(vcf_gz_path, index_path) -> Result<Self>` | Open with explicit index |
| `query(region) -> Result<Vec<Variant>>` | Fetch variants overlapping a region |

### CRAM (`cram.rs`, `cram` feature)

| Type/Function | Description |
|---------------|-------------|
| `CramConfig` | Configuration: `reference_path: Option<PathBuf>` |
| `parse_cram(path, config) -> Result<Vec<SamRecord>>` | Parse CRAM format |
| `parse_cram_default(path) -> Result<Vec<SamRecord>>` | Parse without external reference |
| `cram_stats(path, config) -> Result<SamStats>` | CRAM summary statistics |

### BCF (`bcf.rs`, `bcf` feature)

| Type/Function | Description |
|---------------|-------------|
| `parse_bcf(path) -> Result<Vec<Variant>>` | Parse BCF2 binary VCF |
| `bcf_stats(path) -> Result<VcfStats>` | BCF summary statistics |

### BCF writer (`bcf_write.rs`, `bcf` feature)

| Function | Description |
|----------|-------------|
| `write_bcf(header, variants, path) -> Result<()>` | Write variants in BCF2 format to file |
| `write_bcf_bytes(header, variants) -> Result<Vec<u8>>` | Write to byte vector (BGZF-compressed) |

### Variant calling (`variant_call.rs`, `variant-calling` feature)

| Type | Description |
|------|-------------|
| `VariantCallConfig` | Configuration: min_depth, min_base_quality, min_mapq, priors, strand bias threshold |
| `VariantCallResult` | Called variant with genotype, quality, depth, allele counts |

| Function | Description |
|----------|-------------|
| `call_variants(pileup, config) -> Vec<VariantCallResult>` | Bayesian genotype calling from pileup |
| `pileup_to_vcf(results, header) -> String` | Format calls as VCF text |

### BLAST tabular (`blast.rs`, `blast` feature)

| Type/Function | Description |
|---------------|-------------|
| `BlastRecord` | 12-column tabular: `query_id`, `subject_id`, `pct_identity`, `evalue`, etc. |
| `parse_blast(path) -> Result<Vec<BlastRecord>>` | Parse `-outfmt 6`/`7` output |
| `BlastStats` | `hit_count`, `unique_queries`, `unique_subjects`, `avg_identity`, `avg_evalue` |
| `blast_stats(path) -> Result<BlastStats>` | Summary statistics |

### BLAST XML (`blast_xml.rs`, `blast` feature)

| Type | Description |
|------|-------------|
| `BlastXmlResult` | Top-level result: program, version, db, query, iterations |
| `BlastXmlIteration` | Search iteration with hits |
| `BlastXmlHit` | Database hit with HSPs |
| `BlastXmlHsp` | High-scoring segment pair: scores, positions, sequences |

| Function | Description |
|----------|-------------|
| `parse_blast_xml(input) -> Result<BlastXmlResult>` | Parse BLAST `-outfmt 5` XML output |

### MAF (`maf.rs`, `maf` feature)

| Type/Function | Description |
|---------------|-------------|
| `MafBlock` | `score`, `sequences: Vec<MafSequence>` |
| `MafSequence` | `src`, `start`, `size`, `strand`, `src_size`, `text` |
| `parse_maf(path) -> Result<Vec<MafBlock>>` | Parse MAF alignment blocks |
| `MafStats` | `block_count`, `total_aligned_bases`, `species` |
| `maf_stats(path) -> Result<MafStats>` | Summary statistics |

### GenBank (`genbank.rs`, `genbank` feature)

| Type/Function | Description |
|---------------|-------------|
| `GenbankRecord` | `locus`, `definition`, `accession`, `version`, `organism`, `sequence`, `features` |
| `GenbankFeature` | `feature_type`, `location`, `qualifiers` |
| `parse_genbank(path) -> Result<Vec<GenbankRecord>>` | Parse multi-record GenBank flat file |
| `GenbankStats` | `record_count`, `total_bases`, `feature_counts` |
| `genbank_stats(path) -> Result<GenbankStats>` | Summary statistics |

### EMBL (`embl.rs`, `genbank` feature)

| Type | Description |
|------|-------------|
| `EmblRecord` | `id`, `accession`, `description`, `sequence`, `features` |

| Function | Description |
|----------|-------------|
| `parse_embl(input) -> Result<Vec<EmblRecord>>` | Parse EMBL/ENA format records |
| `write_embl(records) -> String` | Write EMBL format |

### bigWig/bigBed (`bigwig.rs`, `bigwig` feature)

| Type/Function | Description |
|---------------|-------------|
| `BigWigHeader` | `version`, `zoom_levels`, `chrom_count`, `total_summary` |
| `BigWigSummary` | `bases_covered`, `min_val`, `max_val`, `sum`, `sum_squares` |
| `BigWigInterval` | `chrom`, `start`, `end`, `value` |
| `BigBedRecord` | `chrom`, `start`, `end`, `rest` |
| `read_bigwig_header(path) -> Result<BigWigHeader>` | Read bigWig header + summary |
| `read_bigwig_intervals(path, chrom, start, end) -> Result<Vec<BigWigInterval>>` | Query intervals |
| `read_bigbed_header(path) -> Result<BigWigHeader>` | Read bigBed header |
| `read_bigbed_records(path, chrom, start, end) -> Result<Vec<BigBedRecord>>` | Query records |

### Parquet (`parquet.rs`, `parquet` feature)

| Type/Function | Description |
|---------------|-------------|
| `ParquetInfo` | `num_rows`, `num_columns`, `column_names`, `num_row_groups`, `created_by` |
| `parquet_info(path) -> Result<ParquetInfo>` | Extract Parquet metadata |
| `write_variants_parquet(variants, path) -> Result<()>` | Write variants to Parquet |
| `read_variants_parquet(path) -> Result<Vec<Variant>>` | Read variants from Parquet |
| `read_variants_parquet_region(path, chrom, start, end) -> Result<Vec<Variant>>` | Region query with pushdown |
| `write_intervals_parquet(intervals, path) -> Result<()>` | Write intervals to Parquet |
| `read_intervals_parquet(path) -> Result<Vec<GenomicInterval>>` | Read intervals |
| `ExpressionMatrix` | `gene_ids`, `sample_ids`, `values` |
| `write_expression_parquet(matrix, path) -> Result<()>` | Write expression matrix |
| `read_expression_parquet(path) -> Result<ExpressionMatrix>` | Read expression matrix |

### Stockholm (`stockholm.rs`, `genbank` feature)

| Type | Description |
|------|-------------|
| `StockholmAlignment` | Aligned sequences, GC/GS/GR annotations |

| Function | Description |
|----------|-------------|
| `parse_stockholm(input) -> Result<Vec<StockholmAlignment>>` | Parse Stockholm 1.0 format |
| `write_stockholm(alignments) -> String` | Write Stockholm format |

### Clustal (`clustal.rs`, `genbank` feature)

| Type | Description |
|------|-------------|
| `ClustalAlignment` | Aligned sequences with conservation string |

| Function | Description |
|----------|-------------|
| `parse_clustal(input) -> Result<ClustalAlignment>` | Parse ClustalW/Omega format |
| `write_clustal(alignment) -> String` | Write Clustal format |

### Phylip (`phylip.rs`, `genbank` feature)

| Type | Description |
|------|-------------|
| `PhylipAlignment` | Aligned sequences with `n_taxa`, `n_sites` |

| Function | Description |
|----------|-------------|
| `parse_phylip(input) -> Result<PhylipAlignment>` | Parse PHYLIP interleaved format |
| `parse_phylip_sequential(input) -> Result<PhylipAlignment>` | Parse PHYLIP sequential format |
| `write_phylip(alignment) -> String` | Write PHYLIP interleaved format |

### PIR/NBRF (`pir.rs`, `genbank` feature)

| Type | Description |
|------|-------------|
| `PirRecord` | `entry_type`, `name`, `description`, `sequence` |

| Function | Description |
|----------|-------------|
| `parse_pir(input) -> Result<Vec<PirRecord>>` | Parse PIR/NBRF format |
| `write_pir(records) -> String` | Write PIR format |

### ABI chromatogram (`abi.rs`, `genbank` feature)

| Type | Description |
|------|-------------|
| `AbiRecord` | `sequence`, `quality`, `traces`, `sample_name`, `peak_positions` |
| `AbiTraces` | Raw fluorescence: `a`, `c`, `g`, `t` channels |

| Function | Description |
|----------|-------------|
| `parse_abi(data) -> Result<AbiRecord>` | Parse ABI (.ab1) binary trace file |

### bedGraph/Wiggle (`bedgraph.rs`, `bed` feature)

| Type | Description |
|------|-------------|
| `BedGraphRecord` | `chrom`, `start`, `end`, `value` |

| Function | Description |
|----------|-------------|
| `parse_bedgraph_str(input) -> Result<Vec<BedGraphRecord>>` | Parse bedGraph format from string |
| `parse_bedgraph(path) -> Result<Vec<BedGraphRecord>>` | Parse bedGraph file |
| `parse_wiggle_str(input) -> Result<Vec<BedGraphRecord>>` | Parse WIG format (variableStep/fixedStep) |
| `write_bedgraph(records) -> String` | Write bedGraph format |

### GFA (`gfa.rs`, `genbank` feature)

| Type | Description |
|------|-------------|
| `GfaSegment` | Segment: `name`, `sequence`, `length` |
| `GfaLink` | Link: `from_segment`, `from_orient`, `to_segment`, `to_orient`, `overlap` |
| `GfaPath` | Path: `name`, `segment_names`, `overlaps` |
| `GfaGraph` | Complete GFA graph: segments, links, paths |

| Function | Description |
|----------|-------------|
| `parse_gfa(input) -> Result<GfaGraph>` | Parse GFA v1 format |
| `write_gfa(graph) -> String` | Write GFA v1 format |

### Fetch helpers (`fetch.rs`, `fetch` feature)

URL builders and response parsers for bioinformatics APIs. No HTTP client -- pure URL construction and text parsing.

| Type | Description |
|------|-------------|
| `EntrezUrl` | NCBI Entrez E-utilities URL builder (esearch, efetch, esummary) |
| `UniProtUrl` | UniProt REST API URL builder |
| `KeggUrl` | KEGG API URL builder |
| `HtsgetUrl` | htsget protocol URL builder |
| `RefgetUrl` | refget protocol URL builder |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `csv` | Yes | CSV parsing (`csv`, `serde`, `serde_json` deps) |
| `vcf` | No | VCF variant parsing (requires `cyanea-omics`) |
| `bed` | No | BED/BEDPE/bedGraph interval parsing (requires `cyanea-omics`) |
| `gff` | No | GFF3 gene structure parsing (requires `cyanea-omics`) |
| `gtf` | No | GTF (GFF2) parsing (requires `cyanea-omics`) |
| `sam` | No | SAM text alignment parsing, pileup generation |
| `bam` | No | BAM binary parsing, BAM ops, indexed BAM (implies `sam`, requires `flate2`) |
| `cram` | No | CRAM format (implies `sam`, requires noodles) |
| `bcf` | No | BCF binary VCF parsing and writing (implies `vcf`, requires `flate2`) |
| `blast` | No | BLAST tabular and XML output parsing |
| `maf` | No | MAF (Multiple Alignment Format) parsing |
| `genbank` | No | GenBank, EMBL, Stockholm, Clustal, Phylip, PIR, ABI, GFA |
| `bigwig` | No | bigWig/bigBed binary format (requires `flate2`) |
| `parquet` | No | Apache Parquet columnar format (implies `vcf` + `bed`, arrow/parquet deps) |
| `variant-calling` | No | Bayesian variant caller (implies `sam` + `vcf`) |
| `fetch` | No | URL builders for NCBI/UniProt/KEGG/htsget/refget |
| `parallel` | No | Rayon parallelism |
| `wasm` | No | WASM target |

## Dependencies

- `cyanea-core` -- error types
- `cyanea-omics` -- genomic types (`Variant`, `GenomicInterval`, `Gene`), optional
- `csv`, `serde`, `serde_json` -- CSV parsing, optional
- `flate2` -- BAM/BCF/bigWig BGZF/zlib decompression, optional
- `noodles-*` -- CRAM, indexed BAM/VCF parsing, optional
- `parquet`, `arrow` -- Apache Parquet/Arrow, optional
- `cyanea-stats` -- Fisher exact test for variant calling, optional

## Tests

357 unit + 7 doc tests with all features enabled. Breakdown by module: CSV (3), VCF (25), VCF ops (32), BED (9), BEDPE (10), GFF3 (9), GTF (6), SAM (32), pileup (36), BAM (35), indexed BAM/VCF (8), BGZF (4), CRAM (7), Parquet (14), BCF (12), BCF write (included in BCF), variant calling (21), BLAST tabular (5), BLAST XML (5), MAF (5), GenBank (6), bigWig (6), Stockholm (6), Clustal (5), Phylip (5), EMBL (5), PIR (6), ABI (5), bedGraph (5), GFA (5), fetch (12).
