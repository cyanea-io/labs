# cyanea-io

Unified file format parsing for bioinformatics and tabular data. Each parser is behind a feature flag to keep the dependency tree minimal.

## Status: Complete

CSV, VCF, BED, and GFF3 parsers are fully implemented. Each format includes both raw record parsing and summary statistics.

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

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `csv` | Yes | CSV parsing (`csv`, `serde`, `serde_json` deps) |
| `vcf` | No | VCF variant parsing (requires `cyanea-omics`) |
| `bed` | No | BED interval parsing (requires `cyanea-omics`) |
| `gff` | No | GFF3 gene structure parsing (requires `cyanea-omics`) |
| `wasm` | No | WASM target |

## Dependencies

- `cyanea-core` -- error types
- `cyanea-omics` -- genomic types (`Variant`, `GenomicInterval`, `Gene`), optional
- `csv`, `serde`, `serde_json` -- CSV parsing, optional

## Tests

26 tests across 4 source files: CSV (4), VCF (5), BED (9), GFF3 (8).

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 33 | Module declarations, re-exports |
| `csv.rs` | 134 | CSV info extraction |
| `vcf.rs` | 262 | VCF variant parsing |
| `bed.rs` | 278 | BED interval parsing |
| `gff.rs` | 618 | GFF3 hierarchical gene parsing |
