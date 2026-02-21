# cyanea-io

> Unified file format parsing for bioinformatics. Each parser is behind a feature flag to keep the dependency tree minimal.

## What's Inside

- **CSV** -- metadata extraction, preview
- **VCF** -- variant parsing into `cyanea-omics::Variant`, streaming stats
- **BED** -- BED3-BED6 records, interval conversion, statistics
- **BEDPE** -- paired-end intervals with inter/intra-chromosomal stats
- **GFF3** -- hierarchical Gene/Transcript/Exon assembly with coordinate conversion
- **GTF** -- GFF2 format with `gene_id`/`transcript_id` hierarchy
- **SAM** -- full record parsing, flag helpers, paired-end stats, pileup generation
- **BAM** -- BGZF-compressed binary alignment parsing
- **CRAM** -- reference-based alignment format via noodles
- **BCF** -- BCF2.1 binary VCF parsing
- **BLAST** -- tabular output (`-outfmt 6`/`7`) parsing
- **MAF** -- Multiple Alignment Format (UCSC/LAST/minimap2)
- **GenBank** -- flat file parsing (multi-record, features table)
- **bigWig/bigBed** -- Kent binary formats with B+ tree and R-tree index
- **Parquet** -- columnar storage for variants, intervals, and expression matrices with predicate pushdown

## Quick Start

```toml
[dependencies]
cyanea-io = { version = "0.1", features = ["vcf", "bed", "sam"] }
```

```rust
use cyanea_io::{vcf::vcf_stats, bed::parse_bed, sam::parse_sam};

let stats = vcf_stats("variants.vcf").unwrap();
println!("{} variants, {} SNVs", stats.variant_count, stats.snv_count);

let regions = parse_bed("regions.bed").unwrap();
let alignments = parse_sam("reads.sam").unwrap();
```

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `csv` | Yes | CSV parsing (csv, serde, serde_json) |
| `vcf` | No | VCF variant parsing (requires cyanea-omics) |
| `bed` | No | BED interval parsing (requires cyanea-omics) |
| `gff` | No | GFF3 gene structure parsing (requires cyanea-omics) |
| `gtf` | No | GTF (GFF2) parsing (requires cyanea-omics) |
| `sam` | No | SAM text alignment parsing |
| `bam` | No | BAM binary parsing (implies `sam`, adds flate2) |
| `cram` | No | CRAM format (implies `sam`, adds noodles) |
| `bcf` | No | BCF binary VCF (implies `vcf`, adds flate2) |
| `blast` | No | BLAST tabular output parsing |
| `maf` | No | MAF alignment format parsing |
| `genbank` | No | GenBank flat file parsing |
| `bigwig` | No | bigWig/bigBed binary format (adds flate2) |
| `parquet` | No | Apache Parquet (implies `vcf` + `bed`, adds arrow/parquet) |
| `variant-calling` | No | Variant calling (implies `sam` + `vcf`) |
| `parallel` | No | Rayon parallelism |
| `wasm` | No | WASM target marker |

## Modules

| Module | Feature | Description |
|--------|---------|-------------|
| `csv` | `csv` | CSV metadata and preview |
| `vcf` | `vcf` | VCF parsing and statistics |
| `vcf_header` | `vcf` | Structured VCF 4.3 header construction/parsing |
| `vcf_ops` | `vcf` | Normalization, multi-allelic split/join, filtering, set ops |
| `indexed_vcf` | `vcf` | Random-access VCF via tabix index (noodles) |
| `bed` | `bed` | BED3-BED6 record parsing |
| `bedpe` | `bed` | BEDPE paired-end intervals |
| `bedgraph` | `bed` | bedGraph/Wiggle signal tracks |
| `gff` | `gff` | GFF3 hierarchical parsing |
| `gtf` | `gtf` | GTF (GFF2) parsing |
| `sam` | `sam` | SAM records, flags, paired stats |
| `pileup` | `sam` | Pileup generation, mpileup output |
| `bgzf` | `bam` | BGZF block decompression, virtual offsets |
| `bam` | `bam` | BAM binary alignment parsing |
| `bam_ops` | `bam` | Sort, merge, mark duplicates, flagstat, depth |
| `indexed_bam` | `bam` | Random-access BAM via BAI/CSI index (noodles) |
| `cram` | `cram` | CRAM format via noodles |
| `bcf` | `bcf` | BCF2.1 binary VCF reader |
| `bcf_write` | `bcf` | BCF2 binary VCF writer |
| `variant_call` | `variant-calling` | Bayesian genotype caller from pileup |
| `blast` | `blast` | BLAST tabular output (-outfmt 6/7) |
| `blast_xml` | `blast` | BLAST XML output (-outfmt 5) |
| `maf` | `maf` | Multiple Alignment Format |
| `genbank` | `genbank` | GenBank flat file parsing |
| `embl` | `genbank` | EMBL/ENA format parsing |
| `stockholm` | `genbank` | Stockholm MSA format (Pfam/Rfam/HMMER) |
| `clustal` | `genbank` | ClustalW/Omega format |
| `phylip` | `genbank` | PHYLIP interleaved/sequential format |
| `pir` | `genbank` | PIR/NBRF protein format |
| `abi` | `genbank` | ABI Sanger chromatogram binary |
| `gfa` | `genbank` | GFA v1 sequence graph format |
| `bigwig` | `bigwig` | bigWig/bigBed Kent binary format |
| `parquet` | `parquet` | Apache Parquet columnar format |
| `fetch` | `fetch` | URL builders for NCBI/UniProt/KEGG/htsget/refget |

## See Also

- [API Reference](docs/API.md)
- [Usage Guide](docs/GUIDE.md)
- [Internal Architecture](docs/ARCHITECTURE.md)
- [Workspace Architecture](../docs/ARCHITECTURE.md)
- [Build Guide](../docs/BUILDING.md)
