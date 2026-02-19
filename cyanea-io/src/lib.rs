//! File format parsing for the Cyanea bioinformatics ecosystem.
//!
//! Supported formats:
//! - **CSV/TSV** — via the `csv` feature (enabled by default)
//! - **VCF** — Variant Call Format, via the `vcf` feature
//! - **BED** — Browser Extensible Data, via the `bed` feature
//! - **GFF3** — General Feature Format, via the `gff` feature
//! - **SAM** — Sequence Alignment/Map (text format), via the `sam` feature

#[cfg(feature = "csv")]
pub mod csv;

#[cfg(feature = "vcf")]
pub mod vcf;

#[cfg(feature = "bed")]
pub mod bed;

#[cfg(feature = "bed")]
pub mod bedpe;

#[cfg(feature = "gff")]
pub mod gff;

#[cfg(feature = "sam")]
pub mod sam;

#[cfg(feature = "sam")]
pub mod pileup;

#[cfg(any(feature = "bam", feature = "bcf"))]
pub mod bgzf;

#[cfg(feature = "bam")]
pub mod bam;

#[cfg(feature = "bam")]
pub mod bam_ops;

#[cfg(feature = "parquet")]
pub mod parquet;

#[cfg(feature = "cram")]
pub mod cram;

#[cfg(feature = "blast")]
pub mod blast;

#[cfg(feature = "maf")]
pub mod maf;

#[cfg(feature = "gtf")]
pub mod gtf;

#[cfg(feature = "genbank")]
pub mod genbank;

#[cfg(feature = "bigwig")]
pub mod bigwig;

#[cfg(feature = "vcf")]
pub mod vcf_header;

#[cfg(feature = "vcf")]
pub mod vcf_ops;

#[cfg(feature = "bcf")]
pub mod bcf;

#[cfg(feature = "bcf")]
pub mod bcf_write;

#[cfg(feature = "indexed-bam")]
pub mod indexed_bam;

#[cfg(feature = "indexed-vcf")]
pub mod indexed_vcf;

#[cfg(feature = "variant-calling")]
pub mod variant_call;

#[cfg(feature = "stockholm")]
pub mod stockholm;

#[cfg(feature = "clustal")]
pub mod clustal;

#[cfg(feature = "phylip")]
pub mod phylip;

#[cfg(feature = "embl")]
pub mod embl;

#[cfg(feature = "pir")]
pub mod pir;

#[cfg(feature = "bedgraph")]
pub mod bedgraph;

#[cfg(feature = "gfa")]
pub mod gfa;

#[cfg(feature = "blast-xml")]
pub mod blast_xml;

#[cfg(feature = "abi")]
pub mod abi;

#[cfg(feature = "fetch")]
pub mod fetch;

// Re-exports for convenience.

#[cfg(feature = "csv")]
pub use csv::{csv_preview, parse_csv_info, CsvInfo};

#[cfg(feature = "vcf")]
pub use vcf::{parse_vcf, parse_vcf_str, vcf_stats, write_vcf, write_vcf_string, VcfStats};

#[cfg(feature = "bed")]
pub use bed::{bed_stats, parse_bed, parse_bed_intervals, parse_bed_str, BedRecord, BedStats};

#[cfg(feature = "bed")]
pub use bedpe::{bedpe_stats, parse_bedpe, BedpeRecord, BedpeStats};

#[cfg(feature = "gff")]
pub use gff::{gff3_stats, parse_gff3, parse_gff3_str, GffStats};

#[cfg(feature = "sam")]
pub use sam::{
    filter_proper_pairs, pair_sam_records, paired_sam_stats, parse_sam, parse_sam_str, sam_stats,
    sam_stats_from_path, PairedSamStats, SamPair, SamRecord, SamStats,
};

#[cfg(feature = "sam")]
pub use pileup::{
    depth_stats, pileup, pileup_region, pileup_to_mpileup, DepthStats, InsertionEvidence, Pileup,
    PileupColumn,
};

#[cfg(all(feature = "sam", feature = "vcf"))]
#[allow(deprecated)]
pub use pileup::call_snps;

#[cfg(feature = "bam")]
pub use bam::{parse_bam, bam_stats, BamReference};

#[cfg(feature = "bam")]
pub use bam_ops::{
    alignment_stats, coordinate_sort, flagstat, fixmate, idxstats, mark_duplicates,
    merge_bam_records, queryname_sort, AlignmentStats, DuplicateStats, FlagStats, IdxStats,
    SortOrder,
};

#[cfg(feature = "parquet")]
pub use parquet::{
    parquet_info, parquet_interval_stats, parquet_variant_stats, read_expression_parquet,
    read_intervals_parquet, read_intervals_parquet_region, read_variants_parquet,
    read_variants_parquet_region, write_expression_parquet, write_intervals_parquet,
    write_variants_parquet, ExpressionMatrix, ParquetInfo,
};

#[cfg(feature = "cram")]
pub use cram::{cram_stats, cram_stats_default, parse_cram, parse_cram_default, CramConfig};

#[cfg(feature = "blast")]
pub use blast::{blast_stats, parse_blast, BlastRecord, BlastStats};

#[cfg(feature = "maf")]
pub use maf::{maf_stats, parse_maf, MafBlock, MafSequence, MafStats};

#[cfg(feature = "gtf")]
pub use gtf::{gtf_stats, parse_gtf, GtfStats};

#[cfg(feature = "genbank")]
pub use genbank::{genbank_stats, parse_genbank, GenbankFeature, GenbankRecord, GenbankStats};

#[cfg(feature = "bigwig")]
pub use bigwig::{
    read_bigbed_header, read_bigbed_records, read_bigwig_header, read_bigwig_intervals,
    BigBedRecord, BigWigHeader, BigWigInterval, BigWigSummary,
};

#[cfg(feature = "vcf")]
pub use vcf_header::{ContigLine, FieldDef, FilterDef, VcfHeader};

#[cfg(feature = "vcf")]
pub use vcf_ops::{
    annotate_variant_ids, annotate_variants, detailed_vcf_stats, eval_filter, filter_variants,
    join_biallelic, merge_variants, normalize_variant, parse_filter, split_multiallelic,
    variant_complement, variant_concordance, variant_intersection, AnnotationSource, CmpOp,
    ConcordanceStats, DetailedVcfStats, FilterExpr,
};

#[cfg(feature = "bcf")]
pub use bcf::{bcf_stats, parse_bcf};

#[cfg(feature = "bcf")]
pub use bcf_write::{write_bcf, write_bcf_bytes};

#[cfg(feature = "indexed-bam")]
pub use indexed_bam::{create_bai_index, fetch_bam, IndexedBamReader};

#[cfg(feature = "indexed-vcf")]
pub use indexed_vcf::{fetch_vcf, IndexedVcfReader};

#[cfg(feature = "variant-calling")]
pub use variant_call::{
    call_variants, call_variants_all, variant_call_stats, CalledVariant, Genotype,
    VariantCallConfig, VariantCallStats,
};

#[cfg(feature = "variant-calling")]
pub use vcf::{write_called_vcf, write_called_vcf_string};

#[cfg(feature = "stockholm")]
pub use stockholm::{parse_stockholm, write_stockholm, StockholmAlignment};

#[cfg(feature = "clustal")]
pub use clustal::{parse_clustal, write_clustal, ClustalAlignment};

#[cfg(feature = "phylip")]
pub use phylip::{
    parse_phylip, parse_phylip_sequential, write_phylip, PhylipAlignment,
};

#[cfg(feature = "embl")]
pub use embl::{parse_embl, write_embl, EmblRecord};

#[cfg(feature = "pir")]
pub use pir::{parse_pir, write_pir, PirRecord};

#[cfg(feature = "bedgraph")]
pub use bedgraph::{
    parse_bedgraph_str, parse_wiggle_str, write_bedgraph_string, BedGraphRecord,
};

#[cfg(feature = "gfa")]
pub use gfa::{
    gfa_stats, parse_gfa_str, write_gfa_string, GfaGraph, GfaLink, GfaPath, GfaSegment, GfaStats,
};

#[cfg(feature = "blast-xml")]
pub use blast_xml::{
    parse_blast_xml_str, BlastXmlHit, BlastXmlHsp, BlastXmlIteration, BlastXmlResult,
};

#[cfg(feature = "abi")]
pub use abi::{parse_abi_bytes, AbiRecord, AbiTraces};

#[cfg(feature = "fetch")]
pub use fetch::{
    parse_efetch_fasta, parse_esearch_response, parse_refget_metadata, EntrezUrl, HtsgetUrl,
    KeggUrl, RefgetAlias, RefgetMetadata, RefgetUrl, UniprotUrl,
};
