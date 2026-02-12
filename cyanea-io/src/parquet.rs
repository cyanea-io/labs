//! Parquet columnar format for genomic tabular data.
//!
//! Provides read/write support for storing variants (VCF-style) and
//! genomic intervals (BED-style) in Apache Parquet columnar format.
//! Parquet provides efficient compression and columnar access patterns
//! ideal for large-scale genomic datasets.

use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{Array as ArrowArray, ArrayRef, Float64Array, StringArray, UInt64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;

use cyanea_core::{CyaneaError, Result};
use cyanea_omics::genomic::{GenomicInterval, Strand};
use cyanea_omics::variant::{Variant, VariantFilter};

use crate::bed::BedStats;
use crate::vcf::VcfStats;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn open_file(path: &Path) -> Result<File> {
    File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })
}

fn create_file(path: &Path) -> Result<File> {
    File::create(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })
}

fn pq_err(e: impl std::fmt::Display, path: &Path) -> CyaneaError {
    CyaneaError::Parse(format!("{}: {e}", path.display()))
}

// ---------------------------------------------------------------------------
// Parquet info
// ---------------------------------------------------------------------------

/// Metadata about a Parquet file.
#[derive(Debug, Clone)]
pub struct ParquetInfo {
    /// Total number of rows across all row groups.
    pub num_rows: usize,
    /// Number of columns in the schema.
    pub num_columns: usize,
    /// Column names from the schema.
    pub column_names: Vec<String>,
    /// Number of row groups in the file.
    pub num_row_groups: usize,
    /// Creator string from file metadata, if present.
    pub created_by: Option<String>,
}

/// Extract metadata from a Parquet file.
pub fn parquet_info(path: impl AsRef<Path>) -> Result<ParquetInfo> {
    let path = path.as_ref();
    let file = open_file(path)?;
    let builder =
        ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| pq_err(e, path))?;
    let metadata = builder.metadata();
    let file_meta = metadata.file_metadata();
    let schema_desc = file_meta.schema_descr();

    Ok(ParquetInfo {
        num_rows: file_meta.num_rows() as usize,
        num_columns: schema_desc.num_columns(),
        column_names: (0..schema_desc.num_columns())
            .map(|i| schema_desc.column(i).name().to_string())
            .collect(),
        num_row_groups: metadata.num_row_groups(),
        created_by: file_meta.created_by().map(|s| s.to_string()),
    })
}

// ---------------------------------------------------------------------------
// Variant Parquet (VCF-as-Parquet)
// ---------------------------------------------------------------------------

fn variant_schema() -> Schema {
    Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("position", DataType::UInt64, false),
        Field::new("id", DataType::Utf8, true),
        Field::new("ref_allele", DataType::Utf8, false),
        Field::new("alt_alleles", DataType::Utf8, false),
        Field::new("quality", DataType::Float64, true),
        Field::new("filter", DataType::Utf8, false),
    ])
}

fn filter_to_owned(f: &VariantFilter) -> String {
    match f {
        VariantFilter::Pass => "PASS".to_string(),
        VariantFilter::Missing => ".".to_string(),
        VariantFilter::Fail(filters) => filters.join(";"),
    }
}

fn string_to_filter(s: &str) -> VariantFilter {
    match s {
        "PASS" => VariantFilter::Pass,
        "." => VariantFilter::Missing,
        other => VariantFilter::Fail(other.split(';').map(|f| f.to_string()).collect()),
    }
}

/// Write a slice of [`Variant`] records to a Parquet file.
pub fn write_variants_parquet(
    variants: &[Variant],
    path: impl AsRef<Path>,
) -> Result<()> {
    let path = path.as_ref();
    let schema = Arc::new(variant_schema());

    let chroms: Vec<&str> = variants.iter().map(|v| v.chrom.as_str()).collect();
    let positions: Vec<u64> = variants.iter().map(|v| v.position).collect();
    let ids: Vec<Option<&str>> = variants.iter().map(|v| v.id.as_deref()).collect();
    let refs: Vec<&str> = variants
        .iter()
        .map(|v| std::str::from_utf8(&v.ref_allele).unwrap_or(""))
        .collect();
    let alts: Vec<String> = variants
        .iter()
        .map(|v| {
            v.alt_alleles
                .iter()
                .map(|a| String::from_utf8_lossy(a).to_string())
                .collect::<Vec<_>>()
                .join(",")
        })
        .collect();
    let alt_refs: Vec<&str> = alts.iter().map(|s| s.as_str()).collect();
    let qualities: Vec<Option<f64>> = variants.iter().map(|v| v.quality).collect();
    let filters: Vec<String> = variants.iter().map(|v| filter_to_owned(&v.filter)).collect();
    let filter_refs: Vec<&str> = filters.iter().map(|s| s.as_str()).collect();

    let columns: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from(chroms)),
        Arc::new(UInt64Array::from(positions)),
        Arc::new(StringArray::from(ids)),
        Arc::new(StringArray::from(refs)),
        Arc::new(StringArray::from(alt_refs)),
        Arc::new(Float64Array::from(qualities)),
        Arc::new(StringArray::from(filter_refs)),
    ];

    let batch =
        RecordBatch::try_new(schema.clone(), columns).map_err(|e| pq_err(e, path))?;

    let file = create_file(path)?;
    let mut writer =
        ArrowWriter::try_new(file, schema, None).map_err(|e| pq_err(e, path))?;
    writer.write(&batch).map_err(|e| pq_err(e, path))?;
    writer.close().map_err(|e| pq_err(e, path))?;

    Ok(())
}

/// Read [`Variant`] records from a Parquet file.
pub fn read_variants_parquet(path: impl AsRef<Path>) -> Result<Vec<Variant>> {
    let path = path.as_ref();
    let file = open_file(path)?;
    let builder =
        ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| pq_err(e, path))?;
    let reader = builder.build().map_err(|e| pq_err(e, path))?;

    let mut variants = Vec::new();
    for batch_result in reader {
        let batch = batch_result.map_err(|e| pq_err(e, path))?;
        let n = batch.num_rows();

        let chrom_col = batch
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| pq_err("expected Utf8 for chrom", path))?;
        let pos_col = batch
            .column(1)
            .as_any()
            .downcast_ref::<UInt64Array>()
            .ok_or_else(|| pq_err("expected UInt64 for position", path))?;
        let id_col = batch
            .column(2)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| pq_err("expected Utf8 for id", path))?;
        let ref_col = batch
            .column(3)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| pq_err("expected Utf8 for ref_allele", path))?;
        let alt_col = batch
            .column(4)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| pq_err("expected Utf8 for alt_alleles", path))?;
        let qual_col = batch
            .column(5)
            .as_any()
            .downcast_ref::<Float64Array>()
            .ok_or_else(|| pq_err("expected Float64 for quality", path))?;
        let filt_col = batch
            .column(6)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| pq_err("expected Utf8 for filter", path))?;

        for i in 0..n {
            let chrom = chrom_col.value(i).to_string();
            let position = pos_col.value(i);
            let id = if id_col.is_null(i) {
                None
            } else {
                Some(id_col.value(i).to_string())
            };
            let ref_allele = ref_col.value(i).as_bytes().to_vec();
            let alt_alleles: Vec<Vec<u8>> = alt_col
                .value(i)
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|a| a.as_bytes().to_vec())
                .collect();
            let quality = if qual_col.is_null(i) {
                None
            } else {
                Some(qual_col.value(i))
            };
            let filter = string_to_filter(filt_col.value(i));

            variants.push(Variant {
                chrom,
                position,
                id,
                ref_allele,
                alt_alleles,
                quality,
                filter,
            });
        }
    }

    Ok(variants)
}

/// Compute variant statistics from a Parquet file.
///
/// Returns the same [`VcfStats`] type used by the VCF module.
pub fn parquet_variant_stats(path: impl AsRef<Path>) -> Result<VcfStats> {
    let variants = read_variants_parquet(path)?;

    let mut snv_count: u64 = 0;
    let mut indel_count: u64 = 0;
    let mut pass_count: u64 = 0;
    let mut chroms = Vec::new();
    let mut seen = std::collections::HashSet::new();

    for v in &variants {
        if v.is_snv() {
            snv_count += 1;
        }
        if v.is_indel() {
            indel_count += 1;
        }
        if v.filter == VariantFilter::Pass {
            pass_count += 1;
        }
        if seen.insert(v.chrom.clone()) {
            chroms.push(v.chrom.clone());
        }
    }

    Ok(VcfStats {
        variant_count: variants.len() as u64,
        snv_count,
        indel_count,
        pass_count,
        chromosomes: chroms,
    })
}

// ---------------------------------------------------------------------------
// Interval Parquet (BED-as-Parquet)
// ---------------------------------------------------------------------------

fn interval_schema() -> Schema {
    Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::UInt64, false),
        Field::new("end", DataType::UInt64, false),
        Field::new("strand", DataType::Utf8, false),
    ])
}

/// Write a slice of [`GenomicInterval`] records to a Parquet file.
pub fn write_intervals_parquet(
    intervals: &[GenomicInterval],
    path: impl AsRef<Path>,
) -> Result<()> {
    let path = path.as_ref();
    let schema = Arc::new(interval_schema());

    let chroms: Vec<&str> = intervals.iter().map(|iv| iv.chrom.as_str()).collect();
    let starts: Vec<u64> = intervals.iter().map(|iv| iv.start).collect();
    let ends: Vec<u64> = intervals.iter().map(|iv| iv.end).collect();
    let strands: Vec<String> = intervals.iter().map(|iv| iv.strand.to_string()).collect();
    let strand_refs: Vec<&str> = strands.iter().map(|s| s.as_str()).collect();

    let columns: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from(chroms)),
        Arc::new(UInt64Array::from(starts)),
        Arc::new(UInt64Array::from(ends)),
        Arc::new(StringArray::from(strand_refs)),
    ];

    let batch =
        RecordBatch::try_new(schema.clone(), columns).map_err(|e| pq_err(e, path))?;

    let file = create_file(path)?;
    let mut writer =
        ArrowWriter::try_new(file, schema, None).map_err(|e| pq_err(e, path))?;
    writer.write(&batch).map_err(|e| pq_err(e, path))?;
    writer.close().map_err(|e| pq_err(e, path))?;

    Ok(())
}

/// Read [`GenomicInterval`] records from a Parquet file.
pub fn read_intervals_parquet(path: impl AsRef<Path>) -> Result<Vec<GenomicInterval>> {
    let path = path.as_ref();
    let file = open_file(path)?;
    let builder =
        ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| pq_err(e, path))?;
    let reader = builder.build().map_err(|e| pq_err(e, path))?;

    let mut intervals = Vec::new();
    for batch_result in reader {
        let batch = batch_result.map_err(|e| pq_err(e, path))?;
        let n = batch.num_rows();

        let chrom_col = batch
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| pq_err("expected Utf8 for chrom", path))?;
        let start_col = batch
            .column(1)
            .as_any()
            .downcast_ref::<UInt64Array>()
            .ok_or_else(|| pq_err("expected UInt64 for start", path))?;
        let end_col = batch
            .column(2)
            .as_any()
            .downcast_ref::<UInt64Array>()
            .ok_or_else(|| pq_err("expected UInt64 for end", path))?;
        let strand_col = batch
            .column(3)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| pq_err("expected Utf8 for strand", path))?;

        for i in 0..n {
            let chrom = chrom_col.value(i).to_string();
            let start = start_col.value(i);
            let end = end_col.value(i);
            let strand = match strand_col.value(i) {
                "+" => Strand::Forward,
                "-" => Strand::Reverse,
                _ => Strand::Unknown,
            };
            intervals.push(GenomicInterval::with_strand(chrom, start, end, strand)?);
        }
    }

    Ok(intervals)
}

/// Compute interval statistics from a Parquet file.
///
/// Returns the same [`BedStats`] type used by the BED module.
pub fn parquet_interval_stats(path: impl AsRef<Path>) -> Result<BedStats> {
    let intervals = read_intervals_parquet(path)?;

    let mut total_bases: u64 = 0;
    let mut chroms = Vec::new();
    let mut seen = std::collections::HashSet::new();

    for iv in &intervals {
        total_bases += iv.len();
        if seen.insert(iv.chrom.clone()) {
            chroms.push(iv.chrom.clone());
        }
    }

    Ok(BedStats {
        record_count: intervals.len() as u64,
        total_bases,
        chromosomes: chroms,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    fn temp_parquet() -> (NamedTempFile, std::path::PathBuf) {
        let f = NamedTempFile::new().unwrap();
        let p = f.path().with_extension("parquet");
        (f, p)
    }

    fn sample_variants() -> Vec<Variant> {
        vec![
            Variant {
                chrom: "chr1".into(),
                position: 100,
                id: Some("rs123".into()),
                ref_allele: b"A".to_vec(),
                alt_alleles: vec![b"G".to_vec()],
                quality: Some(30.0),
                filter: VariantFilter::Pass,
            },
            Variant {
                chrom: "chr1".into(),
                position: 200,
                id: None,
                ref_allele: b"AC".to_vec(),
                alt_alleles: vec![b"A".to_vec()],
                quality: None,
                filter: VariantFilter::Missing,
            },
            Variant {
                chrom: "chr2".into(),
                position: 300,
                id: None,
                ref_allele: b"T".to_vec(),
                alt_alleles: vec![b"TA".to_vec(), b"TG".to_vec()],
                quality: Some(50.5),
                filter: VariantFilter::Fail(vec!["LowQual".into()]),
            },
        ]
    }

    fn sample_intervals() -> Vec<GenomicInterval> {
        vec![
            GenomicInterval::with_strand("chr1", 100, 200, Strand::Forward).unwrap(),
            GenomicInterval::with_strand("chr1", 300, 400, Strand::Reverse).unwrap(),
            GenomicInterval::with_strand("chr2", 500, 600, Strand::Unknown).unwrap(),
        ]
    }

    #[test]
    fn variant_roundtrip() {
        let (_tmp, path) = temp_parquet();
        let variants = sample_variants();
        write_variants_parquet(&variants, &path).unwrap();
        let loaded = read_variants_parquet(&path).unwrap();
        assert_eq!(loaded.len(), 3);

        assert_eq!(loaded[0].chrom, "chr1");
        assert_eq!(loaded[0].position, 100);
        assert_eq!(loaded[0].id, Some("rs123".into()));
        assert_eq!(loaded[0].ref_allele, b"A");
        assert_eq!(loaded[0].alt_alleles, vec![b"G".to_vec()]);
        assert_eq!(loaded[0].quality, Some(30.0));
        assert_eq!(loaded[0].filter, VariantFilter::Pass);

        assert_eq!(loaded[1].id, None);
        assert_eq!(loaded[1].quality, None);
        assert_eq!(loaded[1].filter, VariantFilter::Missing);

        assert_eq!(loaded[2].alt_alleles.len(), 2);
        assert_eq!(
            loaded[2].filter,
            VariantFilter::Fail(vec!["LowQual".into()])
        );
    }

    #[test]
    fn interval_roundtrip() {
        let (_tmp, path) = temp_parquet();
        let intervals = sample_intervals();
        write_intervals_parquet(&intervals, &path).unwrap();
        let loaded = read_intervals_parquet(&path).unwrap();
        assert_eq!(loaded.len(), 3);
        assert_eq!(loaded[0].chrom, "chr1");
        assert_eq!(loaded[0].start, 100);
        assert_eq!(loaded[0].end, 200);
        assert_eq!(loaded[0].strand, Strand::Forward);
        assert_eq!(loaded[1].strand, Strand::Reverse);
        assert_eq!(loaded[2].strand, Strand::Unknown);
    }

    #[test]
    fn variant_stats() {
        let (_tmp, path) = temp_parquet();
        write_variants_parquet(&sample_variants(), &path).unwrap();
        let stats = parquet_variant_stats(&path).unwrap();
        assert_eq!(stats.variant_count, 3);
        assert_eq!(stats.snv_count, 1);
        assert_eq!(stats.indel_count, 2);
        assert_eq!(stats.pass_count, 1);
        assert_eq!(stats.chromosomes, vec!["chr1", "chr2"]);
    }

    #[test]
    fn interval_stats() {
        let (_tmp, path) = temp_parquet();
        write_intervals_parquet(&sample_intervals(), &path).unwrap();
        let stats = parquet_interval_stats(&path).unwrap();
        assert_eq!(stats.record_count, 3);
        assert_eq!(stats.total_bases, 300);
        assert_eq!(stats.chromosomes, vec!["chr1", "chr2"]);
    }

    #[test]
    fn parquet_info_metadata() {
        let (_tmp, path) = temp_parquet();
        write_variants_parquet(&sample_variants(), &path).unwrap();
        let info = parquet_info(&path).unwrap();
        assert_eq!(info.num_rows, 3);
        assert_eq!(info.num_columns, 7);
        assert_eq!(info.column_names[0], "chrom");
        assert_eq!(info.column_names[6], "filter");
        assert_eq!(info.num_row_groups, 1);
    }

    #[test]
    fn empty_file() {
        let (_tmp, path) = temp_parquet();
        write_variants_parquet(&[], &path).unwrap();
        let loaded = read_variants_parquet(&path).unwrap();
        assert!(loaded.is_empty());

        let info = parquet_info(&path).unwrap();
        assert_eq!(info.num_rows, 0);
    }

    #[test]
    fn nullable_fields() {
        let (_tmp, path) = temp_parquet();
        let variants = vec![
            Variant {
                chrom: "chr1".into(),
                position: 1,
                id: None,
                ref_allele: b"A".to_vec(),
                alt_alleles: vec![b"G".to_vec()],
                quality: None,
                filter: VariantFilter::Missing,
            },
            Variant {
                chrom: "chr1".into(),
                position: 2,
                id: Some("rs1".into()),
                ref_allele: b"C".to_vec(),
                alt_alleles: vec![b"T".to_vec()],
                quality: Some(99.0),
                filter: VariantFilter::Pass,
            },
        ];
        write_variants_parquet(&variants, &path).unwrap();
        let loaded = read_variants_parquet(&path).unwrap();
        assert_eq!(loaded[0].id, None);
        assert_eq!(loaded[0].quality, None);
        assert_eq!(loaded[1].id, Some("rs1".into()));
        assert_eq!(loaded[1].quality, Some(99.0));
    }

    #[test]
    fn nonexistent_file_error() {
        let result = read_variants_parquet("/nonexistent/file.parquet");
        assert!(result.is_err());
    }
}
