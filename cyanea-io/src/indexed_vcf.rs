//! Indexed VCF reader for random-access region queries.
//!
//! Supports tabix-indexed bgzipped VCF files (`.vcf.gz` + `.vcf.gz.tbi`).
//! Uses noodles for index parsing and BGZF seeking, then converts records
//! to our [`Variant`] type.

use std::path::Path;

use cyanea_core::{CyaneaError, Result};
use cyanea_omics::variant::{Variant, VariantFilter};

/// An indexed VCF reader for random-access region queries.
pub struct IndexedVcfReader {
    reader: noodles_vcf::io::IndexedReader<noodles_bgzf::Reader<std::fs::File>>,
    header: noodles_vcf::Header,
}

impl IndexedVcfReader {
    /// Open a bgzipped VCF file with its tabix index.
    ///
    /// Looks for the index at `<vcf_gz_path>.tbi` automatically.
    pub fn open(vcf_gz_path: impl AsRef<Path>) -> Result<Self> {
        let vcf_gz_path = vcf_gz_path.as_ref();

        let reader = noodles_vcf::io::indexed_reader::Builder::default()
            .build_from_path(vcf_gz_path)
            .map_err(|e| {
                CyaneaError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("{}: {}", vcf_gz_path.display(), e),
                ))
            })?;

        Self::from_reader(reader)
    }

    /// Open a bgzipped VCF file with an explicitly specified index path.
    pub fn open_with_index(
        vcf_gz_path: impl AsRef<Path>,
        index_path: impl AsRef<Path>,
    ) -> Result<Self> {
        let vcf_gz_path = vcf_gz_path.as_ref();
        let index_path = index_path.as_ref();

        let index = noodles_tabix::read(index_path).map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("{}: {}", index_path.display(), e),
            ))
        })?;

        let file = std::fs::File::open(vcf_gz_path).map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: {}", vcf_gz_path.display(), e),
            ))
        })?;

        // IndexedReader::new wraps the file in a BGZF reader internally
        let reader = noodles_vcf::io::IndexedReader::new(file, index);

        Self::from_reader(reader)
    }

    fn from_reader(
        mut reader: noodles_vcf::io::IndexedReader<noodles_bgzf::Reader<std::fs::File>>,
    ) -> Result<Self> {
        let header = reader
            .read_header()
            .map_err(|e| CyaneaError::Parse(format!("failed to read VCF header: {e}")))?;

        Ok(Self { reader, header })
    }

    /// Fetch all variants overlapping a genomic region.
    ///
    /// Coordinates are 0-based, half-open `[start, end)`.
    /// VCF internally uses 1-based positions; conversion is handled automatically.
    pub fn fetch(&mut self, chrom: &str, start: u64, end: u64) -> Result<Vec<Variant>> {
        let start_pos =
            noodles_core::Position::try_from((start + 1) as usize).map_err(|_| {
                CyaneaError::InvalidInput(format!("invalid start position: {}", start))
            })?;
        let end_pos = noodles_core::Position::try_from(end as usize).map_err(|_| {
            CyaneaError::InvalidInput(format!("invalid end position: {}", end))
        })?;

        let region = noodles_core::Region::new(chrom, start_pos..=end_pos);

        let query = self
            .reader
            .query(&self.header, &region)
            .map_err(|e| CyaneaError::Parse(format!("VCF query failed: {e}")))?;

        let mut variants = Vec::new();
        for result in query {
            let noodles_rec = result
                .map_err(|e| CyaneaError::Parse(format!("error reading VCF record: {e}")))?;

            let variant = convert_noodles_vcf_record(&noodles_rec, &self.header)?;
            variants.push(variant);
        }

        Ok(variants)
    }
}

/// Convenience function: open a bgzipped VCF file and fetch variants in a region.
pub fn fetch_vcf(
    path: impl AsRef<Path>,
    chrom: &str,
    start: u64,
    end: u64,
) -> Result<Vec<Variant>> {
    let mut reader = IndexedVcfReader::open(path)?;
    reader.fetch(chrom, start, end)
}

/// Convert a noodles VCF record to our Variant type.
fn convert_noodles_vcf_record(
    record: &noodles_vcf::Record,
    header: &noodles_vcf::Header,
) -> Result<Variant> {
    use noodles_vcf::variant::record::AlternateBases as _;
    use noodles_vcf::variant::record::Filters as _;
    use noodles_vcf::variant::record::Ids as _;

    let chrom = record.reference_sequence_name().to_string();

    let position = record
        .variant_start()
        .and_then(|r| r.ok())
        .map(|p| usize::from(p) as u64)
        .unwrap_or(0);

    let ids = record.ids();
    let id = if ids.is_empty() {
        None
    } else {
        let id_strs: Vec<&str> = ids.iter().collect();
        Some(id_strs.join(";"))
    };

    let ref_allele = record.reference_bases().as_bytes().to_vec();

    let alt_alleles: Vec<Vec<u8>> = record
        .alternate_bases()
        .iter()
        .filter_map(|alt_result: std::io::Result<&str>| {
            let s = alt_result.ok()?;
            if s == "." {
                None
            } else {
                Some(s.as_bytes().to_vec())
            }
        })
        .collect();

    let quality = record.quality_score().and_then(|r| r.ok()).map(|q| q as f64);

    let filters = record.filters();
    let filter = if filters.is_empty() {
        VariantFilter::Missing
    } else {
        let filter_names: Vec<String> = filters
            .iter(header)
            .filter_map(|r: std::io::Result<&str>| r.ok().map(String::from))
            .collect();
        if filter_names.len() == 1 && filter_names[0] == "PASS" {
            VariantFilter::Pass
        } else {
            VariantFilter::Fail(filter_names)
        }
    };

    Ok(Variant {
        chrom,
        position,
        id,
        ref_allele,
        alt_alleles,
        quality,
        filter,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Write a bgzipped VCF file and its tabix index.
    fn write_indexed_vcf(
        vcf_content: &str,
    ) -> (tempfile::NamedTempFile, tempfile::NamedTempFile) {
        // Write bgzipped VCF
        let mut vcf_gz_file = tempfile::NamedTempFile::with_suffix(".vcf.gz").unwrap();
        {
            let mut bgzf_writer =
                noodles_bgzf::Writer::new(std::io::BufWriter::new(&mut vcf_gz_file));
            bgzf_writer.write_all(vcf_content.as_bytes()).unwrap();
            bgzf_writer.try_finish().unwrap();
        }
        vcf_gz_file.flush().unwrap();

        // Create tabix index from the bgzipped VCF path
        let index = noodles_vcf::index(vcf_gz_file.path()).unwrap();
        let tbi_file = tempfile::NamedTempFile::with_suffix(".tbi").unwrap();
        noodles_tabix::write(tbi_file.path(), &index).unwrap();

        (vcf_gz_file, tbi_file)
    }

    #[test]
    fn indexed_vcf_fetch_region() {
        let vcf = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30.0\tPASS\t.
chr1\t200\t.\tC\tT\t40.0\tPASS\t.
chr1\t500\t.\tG\tA\t50.0\tPASS\t.
chr2\t100\t.\tT\tC\t25.0\tPASS\t.
";
        let (vcf_gz, tbi) = write_indexed_vcf(vcf);
        let mut reader =
            IndexedVcfReader::open_with_index(vcf_gz.path(), tbi.path()).unwrap();

        let results = reader.fetch("chr1", 50, 250).unwrap();
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].position, 100);
        assert_eq!(results[1].position, 200);
    }

    #[test]
    fn indexed_vcf_empty_region() {
        let vcf = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=248956422>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30.0\tPASS\t.
";
        let (vcf_gz, tbi) = write_indexed_vcf(vcf);
        let mut reader =
            IndexedVcfReader::open_with_index(vcf_gz.path(), tbi.path()).unwrap();

        let results = reader.fetch("chr1", 10000, 20000).unwrap();
        assert!(results.is_empty());
    }

    #[test]
    fn indexed_vcf_fields_preserved() {
        let vcf = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=248956422>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\trs123\tA\tG\t30.0\tPASS\t.
";
        let (vcf_gz, tbi) = write_indexed_vcf(vcf);
        let mut reader =
            IndexedVcfReader::open_with_index(vcf_gz.path(), tbi.path()).unwrap();

        let results = reader.fetch("chr1", 0, 200).unwrap();
        assert_eq!(results.len(), 1);
        let v = &results[0];
        assert_eq!(v.chrom, "chr1");
        assert_eq!(v.position, 100);
        assert_eq!(v.id.as_deref(), Some("rs123"));
        assert_eq!(v.ref_allele, b"A");
        assert_eq!(v.alt_alleles, vec![b"G".to_vec()]);
    }

    #[test]
    fn indexed_vcf_multi_chrom() {
        let vcf = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30.0\tPASS\t.
chr2\t100\t.\tT\tC\t25.0\tPASS\t.
";
        let (vcf_gz, tbi) = write_indexed_vcf(vcf);
        let mut reader =
            IndexedVcfReader::open_with_index(vcf_gz.path(), tbi.path()).unwrap();

        let chr1_results = reader.fetch("chr1", 0, 200).unwrap();
        assert_eq!(chr1_results.len(), 1);
        assert_eq!(chr1_results[0].chrom, "chr1");

        let chr2_results = reader.fetch("chr2", 0, 200).unwrap();
        assert_eq!(chr2_results.len(), 1);
        assert_eq!(chr2_results[0].chrom, "chr2");
    }
}
