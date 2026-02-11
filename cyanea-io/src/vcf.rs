//! VCF (Variant Call Format) parser.
//!
//! Parses VCF files into [`Variant`] records. Supports VCF 4.x with the
//! standard 8 mandatory columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO).
//! INFO and genotype columns are currently skipped.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use cyanea_core::{CyaneaError, Result};
use cyanea_omics::variant::{Variant, VariantFilter};

/// Parse a VCF file and return all variant records.
///
/// Header lines (starting with `#`) are skipped. Each data line is parsed
/// into a [`Variant`] using VCF's 1-based coordinate system.
pub fn parse_vcf(path: impl AsRef<Path>) -> Result<Vec<Variant>> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);

    // Read all data lines (skip headers)
    let mut data_lines: Vec<(usize, String)> = Vec::new();
    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| {
            CyaneaError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: line {}: {}", path.display(), line_num + 1, e),
            ))
        })?;
        let trimmed = line.trim().to_string();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        data_lines.push((line_num + 1, trimmed));
    }

    // Parse data lines (optionally in parallel)
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        data_lines
            .par_iter()
            .map(|(line_num, line)| parse_vcf_line(line, *line_num, path))
            .collect()
    }
    #[cfg(not(feature = "parallel"))]
    data_lines
        .iter()
        .map(|(line_num, line)| parse_vcf_line(line, *line_num, path))
        .collect()
}

/// Summary statistics for a VCF file.
#[derive(Debug, Clone)]
pub struct VcfStats {
    pub variant_count: u64,
    pub snv_count: u64,
    pub indel_count: u64,
    pub pass_count: u64,
    pub chromosomes: Vec<String>,
}

/// Parse a VCF file and return summary statistics without storing all records.
pub fn vcf_stats(path: impl AsRef<Path>) -> Result<VcfStats> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);

    let mut variant_count: u64 = 0;
    let mut snv_count: u64 = 0;
    let mut indel_count: u64 = 0;
    let mut pass_count: u64 = 0;
    let mut chroms = Vec::new();
    let mut seen_chroms = std::collections::HashSet::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(CyaneaError::Io)?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let variant = parse_vcf_line(line, line_num + 1, path)?;
        variant_count += 1;

        if variant.is_snv() {
            snv_count += 1;
        }
        if variant.is_indel() {
            indel_count += 1;
        }
        if variant.filter == VariantFilter::Pass {
            pass_count += 1;
        }
        if seen_chroms.insert(variant.chrom.clone()) {
            chroms.push(variant.chrom);
        }
    }

    Ok(VcfStats {
        variant_count,
        snv_count,
        indel_count,
        pass_count,
        chromosomes: chroms,
    })
}

/// Parse a single VCF data line into a Variant.
fn parse_vcf_line(line: &str, line_num: usize, path: &Path) -> Result<Variant> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        return Err(CyaneaError::Parse(format!(
            "{}: line {}: expected at least 8 tab-separated columns, found {}",
            path.display(),
            line_num,
            fields.len()
        )));
    }

    let chrom = fields[0].to_string();

    let position: u64 = fields[1].parse().map_err(|_| {
        CyaneaError::Parse(format!(
            "{}: line {}: invalid position '{}'",
            path.display(),
            line_num,
            fields[1]
        ))
    })?;

    let id = match fields[2] {
        "." => None,
        s => Some(s.to_string()),
    };

    let ref_allele = fields[3].as_bytes().to_vec();

    let alt_alleles: Vec<Vec<u8>> = fields[4]
        .split(',')
        .filter(|a| *a != ".")
        .map(|a| a.as_bytes().to_vec())
        .collect();

    let quality: Option<f64> = match fields[5] {
        "." => None,
        s => Some(s.parse().map_err(|_| {
            CyaneaError::Parse(format!(
                "{}: line {}: invalid quality '{}'",
                path.display(),
                line_num,
                s
            ))
        })?),
    };

    let filter = match fields[6] {
        "PASS" => VariantFilter::Pass,
        "." => VariantFilter::Missing,
        s => VariantFilter::Fail(s.split(';').map(|f| f.to_string()).collect()),
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
    use tempfile::NamedTempFile;

    fn write_vcf(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(".vcf").unwrap();
        write!(file, "{}", content).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_parse_basic_vcf() {
        let file = write_vcf(
            "##fileformat=VCFv4.3\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t100\trs123\tA\tG\t30.0\tPASS\t.\n\
             chr1\t200\t.\tAC\tA\t.\t.\t.\n\
             chr2\t300\t.\tT\tTA,TG\t50.5\tLowQual\t.\n",
        );

        let variants = parse_vcf(file.path()).unwrap();
        assert_eq!(variants.len(), 3);

        // SNV
        assert_eq!(variants[0].chrom, "chr1");
        assert_eq!(variants[0].position, 100);
        assert_eq!(variants[0].id, Some("rs123".to_string()));
        assert_eq!(variants[0].ref_allele, b"A");
        assert_eq!(variants[0].alt_alleles, vec![b"G".to_vec()]);
        assert_eq!(variants[0].quality, Some(30.0));
        assert_eq!(variants[0].filter, VariantFilter::Pass);
        assert!(variants[0].is_snv());

        // Deletion
        assert_eq!(variants[1].ref_allele, b"AC");
        assert_eq!(variants[1].alt_alleles, vec![b"A".to_vec()]);
        assert!(variants[1].is_indel());
        assert_eq!(variants[1].quality, None);
        assert_eq!(variants[1].filter, VariantFilter::Missing);

        // Multi-allelic
        assert_eq!(variants[2].alt_alleles.len(), 2);
        assert_eq!(variants[2].filter, VariantFilter::Fail(vec!["LowQual".to_string()]));
    }

    #[test]
    fn test_vcf_stats() {
        let file = write_vcf(
            "##fileformat=VCFv4.3\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t100\t.\tA\tG\t30\tPASS\t.\n\
             chr1\t200\t.\tAC\tA\t.\tPASS\t.\n\
             chr2\t300\t.\tT\tC\t50\tLowQual\t.\n",
        );

        let stats = vcf_stats(file.path()).unwrap();
        assert_eq!(stats.variant_count, 3);
        assert_eq!(stats.snv_count, 2);
        assert_eq!(stats.indel_count, 1);
        assert_eq!(stats.pass_count, 2);
        assert_eq!(stats.chromosomes, vec!["chr1", "chr2"]);
    }

    #[test]
    fn test_vcf_empty_file() {
        let file = write_vcf("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        let variants = parse_vcf(file.path()).unwrap();
        assert!(variants.is_empty());
    }

    #[test]
    fn test_vcf_missing_columns() {
        let file = write_vcf("#CHROM\tPOS\n\nchr1\t100\n");
        let result = parse_vcf(file.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_vcf_file_not_found() {
        let result = parse_vcf("/nonexistent/file.vcf");
        assert!(result.is_err());
    }
}
