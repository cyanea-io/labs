//! Indexed FASTA reading with random access via `.fai` index files.
//!
//! Provides [`FastaIndex`] for parsing, building, and writing samtools-compatible
//! `.fai` index files, and [`IndexedFastaReader`] for seeking directly to named
//! regions within a FASTA file without scanning the entire file.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::Path;

use cyanea_core::{CyaneaError, Result};

/// A single entry in a FASTA index (.fai file).
#[derive(Debug, Clone)]
pub struct FastaIndexEntry {
    /// Sequence name.
    pub name: String,
    /// Total length of the sequence in bases.
    pub length: usize,
    /// Byte offset of the first base in the FASTA file.
    pub offset: u64,
    /// Number of bases per line.
    pub line_bases: usize,
    /// Number of bytes per line (including newline).
    pub line_width: usize,
}

/// FASTA index parsed from a `.fai` file.
#[derive(Debug, Clone)]
pub struct FastaIndex {
    entries: Vec<FastaIndexEntry>,
    name_to_idx: HashMap<String, usize>,
}

impl FastaIndex {
    /// Parse a `.fai` index file (tab-separated: name, length, offset, line_bases, line_width).
    pub fn from_file(fai_path: &Path) -> Result<Self> {
        let file = File::open(fai_path)?;
        let reader = BufReader::new(file);
        let mut entries = Vec::new();
        let mut name_to_idx = HashMap::new();

        for (line_num, line) in reader.lines().enumerate() {
            let line = line?;
            let line = line.trim_end();
            if line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 {
                return Err(CyaneaError::Parse(format!(
                    "line {}: expected 5 tab-separated fields, got {}",
                    line_num + 1,
                    fields.len()
                )));
            }
            let name = fields[0].to_string();
            let length = fields[1].parse::<usize>().map_err(|e| {
                CyaneaError::Parse(format!("line {}: invalid length: {}", line_num + 1, e))
            })?;
            let offset = fields[2].parse::<u64>().map_err(|e| {
                CyaneaError::Parse(format!("line {}: invalid offset: {}", line_num + 1, e))
            })?;
            let line_bases = fields[3].parse::<usize>().map_err(|e| {
                CyaneaError::Parse(format!("line {}: invalid line_bases: {}", line_num + 1, e))
            })?;
            let line_width = fields[4].parse::<usize>().map_err(|e| {
                CyaneaError::Parse(format!("line {}: invalid line_width: {}", line_num + 1, e))
            })?;

            name_to_idx.insert(name.clone(), entries.len());
            entries.push(FastaIndexEntry {
                name,
                length,
                offset,
                line_bases,
                line_width,
            });
        }

        Ok(Self {
            entries,
            name_to_idx,
        })
    }

    /// Build an index by scanning a FASTA file.
    pub fn build(fasta_path: &Path) -> Result<Self> {
        let file = File::open(fasta_path)?;
        let mut reader = BufReader::new(file);
        let mut entries = Vec::new();
        let mut name_to_idx = HashMap::new();
        let mut byte_pos: u64 = 0;

        // State for the current sequence being scanned.
        let mut current_name: Option<String> = None;
        let mut seq_offset: u64 = 0;
        let mut seq_length: usize = 0;
        let mut line_bases: usize = 0;
        let mut line_width: usize = 0;
        let mut first_seq_line = true;

        let mut line_buf = String::new();
        loop {
            line_buf.clear();
            let bytes_read = reader.read_line(&mut line_buf)?;
            if bytes_read == 0 {
                break;
            }

            if line_buf.starts_with('>') {
                // Flush previous record.
                if let Some(name) = current_name.take() {
                    name_to_idx.insert(name.clone(), entries.len());
                    entries.push(FastaIndexEntry {
                        name,
                        length: seq_length,
                        offset: seq_offset,
                        line_bases,
                        line_width,
                    });
                }

                let header = line_buf[1..].trim_end_matches(&['\n', '\r'][..]);
                let name = header
                    .split_whitespace()
                    .next()
                    .unwrap_or(header)
                    .to_string();
                current_name = Some(name);
                seq_length = 0;
                line_bases = 0;
                line_width = 0;
                first_seq_line = true;
                seq_offset = byte_pos + bytes_read as u64;
            } else if current_name.is_some() {
                let bases = line_buf.trim_end_matches(&['\n', '\r'][..]).len();
                seq_length += bases;
                if first_seq_line {
                    line_bases = bases;
                    line_width = bytes_read;
                    first_seq_line = false;
                }
            }

            byte_pos += bytes_read as u64;
        }

        // Flush final record.
        if let Some(name) = current_name.take() {
            name_to_idx.insert(name.clone(), entries.len());
            entries.push(FastaIndexEntry {
                name,
                length: seq_length,
                offset: seq_offset,
                line_bases,
                line_width,
            });
        }

        Ok(Self {
            entries,
            name_to_idx,
        })
    }

    /// Write the index to a `.fai` file.
    pub fn write(&self, path: &Path) -> Result<()> {
        let mut file = File::create(path)?;
        for entry in &self.entries {
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}",
                entry.name, entry.length, entry.offset, entry.line_bases, entry.line_width
            )?;
        }
        Ok(())
    }

    /// Look up an entry by sequence name.
    pub fn get(&self, name: &str) -> Option<&FastaIndexEntry> {
        self.name_to_idx.get(name).map(|&idx| &self.entries[idx])
    }

    /// List all sequence names in file order.
    pub fn names(&self) -> Vec<&str> {
        self.entries.iter().map(|e| e.name.as_str()).collect()
    }

    /// Number of indexed sequences.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
}

/// Reader for random access into an indexed FASTA file.
pub struct IndexedFastaReader {
    index: FastaIndex,
    reader: BufReader<File>,
}

impl IndexedFastaReader {
    /// Open a FASTA file with a pre-built index.
    pub fn new(fasta_path: &Path, index: FastaIndex) -> Result<Self> {
        let file = File::open(fasta_path)?;
        Ok(Self {
            index,
            reader: BufReader::new(file),
        })
    }

    /// Open a FASTA file and automatically load its `.fai` sidecar index.
    pub fn open(fasta_path: &Path) -> Result<Self> {
        let fai_path = fasta_path.with_extension({
            let mut ext = fasta_path
                .extension()
                .unwrap_or_default()
                .to_os_string();
            ext.push(".fai");
            ext
        });
        let index = FastaIndex::from_file(&fai_path)?;
        Self::new(fasta_path, index)
    }

    /// Fetch bases `[start, end)` from the named sequence.
    pub fn fetch(&mut self, name: &str, start: usize, end: usize) -> Result<Vec<u8>> {
        let entry = self.index.get(name).ok_or_else(|| {
            CyaneaError::InvalidInput(format!("sequence not found: {}", name))
        })?;

        if start > end || end > entry.length {
            return Err(CyaneaError::InvalidInput(format!(
                "region [{}..{}) out of bounds for {} (length {})",
                start, end, name, entry.length
            )));
        }

        let len = end - start;
        if len == 0 {
            return Ok(Vec::new());
        }

        let byte_start = entry.offset
            + (start / entry.line_bases) as u64 * entry.line_width as u64
            + (start % entry.line_bases) as u64;

        self.reader.seek(SeekFrom::Start(byte_start))?;

        let mut result = Vec::with_capacity(len);
        let mut remaining = len;
        let mut buf = [0u8; 8192];

        while remaining > 0 {
            let to_read = remaining.min(buf.len());
            let n = self.reader.read(&mut buf[..to_read])?;
            if n == 0 {
                return Err(CyaneaError::Parse("unexpected end of FASTA file".into()));
            }
            for &b in &buf[..n] {
                if b != b'\n' && b != b'\r' {
                    result.push(b);
                    remaining -= 1;
                    if remaining == 0 {
                        break;
                    }
                }
            }
        }

        Ok(result)
    }

    /// Fetch the entire sequence by name.
    pub fn fetch_all(&mut self, name: &str) -> Result<Vec<u8>> {
        let length = self
            .index
            .get(name)
            .ok_or_else(|| {
                CyaneaError::InvalidInput(format!("sequence not found: {}", name))
            })?
            .length;
        self.fetch(name, 0, length)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    fn write_fasta(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f.flush().unwrap();
        f
    }

    const MULTI_FASTA: &str = "\
>chr1
ACGTACGT
AAAACCCC
>chr2
TTTTGGGG
";

    #[test]
    fn build_index_from_fasta() {
        let f = write_fasta(MULTI_FASTA);
        let idx = FastaIndex::build(f.path()).unwrap();
        assert_eq!(idx.len(), 2);
        assert_eq!(idx.names(), vec!["chr1", "chr2"]);

        let e1 = idx.get("chr1").unwrap();
        assert_eq!(e1.length, 16);
        assert_eq!(e1.line_bases, 8);
        assert_eq!(e1.line_width, 9); // 8 bases + newline

        let e2 = idx.get("chr2").unwrap();
        assert_eq!(e2.length, 8);
    }

    #[test]
    fn write_and_read_roundtrip() {
        let fa = write_fasta(MULTI_FASTA);
        let idx = FastaIndex::build(fa.path()).unwrap();

        let fai_file = NamedTempFile::new().unwrap();
        idx.write(fai_file.path()).unwrap();

        let idx2 = FastaIndex::from_file(fai_file.path()).unwrap();
        assert_eq!(idx2.len(), idx.len());
        for name in idx.names() {
            let a = idx.get(name).unwrap();
            let b = idx2.get(name).unwrap();
            assert_eq!(a.length, b.length);
            assert_eq!(a.offset, b.offset);
            assert_eq!(a.line_bases, b.line_bases);
            assert_eq!(a.line_width, b.line_width);
        }
    }

    #[test]
    fn parse_fai_format() {
        let fai_content = "chr1\t16\t6\t8\t9\nchr2\t8\t31\t8\t9\n";
        let f = write_fasta(fai_content);
        let idx = FastaIndex::from_file(f.path()).unwrap();
        assert_eq!(idx.len(), 2);
        let e = idx.get("chr1").unwrap();
        assert_eq!(e.length, 16);
        assert_eq!(e.offset, 6);
    }

    #[test]
    fn fetch_subsequence() {
        let fa = write_fasta(MULTI_FASTA);
        let idx = FastaIndex::build(fa.path()).unwrap();
        let mut reader = IndexedFastaReader::new(fa.path(), idx).unwrap();

        let sub = reader.fetch("chr1", 2, 6).unwrap();
        assert_eq!(sub, b"GTAC");
    }

    #[test]
    fn fetch_across_line_boundary() {
        let fa = write_fasta(MULTI_FASTA);
        let idx = FastaIndex::build(fa.path()).unwrap();
        let mut reader = IndexedFastaReader::new(fa.path(), idx).unwrap();

        // Spans the newline between the two lines of chr1.
        let sub = reader.fetch("chr1", 6, 10).unwrap();
        assert_eq!(sub, b"GTAA");
    }

    #[test]
    fn fetch_entire_sequence() {
        let fa = write_fasta(MULTI_FASTA);
        let idx = FastaIndex::build(fa.path()).unwrap();
        let mut reader = IndexedFastaReader::new(fa.path(), idx).unwrap();

        let all = reader.fetch_all("chr2").unwrap();
        assert_eq!(all, b"TTTTGGGG");
    }

    #[test]
    fn error_unknown_sequence() {
        let fa = write_fasta(MULTI_FASTA);
        let idx = FastaIndex::build(fa.path()).unwrap();
        let mut reader = IndexedFastaReader::new(fa.path(), idx).unwrap();

        let err = reader.fetch("chrX", 0, 1).unwrap_err();
        assert!(err.to_string().contains("sequence not found"));
    }

    #[test]
    fn error_out_of_bounds() {
        let fa = write_fasta(MULTI_FASTA);
        let idx = FastaIndex::build(fa.path()).unwrap();
        let mut reader = IndexedFastaReader::new(fa.path(), idx).unwrap();

        let err = reader.fetch("chr2", 0, 100).unwrap_err();
        assert!(err.to_string().contains("out of bounds"));
    }

    #[test]
    fn multi_sequence_independent_fetch() {
        let fa = write_fasta(MULTI_FASTA);
        let idx = FastaIndex::build(fa.path()).unwrap();
        let mut reader = IndexedFastaReader::new(fa.path(), idx).unwrap();

        let s1 = reader.fetch_all("chr1").unwrap();
        let s2 = reader.fetch_all("chr2").unwrap();
        assert_eq!(s1, b"ACGTACGTAAAACCCC");
        assert_eq!(s2, b"TTTTGGGG");

        // Re-fetch chr1 after chr2 to confirm seeking works both ways.
        let s1_again = reader.fetch("chr1", 0, 4).unwrap();
        assert_eq!(s1_again, b"ACGT");
    }
}
