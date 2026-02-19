//! PHYLIP alignment format parser and writer.
//!
//! Supports both interleaved and sequential PHYLIP formats. The format
//! starts with a dimension line giving the number of taxa and sites,
//! followed by sequence data in either interleaved or sequential layout.
//!
//! Uses the "relaxed" naming convention where names are whitespace-delimited
//! tokens of arbitrary length (not the strict 10-character fixed-width format).

use cyanea_core::{CyaneaError, Result};

/// A parsed PHYLIP multiple sequence alignment.
#[derive(Debug, Clone)]
pub struct PhylipAlignment {
    /// Aligned sequences as `(name, aligned_sequence)` pairs.
    pub sequences: Vec<(String, String)>,
    /// Number of taxa declared in the header.
    pub n_taxa: usize,
    /// Number of sites (alignment columns) declared in the header.
    pub n_sites: usize,
}

/// Parse a PHYLIP interleaved alignment from a string.
///
/// The first line must contain `n_taxa n_sites`. The first block contains
/// sequence names followed by sequence data; subsequent blocks contain only
/// sequence data. Blocks are separated by blank lines.
///
/// # Examples
///
/// ```
/// # use cyanea_io::phylip::parse_phylip;
/// let input = " 2 8\nseq1  ACGT\nseq2  TGCA\n\nACGT\nTGCA\n";
/// let aln = parse_phylip(input).unwrap();
/// assert_eq!(aln.sequences[0].1, "ACGTACGT");
/// ```
pub fn parse_phylip(input: &str) -> Result<PhylipAlignment> {
    let mut lines = input.lines().peekable();

    // Parse dimension line
    let (n_taxa, n_sites) = parse_dimensions(&mut lines)?;

    // Skip blank lines after header
    skip_blank_lines(&mut lines);

    // First block: names + sequences
    let mut seq_names: Vec<String> = Vec::with_capacity(n_taxa);
    let mut seq_data: Vec<String> = Vec::with_capacity(n_taxa);

    for _ in 0..n_taxa {
        let line = lines.next().ok_or_else(|| {
            CyaneaError::Parse("unexpected end of input in first PHYLIP block".to_string())
        })?;
        let line = line.trim();
        if line.is_empty() {
            return Err(CyaneaError::Parse(
                "unexpected blank line in first PHYLIP block".to_string(),
            ));
        }

        let (name, seq) = split_name_seq(line)?;
        seq_names.push(name);
        seq_data.push(seq);
    }

    if seq_names.len() != n_taxa {
        return Err(CyaneaError::Parse(format!(
            "expected {} taxa in first block, found {}",
            n_taxa,
            seq_names.len()
        )));
    }

    // Subsequent blocks: sequences only (interleaved)
    loop {
        // Skip blank lines between blocks
        let had_content = skip_blank_lines(&mut lines);
        if lines.peek().is_none() {
            break;
        }
        if !had_content && lines.peek().is_none() {
            break;
        }

        for i in 0..n_taxa {
            let line = match lines.next() {
                Some(l) => l,
                None => break,
            };
            let trimmed = line.trim();
            if trimmed.is_empty() {
                break;
            }
            let seq: String = trimmed.chars().filter(|c| !c.is_whitespace()).collect();
            seq_data[i].push_str(&seq);
        }
    }

    // Validate site counts
    for (i, seq) in seq_data.iter().enumerate() {
        if seq.len() != n_sites {
            return Err(CyaneaError::Parse(format!(
                "taxon '{}': expected {} sites, found {}",
                seq_names[i], n_sites, seq.len()
            )));
        }
    }

    let sequences: Vec<(String, String)> = seq_names
        .into_iter()
        .zip(seq_data)
        .collect();

    Ok(PhylipAlignment {
        sequences,
        n_taxa,
        n_sites,
    })
}

/// Parse a PHYLIP sequential alignment from a string.
///
/// The first line must contain `n_taxa n_sites`. Each taxon entry begins
/// with a name line followed by one or more sequence lines until the
/// expected number of sites is reached.
///
/// # Examples
///
/// ```
/// # use cyanea_io::phylip::parse_phylip_sequential;
/// let input = " 2 8\nseq1\nACGTACGT\nseq2\nTGCATGCA\n";
/// let aln = parse_phylip_sequential(input).unwrap();
/// assert_eq!(aln.sequences.len(), 2);
/// ```
pub fn parse_phylip_sequential(input: &str) -> Result<PhylipAlignment> {
    let mut lines = input.lines().peekable();

    let (n_taxa, n_sites) = parse_dimensions(&mut lines)?;

    skip_blank_lines(&mut lines);

    let mut sequences: Vec<(String, String)> = Vec::with_capacity(n_taxa);

    for _ in 0..n_taxa {
        // Skip blank lines between taxa
        skip_blank_lines(&mut lines);

        // First line for each taxon: name [optional_seq_fragment]
        let first_line = lines.next().ok_or_else(|| {
            CyaneaError::Parse("unexpected end of input in sequential PHYLIP".to_string())
        })?;
        let first_line = first_line.trim();

        let (name, initial_seq) = split_name_seq(first_line)?;
        let mut seq = initial_seq;

        // Read more lines until we have n_sites characters
        while seq.len() < n_sites {
            let line = lines.next().ok_or_else(|| {
                CyaneaError::Parse(format!(
                    "unexpected end of input for taxon '{}': got {} of {} sites",
                    name,
                    seq.len(),
                    n_sites
                ))
            })?;
            let fragment: String = line.trim().chars().filter(|c| !c.is_whitespace()).collect();
            seq.push_str(&fragment);
        }

        if seq.len() != n_sites {
            return Err(CyaneaError::Parse(format!(
                "taxon '{}': expected {} sites, found {}",
                name, n_sites, seq.len()
            )));
        }

        sequences.push((name, seq));
    }

    if sequences.len() != n_taxa {
        return Err(CyaneaError::Parse(format!(
            "expected {} taxa, found {}",
            n_taxa,
            sequences.len()
        )));
    }

    Ok(PhylipAlignment {
        sequences,
        n_taxa,
        n_sites,
    })
}

/// Write a PHYLIP alignment in interleaved format.
///
/// Uses relaxed format with whitespace-delimited names.
///
/// # Examples
///
/// ```
/// # use cyanea_io::phylip::{PhylipAlignment, write_phylip};
/// let aln = PhylipAlignment {
///     sequences: vec![("s1".to_string(), "ACGT".to_string())],
///     n_taxa: 1,
///     n_sites: 4,
/// };
/// let output = write_phylip(&aln);
/// assert!(output.starts_with(" 1 4"));
/// ```
pub fn write_phylip(aln: &PhylipAlignment) -> String {
    const BLOCK_SIZE: usize = 60;

    let mut out = String::new();
    out.push_str(&format!(" {} {}\n", aln.n_taxa, aln.n_sites));

    if aln.sequences.is_empty() {
        return out;
    }

    let max_name_len = aln.sequences.iter().map(|(n, _)| n.len()).max().unwrap_or(0);
    let pad = max_name_len + 2;

    let total_len = aln.n_sites;
    let mut offset = 0;

    while offset < total_len {
        let end = std::cmp::min(offset + BLOCK_SIZE, total_len);

        for (name, seq) in &aln.sequences {
            if offset == 0 {
                // First block: include names
                let fragment = if end <= seq.len() {
                    &seq[offset..end]
                } else {
                    &seq[offset..]
                };
                out.push_str(&format!("{:<width$}{}\n", name, fragment, width = pad));
            } else {
                // Subsequent blocks: sequences only
                let fragment = if end <= seq.len() {
                    &seq[offset..end]
                } else if offset < seq.len() {
                    &seq[offset..]
                } else {
                    ""
                };
                out.push_str(&format!("{}\n", fragment));
            }
        }

        if offset + BLOCK_SIZE < total_len {
            out.push('\n');
        }

        offset = end;
    }

    out
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn parse_dimensions<'a, I: Iterator<Item = &'a str>>(
    lines: &mut std::iter::Peekable<I>,
) -> Result<(usize, usize)> {
    let dim_line = loop {
        match lines.next() {
            Some(l) => {
                let trimmed = l.trim();
                if !trimmed.is_empty() {
                    break trimmed.to_string();
                }
            }
            None => {
                return Err(CyaneaError::Parse(
                    "empty input: no PHYLIP dimension line".to_string(),
                ))
            }
        }
    };

    let parts: Vec<&str> = dim_line.split_whitespace().collect();
    if parts.len() < 2 {
        return Err(CyaneaError::Parse(format!(
            "invalid PHYLIP dimension line: '{}'",
            dim_line
        )));
    }

    let n_taxa: usize = parts[0].parse().map_err(|_| {
        CyaneaError::Parse(format!("invalid n_taxa '{}' in dimension line", parts[0]))
    })?;

    let n_sites: usize = parts[1].parse().map_err(|_| {
        CyaneaError::Parse(format!("invalid n_sites '{}' in dimension line", parts[1]))
    })?;

    Ok((n_taxa, n_sites))
}

fn skip_blank_lines<'a, I: Iterator<Item = &'a str>>(
    lines: &mut std::iter::Peekable<I>,
) -> bool {
    let mut skipped = false;
    while let Some(line) = lines.peek() {
        if line.trim().is_empty() {
            skipped = true;
            lines.next();
        } else {
            break;
        }
    }
    skipped
}

fn split_name_seq(line: &str) -> Result<(String, String)> {
    let parts: Vec<&str> = line.splitn(2, char::is_whitespace).collect();
    if parts.is_empty() || parts[0].is_empty() {
        return Err(CyaneaError::Parse(format!(
            "could not parse name from line: '{}'",
            line
        )));
    }
    let name = parts[0].to_string();
    let seq = if parts.len() > 1 {
        parts[1].chars().filter(|c| !c.is_whitespace()).collect()
    } else {
        String::new()
    };
    Ok((name, seq))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phylip_round_trip() {
        let aln = PhylipAlignment {
            sequences: vec![
                ("seq1".to_string(), "ACGTACGT".to_string()),
                ("seq2".to_string(), "TGCATGCA".to_string()),
            ],
            n_taxa: 2,
            n_sites: 8,
        };

        let written = write_phylip(&aln);
        let parsed = parse_phylip(&written).unwrap();
        assert_eq!(parsed.n_taxa, 2);
        assert_eq!(parsed.n_sites, 8);
        assert_eq!(parsed.sequences[0].0, "seq1");
        assert_eq!(parsed.sequences[0].1, "ACGTACGT");
        assert_eq!(parsed.sequences[1].0, "seq2");
        assert_eq!(parsed.sequences[1].1, "TGCATGCA");
    }

    #[test]
    fn phylip_sequential_parse() {
        let input = " 2 10\nAlpha\nACGTACGTAC\nBeta\nTGCATGCATG\n";
        let aln = parse_phylip_sequential(input).unwrap();
        assert_eq!(aln.n_taxa, 2);
        assert_eq!(aln.n_sites, 10);
        assert_eq!(aln.sequences[0].0, "Alpha");
        assert_eq!(aln.sequences[0].1, "ACGTACGTAC");
        assert_eq!(aln.sequences[1].0, "Beta");
        assert_eq!(aln.sequences[1].1, "TGCATGCATG");
    }

    #[test]
    fn phylip_interleaved_parse() {
        let input = " 2 8\nseq1  ACGT\nseq2  TGCA\n\nACGT\nTGCA\n";
        let aln = parse_phylip(input).unwrap();
        assert_eq!(aln.n_taxa, 2);
        assert_eq!(aln.n_sites, 8);
        assert_eq!(aln.sequences[0].1, "ACGTACGT");
        assert_eq!(aln.sequences[1].1, "TGCATGCA");
    }

    #[test]
    fn phylip_dimension_mismatch_error() {
        let input = " 2 10\nseq1  ACGT\nseq2  TGCA\n";
        let result = parse_phylip(input);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("expected 10 sites"));
    }

    #[test]
    fn phylip_whitespace_handling() {
        // Extra whitespace in sequences should be stripped
        let input = " 3 8\nalpha     ACGT ACGT\nbeta      TGCA TGCA\ngamma     AAAA CCCC\n";
        let aln = parse_phylip(input).unwrap();
        assert_eq!(aln.n_taxa, 3);
        assert_eq!(aln.sequences[0].1, "ACGTACGT");
        assert_eq!(aln.sequences[1].1, "TGCATGCA");
        assert_eq!(aln.sequences[2].1, "AAAACCCC");
    }
}
