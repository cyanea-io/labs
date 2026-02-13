//! Open Reading Frame (ORF) finder.
//!
//! Scans nucleotide sequences for ORFs across all six reading frames
//! (three forward, three reverse complement). An ORF begins at a start
//! codon (ATG by default) and extends to the next in-frame stop codon
//! (TAA, TAG, TGA) or the end of the sequence.

/// Result of an ORF search.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OrfResult {
    /// Start position in the input sequence (0-indexed).
    pub start: usize,
    /// End position (exclusive) in the input sequence.
    pub end: usize,
    /// Reading frame (0, 1, or 2).
    pub frame: usize,
    /// Strand: `Forward` or `Reverse`.
    pub strand: Strand,
    /// The nucleotide sequence of the ORF.
    pub sequence: Vec<u8>,
}

/// Strand orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

const DEFAULT_STARTS: &[&[u8]] = &[b"ATG"];
const DEFAULT_STOPS: &[&[u8]] = &[b"TAA", b"TAG", b"TGA"];

/// Compute the reverse complement of a DNA sequence.
///
/// Complements each base (A<->T, C<->G) and reverses the result.
/// Non-ACGT characters are left unchanged.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => other,
        })
        .collect()
}

/// Check whether `codon` matches any entry in `codons`.
fn matches_codon(codon: &[u8], codons: &[&[u8]]) -> bool {
    codons.iter().any(|c| c == &codon)
}

/// Find ORFs on the forward strand using configurable start/stop codons.
///
/// An ORF starts at a codon matching `start_codons` and extends to the
/// next in-frame codon matching `stop_codons` (inclusive) or end of
/// sequence. Only ORFs with nucleotide length >= `min_length` are returned.
pub fn find_orfs_with_codons(
    seq: &[u8],
    min_length: usize,
    start_codons: &[&[u8]],
    stop_codons: &[&[u8]],
) -> Vec<OrfResult> {
    let mut results = Vec::new();
    let len = seq.len();

    for frame in 0..3 {
        let mut pos = frame;
        let mut orf_start: Option<usize> = None;

        while pos + 3 <= len {
            let codon = &seq[pos..pos + 3];

            if orf_start.is_none() && matches_codon(codon, start_codons) {
                orf_start = Some(pos);
            } else if orf_start.is_some() && matches_codon(codon, stop_codons) {
                let start = orf_start.unwrap();
                let end = pos + 3;
                if end - start >= min_length {
                    results.push(OrfResult {
                        start,
                        end,
                        frame,
                        strand: Strand::Forward,
                        sequence: seq[start..end].to_vec(),
                    });
                }
                orf_start = None;
            }

            pos += 3;
        }

        // ORF extending to end of sequence (no stop codon found)
        if let Some(start) = orf_start {
            let end = len;
            if end - start >= min_length {
                results.push(OrfResult {
                    start,
                    end,
                    frame,
                    strand: Strand::Forward,
                    sequence: seq[start..end].to_vec(),
                });
            }
        }
    }

    results
}

/// Find ORFs in all three forward reading frames.
///
/// Uses the standard start codon (ATG) and stop codons (TAA, TAG, TGA).
/// `min_length` is in nucleotides. Input should be uppercase DNA.
pub fn find_orfs(seq: &[u8], min_length: usize) -> Vec<OrfResult> {
    find_orfs_with_codons(seq, min_length, DEFAULT_STARTS, DEFAULT_STOPS)
}

/// Find ORFs in all six reading frames (forward + reverse complement).
///
/// For the reverse strand, the reverse complement is computed internally.
/// Coordinates in `OrfResult` for reverse-strand ORFs refer to positions
/// on the original (input) sequence.
pub fn find_orfs_both_strands(seq: &[u8], min_length: usize) -> Vec<OrfResult> {
    let mut results = find_orfs(seq, min_length);

    let rc = reverse_complement(seq);
    let rc_orfs = find_orfs(&rc, min_length);
    let len = seq.len();

    for mut orf in rc_orfs {
        // Map reverse-complement positions back to the original sequence.
        // Position `p` on the RC corresponds to `len - p` on the original.
        let orig_start = len - orf.end;
        let orig_end = len - orf.start;
        orf.start = orig_start;
        orf.end = orig_end;
        orf.strand = Strand::Reverse;
        // sequence stays as the ORF nucleotides on the RC strand
        results.push(orf);
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_orf() {
        // ATGAAATAA: ATG (start) + AAA + TAA (stop) = 9 nt
        let orfs = find_orfs(b"ATGAAATAA", 1);
        assert_eq!(orfs.len(), 1);
        assert_eq!(orfs[0].start, 0);
        assert_eq!(orfs[0].end, 9);
        assert_eq!(orfs[0].frame, 0);
        assert_eq!(orfs[0].strand, Strand::Forward);
        assert_eq!(orfs[0].sequence, b"ATGAAATAA");
    }

    #[test]
    fn multiple_orfs_different_frames() {
        // Frame 0: ATG CCC TAA (pos 0..9)
        // Frame 1: xATG AAA TGA x (pos 1..10)
        let seq = b"AATGAAATGAX";
        let orfs = find_orfs(seq, 1);
        // Frame 0 starting at 0: no ATG at pos 0 (AAT), pos 3 (GAA), pos 6 (ATG) -> ATG at 6, then no stop -> extends to end
        // Frame 1 starting at 1: pos 1 (ATG) -> start, pos 4 (AAA), pos 7 (TGA) -> stop => ORF 1..10
        // Frame 2 starting at 2: pos 2 (TGA) no start, pos 5 (AAT), pos 8 (GAX) no start
        // Also check frame 0: pos 0=AAT, pos 3=GAA, pos 6=ATG -> start at 6, no stop -> ORF 6..11
        assert!(orfs.len() >= 2);
        let frame1: Vec<_> = orfs.iter().filter(|o| o.frame == 1).collect();
        assert_eq!(frame1.len(), 1);
        assert_eq!(frame1[0].start, 1);
        assert_eq!(frame1[0].end, 10);
    }

    #[test]
    fn no_orfs_no_atg() {
        let orfs = find_orfs(b"CCCGGGTTTAAA", 1);
        assert!(orfs.is_empty());
    }

    #[test]
    fn orf_extending_to_end() {
        // ATG followed by non-stop codons, no stop before end
        let orfs = find_orfs(b"ATGAAACCC", 1);
        assert_eq!(orfs.len(), 1);
        assert_eq!(orfs[0].start, 0);
        assert_eq!(orfs[0].end, 9);
        assert_eq!(orfs[0].sequence, b"ATGAAACCC");
    }

    #[test]
    fn overlapping_orfs_same_frame() {
        // Two ORFs in frame 0: ATG AAA TAA ATG CCC TAG
        let seq = b"ATGAAATAAATGCCCTAG";
        let orfs: Vec<_> = find_orfs(seq, 1)
            .into_iter()
            .filter(|o| o.frame == 0)
            .collect();
        assert_eq!(orfs.len(), 2);
        assert_eq!(orfs[0].start, 0);
        assert_eq!(orfs[0].end, 9);
        assert_eq!(orfs[1].start, 9);
        assert_eq!(orfs[1].end, 18);
    }

    #[test]
    fn min_length_filtering() {
        // ORF is 9 nt (ATG AAA TAA)
        let orfs = find_orfs(b"ATGAAATAA", 10);
        assert!(orfs.is_empty());

        let orfs = find_orfs(b"ATGAAATAA", 9);
        assert_eq!(orfs.len(), 1);
    }

    #[test]
    fn both_strands() {
        // Forward: ATG AAA TAA (ORF at 0..9)
        // Reverse complement of ATGAAATAA is TTATTTTCAT — no ATG on RC
        // Use a sequence that has ORFs on both strands.
        // Forward: ...ATG CCC TAA...
        // If we also want a reverse-strand ORF, the RC must contain ATG...stop.
        // RC of TTA GGG CAT = ATG CCC TAA → so the original is TTA GGG CAT
        // Combined: ATG CCC TAA TTA GGG CAT
        let seq = b"ATGCCCTAATTAGGGCAT";
        let orfs = find_orfs_both_strands(seq, 1);
        let fwd: Vec<_> = orfs.iter().filter(|o| o.strand == Strand::Forward).collect();
        let rev: Vec<_> = orfs.iter().filter(|o| o.strand == Strand::Reverse).collect();
        assert!(!fwd.is_empty(), "expected forward ORFs");
        assert!(!rev.is_empty(), "expected reverse ORFs");
    }

    #[test]
    fn reverse_strand_coordinates() {
        // Sequence: TTAGGGCAT (len 9)
        // RC:       ATGCCCTAA — ORF at RC positions 0..9
        // Mapped back: start = 9-9 = 0, end = 9-0 = 9
        let seq = b"TTAGGGCAT";
        let orfs = find_orfs_both_strands(seq, 1);
        let rev: Vec<_> = orfs.iter().filter(|o| o.strand == Strand::Reverse).collect();
        assert_eq!(rev.len(), 1);
        assert_eq!(rev[0].start, 0);
        assert_eq!(rev[0].end, 9);
        assert_eq!(rev[0].sequence, b"ATGCCCTAA");
    }

    #[test]
    fn configurable_start_stop_codons() {
        // Use GTG as an alternative start codon
        let seq = b"GTGAAATAA";
        let orfs = find_orfs(seq, 1);
        assert!(orfs.is_empty(), "standard ATG-only should miss GTG start");

        let orfs = find_orfs_with_codons(seq, 1, &[b"ATG", b"GTG"], DEFAULT_STOPS);
        assert_eq!(orfs.len(), 1);
        assert_eq!(orfs[0].start, 0);
        assert_eq!(orfs[0].end, 9);
    }

    #[test]
    fn empty_sequence() {
        let orfs = find_orfs(b"", 1);
        assert!(orfs.is_empty());
    }

    #[test]
    fn sequence_shorter_than_codon() {
        let orfs = find_orfs(b"AT", 1);
        assert!(orfs.is_empty());

        let orfs = find_orfs(b"A", 1);
        assert!(orfs.is_empty());
    }

    #[test]
    fn reverse_complement_basic() {
        assert_eq!(reverse_complement(b"ATGC"), b"GCAT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b""), b"");
    }
}
