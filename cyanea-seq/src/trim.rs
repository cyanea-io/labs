//! Quality trimming, adapter removal, and read filtering for FASTQ records.
//!
//! Two-level API:
//!
//! 1. **Low-level functions** operate on `&[u8]` quality/sequence slices and return
//!    [`TrimRange`] values that can be composed via [`intersect_ranges`].
//! 2. **High-level functions** operate on [`FastqRecord`] and return
//!    `Option<FastqRecord>` (None = filtered out).
//! 3. **[`TrimPipeline`]** builder chains operations in Trimmomatic-style order
//!    and collects statistics via [`TrimReport`].
//!
//! # Example
//!
//! ```
//! use cyanea_seq::trim::{TrimPipeline, adapters};
//! use cyanea_seq::{FastqRecord, DnaSequence, QualityScores};
//!
//! let pipeline = TrimPipeline::new()
//!     .adapter(adapters::TRUSEQ_PREFIX)
//!     .leading(3)
//!     .trailing(3)
//!     .sliding_window(4, 15.0)
//!     .min_length(4);
//!
//! let seq = DnaSequence::new(b"ACGTACGTACGTACGT").unwrap();
//! let qual = QualityScores::from_raw(vec![30; 16]);
//! let record = FastqRecord::new("read1".into(), None, seq, qual).unwrap();
//!
//! let result = pipeline.process(&record);
//! assert!(result.is_some());
//! ```

use crate::fastq::FastqRecord;
use crate::quality::QualityScores;
use crate::seq::ValidatedSeq;
use crate::alphabet::DnaAlphabet;
use cyanea_core::{Annotated, Sequence};

/// A half-open range `[start, end)` describing which portion of a read to keep.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TrimRange {
    pub start: usize,
    pub end: usize,
}

impl TrimRange {
    /// Length of the range (0 if empty).
    pub fn len(&self) -> usize {
        if self.end > self.start {
            self.end - self.start
        } else {
            0
        }
    }

    /// Whether this range is empty.
    pub fn is_empty(&self) -> bool {
        self.end <= self.start
    }
}

/// Intersect multiple trim ranges, returning the overlap.
///
/// Returns a range covering only positions present in *all* input ranges.
/// If `ranges` is empty, returns an empty range.
pub fn intersect_ranges(ranges: &[TrimRange]) -> TrimRange {
    if ranges.is_empty() {
        return TrimRange { start: 0, end: 0 };
    }
    let start = ranges.iter().map(|r| r.start).max().unwrap();
    let end = ranges.iter().map(|r| r.end).min().unwrap();
    TrimRange {
        start,
        end: end.max(start),
    }
}

// ---------------------------------------------------------------------------
// Low-level trim functions
// ---------------------------------------------------------------------------

/// Trimmomatic-style sliding window trim.
///
/// Scans from the 5' end with a window of `window_size` bases. When the mean
/// quality in the window drops below `threshold`, the read is cut at the start
/// of that window. Returns the range of bases to keep.
pub fn trim_sliding_window(quality: &[u8], window_size: usize, threshold: f64) -> TrimRange {
    let len = quality.len();
    if len == 0 || window_size == 0 {
        return TrimRange { start: 0, end: 0 };
    }

    let ws = window_size.min(len);
    let threshold_sum = threshold * ws as f64;

    // Initial window sum
    let mut sum: u64 = quality[..ws].iter().map(|&q| q as u64).sum();
    if (sum as f64) < threshold_sum {
        return TrimRange { start: 0, end: 0 };
    }

    for i in 1..=(len - ws) {
        // Slide: remove left, add right
        sum -= quality[i - 1] as u64;
        sum += quality[i + ws - 1] as u64;
        if (sum as f64) < threshold_sum {
            return TrimRange { start: 0, end: i + ws - 1 };
        }
    }

    TrimRange { start: 0, end: len }
}

/// Trim low-quality bases from the 5' (leading) end.
///
/// Removes consecutive bases from the start that have quality below `threshold`.
pub fn trim_leading(quality: &[u8], threshold: u8) -> TrimRange {
    let start = quality.iter().position(|&q| q >= threshold).unwrap_or(quality.len());
    TrimRange {
        start,
        end: quality.len(),
    }
}

/// Trim low-quality bases from the 3' (trailing) end.
///
/// Removes consecutive bases from the end that have quality below `threshold`.
pub fn trim_trailing(quality: &[u8], threshold: u8) -> TrimRange {
    let end = quality
        .iter()
        .rposition(|&q| q >= threshold)
        .map(|i| i + 1)
        .unwrap_or(0);
    TrimRange { start: 0, end }
}

/// BWA-style quality trimming from the 3' end.
///
/// Scans right-to-left, accumulating `(threshold - Q[i])`. Resets the running
/// sum to 0 when it goes negative. Cuts at the position where the sum was
/// maximized. This effectively finds the longest suffix with low aggregate
/// quality and removes it.
pub fn trim_quality_3prime(quality: &[u8], threshold: u8) -> TrimRange {
    let len = quality.len();
    if len == 0 {
        return TrimRange { start: 0, end: 0 };
    }

    let mut max_sum: i64 = 0;
    let mut sum: i64 = 0;
    let mut cut_pos = len;

    for i in (0..len).rev() {
        sum += threshold as i64 - quality[i] as i64;
        if sum > max_sum {
            max_sum = sum;
            cut_pos = i;
        }
    }

    TrimRange {
        start: 0,
        end: cut_pos,
    }
}

/// Find the position of a 3' adapter in a sequence.
///
/// Checks overlaps at the 3' end of the read, longest first, requiring at most
/// `max_mismatches` mismatches and a minimum overlap of `max(8, adapter.len()/3)`.
/// Returns the position where the adapter starts (i.e., where to cut).
/// If no adapter is found, returns `seq.len()`.
pub fn find_adapter_3prime(seq: &[u8], adapter: &[u8], max_mismatches: usize) -> usize {
    let slen = seq.len();
    let alen = adapter.len();
    if slen == 0 || alen == 0 {
        return slen;
    }

    let min_overlap = 8.max(alen / 3);

    // Check overlaps longest-first: adapter starts at position `start` in the read
    // Overlap length = min(slen - start, alen)
    for start in 0..slen {
        let overlap = (slen - start).min(alen);
        if overlap < min_overlap {
            break;
        }

        let mismatches = seq[start..start + overlap]
            .iter()
            .zip(&adapter[..overlap])
            .filter(|(&a, &b)| a != b)
            .count();

        if mismatches <= max_mismatches {
            return start;
        }
    }

    slen
}

/// Shannon entropy of base composition (bits, max 2.0 for 4-letter DNA).
///
/// Computes `-sum(p * log2(p))` over the frequencies of A, C, G, T.
/// Bases not in {A, C, G, T} are ignored.
pub fn shannon_entropy(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }

    let mut counts = [0u64; 4]; // A, C, G, T
    for &b in seq {
        match b {
            b'A' | b'a' => counts[0] += 1,
            b'C' | b'c' => counts[1] += 1,
            b'G' | b'g' => counts[2] += 1,
            b'T' | b't' => counts[3] += 1,
            _ => {}
        }
    }

    let total: u64 = counts.iter().sum();
    if total == 0 {
        return 0.0;
    }

    let mut entropy = 0.0;
    for &c in &counts {
        if c > 0 {
            let p = c as f64 / total as f64;
            entropy -= p * p.log2();
        }
    }
    entropy
}

// ---------------------------------------------------------------------------
// High-level functions
// ---------------------------------------------------------------------------

/// Apply a trim range to a FASTQ record, producing a new trimmed record.
///
/// Returns `None` if the range is empty (nothing left after trimming).
/// Uses `from_validated()` internally since subslicing validated data is safe.
pub fn apply_trim(record: &FastqRecord, range: TrimRange) -> Option<FastqRecord> {
    if range.is_empty() || range.start >= record.sequence().len() {
        return None;
    }
    let end = range.end.min(record.sequence().len());
    if end <= range.start {
        return None;
    }

    let seq_bytes = record.sequence().as_bytes()[range.start..end].to_vec();
    let qual_bytes = record.quality().as_slice()[range.start..end].to_vec();

    let sequence = ValidatedSeq::<DnaAlphabet>::from_validated(seq_bytes);
    let quality = QualityScores::from_raw(qual_bytes);

    // FastqRecord::new checks length match, but we guarantee it here.
    FastqRecord::new(
        record.name().to_string(),
        record.description().map(|d| d.to_string()),
        sequence,
        quality,
    )
    .ok()
}

/// Remove a 3' adapter from a record.
///
/// Returns a new record with the adapter (and everything after it) removed.
/// If no adapter is found, returns a clone of the original record.
pub fn trim_adapter(record: &FastqRecord, adapter: &[u8], max_mismatches: usize) -> FastqRecord {
    let cut = find_adapter_3prime(record.sequence().as_bytes(), adapter, max_mismatches);
    if cut >= record.sequence().len() {
        return record.clone();
    }
    let range = TrimRange { start: 0, end: cut };
    apply_trim(record, range).unwrap_or_else(|| record.clone())
}

/// Filter a record by length.
///
/// Returns `None` if the record's length is outside `[min_len, max_len]`.
pub fn filter_by_length<'a>(record: &'a FastqRecord, min_len: usize, max_len: usize) -> Option<&'a FastqRecord> {
    let len = record.sequence().len();
    if len >= min_len && len <= max_len {
        Some(record)
    } else {
        None
    }
}

/// Filter a record by low complexity.
///
/// Returns `None` if the Shannon entropy is below `min_entropy`.
pub fn filter_low_complexity<'a>(record: &'a FastqRecord, min_entropy: f64) -> Option<&'a FastqRecord> {
    if shannon_entropy(record.sequence().as_bytes()) >= min_entropy {
        Some(record)
    } else {
        None
    }
}

/// Filter a record by mean quality score.
///
/// Returns `None` if the mean quality is below `min_quality`.
pub fn filter_by_quality<'a>(record: &'a FastqRecord, min_quality: f64) -> Option<&'a FastqRecord> {
    if record.quality().mean() >= min_quality {
        Some(record)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Adapter constants
// ---------------------------------------------------------------------------

/// Common sequencing adapter sequences.
pub mod adapters {
    /// Illumina TruSeq Universal Adapter.
    pub const TRUSEQ_UNIVERSAL: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    /// Illumina TruSeq Indexed Adapter.
    pub const TRUSEQ_INDEXED: &[u8] = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
    /// Nextera Transposase Read 1.
    pub const NEXTERA_READ1: &[u8] = b"TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG";
    /// Nextera Transposase Read 2.
    pub const NEXTERA_READ2: &[u8] = b"GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG";
    /// Illumina Small RNA 3' Adapter.
    pub const SMALL_RNA_3P: &[u8] = b"TGGAATTCTCGGGTGCCAAGG";
    /// Common 12-base prefix shared by TruSeq adapters.
    pub const TRUSEQ_PREFIX: &[u8] = b"AGATCGGAAGAG";
    /// All standard Illumina adapters for batch searching.
    pub const ALL_ILLUMINA: &[&[u8]] = &[
        TRUSEQ_UNIVERSAL,
        TRUSEQ_INDEXED,
        NEXTERA_READ1,
        NEXTERA_READ2,
        SMALL_RNA_3P,
    ];
}

// ---------------------------------------------------------------------------
// TrimPipeline
// ---------------------------------------------------------------------------

/// Which quality trimming algorithm to use in the pipeline.
#[derive(Debug, Clone)]
enum QualityTrimAlgo {
    None,
    SlidingWindow { window_size: usize, threshold: f64 },
    Bwa { threshold: u8 },
}

/// A configurable read-processing pipeline.
///
/// Operations are applied in a fixed order matching Trimmomatic convention:
/// 1. Adapter trimming
/// 2. Leading quality trim
/// 3. Trailing quality trim
/// 4. Sliding window / BWA quality trim
/// 5. Length filter
/// 6. Quality filter
/// 7. Complexity filter
///
/// # Example
///
/// ```
/// use cyanea_seq::trim::TrimPipeline;
///
/// let pipeline = TrimPipeline::new()
///     .leading(3)
///     .trailing(3)
///     .sliding_window(4, 15.0)
///     .min_length(36)
///     .min_mean_quality(20.0);
/// ```
#[derive(Debug, Clone)]
pub struct TrimPipeline {
    adapters: Vec<Vec<u8>>,
    adapter_max_mismatches: usize,
    leading_threshold: Option<u8>,
    trailing_threshold: Option<u8>,
    quality_trim: QualityTrimAlgo,
    min_length: Option<usize>,
    max_length: Option<usize>,
    min_mean_quality: Option<f64>,
    min_entropy: Option<f64>,
}

impl TrimPipeline {
    /// Create a new empty pipeline (no-op by default).
    pub fn new() -> Self {
        Self {
            adapters: Vec::new(),
            adapter_max_mismatches: 1,
            leading_threshold: None,
            trailing_threshold: None,
            quality_trim: QualityTrimAlgo::None,
            min_length: None,
            max_length: None,
            min_mean_quality: None,
            min_entropy: None,
        }
    }

    /// Add a single adapter sequence to search for and remove.
    pub fn adapter(mut self, adapter: &[u8]) -> Self {
        self.adapters.push(adapter.to_vec());
        self
    }

    /// Add all standard Illumina adapters.
    pub fn illumina_adapters(mut self) -> Self {
        for &a in adapters::ALL_ILLUMINA {
            self.adapters.push(a.to_vec());
        }
        self
    }

    /// Set the maximum number of mismatches allowed for adapter matching.
    pub fn adapter_mismatches(mut self, max: usize) -> Self {
        self.adapter_max_mismatches = max;
        self
    }

    /// Trim bases below `threshold` quality from the 5' end.
    pub fn leading(mut self, threshold: u8) -> Self {
        self.leading_threshold = Some(threshold);
        self
    }

    /// Trim bases below `threshold` quality from the 3' end.
    pub fn trailing(mut self, threshold: u8) -> Self {
        self.trailing_threshold = Some(threshold);
        self
    }

    /// Use Trimmomatic-style sliding window trimming.
    ///
    /// Mutually exclusive with [`bwa_quality`](Self::bwa_quality) — setting
    /// one clears the other.
    pub fn sliding_window(mut self, window_size: usize, threshold: f64) -> Self {
        self.quality_trim = QualityTrimAlgo::SlidingWindow {
            window_size,
            threshold,
        };
        self
    }

    /// Use BWA-style 3' quality trimming.
    ///
    /// Mutually exclusive with [`sliding_window`](Self::sliding_window) — setting
    /// one clears the other.
    pub fn bwa_quality(mut self, threshold: u8) -> Self {
        self.quality_trim = QualityTrimAlgo::Bwa { threshold };
        self
    }

    /// Set the minimum read length (reads shorter than this are discarded).
    pub fn min_length(mut self, len: usize) -> Self {
        self.min_length = Some(len);
        self
    }

    /// Set the maximum read length (reads longer than this are discarded).
    pub fn max_length(mut self, len: usize) -> Self {
        self.max_length = Some(len);
        self
    }

    /// Set the minimum mean quality (reads below this are discarded).
    pub fn min_mean_quality(mut self, quality: f64) -> Self {
        self.min_mean_quality = Some(quality);
        self
    }

    /// Set the minimum Shannon entropy (reads below this are discarded).
    pub fn min_entropy(mut self, entropy: f64) -> Self {
        self.min_entropy = Some(entropy);
        self
    }

    /// Process a single record through the pipeline.
    ///
    /// Returns `Some(trimmed_record)` if the record passes all filters,
    /// or `None` if it was filtered out.
    pub fn process(&self, record: &FastqRecord) -> Option<FastqRecord> {
        let mut current = record.clone();

        // 1. Adapter trimming — try each adapter, take the earliest hit
        if !self.adapters.is_empty() {
            let seq = current.sequence().as_bytes();
            let mut best_cut = seq.len();
            for adapter in &self.adapters {
                let cut = find_adapter_3prime(seq, adapter, self.adapter_max_mismatches);
                if cut < best_cut {
                    best_cut = cut;
                }
            }
            if best_cut < seq.len() {
                let range = TrimRange { start: 0, end: best_cut };
                current = apply_trim(&current, range)?;
            }
        }

        // 2-4. Quality trimming via TrimRange composition
        let quality = current.quality().as_slice();
        let mut ranges = Vec::new();

        if let Some(threshold) = self.leading_threshold {
            ranges.push(trim_leading(quality, threshold));
        }

        if let Some(threshold) = self.trailing_threshold {
            ranges.push(trim_trailing(quality, threshold));
        }

        match &self.quality_trim {
            QualityTrimAlgo::SlidingWindow {
                window_size,
                threshold,
            } => {
                ranges.push(trim_sliding_window(quality, *window_size, *threshold));
            }
            QualityTrimAlgo::Bwa { threshold } => {
                ranges.push(trim_quality_3prime(quality, *threshold));
            }
            QualityTrimAlgo::None => {}
        }

        if !ranges.is_empty() {
            // Start with the full range, then intersect with each trim result
            let full = TrimRange {
                start: 0,
                end: quality.len(),
            };
            ranges.insert(0, full);
            let combined = intersect_ranges(&ranges);
            current = apply_trim(&current, combined)?;
        }

        // 5. Length filter
        let len = current.sequence().len();
        if let Some(min) = self.min_length {
            if len < min {
                return None;
            }
        }
        if let Some(max) = self.max_length {
            if len > max {
                return None;
            }
        }

        // 6. Quality filter
        if let Some(min_q) = self.min_mean_quality {
            if current.quality().mean() < min_q {
                return None;
            }
        }

        // 7. Complexity filter
        if let Some(min_e) = self.min_entropy {
            if shannon_entropy(current.sequence().as_bytes()) < min_e {
                return None;
            }
        }

        Some(current)
    }

    /// Process a batch of records, returning only those that pass.
    pub fn process_batch(&self, records: &[FastqRecord]) -> Vec<FastqRecord> {
        records.iter().filter_map(|r| self.process(r)).collect()
    }

    /// Process a batch and collect detailed statistics.
    pub fn process_batch_with_stats(&self, records: &[FastqRecord]) -> TrimReport {
        let total_input = records.len();
        let mut total_bases_input: u64 = 0;
        let mut total_bases_output: u64 = 0;
        let mut filtered_by_length: usize = 0;
        let mut filtered_by_quality: usize = 0;
        let mut filtered_by_complexity: usize = 0;
        let mut adapters_found: usize = 0;
        let mut kept = Vec::new();

        for record in records {
            total_bases_input += record.sequence().len() as u64;

            // Track adapter detection
            if !self.adapters.is_empty() {
                let seq = record.sequence().as_bytes();
                for adapter in &self.adapters {
                    if find_adapter_3prime(seq, adapter, self.adapter_max_mismatches) < seq.len() {
                        adapters_found += 1;
                        break;
                    }
                }
            }

            // Run the full pipeline and track which filter rejected the record
            match self.process(record) {
                Some(trimmed) => {
                    total_bases_output += trimmed.sequence().len() as u64;
                    kept.push(trimmed);
                }
                None => {
                    // Determine which filter rejected it by running steps incrementally
                    let rejection = self.find_rejection_reason(record);
                    match rejection {
                        Rejection::Length => filtered_by_length += 1,
                        Rejection::Quality => filtered_by_quality += 1,
                        Rejection::Complexity => filtered_by_complexity += 1,
                        Rejection::Trimmed => filtered_by_length += 1,
                    }
                }
            }
        }

        TrimReport {
            kept,
            total_input,
            total_output: total_input - filtered_by_length - filtered_by_quality - filtered_by_complexity,
            filtered_by_length,
            filtered_by_quality,
            filtered_by_complexity,
            adapters_found,
            total_bases_input,
            total_bases_output,
        }
    }
}

impl Default for TrimPipeline {
    fn default() -> Self {
        Self::new()
    }
}

enum Rejection {
    Trimmed,
    Length,
    Quality,
    Complexity,
}

impl TrimPipeline {
    fn find_rejection_reason(&self, record: &FastqRecord) -> Rejection {
        let mut current = record.clone();

        // Adapter trimming
        if !self.adapters.is_empty() {
            let seq = current.sequence().as_bytes();
            let mut best_cut = seq.len();
            for adapter in &self.adapters {
                let cut = find_adapter_3prime(seq, adapter, self.adapter_max_mismatches);
                if cut < best_cut {
                    best_cut = cut;
                }
            }
            if best_cut < seq.len() {
                let range = TrimRange { start: 0, end: best_cut };
                match apply_trim(&current, range) {
                    Some(t) => current = t,
                    None => return Rejection::Trimmed,
                }
            }
        }

        // Quality trimming
        let quality = current.quality().as_slice();
        let mut ranges = Vec::new();

        if let Some(threshold) = self.leading_threshold {
            ranges.push(trim_leading(quality, threshold));
        }
        if let Some(threshold) = self.trailing_threshold {
            ranges.push(trim_trailing(quality, threshold));
        }
        match &self.quality_trim {
            QualityTrimAlgo::SlidingWindow { window_size, threshold } => {
                ranges.push(trim_sliding_window(quality, *window_size, *threshold));
            }
            QualityTrimAlgo::Bwa { threshold } => {
                ranges.push(trim_quality_3prime(quality, *threshold));
            }
            QualityTrimAlgo::None => {}
        }
        if !ranges.is_empty() {
            let full = TrimRange { start: 0, end: quality.len() };
            ranges.insert(0, full);
            let combined = intersect_ranges(&ranges);
            match apply_trim(&current, combined) {
                Some(t) => current = t,
                None => return Rejection::Trimmed,
            }
        }

        // Length filter
        let len = current.sequence().len();
        if let Some(min) = self.min_length {
            if len < min { return Rejection::Length; }
        }
        if let Some(max) = self.max_length {
            if len > max { return Rejection::Length; }
        }

        // Quality filter
        if let Some(min_q) = self.min_mean_quality {
            if current.quality().mean() < min_q { return Rejection::Quality; }
        }

        // Complexity filter
        if let Some(min_e) = self.min_entropy {
            if shannon_entropy(current.sequence().as_bytes()) < min_e {
                return Rejection::Complexity;
            }
        }

        // Shouldn't reach here if process() returned None, but default to trimmed
        Rejection::Trimmed
    }
}

/// Summary statistics from [`TrimPipeline::process_batch_with_stats`].
#[derive(Debug, Clone)]
pub struct TrimReport {
    /// Records that passed all filters.
    pub kept: Vec<FastqRecord>,
    /// Total number of input records.
    pub total_input: usize,
    /// Total number of output records.
    pub total_output: usize,
    /// Records filtered out for being too short or too long.
    pub filtered_by_length: usize,
    /// Records filtered out for low mean quality.
    pub filtered_by_quality: usize,
    /// Records filtered out for low complexity.
    pub filtered_by_complexity: usize,
    /// Number of records where an adapter was detected.
    pub adapters_found: usize,
    /// Total bases in input records.
    pub total_bases_input: u64,
    /// Total bases in output records.
    pub total_bases_output: u64,
}

// ---------------------------------------------------------------------------
// Paired-end trimming
// ---------------------------------------------------------------------------

#[cfg(feature = "std")]
use crate::paired::PairedFastqRecord;

/// How to handle orphan reads where one mate passes and the other doesn't.
#[cfg(feature = "std")]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrphanPolicy {
    /// Drop both reads if either fails.
    DropBoth,
    /// Keep R1 if it passes, even if R2 fails.
    KeepFirst,
    /// Keep R2 if it passes, even if R1 fails.
    KeepSecond,
}

/// Result of processing a single read pair through a trim pipeline.
#[cfg(feature = "std")]
#[derive(Debug, Clone)]
pub enum PairedTrimResult {
    /// Both reads passed all filters.
    BothPassed(FastqRecord, FastqRecord),
    /// Only R1 passed.
    OnlyFirst(FastqRecord),
    /// Only R2 passed.
    OnlySecond(FastqRecord),
    /// Both reads were filtered out.
    Dropped,
}

/// Summary statistics from paired-end trim processing.
#[cfg(feature = "std")]
#[derive(Debug, Clone)]
pub struct PairedTrimReport {
    /// Pairs where both reads passed all filters.
    pub kept: Vec<PairedFastqRecord>,
    /// Total number of input pairs.
    pub total_input: usize,
    /// Pairs where both reads passed.
    pub both_passed: usize,
    /// Pairs where only R1 passed.
    pub r1_only_passed: usize,
    /// Pairs where only R2 passed.
    pub r2_only_passed: usize,
    /// Pairs where both reads failed.
    pub both_failed: usize,
    /// Total bases across all input reads (R1 + R2).
    pub total_bases_input: u64,
    /// Total bases across kept output reads (R1 + R2).
    pub total_bases_output: u64,
}

#[cfg(feature = "std")]
impl PairedTrimReport {
    /// Number of orphan reads (one mate passed, the other didn't).
    pub fn orphans(&self) -> usize {
        self.r1_only_passed + self.r2_only_passed
    }

    /// Fraction of input pairs where both reads survived.
    pub fn survival_rate(&self) -> f64 {
        if self.total_input == 0 {
            return 0.0;
        }
        self.both_passed as f64 / self.total_input as f64
    }
}

#[cfg(feature = "std")]
impl TrimPipeline {
    /// Process a single read pair through the pipeline.
    ///
    /// Applies [`process`](Self::process) to each read independently, then
    /// applies the orphan policy to decide what to keep.
    pub fn process_paired(
        &self,
        r1: &FastqRecord,
        r2: &FastqRecord,
        policy: OrphanPolicy,
    ) -> PairedTrimResult {
        let r1_result = self.process(r1);
        let r2_result = self.process(r2);

        match (r1_result, r2_result) {
            (Some(r1), Some(r2)) => PairedTrimResult::BothPassed(r1, r2),
            (Some(r1), None) => match policy {
                OrphanPolicy::KeepFirst => PairedTrimResult::OnlyFirst(r1),
                _ => PairedTrimResult::Dropped,
            },
            (None, Some(r2)) => match policy {
                OrphanPolicy::KeepSecond => PairedTrimResult::OnlySecond(r2),
                _ => PairedTrimResult::Dropped,
            },
            (None, None) => PairedTrimResult::Dropped,
        }
    }

    /// Process a batch of pairs, keeping only those where both reads pass.
    ///
    /// Uses [`OrphanPolicy::DropBoth`] — pairs with a single surviving read
    /// are discarded.
    pub fn process_paired_batch(
        &self,
        pairs: &[PairedFastqRecord],
    ) -> Vec<PairedFastqRecord> {
        pairs
            .iter()
            .filter_map(|pair| {
                let r1 = self.process(pair.r1())?;
                let r2 = self.process(pair.r2())?;
                Some(PairedFastqRecord::new_unchecked(r1, r2))
            })
            .collect()
    }

    /// Process a batch of pairs and collect detailed statistics.
    pub fn process_paired_batch_with_stats(
        &self,
        pairs: &[PairedFastqRecord],
    ) -> PairedTrimReport {
        let total_input = pairs.len();
        let mut total_bases_input: u64 = 0;
        let mut total_bases_output: u64 = 0;
        let mut both_passed: usize = 0;
        let mut r1_only_passed: usize = 0;
        let mut r2_only_passed: usize = 0;
        let mut both_failed: usize = 0;
        let mut kept = Vec::new();

        for pair in pairs {
            total_bases_input += pair.r1().sequence().len() as u64;
            total_bases_input += pair.r2().sequence().len() as u64;

            let r1_result = self.process(pair.r1());
            let r2_result = self.process(pair.r2());

            match (r1_result, r2_result) {
                (Some(r1), Some(r2)) => {
                    total_bases_output += r1.sequence().len() as u64;
                    total_bases_output += r2.sequence().len() as u64;
                    both_passed += 1;
                    kept.push(PairedFastqRecord::new_unchecked(r1, r2));
                }
                (Some(_), None) => {
                    r1_only_passed += 1;
                }
                (None, Some(_)) => {
                    r2_only_passed += 1;
                }
                (None, None) => {
                    both_failed += 1;
                }
            }
        }

        PairedTrimReport {
            kept,
            total_input,
            both_passed,
            r1_only_passed,
            r2_only_passed,
            both_failed,
            total_bases_input,
            total_bases_output,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::DnaSequence;

    /// Helper to create a FastqRecord for testing.
    fn make_record(seq: &[u8], quals: &[u8]) -> FastqRecord {
        let sequence = DnaSequence::new(seq).unwrap();
        let quality = QualityScores::from_raw(quals.to_vec());
        FastqRecord::new("test".into(), None, sequence, quality).unwrap()
    }

    // --- Sliding window ---

    #[test]
    fn sliding_window_all_high() {
        let q = &[30, 30, 30, 30, 30, 30, 30, 30];
        let r = trim_sliding_window(q, 4, 15.0);
        assert_eq!(r, TrimRange { start: 0, end: 8 });
    }

    #[test]
    fn sliding_window_mid_drop() {
        // Quality drops in the middle
        let q = &[30, 30, 30, 30, 5, 5, 5, 5];
        let r = trim_sliding_window(q, 4, 15.0);
        // Window [30,5,5,5] at pos 3 has mean 11.25 < 15 → cut at pos 6
        // Window [30,30,5,5] at pos 2 has mean 17.5 >= 15 → ok
        // Window [30,5,5,5] at pos 3 has mean 11.25 < 15 → cut at end of prev window = 6
        assert!(r.end <= 7);
        assert!(r.end >= 4);
    }

    #[test]
    fn sliding_window_immediate_drop() {
        let q = &[2, 2, 2, 2];
        let r = trim_sliding_window(q, 4, 15.0);
        assert_eq!(r, TrimRange { start: 0, end: 0 });
    }

    #[test]
    fn sliding_window_window_1() {
        // Window size 1 = per-base threshold
        let q = &[30, 30, 5, 30];
        let r = trim_sliding_window(q, 1, 15.0);
        assert_eq!(r, TrimRange { start: 0, end: 2 });
    }

    #[test]
    fn sliding_window_last_window_drop() {
        let q = &[30, 30, 30, 30, 30, 2, 2, 2];
        let r = trim_sliding_window(q, 4, 15.0);
        // Window at pos 4 = [30,2,2,2] mean=9 < 15 → cut at pos 7
        assert!(r.end < 8);
    }

    // --- Leading / trailing ---

    #[test]
    fn leading_no_trim() {
        let q = &[30, 30, 30, 30];
        assert_eq!(trim_leading(q, 20), TrimRange { start: 0, end: 4 });
    }

    #[test]
    fn leading_partial() {
        let q = &[2, 5, 30, 30];
        assert_eq!(trim_leading(q, 20), TrimRange { start: 2, end: 4 });
    }

    #[test]
    fn trailing_no_trim() {
        let q = &[30, 30, 30, 30];
        assert_eq!(trim_trailing(q, 20), TrimRange { start: 0, end: 4 });
    }

    #[test]
    fn trailing_partial() {
        let q = &[30, 30, 5, 2];
        assert_eq!(trim_trailing(q, 20), TrimRange { start: 0, end: 2 });
    }

    #[test]
    fn leading_trailing_combined() {
        let q = &[2, 5, 30, 30, 5, 2];
        let r1 = trim_leading(q, 20);
        let r2 = trim_trailing(q, 20);
        let combined = intersect_ranges(&[r1, r2]);
        assert_eq!(combined, TrimRange { start: 2, end: 4 });
    }

    // --- BWA quality trim ---

    #[test]
    fn bwa_clean_read() {
        let q = &[30, 30, 30, 30, 30];
        let r = trim_quality_3prime(q, 20);
        assert_eq!(r.end, 5);
    }

    #[test]
    fn bwa_3prime_ramp_down() {
        let q = &[30, 30, 30, 10, 5, 2];
        let r = trim_quality_3prime(q, 20);
        assert!(r.end <= 3);
    }

    #[test]
    fn bwa_all_low() {
        let q = &[2, 2, 2, 2];
        let r = trim_quality_3prime(q, 20);
        assert_eq!(r.end, 0);
    }

    #[test]
    fn bwa_isolated_low_base() {
        // One low base among high — BWA should keep most of the read
        let q = &[30, 30, 5, 30, 30];
        let r = trim_quality_3prime(q, 20);
        assert_eq!(r.end, 5);
    }

    // --- Adapter detection ---

    #[test]
    fn adapter_exact_match() {
        // Read ends with the adapter prefix
        let seq = b"ACGTACGTACGTAACCAGATCGGAAGAG";
        let cut = find_adapter_3prime(seq, adapters::TRUSEQ_PREFIX, 0);
        assert_eq!(cut, 16);
    }

    #[test]
    fn adapter_partial_3prime() {
        // Only 8 bases of adapter at the end of the read
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let seq = b"ACGTACGTACGTACGTAACCAGATCGGAAGAG";
        // Adapter prefix at position 21
        let cut = find_adapter_3prime(seq, adapter, 0);
        assert!(cut < seq.len());
    }

    #[test]
    fn adapter_one_mismatch() {
        // Adapter with one mismatch
        let seq = b"ACGTACGTACGTAACCAGATCGGAATAG";
        // Mismatch at position 25 (G→T in "AAGAG" → "AATAG")
        let cut = find_adapter_3prime(seq, adapters::TRUSEQ_PREFIX, 1);
        assert_eq!(cut, 16);
    }

    #[test]
    fn adapter_too_many_mismatches() {
        let seq = b"ACGTACGTACGTAACCNNNNNNNNNNN";
        let cut = find_adapter_3prime(seq, adapters::TRUSEQ_PREFIX, 1);
        assert_eq!(cut, seq.len());
    }

    #[test]
    fn adapter_no_adapter() {
        let seq = b"ACGTACGTACGTACGT";
        let cut = find_adapter_3prime(seq, adapters::TRUSEQ_PREFIX, 1);
        assert_eq!(cut, seq.len());
    }

    // --- Shannon entropy ---

    #[test]
    fn entropy_homopolymer() {
        let e = shannon_entropy(b"AAAAAAAAAA");
        assert!((e - 0.0).abs() < 1e-10);
    }

    #[test]
    fn entropy_equiprobable() {
        let e = shannon_entropy(b"ACGTACGTACGTACGT");
        assert!((e - 2.0).abs() < 1e-10);
    }

    #[test]
    fn entropy_dinucleotide() {
        let e = shannon_entropy(b"ACACACACACACACAC");
        assert!((e - 1.0).abs() < 1e-10);
    }

    #[test]
    fn entropy_empty() {
        assert_eq!(shannon_entropy(b""), 0.0);
    }

    // --- Filter functions ---

    #[test]
    fn filter_length_pass() {
        let r = make_record(b"ACGTACGT", &[30; 8]);
        assert!(filter_by_length(&r, 4, 100).is_some());
    }

    #[test]
    fn filter_length_too_short() {
        let r = make_record(b"ACG", &[30; 3]);
        assert!(filter_by_length(&r, 4, 100).is_none());
    }

    #[test]
    fn filter_length_too_long() {
        let r = make_record(b"ACGTACGT", &[30; 8]);
        assert!(filter_by_length(&r, 1, 4).is_none());
    }

    #[test]
    fn filter_quality_pass() {
        let r = make_record(b"ACGT", &[30, 30, 30, 30]);
        assert!(filter_by_quality(&r, 20.0).is_some());
    }

    #[test]
    fn filter_quality_fail() {
        let r = make_record(b"ACGT", &[5, 5, 5, 5]);
        assert!(filter_by_quality(&r, 20.0).is_none());
    }

    #[test]
    fn filter_complexity_pass() {
        let r = make_record(b"ACGTACGTACGTACGT", &[30; 16]);
        assert!(filter_low_complexity(&r, 1.5).is_some());
    }

    #[test]
    fn filter_complexity_fail() {
        let r = make_record(b"AAAAAAAAAAAAAAAA", &[30; 16]);
        assert!(filter_low_complexity(&r, 1.0).is_none());
    }

    // --- TrimRange intersection ---

    #[test]
    fn intersect_overlapping() {
        let r = intersect_ranges(&[
            TrimRange { start: 0, end: 8 },
            TrimRange { start: 2, end: 10 },
        ]);
        assert_eq!(r, TrimRange { start: 2, end: 8 });
    }

    #[test]
    fn intersect_non_overlapping() {
        let r = intersect_ranges(&[
            TrimRange { start: 0, end: 3 },
            TrimRange { start: 5, end: 10 },
        ]);
        assert!(r.is_empty());
    }

    // --- apply_trim ---

    #[test]
    fn apply_trim_basic() {
        let r = make_record(b"ACGTACGT", &[10, 20, 30, 40, 30, 20, 10, 5]);
        let trimmed = apply_trim(&r, TrimRange { start: 2, end: 6 }).unwrap();
        assert_eq!(trimmed.sequence().as_bytes(), b"GTAC");
        assert_eq!(trimmed.quality().as_slice(), &[30, 40, 30, 20]);
    }

    #[test]
    fn apply_trim_empty_range() {
        let r = make_record(b"ACGT", &[30; 4]);
        assert!(apply_trim(&r, TrimRange { start: 5, end: 2 }).is_none());
    }

    // --- Pipeline integration ---

    #[test]
    fn pipeline_noop() {
        let pipeline = TrimPipeline::new();
        let r = make_record(b"ACGTACGT", &[30; 8]);
        let result = pipeline.process(&r).unwrap();
        assert_eq!(result.sequence().as_bytes(), b"ACGTACGT");
    }

    #[test]
    fn pipeline_full() {
        let pipeline = TrimPipeline::new()
            .leading(20)
            .trailing(20)
            .sliding_window(4, 15.0)
            .min_length(2);

        let r = make_record(b"ACGTACGT", &[5, 30, 30, 30, 30, 30, 30, 5]);
        let result = pipeline.process(&r).unwrap();
        // Leading trims first base (Q5 < 20), trailing trims last (Q5 < 20)
        assert_eq!(result.sequence().as_bytes(), b"CGTACG");
    }

    #[test]
    fn pipeline_everything_filtered() {
        let pipeline = TrimPipeline::new().min_mean_quality(40.0);
        let r = make_record(b"ACGT", &[10, 10, 10, 10]);
        assert!(pipeline.process(&r).is_none());
    }

    #[test]
    fn pipeline_batch_stats() {
        let pipeline = TrimPipeline::new().min_mean_quality(20.0);
        let records = vec![
            make_record(b"ACGT", &[30, 30, 30, 30]), // passes
            make_record(b"ACGT", &[5, 5, 5, 5]),     // filtered by quality
            make_record(b"ACGT", &[25, 25, 25, 25]), // passes
        ];
        let report = pipeline.process_batch_with_stats(&records);
        assert_eq!(report.total_input, 3);
        assert_eq!(report.kept.len(), 2);
        assert_eq!(report.filtered_by_quality, 1);
    }

    // --- Edge cases ---

    #[test]
    fn empty_sequence() {
        let seq = DnaSequence::new(b"").unwrap();
        let qual = QualityScores::from_raw(vec![]);
        let r = FastqRecord::new("empty".into(), None, seq, qual).unwrap();

        assert_eq!(trim_sliding_window(&[], 4, 15.0), TrimRange { start: 0, end: 0 });
        assert!(TrimPipeline::new().process(&r).is_some());
    }

    #[test]
    fn single_base() {
        let r = make_record(b"A", &[30]);
        let pipeline = TrimPipeline::new().min_length(1);
        assert!(pipeline.process(&r).is_some());
    }

    #[test]
    fn uniform_quality() {
        let q = &[20, 20, 20, 20, 20];
        assert_eq!(trim_sliding_window(q, 4, 20.0), TrimRange { start: 0, end: 5 });
        assert_eq!(trim_leading(q, 20), TrimRange { start: 0, end: 5 });
        assert_eq!(trim_trailing(q, 20), TrimRange { start: 0, end: 5 });
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use crate::types::DnaSequence;
    use proptest::prelude::*;

    fn dna_and_quality(max_len: usize) -> impl Strategy<Value = (Vec<u8>, Vec<u8>)> {
        (1..=max_len).prop_flat_map(|len| {
            let seq = proptest::collection::vec(
                prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')],
                len,
            );
            let qual = proptest::collection::vec(0..=41u8, len);
            (seq, qual)
        })
    }

    proptest! {
        #[test]
        fn trimmed_never_longer(
            (seq, qual) in dna_and_quality(200)
        ) {
            let record = {
                let s = DnaSequence::new(&seq).unwrap();
                let q = QualityScores::from_raw(qual.clone());
                FastqRecord::new("test".into(), None, s, q).unwrap()
            };
            let pipeline = TrimPipeline::new()
                .leading(10)
                .trailing(10)
                .sliding_window(4, 15.0);
            if let Some(trimmed) = pipeline.process(&record) {
                prop_assert!(trimmed.sequence().len() <= record.sequence().len());
            }
        }

        #[test]
        fn intersect_valid_subrange(
            s1 in 0..50usize,
            e1 in 0..50usize,
            s2 in 0..50usize,
            e2 in 0..50usize,
        ) {
            let r1 = TrimRange { start: s1, end: e1 };
            let r2 = TrimRange { start: s2, end: e2 };
            let result = intersect_ranges(&[r1, r2]);
            prop_assert!(result.start <= result.end);
            if !r1.is_empty() && !r2.is_empty() {
                prop_assert!(result.start >= s1.max(s2));
                prop_assert!(result.end <= e1.min(e2).max(result.start));
            }
        }
    }
}

#[cfg(test)]
#[cfg(feature = "std")]
mod paired_tests {
    use super::*;
    use crate::paired::PairedFastqRecord;
    use crate::types::DnaSequence;

    fn make_record(seq: &[u8], quals: &[u8]) -> FastqRecord {
        let sequence = DnaSequence::new(seq).unwrap();
        let quality = QualityScores::from_raw(quals.to_vec());
        FastqRecord::new("test".into(), None, sequence, quality).unwrap()
    }

    fn make_pair(
        seq1: &[u8],
        quals1: &[u8],
        seq2: &[u8],
        quals2: &[u8],
    ) -> PairedFastqRecord {
        PairedFastqRecord::new_unchecked(make_record(seq1, quals1), make_record(seq2, quals2))
    }

    #[test]
    fn paired_both_pass() {
        let pipeline = TrimPipeline::new().min_mean_quality(10.0);
        let r1 = make_record(b"ACGT", &[30; 4]);
        let r2 = make_record(b"TGCA", &[30; 4]);
        match pipeline.process_paired(&r1, &r2, OrphanPolicy::DropBoth) {
            PairedTrimResult::BothPassed(_, _) => {}
            other => panic!("expected BothPassed, got {:?}", other),
        }
    }

    #[test]
    fn paired_r1_fails_drop_both() {
        let pipeline = TrimPipeline::new().min_mean_quality(20.0);
        let r1 = make_record(b"ACGT", &[5; 4]);
        let r2 = make_record(b"TGCA", &[30; 4]);
        match pipeline.process_paired(&r1, &r2, OrphanPolicy::DropBoth) {
            PairedTrimResult::Dropped => {}
            other => panic!("expected Dropped, got {:?}", other),
        }
    }

    #[test]
    fn paired_r1_fails_keep_second() {
        let pipeline = TrimPipeline::new().min_mean_quality(20.0);
        let r1 = make_record(b"ACGT", &[5; 4]);
        let r2 = make_record(b"TGCA", &[30; 4]);
        match pipeline.process_paired(&r1, &r2, OrphanPolicy::KeepSecond) {
            PairedTrimResult::OnlySecond(_) => {}
            other => panic!("expected OnlySecond, got {:?}", other),
        }
    }

    #[test]
    fn paired_r2_fails_keep_first() {
        let pipeline = TrimPipeline::new().min_mean_quality(20.0);
        let r1 = make_record(b"ACGT", &[30; 4]);
        let r2 = make_record(b"TGCA", &[5; 4]);
        match pipeline.process_paired(&r1, &r2, OrphanPolicy::KeepFirst) {
            PairedTrimResult::OnlyFirst(_) => {}
            other => panic!("expected OnlyFirst, got {:?}", other),
        }
    }

    #[test]
    fn paired_r2_fails_drop_both() {
        let pipeline = TrimPipeline::new().min_mean_quality(20.0);
        let r1 = make_record(b"ACGT", &[30; 4]);
        let r2 = make_record(b"TGCA", &[5; 4]);
        match pipeline.process_paired(&r1, &r2, OrphanPolicy::DropBoth) {
            PairedTrimResult::Dropped => {}
            other => panic!("expected Dropped, got {:?}", other),
        }
    }

    #[test]
    fn paired_both_fail() {
        let pipeline = TrimPipeline::new().min_mean_quality(20.0);
        let r1 = make_record(b"ACGT", &[5; 4]);
        let r2 = make_record(b"TGCA", &[5; 4]);
        match pipeline.process_paired(&r1, &r2, OrphanPolicy::KeepFirst) {
            PairedTrimResult::Dropped => {}
            other => panic!("expected Dropped, got {:?}", other),
        }
    }

    #[test]
    fn paired_batch_drop_both() {
        let pipeline = TrimPipeline::new().min_mean_quality(20.0);
        let pairs = vec![
            make_pair(b"ACGT", &[30; 4], b"TGCA", &[30; 4]),
            make_pair(b"ACGT", &[5; 4], b"TGCA", &[30; 4]),
        ];
        let kept = pipeline.process_paired_batch(&pairs);
        assert_eq!(kept.len(), 1);
    }

    #[test]
    fn paired_batch_stats() {
        let pipeline = TrimPipeline::new().min_mean_quality(20.0);
        let pairs = vec![
            make_pair(b"ACGT", &[30; 4], b"TGCA", &[30; 4]), // both pass
            make_pair(b"ACGT", &[5; 4], b"TGCA", &[30; 4]),  // r1 fails
            make_pair(b"ACGT", &[30; 4], b"TGCA", &[5; 4]),  // r2 fails
            make_pair(b"ACGT", &[5; 4], b"TGCA", &[5; 4]),   // both fail
        ];
        let report = pipeline.process_paired_batch_with_stats(&pairs);
        assert_eq!(report.total_input, 4);
        assert_eq!(report.both_passed, 1);
        assert_eq!(report.r1_only_passed, 1);
        assert_eq!(report.r2_only_passed, 1);
        assert_eq!(report.both_failed, 1);
        assert_eq!(report.orphans(), 2);
        assert!((report.survival_rate() - 0.25).abs() < 1e-10);
        assert_eq!(report.kept.len(), 1);
    }

    #[test]
    fn paired_batch_stats_bases() {
        let pipeline = TrimPipeline::new();
        let pairs = vec![make_pair(b"ACGTACGT", &[30; 8], b"TGCA", &[30; 4])];
        let report = pipeline.process_paired_batch_with_stats(&pairs);
        assert_eq!(report.total_bases_input, 12);
        assert_eq!(report.total_bases_output, 12);
    }

    #[test]
    fn paired_batch_stats_empty() {
        let pipeline = TrimPipeline::new();
        let pairs: Vec<PairedFastqRecord> = vec![];
        let report = pipeline.process_paired_batch_with_stats(&pairs);
        assert_eq!(report.total_input, 0);
        assert_eq!(report.both_passed, 0);
        assert!((report.survival_rate() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn paired_full_pipeline() {
        let pipeline = TrimPipeline::new()
            .adapter(b"AGATCGGAAGAG")
            .leading(3)
            .trailing(3)
            .sliding_window(4, 15.0)
            .min_length(4);

        let r1 = make_record(b"ACGTACGTACGTACGT", &[30; 16]);
        let r2 = make_record(b"TGCATGCATGCATGCA", &[30; 16]);
        match pipeline.process_paired(&r1, &r2, OrphanPolicy::DropBoth) {
            PairedTrimResult::BothPassed(_, _) => {}
            other => panic!("expected BothPassed, got {:?}", other),
        }
    }
}

#[cfg(test)]
#[cfg(feature = "std")]
mod paired_proptests {
    use super::*;
    use crate::types::DnaSequence;
    use proptest::prelude::*;

    fn dna_and_quality(max_len: usize) -> impl Strategy<Value = (Vec<u8>, Vec<u8>)> {
        (1..=max_len).prop_flat_map(|len| {
            let seq = proptest::collection::vec(
                prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')],
                len,
            );
            let qual = proptest::collection::vec(0..=41u8, len);
            (seq, qual)
        })
    }

    proptest! {
        #[test]
        fn paired_trimmed_never_longer(
            (seq1, qual1) in dna_and_quality(100),
            (seq2, qual2) in dna_and_quality(100),
        ) {
            let r1 = {
                let s = DnaSequence::new(&seq1).unwrap();
                let q = QualityScores::from_raw(qual1);
                FastqRecord::new("test".into(), None, s, q).unwrap()
            };
            let r2 = {
                let s = DnaSequence::new(&seq2).unwrap();
                let q = QualityScores::from_raw(qual2);
                FastqRecord::new("test".into(), None, s, q).unwrap()
            };

            let pipeline = TrimPipeline::new()
                .leading(10)
                .trailing(10)
                .sliding_window(4, 15.0);

            match pipeline.process_paired(&r1, &r2, OrphanPolicy::KeepFirst) {
                PairedTrimResult::BothPassed(tr1, tr2) => {
                    prop_assert!(tr1.sequence().len() <= r1.sequence().len());
                    prop_assert!(tr2.sequence().len() <= r2.sequence().len());
                }
                PairedTrimResult::OnlyFirst(tr1) => {
                    prop_assert!(tr1.sequence().len() <= r1.sequence().len());
                }
                _ => {}
            }
        }
    }
}
