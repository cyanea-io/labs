# cyanea-seq

Sequence types, validation, manipulation, and file parsing for DNA, RNA, and protein sequences.

## Status: Complete

All planned functionality is implemented and tested. FASTA/FASTQ parsing uses `needletail` for streaming performance.

## Public API

### Alphabet system (`alphabet.rs`)

| Type | Description |
|------|-------------|
| `trait Alphabet` | Defines valid bytes for a sequence type (`NAME`, `VALID_BYTES`, `is_valid`) |
| `DnaAlphabet` | IUPAC DNA: `ACGTNRYSWKMBDHV` |
| `RnaAlphabet` | IUPAC RNA: `ACGUNRYSWKMBDHV` |
| `ProteinAlphabet` | 20 standard amino acids + `XBZJUO*` |

### Sequence types (`types.rs`, `seq.rs`)

All sequence types are `ValidatedSeq<A: Alphabet>` -- validated, uppercase byte sequences.

| Type | Description |
|------|-------------|
| `DnaSequence` | Validated DNA sequence |
| `RnaSequence` | Validated RNA sequence |
| `ProteinSequence` | Validated protein sequence |

**DnaSequence methods:**

| Method | Description |
|--------|-------------|
| `new(data: &[u8]) -> Result<Self>` | Validate and construct |
| `reverse_complement() -> DnaSequence` | Reverse complement |
| `transcribe() -> RnaSequence` | DNA to RNA (T -> U) |
| `translate() -> Result<ProteinSequence>` | Transcribe then translate (NCBI Table 1) |
| `gc_content() -> f64` | GC fraction in [0.0, 1.0] |
| `kmers(k) -> Result<KmerIter>` | Extract k-mers |

**RnaSequence methods:**

| Method | Description |
|--------|-------------|
| `reverse_complement() -> RnaSequence` | Reverse complement |
| `reverse_transcribe() -> DnaSequence` | RNA to DNA (U -> T) |
| `translate() -> Result<ProteinSequence>` | Translate (NCBI Table 1) |
| `translate_frames() -> [Result<ProteinSequence>; 3]` | All three reading frames |
| `kmers(k) -> Result<KmerIter>` | Extract k-mers |

**ProteinSequence methods:**

| Method | Description |
|--------|-------------|
| `molecular_weight() -> f64` | Molecular weight in Daltons |
| `kmers(k) -> Result<KmerIter>` | Extract k-mers |

### Codon translation (`codon.rs`)

| Function | Description |
|----------|-------------|
| `translate_codon(codon: &[u8]) -> Option<u8>` | Single codon to amino acid (NCBI standard code) |
| `translate_sequence(seq: &[u8]) -> Vec<u8>` | Translate full coding sequence |

### K-mers (`kmer.rs`)

| Type | Description |
|------|-------------|
| `KmerIter<'a>` | Iterator over k-mers; implements `Iterator`, `ExactSizeIterator`, `DoubleEndedIterator` |

### Quality scores (`quality.rs`)

| Type/Function | Description |
|---------------|-------------|
| `PhredEncoding` | Enum: `Phred33`, `Phred64` |
| `QualityScores` | Decoded Phred scores with statistics |
| `QualityScores::mean()` | Mean quality score |
| `QualityScores::fraction_above(threshold)` | Fraction of bases above a quality threshold |
| `QualityScores::error_probability(phred)` | Convert Phred score to error probability |

### FASTA parsing (`fasta.rs`)

| Type/Function | Description |
|---------------|-------------|
| `FastaStats` | `sequence_count`, `total_bases`, `gc_content`, `avg_length` |
| `parse_fasta_stats(path) -> Result<FastaStats>` | Streaming FASTA statistics |

### FASTQ parsing (`fastq.rs`)

| Type/Function | Description |
|---------------|-------------|
| `FastqRecord` | Single record: `id`, `seq`, `quality` |
| `FastqStats` | `sequence_count`, `total_bases`, `gc_content`, `avg_length`, `mean_quality`, `q20_fraction`, `q30_fraction` |
| `parse_fastq_file(path) -> Result<Vec<FastqRecord>>` | Parse all records |
| `parse_fastq_stats(path) -> Result<FastqStats>` | Streaming statistics |

### 2-bit encoding (`twobit.rs`)

Compact 2-bit-per-base representation for unambiguous DNA (A, C, G, T only). Achieves 4x compression over ASCII. Encoding: A=00, C=01, G=10, T=11.

| Type | Description |
|------|-------------|
| `TwoBitSequence` | DNA sequence stored in 2-bit packed representation (4 bases per byte) |

**TwoBitSequence methods:**

| Method | Description |
|--------|-------------|
| `encode(seq: &[u8]) -> Result<Self>` | Encode ASCII DNA into 2-bit packed form (case-insensitive, ACGT only) |
| `decode(&self) -> Vec<u8>` | Decode back to uppercase ASCII DNA bytes |
| `get(&self, index: usize) -> Option<u8>` | Get the base at a specific position |
| `len(&self) -> usize` | Number of bases in the sequence |
| `is_empty(&self) -> bool` | Whether the sequence is empty |
| `kmer(&self, pos: usize, k: usize) -> Option<u64>` | Extract a k-mer as a 2-bit integer encoding (max k=32) |
| `complement(&self) -> Self` | Bitwise complement (A<->T, C<->G) via XOR |

### Suffix array (`suffix.rs`)

O(n) suffix array construction via the SA-IS algorithm (Nong, Zhang & Chan, 2009) with O(m log n) binary search for pattern matching.

| Type | Description |
|------|-------------|
| `SuffixArray` | Suffix array built from a text, storing sorted suffix positions |

**SuffixArray methods:**

| Method | Description |
|--------|-------------|
| `build(text: &[u8]) -> Self` | Build a suffix array using SA-IS (appends sentinel, O(n) construction) |
| `search(&self, text: &[u8], pattern: &[u8]) -> Vec<usize>` | Find all occurrences of pattern in text (O(m log n), sorted positions) |
| `len(&self) -> usize` | Number of entries in the suffix array (text length + 1 for sentinel) |

### FM-Index (`fm_index.rs`)

FM-Index for O(m) exact pattern matching on DNA sequences via backward search. Built on the Burrows-Wheeler Transform (BWT), occurrence table, and C table. Supports the DNA alphabet (A, C, G, T) plus sentinel '$'.

| Type | Description |
|------|-------------|
| `FmIndex` | FM-Index for DNA sequences with BWT, suffix array, occurrence table, and C table |

**FmIndex methods:**

| Method | Description |
|--------|-------------|
| `build(text: &[u8]) -> Self` | Build an FM-Index from a DNA sequence (ACGT only, sentinel appended internally) |
| `search(&self, pattern: &[u8]) -> Vec<usize>` | Find all occurrence positions of pattern (O(m) search + O(k) lookup, sorted) |
| `count(&self, pattern: &[u8]) -> usize` | Count occurrences without returning positions (more efficient than `search`) |

### Pattern matching (`pattern.rs`)

7 exact and approximate string matching algorithms operating on `&[u8]` for generality. Bitparallel methods (Shift-And, BNDM, Myers) limited to patterns of length 64 or less.

**Exact matching:**

| Function | Description |
|----------|-------------|
| `horspool(text: &[u8], pattern: &[u8]) -> Vec<usize>` | Boyer-Moore-Horspool exact matching, O(n/m) average case |
| `kmp(text: &[u8], pattern: &[u8]) -> Vec<usize>` | Knuth-Morris-Pratt exact matching, O(n+m) |
| `shift_and(text: &[u8], pattern: &[u8]) -> Vec<usize>` | Shift-And bitparallel exact matching (pattern <= 64) |
| `bndm(text: &[u8], pattern: &[u8]) -> Vec<usize>` | BNDM bitparallel exact matching (pattern <= 64) |
| `bom(text: &[u8], pattern: &[u8]) -> Vec<usize>` | Backward Oracle Matching exact matching |

**Approximate matching:**

| Function | Description |
|----------|-------------|
| `myers_bitparallel(text: &[u8], pattern: &[u8], max_dist: usize) -> Vec<(usize, usize)>` | Myers bit-parallel approximate matching (pattern <= 64), returns (end_position, edit_distance) |
| `ukkonen(text: &[u8], pattern: &[u8], max_dist: usize) -> Vec<(usize, usize)>` | Ukkonen cut-off approximate matching, returns (end_position, edit_distance) |

### PSSM / Motif scanning (`pssm.rs`)

Position-Specific Scoring Matrix with const generic alphabet size for DNA and protein motif scanning.

| Type | Description |
|------|-------------|
| `Pssm<const N: usize>` | Position-Specific Scoring Matrix with const generic alphabet size |
| `PssmDna` | Type alias for `Pssm<4>` (DNA alphabet) |
| `PssmProtein` | Type alias for `Pssm<20>` (protein alphabet) |

**Pssm methods:**

| Method | Description |
|--------|-------------|
| `from_counts(counts: &[Vec<f64>], pseudocount: f64, background: &[f64], mapping: fn(u8) -> Option<usize>) -> Result<Self>` | Build from count matrix with pseudocounts and background frequencies |
| `score(seq: &[u8], mapping: fn(u8) -> Option<usize>) -> Result<f64>` | Score a sequence window of length equal to the PSSM |
| `scan(seq: &[u8], threshold: f64, mapping: fn(u8) -> Option<usize>) -> Result<Vec<(usize, f64)>>` | Scan for motif hits above threshold, returns (position, score) |
| `information_content() -> Vec<f64>` | Per-position information content in bits |

**Helper functions:**

| Function | Description |
|----------|-------------|
| `dna_mapping(b: u8) -> Option<usize>` | Map DNA base to column index (A=0, C=1, G=2, T=3) |
| `protein_mapping(b: u8) -> Option<usize>` | Map amino acid to column index (20 standard residues) |

### ORF finder (`orf.rs`)

Open reading frame detection across forward and reverse complement strands with configurable start/stop codons.

| Type | Description |
|------|-------------|
| `OrfResult` | ORF result: `start`, `end`, `frame`, `strand`, `sequence` |
| `Strand` | Enum: `Forward`, `Reverse` |

| Function | Description |
|----------|-------------|
| `find_orfs(seq: &[u8], min_length: usize) -> Vec<OrfResult>` | Find ORFs in all 3 forward reading frames |
| `find_orfs_both_strands(seq: &[u8], min_length: usize) -> Vec<OrfResult>` | Find ORFs in all 6 reading frames (forward + reverse complement) |
| `find_orfs_with_codons(seq: &[u8], min_length: usize, start_codons: &[&[u8]], stop_codons: &[&[u8]]) -> Vec<OrfResult>` | Find ORFs with configurable start and stop codons |

### FASTA indexed reader (`fasta_index.rs`)

Random access FASTA reading using `.fai` index files. Feature-gated behind `std`.

| Type | Description |
|------|-------------|
| `FastaIndex` | Parsed `.fai` index with entries keyed by sequence name |
| `FastaIndexEntry` | Entry for one sequence: `name`, `length`, `offset`, `line_bases`, `line_width` |
| `IndexedFastaReader` | Random access FASTA reader using an index |

**FastaIndex methods:**

| Method | Description |
|--------|-------------|
| `from_file(path: &Path) -> Result<Self>` | Parse an existing `.fai` index file |
| `build(fasta_path: &Path) -> Result<Self>` | Build an index from a FASTA file |
| `write(path: &Path) -> Result<()>` | Write the index to a `.fai` file |

**IndexedFastaReader methods:**

| Method | Description |
|--------|-------------|
| `open(fasta_path: &Path, fai_path: &Path) -> Result<Self>` | Open a FASTA file with its index |
| `fetch(name: &str, start: usize, end: usize) -> Result<Vec<u8>>` | Fetch a region by sequence name and coordinates |
| `fetch_all(name: &str) -> Result<Vec<u8>>` | Fetch the full sequence by name |

### FMD-Index (`fmd_index.rs`)

Bidirectional FM-Index for strand-aware search on DNA sequences. Indexes `text#reverse_complement(text)$` to support both forward and backward extension of search intervals, enabling super-maximal exact match (SMEM) enumeration.

| Type | Description |
|------|-------------|
| `FmdIndex` | Bidirectional FM-Index over `text#revcomp(text)$` |
| `BiInterval` | Bidirectional interval: `lower`, `size`, `lower_rev` |

**FmdIndex methods:**

| Method | Description |
|--------|-------------|
| `new(seq: &[u8]) -> Self` | Build from a DNA sequence |
| `init_interval(c: u8) -> BiInterval` | Initialize a bidirectional interval with a single character |
| `extend_backward(&self, interval: &BiInterval, c: u8) -> BiInterval` | Extend search backward (prepend character) |
| `extend_forward(&self, interval: &BiInterval, c: u8) -> BiInterval` | Extend search forward (append character) |
| `locate(&self, interval: &BiInterval) -> Vec<usize>` | Get text positions from a bidirectional interval |
| `backward_search(&self, pattern: &[u8]) -> BiInterval` | Full backward search for a pattern |
| `smems(&self, query: &[u8], min_len: usize) -> Vec<(usize, usize, BiInterval)>` | Enumerate super-maximal exact matches (query_start, query_end, interval) |

### BWT utilities (`bwt.rs`)

Standalone Burrows-Wheeler Transform construction and inversion.

| Type | Description |
|------|-------------|
| `Bwt` | Burrows-Wheeler Transform with primary index |

**Bwt methods:**

| Method | Description |
|--------|-------------|
| `build(text: &[u8]) -> Self` | Build BWT from text (SA-via-sort, sentinel appended internally) |
| `as_bytes(&self) -> &[u8]` | The BWT string |
| `primary_index(&self) -> usize` | Position of sentinel in the BWT |
| `len(&self) -> usize` | Length of BWT (text length + 1) |
| `invert(&self) -> Vec<u8>` | Reconstruct original text via LF-mapping |

### Paired-end FASTQ (`paired.rs`)

Paired-end read support for Illumina short-read workflows. Handles separate R1/R2 files and interleaved files with mate validation. Feature-gated behind `std`.

**Core types:**

| Type | Description |
|------|-------------|
| `PairedFastqRecord` | Paired-end record containing R1 (forward) and R2 (reverse) `FastqRecord` |
| `MateValidation` | Enum: `Strict` (requires `/1` `/2`), `Relaxed` (any suffix), `None` (trust order) |
| `PairedFastqStats` | Per-file `FastqStats` for R1 and R2, plus `pair_count` |

**PairedFastqRecord methods:**

| Method | Description |
|--------|-------------|
| `new(r1, r2, validation) -> Result<Self>` | Create with mate validation |
| `new_unchecked(r1, r2) -> Self` | Create without validation |
| `r1() -> &FastqRecord` | Forward read |
| `r2() -> &FastqRecord` | Reverse read |
| `into_reads() -> (FastqRecord, FastqRecord)` | Consume and return both reads |
| `pair_name() -> &str` | Shared name (R1 name with read-number suffix stripped) |

**Helper functions:**

| Function | Description |
|----------|-------------|
| `strip_read_suffix(name) -> &str` | Strip `/1`, `/2`, `_1`, `_2` suffixes |
| `validate_mate_pair(r1, r2) -> Result<()>` | Relaxed validation (name prefixes match) |
| `validate_mate_pair_strict(r1, r2) -> Result<()>` | Strict validation (requires `/1` and `/2` suffixes) |

**Parsing:**

| Function | Description |
|----------|-------------|
| `parse_paired_fastq_files(r1_path, r2_path, validation) -> Result<Vec<PairedFastqRecord>>` | Parse separate R1/R2 files in lockstep |
| `parse_interleaved_fastq(path, validation) -> Result<Vec<PairedFastqRecord>>` | Parse interleaved file (alternating R1/R2) |
| `parse_paired_fastq_stats(r1_path, r2_path) -> Result<PairedFastqStats>` | Streaming statistics without storing records |

**Writing:**

| Function | Description |
|----------|-------------|
| `write_paired_fastq(pairs, r1_path, r2_path, encoding) -> Result<()>` | Write to separate R1/R2 files |
| `write_interleaved_fastq(pairs, path, encoding) -> Result<()>` | Write to interleaved file |

**Streaming interleave/deinterleave:**

| Function | Description |
|----------|-------------|
| `interleave_fastq_files(r1_path, r2_path, output_path, validation) -> Result<u64>` | Interleave two files into one (returns pair count) |
| `deinterleave_fastq_file(input_path, r1_path, r2_path, validation) -> Result<u64>` | Split interleaved file into R1/R2 (returns pair count) |

### Quality trimming & filtering (`trim.rs`)

Quality trimming, adapter removal, and read filtering for FASTQ records. Two-level API: low-level functions on `&[u8]` slices return composable `TrimRange` values, high-level functions operate on `FastqRecord`. `TrimPipeline` builder chains operations in Trimmomatic-style order with batch statistics.

**Core type:**

| Type | Description |
|------|-------------|
| `TrimRange` | Half-open range `[start, end)` describing which bases to keep |
| `TrimPipeline` | Configurable read-processing pipeline builder |
| `TrimReport` | Batch processing statistics (kept/filtered counts, adapter hits, base counts) |

**Low-level trim functions** (operate on `&[u8]` quality/sequence slices, return `TrimRange`):

| Function | Description |
|----------|-------------|
| `trim_sliding_window(quality, window_size, threshold) -> TrimRange` | Trimmomatic SLIDINGWINDOW — O(n) running sum, cut when window mean drops below threshold |
| `trim_leading(quality, threshold) -> TrimRange` | Remove 5' bases below quality threshold |
| `trim_trailing(quality, threshold) -> TrimRange` | Remove 3' bases below quality threshold |
| `trim_quality_3prime(quality, threshold) -> TrimRange` | BWA-style — maximize suffix sum of (threshold - Q) from 3' end |
| `find_adapter_3prime(seq, adapter, max_mismatches) -> usize` | Sliding overlap at 3' end, longest-first, returns cut position |
| `shannon_entropy(seq) -> f64` | Base composition entropy (max 2.0 for DNA) |
| `intersect_ranges(ranges) -> TrimRange` | Compose multiple trim ranges via intersection |

**High-level functions** (operate on `&FastqRecord`):

| Function | Description |
|----------|-------------|
| `apply_trim(record, range) -> Option<FastqRecord>` | New trimmed record, or None if range is empty |
| `trim_adapter(record, adapter, max_mismatches) -> FastqRecord` | Record with adapter removed |
| `filter_by_length(record, min, max) -> Option<&FastqRecord>` | None if outside range |
| `filter_low_complexity(record, min_entropy) -> Option<&FastqRecord>` | None if entropy below threshold |
| `filter_by_quality(record, min_quality) -> Option<&FastqRecord>` | None if mean Q below threshold |

**Adapter constants** (`trim::adapters`):

| Constant | Description |
|----------|-------------|
| `TRUSEQ_UNIVERSAL` | Illumina TruSeq Universal Adapter (33 bp) |
| `TRUSEQ_INDEXED` | Illumina TruSeq Indexed Adapter (33 bp) |
| `NEXTERA_READ1` | Nextera Transposase Read 1 (33 bp) |
| `NEXTERA_READ2` | Nextera Transposase Read 2 (34 bp) |
| `SMALL_RNA_3P` | Illumina Small RNA 3' Adapter (21 bp) |
| `TRUSEQ_PREFIX` | Common 12-base prefix shared by TruSeq adapters |
| `ALL_ILLUMINA` | All standard Illumina adapters for batch searching |

**TrimPipeline methods:**

| Method | Description |
|--------|-------------|
| `new() -> Self` | Create an empty pipeline (no-op by default) |
| `adapter(self, adapter: &[u8]) -> Self` | Add a single adapter sequence |
| `illumina_adapters(self) -> Self` | Add all standard Illumina adapters |
| `adapter_mismatches(self, max: usize) -> Self` | Set max mismatches for adapter matching (default 1) |
| `leading(self, threshold: u8) -> Self` | Trim 5' bases below threshold |
| `trailing(self, threshold: u8) -> Self` | Trim 3' bases below threshold |
| `sliding_window(self, window_size: usize, threshold: f64) -> Self` | Trimmomatic-style sliding window (mutually exclusive with `bwa_quality`) |
| `bwa_quality(self, threshold: u8) -> Self` | BWA-style 3' trim (mutually exclusive with `sliding_window`) |
| `min_length(self, len: usize) -> Self` | Discard reads shorter than `len` |
| `max_length(self, len: usize) -> Self` | Discard reads longer than `len` |
| `min_mean_quality(self, quality: f64) -> Self` | Discard reads with mean Q below threshold |
| `min_entropy(self, entropy: f64) -> Self` | Discard low-complexity reads below entropy threshold |
| `process(&self, record: &FastqRecord) -> Option<FastqRecord>` | Process a single record |
| `process_batch(&self, records: &[FastqRecord]) -> Vec<FastqRecord>` | Process a batch, returning only passing records |
| `process_batch_with_stats(&self, records: &[FastqRecord]) -> TrimReport` | Process a batch with detailed statistics |

Pipeline operation order (Trimmomatic convention): adapter trim → leading → trailing → sliding window/BWA → length filter → quality filter → complexity filter.

**Paired-end trimming** (feature-gated behind `std`):

| Type | Description |
|------|-------------|
| `OrphanPolicy` | Enum: `DropBoth`, `KeepFirst` (keep R1 orphans), `KeepSecond` (keep R2 orphans) |
| `PairedTrimResult` | Enum: `BothPassed(r1, r2)`, `OnlyFirst(r1)`, `OnlySecond(r2)`, `Dropped` |
| `PairedTrimReport` | Batch stats: `kept`, `total_input`, `both_passed`, `r1_only_passed`, `r2_only_passed`, `both_failed`, base counts |

| Method | Description |
|--------|-------------|
| `process_paired(&self, r1, r2, policy) -> PairedTrimResult` | Process a single pair with orphan policy |
| `process_paired_batch(&self, pairs) -> Vec<PairedFastqRecord>` | Process batch, DropBoth policy |
| `process_paired_batch_with_stats(&self, pairs) -> PairedTrimReport` | Process batch with detailed statistics |
| `PairedTrimReport::orphans() -> usize` | Count of pairs where only one mate survived |
| `PairedTrimReport::survival_rate() -> f64` | Fraction of pairs where both reads passed |

### MinHash sketching (`minhash.rs`)

Bottom-k MinHash and scaled FracMinHash sketching for rapid genome comparison. Estimates Jaccard similarity, containment, and average nucleotide identity (ANI) between DNA sequences without full alignment. Uses canonical k-mers (min of forward and reverse complement hash) for strand-agnostic sketching.

| Type | Description |
|------|-------------|
| `MinHash` | Bottom-k MinHash sketch: keeps the `sketch_size` smallest canonical k-mer hash values |
| `FracMinHash` | Scaled (fractional) MinHash sketch: keeps all hash values below `u64::MAX / scale` |

**MinHash methods:**

| Method | Description |
|--------|-------------|
| `new(k: usize, sketch_size: usize) -> Result<Self>` | Create an empty bottom-k sketch |
| `from_sequence(seq: &[u8], k: usize, sketch_size: usize) -> Result<Self>` | Build a sketch from a DNA sequence |
| `add_sequence(&mut self, seq: &[u8])` | Add k-mers from a DNA sequence to the sketch |
| `jaccard(&self, other: &MinHash) -> Result<f64>` | Estimate Jaccard similarity (merge-based estimator) |
| `containment(&self, other: &MinHash) -> Result<f64>` | Estimate containment C(A,B) = \|A intersect B\| / \|A\| |
| `ani(&self, other: &MinHash) -> Result<f64>` | Estimate ANI from Jaccard (Mash formula) |
| `len(&self) -> usize` | Number of hash values in the sketch |
| `is_empty(&self) -> bool` | Whether the sketch is empty |
| `k(&self) -> usize` | The k-mer size |
| `sketch_size(&self) -> usize` | The target sketch size (bottom-k parameter) |
| `hashes(&self) -> &[u64]` | Reference to the sorted hash values |

**FracMinHash methods:**

| Method | Description |
|--------|-------------|
| `new(k: usize, scale: u64) -> Result<Self>` | Create an empty scaled sketch |
| `from_sequence(seq: &[u8], k: usize, scale: u64) -> Result<Self>` | Build a scaled sketch from a DNA sequence |
| `add_sequence(&mut self, seq: &[u8])` | Add k-mers from a DNA sequence to the sketch |
| `jaccard(&self, other: &FracMinHash) -> Result<f64>` | Estimate Jaccard similarity |
| `containment(&self, other: &FracMinHash) -> Result<f64>` | Estimate containment C(A,B) = \|A intersect B\| / \|A\| |
| `ani(&self, other: &FracMinHash) -> Result<f64>` | Estimate ANI from Jaccard (Mash formula) |
| `len(&self) -> usize` | Number of hash values in the sketch |
| `is_empty(&self) -> bool` | Whether the sketch is empty |
| `k(&self) -> usize` | The k-mer size |
| `scale(&self) -> u64` | The scale factor |
| `hashes(&self) -> &[u64]` | Reference to the sorted hash values |

### RNA secondary structure prediction (`rna_structure.rs`)

RNA secondary structure prediction using dynamic programming. Operates on RNA/DNA sequences (T→U normalized). Implements three classical algorithms plus structure comparison.

**Types:**

| Type | Description |
|------|-------------|
| `RnaSecondaryStructure` | Pair table representation: `pairs[i] = Some(j)` if paired, `None` if unpaired |
| `NussinovResult` | Nussinov prediction: `structure` + `max_pairs` |
| `MfeResult` | Zuker MFE prediction: `structure` + `energy` (kcal/mol) |
| `PartitionResult` | McCaskill partition function: `pair_probabilities` (n×n) + `ensemble_energy` |

**RnaSecondaryStructure methods:**

| Method | Description |
|--------|-------------|
| `from_dot_bracket(s: &str) -> Result<Self>` | Parse dot-bracket notation (`(`, `)`, `.`) |
| `to_dot_bracket(&self) -> String` | Convert to dot-bracket string |
| `base_pairs(&self) -> Vec<(usize, usize)>` | Sorted list of pairs (i < j) |
| `is_paired(&self, i: usize) -> bool` | Whether position i is paired |
| `partner(&self, i: usize) -> Option<usize>` | Pairing partner of position i |
| `num_pairs(&self) -> usize` | Number of base pairs |

**PartitionResult methods:**

| Method | Description |
|--------|-------------|
| `pair_probability(&self, i: usize, j: usize) -> f64` | Probability that i and j are paired |
| `unpaired_probability(&self, i: usize) -> f64` | Probability that i is unpaired |

**Prediction functions:**

| Function | Description |
|----------|-------------|
| `nussinov(seq: &[u8], min_loop_size: usize) -> Result<NussinovResult>` | Maximize base pair count (O(n³) DP) |
| `zuker_mfe(seq: &[u8]) -> Result<MfeResult>` | Minimum free energy with Turner nearest-neighbor parameters |
| `mccaskill(seq: &[u8], temperature: f64) -> Result<PartitionResult>` | Base pair probabilities via inside-outside algorithm |

**Comparison functions:**

| Function | Description |
|----------|-------------|
| `base_pair_distance(a, b: &RnaSecondaryStructure) -> Result<usize>` | Symmetric difference of base pair sets |
| `mountain_distance(a, b: &RnaSecondaryStructure) -> Result<f64>` | L1 distance of mountain representations |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |
| `minhash` | No | MinHash/FracMinHash sketching |

## Dependencies

- `cyanea-core` -- traits, error types
- `needletail` -- high-performance FASTA/FASTQ parsing

## Tests

398 tests across 22 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 140 | Module declarations, re-exports |
| `alphabet.rs` | 94 | Alphabet trait and implementations |
| `types.rs` | 283 | ValidatedSeq generic type |
| `seq.rs` | 197 | Sequence-specific methods |
| `codon.rs` | 159 | Codon translation table |
| `kmer.rs` | 113 | K-mer iterator |
| `quality.rs` | 164 | Phred quality scores |
| `fasta.rs` | 98 | FASTA parsing |
| `fastq.rs` | 256 | FASTQ parsing |
| `twobit.rs` | 373 | 2-bit packed DNA encoding |
| `suffix.rs` | 635 | Suffix array (SA-IS algorithm) |
| `fm_index.rs` | 355 | FM-Index (BWT backward search) |
| `minhash.rs` | 833 | MinHash and FracMinHash sketching |
| `pattern.rs` | 631 | Exact and approximate pattern matching (7 algorithms) |
| `pssm.rs` | 370 | Position-Specific Scoring Matrix / motif scanning |
| `orf.rs` | 288 | Open reading frame finder |
| `fasta_index.rs` | 421 | FASTA indexed reader (.fai) |
| `fmd_index.rs` | 732 | Bidirectional FM-Index (FMD-Index) |
| `paired.rs` | 1084 | Paired-end FASTQ: types, parsing, writing, interleave/deinterleave |
| `trim.rs` | 1535 | Quality trimming, adapter removal, filtering, TrimPipeline, paired trimming |
| `bwt.rs` | 225 | Standalone BWT construction and inversion |
| `rna_structure.rs` | 1110 | RNA secondary structure: Nussinov, Zuker MFE, McCaskill, dot-bracket, comparison |
