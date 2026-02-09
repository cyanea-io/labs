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

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |

## Dependencies

- `cyanea-core` -- traits, error types
- `needletail` -- high-performance FASTA/FASTQ parsing

## Tests

43 unit tests + 1 doc test across 8 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 68 | Module declarations, re-exports |
| `alphabet.rs` | 94 | Alphabet trait and implementations |
| `types.rs` | 283 | ValidatedSeq generic type |
| `seq.rs` | 197 | Sequence-specific methods |
| `codon.rs` | 159 | Codon translation table |
| `kmer.rs` | 113 | K-mer iterator |
| `quality.rs` | 164 | Phred quality scores |
| `fasta.rs` | 98 | FASTA parsing |
| `fastq.rs` | 256 | FASTQ parsing |
