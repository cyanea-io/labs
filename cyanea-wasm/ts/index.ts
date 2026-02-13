// ---------------------------------------------------------------------------
// @cyanea/bio -- Typed wrapper around cyanea-wasm WASM bindings
// ---------------------------------------------------------------------------
// Imports the raw WASM functions (which return JSON strings), parses the
// JSON envelope, and exposes a clean, fully typed TypeScript API.
//
// Usage:
//   import init from "../pkg/cyanea_wasm.js";
//   import { Seq, Align, Stats, ML, Chem, StructBio, Phylo, Core } from "./index.js";
//
//   await init();
//   const stats = Seq.parseFasta(">seq1\nACGT\n");
//   console.log(stats.gc_content);
// ---------------------------------------------------------------------------

import type {
  FastaStats,
  FastqRecord,
  CigarOp,
  CigarStats,
  AlignmentResult,
  AlignmentMode,
  SubstitutionMatrix,
  SeqPair,
  DescriptiveStats,
  TestResult,
  KmerCounts,
  UmapResult,
  DistanceMetric,
  MolecularProperties,
  Fingerprint,
  SubstructureResult,
  StructureInfo,
  SecondaryStructure,
  TreeInfo,
  RFDistance,
  DistanceModel,
  Alphabet,
  WasmResult,
} from "./types.js";

export type {
  FastaStats,
  FastqRecord,
  CigarOp,
  CigarStats,
  AlignmentResult,
  AlignmentMode,
  SubstitutionMatrix,
  SeqPair,
  DescriptiveStats,
  TestResult,
  KmerCounts,
  UmapResult,
  DistanceMetric,
  MolecularProperties,
  Fingerprint,
  SubstructureResult,
  StructureInfo,
  SecondaryStructure,
  TreeInfo,
  RFDistance,
  DistanceModel,
  Alphabet,
} from "./types.js";

// Re-export the envelope types for advanced users who want raw access.
export type { WasmOk, WasmErr, WasmResult } from "./types.js";

// ── Raw WASM imports ───────────────────────────────────────────────────────
// These are the raw wasm-bindgen functions that return JSON strings.
// The import path assumes wasm-pack output in ../pkg/.
// All imports are aliased with a _raw prefix to avoid shadowing by the
// namespace wrapper functions that share similar names.

import {
  parse_fasta as _raw_parse_fasta,
  parse_fastq as _raw_parse_fastq,
  gc_content_json as _raw_gc_content_json,
  reverse_complement as _raw_reverse_complement,
  transcribe as _raw_transcribe,
  translate as _raw_translate,
  validate as _raw_validate,
  align_dna as _raw_align_dna,
  align_dna_custom as _raw_align_dna_custom,
  align_protein as _raw_align_protein,
  align_batch as _raw_align_batch,
  parse_cigar as _raw_parse_cigar,
  validate_cigar as _raw_validate_cigar,
  cigar_stats as _raw_cigar_stats,
  cigar_to_alignment as _raw_cigar_to_alignment,
  alignment_to_cigar as _raw_alignment_to_cigar,
  generate_md_tag as _raw_generate_md_tag,
  merge_cigar as _raw_merge_cigar,
  reverse_cigar as _raw_reverse_cigar,
  collapse_cigar as _raw_collapse_cigar,
  hard_clip_to_soft as _raw_hard_clip_to_soft,
  split_cigar as _raw_split_cigar,
  describe as _raw_describe,
  pearson as _raw_pearson,
  spearman as _raw_spearman,
  t_test as _raw_t_test,
  t_test_two_sample as _raw_t_test_two_sample,
  mann_whitney_u as _raw_mann_whitney_u,
  bonferroni as _raw_bonferroni,
  benjamini_hochberg as _raw_benjamini_hochberg,
  kmer_count as _raw_kmer_count,
  euclidean_distance as _raw_euclidean_distance,
  manhattan_distance as _raw_manhattan_distance,
  hamming_distance as _raw_hamming_distance,
  cosine_similarity as _raw_cosine_similarity,
  umap as _raw_umap,
  smiles_properties as _raw_smiles_properties,
  canonical as _raw_canonical,
  smiles_fingerprint as _raw_smiles_fingerprint,
  tanimoto as _raw_tanimoto,
  smiles_substructure as _raw_smiles_substructure,
  pdb_info as _raw_pdb_info,
  pdb_secondary_structure as _raw_pdb_secondary_structure,
  rmsd as _raw_rmsd,
  newick_info as _raw_newick_info,
  evolutionary_distance as _raw_evolutionary_distance,
  build_upgma as _raw_build_upgma,
  build_nj as _raw_build_nj,
  rf_distance as _raw_rf_distance,
  sha256 as _raw_sha256,
  zstd_compress as _raw_zstd_compress,
  zstd_decompress as _raw_zstd_decompress,
} from "../pkg/cyanea_wasm.js";

// ── Error type ─────────────────────────────────────────────────────────────

/**
 * Error thrown when a WASM function returns an error envelope.
 * The `message` property contains the error string from the Rust side.
 */
export class CyaneaError extends Error {
  override readonly name = "CyaneaError";

  constructor(message: string) {
    super(message);
  }
}

// ── JSON envelope unwrapper ────────────────────────────────────────────────

/**
 * Parse a JSON envelope string and return the success value, or throw
 * a `CyaneaError` if the envelope contains an error.
 */
function unwrap<T>(json: string): T {
  const parsed: WasmResult<T> = JSON.parse(json);
  if ("error" in parsed) {
    throw new CyaneaError(parsed.error);
  }
  return parsed.ok;
}

// ── Seq ────────────────────────────────────────────────────────────────────

export namespace Seq {
  /**
   * Parse FASTA-formatted text and return summary statistics.
   *
   * @param data - FASTA-formatted string (e.g. ">seq1\nACGT\n")
   */
  export function parseFasta(data: string): FastaStats {
    return unwrap<FastaStats>(_raw_parse_fasta(data));
  }

  /**
   * Parse FASTQ-formatted text and return an array of records.
   *
   * @param data - FASTQ-formatted string
   */
  export function parseFastq(data: string): FastqRecord[] {
    return unwrap<FastqRecord[]>(_raw_parse_fastq(data));
  }

  /**
   * Compute GC content as a percentage (0-100) of a nucleotide string.
   *
   * @param seq - Raw nucleotide sequence string (e.g. "ACGT")
   */
  export function gcContent(seq: string): number {
    return unwrap<number>(_raw_gc_content_json(seq));
  }

  /**
   * Compute the reverse complement of a DNA sequence.
   *
   * @param seq - DNA sequence string
   */
  export function reverseComplement(seq: string): string {
    return unwrap<string>(_raw_reverse_complement(seq));
  }

  /**
   * Transcribe a DNA sequence to RNA (T -> U).
   *
   * @param seq - DNA sequence string
   */
  export function transcribe(seq: string): string {
    return unwrap<string>(_raw_transcribe(seq));
  }

  /**
   * Translate a DNA sequence to a protein string (standard codon table).
   *
   * @param seq - DNA sequence string (length must be divisible by 3)
   */
  export function translate(seq: string): string {
    return unwrap<string>(_raw_translate(seq));
  }

  /**
   * Validate a sequence against an alphabet.
   *
   * @param seq - Sequence string to validate
   * @param alphabet - One of "dna", "rna", or "protein"
   * @returns `true` if the sequence is valid
   * @throws CyaneaError if the sequence contains invalid characters
   */
  export function validate(seq: string, alphabet: Alphabet): boolean {
    return unwrap<boolean>(_raw_validate(seq, alphabet));
  }
}

// ── Align ──────────────────────────────────────────────────────────────────

export namespace Align {
  /**
   * Align two DNA sequences with default scoring (+2/-1/-5/-2).
   *
   * @param query - Query DNA sequence
   * @param target - Target DNA sequence
   * @param mode - Alignment mode: "local", "global", or "semiglobal"
   */
  export function alignDna(
    query: string,
    target: string,
    mode: AlignmentMode,
  ): AlignmentResult {
    return unwrap<AlignmentResult>(_raw_align_dna(query, target, mode));
  }

  /**
   * Align two DNA sequences with custom scoring parameters.
   *
   * @param query - Query DNA sequence
   * @param target - Target DNA sequence
   * @param mode - Alignment mode
   * @param matchScore - Score for matching bases (positive)
   * @param mismatchScore - Score for mismatching bases (negative)
   * @param gapOpen - Gap opening penalty (negative)
   * @param gapExtend - Gap extension penalty (negative)
   */
  export function alignDnaCustom(
    query: string,
    target: string,
    mode: AlignmentMode,
    matchScore: number,
    mismatchScore: number,
    gapOpen: number,
    gapExtend: number,
  ): AlignmentResult {
    return unwrap<AlignmentResult>(
      _raw_align_dna_custom(query, target, mode, matchScore, mismatchScore, gapOpen, gapExtend),
    );
  }

  /**
   * Align two protein sequences using a named substitution matrix.
   *
   * @param query - Query protein sequence
   * @param target - Target protein sequence
   * @param mode - Alignment mode
   * @param matrix - Substitution matrix name: "blosum62", "blosum45", "blosum80", or "pam250"
   */
  export function alignProtein(
    query: string,
    target: string,
    mode: AlignmentMode,
    matrix: SubstitutionMatrix,
  ): AlignmentResult {
    return unwrap<AlignmentResult>(_raw_align_protein(query, target, mode, matrix));
  }

  /**
   * Batch-align multiple sequence pairs with custom scoring.
   *
   * @param pairs - Array of query/target pairs
   * @param mode - Alignment mode
   * @param matchScore - Score for matching bases
   * @param mismatchScore - Score for mismatching bases
   * @param gapOpen - Gap opening penalty
   * @param gapExtend - Gap extension penalty
   */
  export function alignBatch(
    pairs: SeqPair[],
    mode: AlignmentMode,
    matchScore: number,
    mismatchScore: number,
    gapOpen: number,
    gapExtend: number,
  ): AlignmentResult[] {
    const pairsJson = JSON.stringify(pairs);
    return unwrap<AlignmentResult[]>(
      _raw_align_batch(pairsJson, mode, matchScore, mismatchScore, gapOpen, gapExtend),
    );
  }

  // -- CIGAR utilities --

  /**
   * Parse a SAM CIGAR string into an array of operations.
   *
   * Accepts the full SAM alphabet (M, I, D, N, S, H, P, =, X) and `*`.
   *
   * @param cigar - CIGAR string (e.g. "10M3I4D2S")
   */
  export function parseCigar(cigar: string): CigarOp[] {
    return unwrap<CigarOp[]>(_raw_parse_cigar(cigar));
  }

  /**
   * Validate a CIGAR string against SAM spec rules.
   *
   * Returns `true` if valid.
   * @throws CyaneaError if the CIGAR is invalid
   */
  export function validateCigar(cigar: string): boolean {
    return unwrap<boolean>(_raw_validate_cigar(cigar));
  }

  /**
   * Compute statistics from a CIGAR string.
   *
   * @param cigar - CIGAR string
   */
  export function cigarStats(cigar: string): CigarStats {
    return unwrap<CigarStats>(_raw_cigar_stats(cigar));
  }

  /**
   * Reconstruct gapped alignment from CIGAR and ungapped sequences.
   *
   * @param cigar - CIGAR string
   * @param query - Ungapped query sequence
   * @param target - Ungapped target/reference sequence
   * @returns Object with `aligned_query` and `aligned_target` byte arrays
   */
  export function cigarToAlignment(
    cigar: string,
    query: string,
    target: string,
  ): { aligned_query: number[]; aligned_target: number[] } {
    return unwrap<{ aligned_query: number[]; aligned_target: number[] }>(
      _raw_cigar_to_alignment(cigar, query, target),
    );
  }

  /**
   * Extract a CIGAR string from a gapped alignment (using =/X distinction).
   *
   * Both sequences must be the same length, with `-` for gaps.
   *
   * @param query - Gapped query sequence
   * @param target - Gapped target sequence
   * @returns CIGAR string
   */
  export function alignmentToCigar(query: string, target: string): string {
    return unwrap<string>(_raw_alignment_to_cigar(query, target));
  }

  /**
   * Generate a SAM MD:Z tag from CIGAR and ungapped sequences.
   *
   * @param cigar - CIGAR string
   * @param query - Ungapped query sequence
   * @param reference - Ungapped reference sequence
   */
  export function generateMdTag(cigar: string, query: string, reference: string): string {
    return unwrap<string>(_raw_generate_md_tag(cigar, query, reference));
  }

  /**
   * Merge adjacent same-type CIGAR operations.
   *
   * @param cigar - CIGAR string
   * @returns Merged CIGAR string
   */
  export function mergeCigar(cigar: string): string {
    return unwrap<string>(_raw_merge_cigar(cigar));
  }

  /**
   * Reverse CIGAR operation order.
   *
   * @param cigar - CIGAR string
   * @returns Reversed CIGAR string
   */
  export function reverseCigar(cigar: string): string {
    return unwrap<string>(_raw_reverse_cigar(cigar));
  }

  /**
   * Collapse =/X operations into M (alignment match).
   *
   * @param cigar - CIGAR string
   * @returns Collapsed CIGAR string
   */
  export function collapseCigar(cigar: string): string {
    return unwrap<string>(_raw_collapse_cigar(cigar));
  }

  /**
   * Convert hard clips (H) to soft clips (S).
   *
   * @param cigar - CIGAR string
   * @returns Converted CIGAR string
   */
  export function hardClipToSoft(cigar: string): string {
    return unwrap<string>(_raw_hard_clip_to_soft(cigar));
  }

  /**
   * Split CIGAR at a reference coordinate.
   *
   * @param cigar - CIGAR string
   * @param refPos - 0-based reference position to split at
   * @returns Object with `left` and `right` CIGAR strings
   */
  export function splitCigar(
    cigar: string,
    refPos: number,
  ): { left: string; right: string } {
    return unwrap<{ left: string; right: string }>(_raw_split_cigar(cigar, refPos));
  }
}

// ── Stats ──────────────────────────────────────────────────────────────────

export namespace Stats {
  /**
   * Compute descriptive statistics for a numeric dataset.
   *
   * @param data - Array of numbers
   */
  export function describe(data: number[]): DescriptiveStats {
    return unwrap<DescriptiveStats>(_raw_describe(JSON.stringify(data)));
  }

  /**
   * Compute Pearson correlation coefficient between two arrays.
   *
   * @param x - First array of numbers
   * @param y - Second array of numbers (same length as x)
   * @returns Pearson r in [-1, 1]
   */
  export function pearson(x: number[], y: number[]): number {
    return unwrap<number>(_raw_pearson(JSON.stringify(x), JSON.stringify(y)));
  }

  /**
   * Compute Spearman rank correlation coefficient between two arrays.
   *
   * @param x - First array of numbers
   * @param y - Second array of numbers (same length as x)
   * @returns Spearman rho in [-1, 1]
   */
  export function spearman(x: number[], y: number[]): number {
    return unwrap<number>(_raw_spearman(JSON.stringify(x), JSON.stringify(y)));
  }

  /**
   * One-sample t-test.
   *
   * @param data - Array of observations
   * @param mu - Hypothesized population mean
   */
  export function tTest(data: number[], mu: number): TestResult {
    return unwrap<TestResult>(_raw_t_test(JSON.stringify(data), mu));
  }

  /**
   * Two-sample t-test (Student's or Welch's).
   *
   * @param x - First sample
   * @param y - Second sample
   * @param equalVar - If true, assume equal variances (Student's t); if false, use Welch's
   */
  export function tTestTwoSample(
    x: number[],
    y: number[],
    equalVar: boolean,
  ): TestResult {
    return unwrap<TestResult>(
      _raw_t_test_two_sample(JSON.stringify(x), JSON.stringify(y), equalVar),
    );
  }

  /**
   * Mann-Whitney U test (non-parametric).
   *
   * @param x - First sample
   * @param y - Second sample
   */
  export function mannWhitneyU(x: number[], y: number[]): TestResult {
    return unwrap<TestResult>(
      _raw_mann_whitney_u(JSON.stringify(x), JSON.stringify(y)),
    );
  }

  /**
   * Bonferroni p-value correction for multiple comparisons.
   *
   * @param pValues - Array of raw p-values
   * @returns Corrected p-values (capped at 1.0)
   */
  export function bonferroni(pValues: number[]): number[] {
    return unwrap<number[]>(_raw_bonferroni(JSON.stringify(pValues)));
  }

  /**
   * Benjamini-Hochberg FDR correction for multiple comparisons.
   *
   * @param pValues - Array of raw p-values
   * @returns Adjusted p-values
   */
  export function benjaminiHochberg(pValues: number[]): number[] {
    return unwrap<number[]>(_raw_benjamini_hochberg(JSON.stringify(pValues)));
  }
}

// ── ML ─────────────────────────────────────────────────────────────────────

export namespace ML {
  /**
   * Count k-mers in a nucleotide or protein sequence.
   *
   * @param seq - Sequence string
   * @param k - K-mer length (must be >= 1)
   */
  export function kmerCount(seq: string, k: number): KmerCounts {
    return unwrap<KmerCounts>(_raw_kmer_count(seq, k));
  }

  /**
   * Euclidean distance between two numeric vectors.
   *
   * @param a - First vector
   * @param b - Second vector (same length as a)
   */
  export function euclideanDistance(a: number[], b: number[]): number {
    return unwrap<number>(
      _raw_euclidean_distance(JSON.stringify(a), JSON.stringify(b)),
    );
  }

  /**
   * Manhattan (L1) distance between two numeric vectors.
   *
   * @param a - First vector
   * @param b - Second vector (same length as a)
   */
  export function manhattanDistance(a: number[], b: number[]): number {
    return unwrap<number>(
      _raw_manhattan_distance(JSON.stringify(a), JSON.stringify(b)),
    );
  }

  /**
   * Hamming distance between two strings (byte-level comparison).
   *
   * @param a - First string
   * @param b - Second string (same length as a)
   * @returns Number of positions where the strings differ
   */
  export function hammingDistance(a: string, b: string): number {
    return unwrap<number>(_raw_hamming_distance(a, b));
  }

  /**
   * Cosine similarity between two numeric vectors.
   *
   * @param a - First vector
   * @param b - Second vector (same length as a)
   * @returns Similarity in [0, 1] (or [-1, 1] if vectors have negative components)
   */
  export function cosineSimilarity(a: number[], b: number[]): number {
    return unwrap<number>(
      _raw_cosine_similarity(JSON.stringify(a), JSON.stringify(b)),
    );
  }

  /**
   * UMAP dimensionality reduction.
   *
   * @param data - Flat row-major matrix of input data
   * @param nFeatures - Number of features (columns) per sample
   * @param nComponents - Output dimensionality (typically 2 or 3)
   * @param nNeighbors - Number of nearest neighbors (default 15)
   * @param minDist - Minimum distance in embedding space (default 0.1)
   * @param nEpochs - Optimization epochs (default 200)
   * @param metric - Distance metric: "euclidean", "manhattan", or "cosine"
   */
  export function umap(
    data: number[],
    nFeatures: number,
    nComponents: number,
    nNeighbors: number,
    minDist: number,
    nEpochs: number,
    metric: DistanceMetric,
  ): UmapResult {
    return unwrap<UmapResult>(
      _raw_umap(JSON.stringify(data), nFeatures, nComponents, nNeighbors, minDist, nEpochs, metric),
    );
  }
}

// ── Chem ───────────────────────────────────────────────────────────────────

export namespace Chem {
  /**
   * Parse a SMILES string and compute molecular properties.
   *
   * @param smiles - SMILES notation string
   */
  export function properties(smiles: string): MolecularProperties {
    return unwrap<MolecularProperties>(_raw_smiles_properties(smiles));
  }

  /**
   * Generate canonical SMILES from an input SMILES string.
   * Canonicalization produces a unique representation for each molecule.
   *
   * @param smiles - Input SMILES string
   * @returns Canonical SMILES string
   */
  export function canonical(smiles: string): string {
    return unwrap<string>(_raw_canonical(smiles));
  }

  /**
   * Compute a Morgan (circular) fingerprint for a molecule.
   *
   * @param smiles - SMILES string
   * @param radius - Fingerprint radius (typically 2)
   * @param nBits - Number of bits in the fingerprint (typically 2048)
   */
  export function fingerprint(
    smiles: string,
    radius: number,
    nBits: number,
  ): Fingerprint {
    return unwrap<Fingerprint>(_raw_smiles_fingerprint(smiles, radius, nBits));
  }

  /**
   * Tanimoto similarity between two molecules (Morgan fingerprint r=2, 2048 bits).
   *
   * @param smiles1 - First molecule SMILES
   * @param smiles2 - Second molecule SMILES
   * @returns Similarity in [0, 1]
   */
  export function tanimoto(smiles1: string, smiles2: string): number {
    return unwrap<number>(_raw_tanimoto(smiles1, smiles2));
  }

  /**
   * Substructure search: check if a molecule contains a pattern.
   *
   * @param molecule - SMILES of the molecule to search in
   * @param pattern - SMILES of the substructure pattern
   */
  export function substructure(
    molecule: string,
    pattern: string,
  ): SubstructureResult {
    return unwrap<SubstructureResult>(_raw_smiles_substructure(molecule, pattern));
  }
}

// ── StructBio ──────────────────────────────────────────────────────────────

export namespace StructBio {
  /**
   * Parse PDB-formatted text and return structure summary information.
   *
   * @param pdbText - PDB file contents as a string
   */
  export function pdbInfo(pdbText: string): StructureInfo {
    return unwrap<StructureInfo>(_raw_pdb_info(pdbText));
  }

  /**
   * Assign secondary structure (DSSP-like) from PDB text.
   * Analyzes the first chain in the structure.
   *
   * @param pdbText - PDB file contents as a string
   */
  export function secondaryStructure(pdbText: string): SecondaryStructure {
    return unwrap<SecondaryStructure>(_raw_pdb_secondary_structure(pdbText));
  }

  /**
   * Compute RMSD between two sets of 3D coordinates.
   *
   * @param coords1 - First set of coordinates as [x, y, z] arrays
   * @param coords2 - Second set of coordinates (same length as coords1)
   * @returns RMSD value in the same units as the input coordinates
   */
  export function rmsd(
    coords1: [number, number, number][],
    coords2: [number, number, number][],
  ): number {
    return unwrap<number>(_raw_rmsd(JSON.stringify(coords1), JSON.stringify(coords2)));
  }
}

// ── Phylo ──────────────────────────────────────────────────────────────────

export namespace Phylo {
  /**
   * Parse a Newick-format tree string and return tree information.
   *
   * @param newick - Newick tree string (e.g. "((A:0.1,B:0.2):0.3,C:0.4);")
   */
  export function newickInfo(newick: string): TreeInfo {
    return unwrap<TreeInfo>(_raw_newick_info(newick));
  }

  /**
   * Compute evolutionary distance between two aligned sequences.
   *
   * @param seq1 - First sequence
   * @param seq2 - Second sequence (same length as seq1)
   * @param model - Distance model: "p" (p-distance), "jc" (Jukes-Cantor), or "k2p" (Kimura 2-parameter)
   */
  export function evolutionaryDistance(
    seq1: string,
    seq2: string,
    model: DistanceModel,
  ): number {
    return unwrap<number>(_raw_evolutionary_distance(seq1, seq2, model));
  }

  /**
   * Build a UPGMA tree from a distance matrix.
   *
   * @param labels - Array of taxon labels
   * @param matrix - Symmetric distance matrix (2D array, labels.length x labels.length)
   * @returns Newick string of the resulting tree
   */
  export function buildUpgma(
    labels: string[],
    matrix: number[][],
  ): string {
    return unwrap<string>(_raw_build_upgma(JSON.stringify(labels), JSON.stringify(matrix)));
  }

  /**
   * Build a Neighbor-Joining tree from a distance matrix.
   *
   * @param labels - Array of taxon labels
   * @param matrix - Symmetric distance matrix (2D array, labels.length x labels.length)
   * @returns Newick string of the resulting tree
   */
  export function buildNj(
    labels: string[],
    matrix: number[][],
  ): string {
    return unwrap<string>(_raw_build_nj(JSON.stringify(labels), JSON.stringify(matrix)));
  }

  /**
   * Compute Robinson-Foulds distance between two trees.
   *
   * @param newick1 - First Newick tree string
   * @param newick2 - Second Newick tree string (must have same leaf set)
   */
  export function rfDistance(
    newick1: string,
    newick2: string,
  ): RFDistance {
    return unwrap<RFDistance>(_raw_rf_distance(newick1, newick2));
  }
}

// ── Core Utilities ─────────────────────────────────────────────────────────

export namespace Core {
  /**
   * Compute SHA-256 hash of a string.
   *
   * @param data - Input string
   * @returns Hex-encoded SHA-256 hash
   */
  export function sha256(data: string): string {
    return unwrap<string>(_raw_sha256(data));
  }

  /**
   * Compress a string with zstd at the given compression level.
   *
   * @param data - Input string
   * @param level - Compression level (1-22, 3 is a good default)
   * @returns Compressed data as a byte array
   */
  export function zstdCompress(data: string, level: number): number[] {
    return unwrap<number[]>(_raw_zstd_compress(data, level));
  }

  /**
   * Decompress zstd-compressed data back to a string.
   *
   * @param data - Compressed byte array (as returned by zstdCompress)
   * @returns Decompressed string
   */
  export function zstdDecompress(data: number[]): string {
    return unwrap<string>(_raw_zstd_decompress(JSON.stringify(data)));
  }
}
