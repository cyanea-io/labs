// ---------------------------------------------------------------------------
// @cyanea/bio -- TypeScript type definitions for the cyanea-wasm WASM bindings
// ---------------------------------------------------------------------------
// Every raw WASM function returns a JSON string with either:
//   { "ok": T }   on success
//   { "error": string }  on failure
//
// The typed wrapper in index.ts unwraps this envelope so callers receive T
// directly or get a thrown CyaneaError.
// ---------------------------------------------------------------------------

// ── JSON Envelope ──────────────────────────────────────────────────────────

/** Success envelope returned by raw WASM functions. */
export interface WasmOk<T> {
  ok: T;
}

/** Error envelope returned by raw WASM functions. */
export interface WasmErr {
  error: string;
}

/** Union of success and error envelopes. */
export type WasmResult<T> = WasmOk<T> | WasmErr;

// ── Seq Module ─────────────────────────────────────────────────────────────

/** Summary statistics from parsing a FASTA string. */
export interface FastaStats {
  sequence_count: number;
  total_bases: number;
  gc_content: number;
  avg_length: number;
}

/** A single FASTQ record. */
export interface FastqRecord {
  name: string;
  sequence: string;
  quality: string;
}

/** A paired FASTQ record (R1 + R2). */
export interface PairedFastqRecord {
  r1: FastqRecord;
  r2: FastqRecord;
}

/** Mate validation mode for paired FASTQ parsing. */
export type MateValidation = "strict" | "relaxed" | "none";

/** How to handle orphan reads when one mate fails trimming. */
export type OrphanPolicy = "drop_both" | "keep_first" | "keep_second";

/** Configuration for read trimming. */
export interface TrimConfig {
  min_quality?: number;
  window_size?: number;
  min_length?: number;
  max_length?: number;
  adapters?: string[];
}

/** Statistics from paired-end trimming. */
export interface PairedTrimStats {
  total_input: number;
  both_passed: number;
  r1_only_passed: number;
  r2_only_passed: number;
  both_failed: number;
  survival_rate: number;
}

/** Result of paired-end trimming. */
export interface PairedTrimResult {
  pairs: PairedFastqRecord[];
  stats: PairedTrimStats;
}

// ── Align Module ───────────────────────────────────────────────────────────

/**
 * CIGAR operation describing how aligned sequences relate.
 *
 * Serde serialization of the Rust enum produces tagged objects like:
 *   { "Match": 4 }      — SAM '=' (sequence match)
 *   { "Mismatch": 1 }   — SAM 'X' (sequence mismatch)
 *   { "Insertion": 2 }  — SAM 'I' (insertion in query)
 *   { "Deletion": 1 }   — SAM 'D' (deletion from query)
 *   { "AlnMatch": 10 }  — SAM 'M' (alignment match, ambiguous)
 *   { "Skip": 100 }     — SAM 'N' (skipped region / intron)
 *   { "SoftClip": 5 }   — SAM 'S' (soft clipping)
 *   { "HardClip": 3 }   — SAM 'H' (hard clipping)
 *   { "Padding": 1 }    — SAM 'P' (silent deletion from padded ref)
 */
export type CigarOp =
  | { Match: number }
  | { Mismatch: number }
  | { Insertion: number }
  | { Deletion: number }
  | { AlnMatch: number }
  | { Skip: number }
  | { SoftClip: number }
  | { HardClip: number }
  | { Padding: number };

/** Statistics computed from a CIGAR string. */
export interface CigarStats {
  /** Compact CIGAR string representation. */
  cigar_string: string;
  /** Reference bases consumed (M, =, X, D, N). */
  reference_consumed: number;
  /** Query bases consumed (M, =, X, I, S). */
  query_consumed: number;
  /** Aligned columns (M, =, X, I, D). */
  alignment_columns: number;
  /** Sequence identity: = / (M + = + X + I + D). */
  identity: number;
  /** Number of gap openings (I or D runs). */
  gap_count: number;
  /** Total bases in gaps (I + D). */
  gap_bases: number;
  /** Soft-clipped bases. */
  soft_clipped: number;
  /** Hard-clipped bases. */
  hard_clipped: number;
}

/** Result of a pairwise sequence alignment. */
export interface AlignmentResult {
  /** Alignment score. */
  score: number;
  /** Aligned query sequence bytes (with 45 = '-' for gaps). */
  aligned_query: number[];
  /** Aligned target sequence bytes (with 45 = '-' for gaps). */
  aligned_target: number[];
  /** Start position in the original query (0-based, inclusive). */
  query_start: number;
  /** End position in the original query (0-based, exclusive). */
  query_end: number;
  /** Start position in the original target (0-based, inclusive). */
  target_start: number;
  /** End position in the original target (0-based, exclusive). */
  target_end: number;
  /** CIGAR operations describing the alignment. */
  cigar: CigarOp[];
}

/** Alignment mode string accepted by alignment functions. */
export type AlignmentMode = "local" | "global" | "semiglobal" | "semi-global";

/** Substitution matrix name accepted by protein alignment. */
export type SubstitutionMatrix = "blosum62" | "blosum45" | "blosum80" | "pam250";

/** A query/target pair for batch alignment. */
export interface SeqPair {
  query: string;
  target: string;
}

// ── Stats Module ───────────────────────────────────────────────────────────

/** Descriptive statistics for a numeric dataset. */
export interface DescriptiveStats {
  count: number;
  mean: number;
  median: number;
  variance: number;
  sample_variance: number;
  std_dev: number;
  sample_std_dev: number;
  min: number;
  max: number;
  range: number;
  q1: number;
  q3: number;
  iqr: number;
  skewness: number;
  kurtosis: number;
}

/** Result of a statistical hypothesis test. */
export interface TestResult {
  statistic: number;
  p_value: number;
  degrees_of_freedom: number | null;
  method: string;
}

// ── ML Module ──────────────────────────────────────────────────────────────

/** K-mer count result. */
export interface KmerCounts {
  k: number;
  total: number;
  distinct: number;
  counts: Record<string, number>;
}

/** UMAP dimensionality reduction result. */
export interface UmapResult {
  /** Flat row-major embedding coordinates. */
  embedding: number[];
  n_samples: number;
  n_components: number;
  n_epochs: number;
}

/** PCA dimensionality reduction result. */
export interface PcaResult {
  components: number[];
  explained_variance: number[];
  explained_variance_ratio: number[];
  transformed: number[];
  mean: number[];
  n_features: number;
  n_components: number;
}

/** t-SNE dimensionality reduction result. */
export interface TsneResult {
  embedding: number[];
  n_samples: number;
  n_components: number;
  kl_divergence: number;
}

/** K-means clustering result. */
export interface KmeansResult {
  centroids: number[];
  labels: number[];
  inertia: number;
  n_iter: number;
  n_features: number;
}

/** Distance metric accepted by UMAP/PCA/t-SNE/K-means. */
export type DistanceMetric = "euclidean" | "manhattan" | "cosine";

// ── Seq Module (MinHash) ──────────────────────────────────────────────────

/** MinHash sketch result. */
export interface MinHashSketch {
  k: number;
  sketch_size: number;
  num_hashes: number;
  hashes: number[];
}

/** MinHash comparison result. */
export interface MinHashComparison {
  jaccard: number;
  containment_a_in_b: number;
  containment_b_in_a: number;
  ani: number;
}

// ── IO Module ─────────────────────────────────────────────────────────────

/** A single pileup column. */
export interface PileupColumn {
  pos: number;
  ref_base: string;
  depth: number;
  base_counts: Record<string, number>;
}

/** Pileup for a single reference sequence. */
export interface Pileup {
  rname: string;
  columns: PileupColumn[];
}

/** Depth statistics for a reference sequence. */
export interface DepthStats {
  rname: string;
  length: number;
  covered: number;
  breadth: number;
  min_depth: number;
  max_depth: number;
  mean_depth: number;
  median_depth: number;
}

// ── Chem Module ────────────────────────────────────────────────────────────

/** Molecular properties computed from a SMILES string. */
export interface MolecularProperties {
  formula: string;
  molecular_weight: number;
  atom_count: number;
  bond_count: number;
  ring_count: number;
  rotatable_bonds: number;
  /** Hydrogen bond donors. */
  hbd: number;
  /** Hydrogen bond acceptors. */
  hba: number;
}

/** Morgan fingerprint result. */
export interface Fingerprint {
  /** Indices of set bits. */
  bits: number[];
  /** Total number of bits in the fingerprint. */
  n_bits: number;
  /** Number of bits that are set. */
  n_on_bits: number;
}

/** Substructure search result. */
export interface SubstructureResult {
  has_match: boolean;
  match_count: number;
}

// ── Struct Bio Module ──────────────────────────────────────────────────────

/** PDB structure summary. */
export interface StructureInfo {
  id: string;
  chain_count: number;
  residue_count: number;
  atom_count: number;
  chains: ChainInfo[];
}

/** Per-chain summary within a PDB structure. */
export interface ChainInfo {
  id: string;
  residue_count: number;
  atom_count: number;
}

/** Secondary structure assignment for a chain. */
export interface SecondaryStructure {
  chain_id: string;
  assignments: SSAssignment[];
}

/** Per-residue secondary structure assignment. */
export interface SSAssignment {
  residue_num: number;
  residue_name: string;
  structure: "Helix" | "Sheet" | "Turn" | "Coil";
}

// ── Phylo Module ───────────────────────────────────────────────────────────

/** Information about a parsed phylogenetic tree. */
export interface TreeInfo {
  leaf_count: number;
  internal_count: number;
  total_nodes: number;
  leaf_names: string[];
  newick: string;
}

/** Robinson-Foulds distance between two trees. */
export interface RFDistance {
  distance: number;
  normalized: number;
}

/** Evolutionary distance model accepted by evolutionary_distance. */
export type DistanceModel = "p" | "jc" | "k2p";

// ── Sequence Alphabet ──────────────────────────────────────────────────────

/** Alphabet name accepted by the validate function. */
export type Alphabet = "dna" | "rna" | "protein";

// ── Worker Message Protocol ────────────────────────────────────────────────

/** Unique request identifier for correlating Worker messages. */
export type RequestId = string;

/** Main thread → Worker: call a function. */
export interface WorkerRequest {
  type: "call";
  id: RequestId;
  method: string;
  args: unknown[];
}

/** Main thread → Worker: cancel an in-progress call. */
export interface WorkerCancel {
  type: "cancel";
  id: RequestId;
}

/** Main thread → Worker: initialize WASM (sent automatically). */
export interface WorkerInit {
  type: "init";
}

/** Worker → Main thread: function returned successfully. */
export interface WorkerResponseOk {
  type: "result";
  id: RequestId;
  ok: unknown;
}

/** Worker → Main thread: function threw an error. */
export interface WorkerResponseErr {
  type: "result";
  id: RequestId;
  error: string;
}

/** Worker → Main thread: progress update for a long-running call. */
export interface WorkerProgress {
  type: "progress";
  id: RequestId;
  current: number;
  total: number;
  message?: string;
}

/** Worker → Main thread: WASM loaded and ready. */
export interface WorkerReady {
  type: "ready";
}

/** Worker → Main thread: WASM initialization failed. */
export interface WorkerInitError {
  type: "init_error";
  error: string;
}

/** Union of all messages sent from the main thread to the Worker. */
export type MainToWorkerMessage = WorkerRequest | WorkerCancel | WorkerInit;

/** Union of all messages sent from the Worker to the main thread. */
export type WorkerToMainMessage =
  | WorkerResponseOk
  | WorkerResponseErr
  | WorkerProgress
  | WorkerReady
  | WorkerInitError;

/** Progress info passed to onProgress callbacks. */
export interface ProgressInfo {
  current: number;
  total: number;
  message?: string;
}

/** Options for Worker calls that support progress and cancellation. */
export interface CallOptions {
  signal?: AbortSignal;
  onProgress?: (info: ProgressInfo) => void;
}
