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

// ── Align Module ───────────────────────────────────────────────────────────

/**
 * CIGAR operation describing how aligned sequences relate.
 *
 * Serde serialization of the Rust enum produces tagged objects like:
 *   { "Match": 4 }
 *   { "Mismatch": 1 }
 *   { "Insertion": 2 }
 *   { "Deletion": 1 }
 */
export type CigarOp =
  | { Match: number }
  | { Mismatch: number }
  | { Insertion: number }
  | { Deletion: number };

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

/** Distance metric accepted by UMAP. */
export type DistanceMetric = "euclidean" | "manhattan" | "cosine";

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
