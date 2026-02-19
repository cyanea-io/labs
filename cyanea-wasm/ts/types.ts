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

// ── Seq Module (RNA / Protein / Reads / Assembly) ───────────────────────────

/** RNA secondary structure prediction result. */
export interface RnaStructure {
  structure: string;
  pairs: [number, number][];
  free_energy: number;
}

/** Protein sequence properties. */
export interface ProteinProperties {
  molecular_weight: number;
  isoelectric_point: number;
  gravy: number;
  instability_index: number;
  extinction_280_reduced: number;
  extinction_280_oxidized: number;
}

/** A simulated sequencing read. */
export interface SimulatedRead {
  name: string;
  sequence: string;
  quality: string;
  position: number;
}

/** Codon usage table. */
export interface CodonUsage {
  codons: Record<string, number>;
  total: number;
}

/** Assembly statistics for a set of contigs. */
export interface AssemblyStats {
  n_contigs: number;
  total_length: number;
  n50: number;
  l50: number;
  largest: number;
  gc_content: number;
}

// ── Align Module (MSA / POA / Banded) ─────────────────────────────────────

/** Multiple sequence alignment result. */
export interface MsaResult {
  alignment: string[];
  score: number;
  n_sequences: number;
  alignment_length: number;
}

/** Partial-order alignment consensus result. */
export interface PoaConsensus {
  consensus: string;
  n_sequences: number;
  n_nodes: number;
}

// ── Stats Module (Survival / Diversity / PopGen / NullModel) ──────────────

/** Kaplan-Meier survival step. */
export interface KmStep {
  time: number;
  survival: number;
  n_at_risk: number;
  n_events: number;
}

/** Kaplan-Meier survival analysis result. */
export interface KmResult {
  steps: KmStep[];
  median_survival: number | null;
  total_events: number;
  total_censored: number;
}

/** Log-rank test result. */
export interface LogRankResult {
  statistic: number;
  p_value: number;
  degrees_of_freedom: number;
}

/** Cox proportional hazards result. */
export interface CoxPhResult {
  coefficients: number[];
  standard_errors: number[];
  hazard_ratios: number[];
  log_likelihood: number;
}

/** Wright-Fisher simulation result. */
export interface WrightFisherResult {
  trajectory: number[];
  fixation_generation: number | null;
  final_frequency: number;
}

/** Tajima's D result. */
export interface TajimaD {
  d: number;
  pi: number;
  s: number;
}

/** Hudson's Fst result. */
export interface FstResult {
  fst: number;
  numerator: number;
  denominator: number;
}

// ── ML Module (Forest / GBDT / HMM / Metrics / CV) ──────────────────────

/** Random forest classification result. */
export interface RandomForestResult {
  predictions: number[];
  oob_error: number | null;
  n_trees: number;
  feature_importances: number[];
}

/** Gradient-boosted tree regression result. */
export interface GbdtRegressionResult {
  predictions: number[];
  n_trees: number;
  feature_importances: number[];
  training_loss: number[];
}

/** Gradient-boosted tree classification result. */
export interface GbdtClassifyResult {
  predictions: number[];
  probabilities: number[];
  n_trees: number;
  feature_importances: number[];
  training_loss: number[];
}

/** HMM Viterbi decoding result. */
export interface HmmViterbiResult {
  path: number[];
  log_probability: number;
  n_states: number;
  n_observations: number;
}

/** Confusion matrix from classification evaluation. */
export interface ConfusionMatrix {
  matrix: number[][];
  labels: number[];
  accuracy: number;
  precision_per_class: number[];
  recall_per_class: number[];
  f1_per_class: number[];
}

/** ROC curve point. */
export interface RocPoint {
  fpr: number;
  tpr: number;
  threshold: number;
}

/** ROC curve result. */
export interface RocCurve {
  points: RocPoint[];
  auc: number;
}

/** Precision-recall curve point. */
export interface PrPoint {
  precision: number;
  recall: number;
  threshold: number;
}

/** Precision-recall curve result. */
export interface PrCurve {
  points: PrPoint[];
  auc: number;
}

/** Cross-validation result. */
export interface CvResult {
  fold_scores: number[];
  mean: number;
  std_dev: number;
}

/** Feature selection result from variance threshold. */
export interface FeatureSelection {
  selected_indices: number[];
  n_selected: number;
  n_original: number;
}

// ── Chem Module (SDF / MACCS) ─────────────────────────────────────────────

/** SDF-parsed molecule. */
export interface SdfMolecule {
  name: string;
  smiles: string;
  n_atoms: number;
  n_bonds: number;
  properties: Record<string, string>;
}

/** MACCS 166-key fingerprint result. */
export interface MaccsFingerprint {
  bits: number[];
  n_bits: number;
  n_on_bits: number;
}

// ── StructBio Module (Contact Map / Ramachandran / mmCIF / Kabsch) ──────

/** Contact map from PDB structure. */
export interface ContactMap {
  chain_id: string;
  size: number;
  n_contacts: number;
  contact_density: number;
  contacts: [number, number][];
}

/** Ramachandran plot entry. */
export interface RamachandranEntry {
  residue_index: number;
  residue_name: string;
  phi: number;
  psi: number;
  region: string;
}

/** mmCIF structure info. */
export interface MmcifInfo {
  id: string;
  n_chains: number;
  n_residues: number;
  n_atoms: number;
}

/** Kabsch superposition result. */
export interface KabschResult {
  rmsd: number;
  rotation: number[];
  translation: number[];
}

// ── Phylo Module (NEXUS / Simulation) ─────────────────────────────────────

/** NEXUS file parse result. */
export interface NexusFile {
  taxa: string[];
  trees: NamedTree[];
}

/** A named tree from a NEXUS file. */
export interface NamedTree {
  name: string;
  newick: string;
}

/** Simulated sequence alignment from tree evolution. */
export interface SimulatedAlignment {
  names: string[];
  sequences: string[];
}

/** Coalescent tree simulation result. */
export interface CoalescentTree {
  newick: string;
  n_samples: number;
  tmrca: number;
}

// ── IO Module (VCF / BED / GFF3 / BLAST XML / bedGraph / GFA / Fetch) ──

/** Parsed VCF variant. */
export interface VcfVariant {
  chrom: string;
  pos: number;
  id: string;
  ref_allele: string;
  alt_allele: string;
  qual: number | null;
  filter: string;
}

/** Parsed BED record. */
export interface BedRecord {
  chrom: string;
  start: number;
  end: number;
  name: string | null;
  score: number | null;
  strand: string;
}

/** Parsed GFF3 gene. */
export interface Gff3Gene {
  id: string;
  chrom: string;
  start: number;
  end: number;
  strand: string;
  gene_type: string;
  name: string | null;
}

/** BLAST XML hit. */
export interface BlastXmlHit {
  id: string;
  definition: string;
  accession: string;
  length: number;
  bit_score: number;
  evalue: number;
  identity: number;
  query_from: number;
  query_to: number;
  hit_from: number;
  hit_to: number;
}

/** BLAST XML parse result. */
export interface BlastXmlResult {
  program: string;
  query: string;
  hits: BlastXmlHit[];
}

/** bedGraph record. */
export interface BedGraphRecord {
  chrom: string;
  start: number;
  end: number;
  value: number;
}

/** GFA assembly graph segment. */
export interface GfaSegment {
  name: string;
  sequence: string;
  length: number;
}

/** GFA assembly graph. */
export interface GfaGraph {
  segments: GfaSegment[];
  n_links: number;
  n_paths: number;
  total_sequence_length: number;
}

// ── Omics Module (Intervals / VEP / CNV / Methylation / Spatial) ──────

/** Genomic interval. */
export interface GenomicInterval {
  chrom: string;
  start: number;
  end: number;
  strand: string;
}

/** Variant effect prediction result. */
export interface VariantEffect {
  consequence: string;
  gene_id: string;
  gene_name: string | null;
  distance: number | null;
}

/** CNV segment from CBS. */
export interface CnvSegment {
  chrom: string;
  start: number;
  end: number;
  mean_value: number;
  n_probes: number;
}

/** CpG island. */
export interface CpgIsland {
  chrom: string;
  start: number;
  end: number;
  cpg_count: number;
  gc_content: number;
  observed_expected: number;
}

/** Spatial autocorrelation result (Moran's I). */
export interface SpatialAutocorrelation {
  statistic: number;
  expected: number;
  variance: number;
  z_score: number;
  p_value: number;
}

/** Geary's C spatial autocorrelation result. */
export interface GearysC {
  statistic: number;
  expected: number;
  variance: number;
  z_score: number;
  p_value: number;
}

/** Liftover result. */
export interface LiftoverResult {
  chrom: string;
  start: number;
  end: number;
  success: boolean;
}

/** Jaccard similarity between interval sets. */
export interface JaccardResult {
  jaccard: number;
  intersection_bases: number;
  union_bases: number;
}

/** Closest interval result. */
export interface ClosestResult {
  query: GenomicInterval;
  closest: GenomicInterval;
  distance: number;
}

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
