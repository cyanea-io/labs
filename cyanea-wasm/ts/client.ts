// ---------------------------------------------------------------------------
// @cyanea/bio -- CyaneaWorker client
// ---------------------------------------------------------------------------
// Provides an async, off-main-thread API that mirrors the synchronous
// namespace API from index.ts.  Each method posts a message to a Worker
// running worker.ts and returns a Promise.
//
// Usage:
//   import { CyaneaWorker } from "@cyanea/bio/client";
//   const worker = new CyaneaWorker();
//   await worker.ready();
//   const stats = await worker.Seq.parseFasta(">seq1\nACGT\n");
//   worker.terminate();
// ---------------------------------------------------------------------------

import type {
  WorkerToMainMessage,
  MainToWorkerMessage,
  FastaStats,
  FastqRecord,
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
  CallOptions,
  ProgressInfo,
} from "./types.js";

// ── Pending call tracking ──────────────────────────────────────────────────

interface PendingCall {
  resolve: (value: unknown) => void;
  reject: (reason: Error) => void;
  onProgress?: (info: ProgressInfo) => void;
}

// ── CyaneaWorker ───────────────────────────────────────────────────────────

export class CyaneaWorker {
  private readonly worker: Worker;
  private readonly pending = new Map<string, PendingCall>();
  private readyPromise: Promise<void>;
  private readyResolve!: () => void;
  private readyReject!: (e: Error) => void;

  /**
   * Create a new CyaneaWorker.
   *
   * @param workerUrl - URL to the worker script. Defaults to resolving
   *   `@cyanea/bio/worker` relative to this module.
   */
  constructor(workerUrl?: URL | string) {
    const url = workerUrl ?? new URL("./worker.js", import.meta.url);
    this.worker = new Worker(url, { type: "module" });

    this.readyPromise = new Promise<void>((resolve, reject) => {
      this.readyResolve = resolve;
      this.readyReject = reject;
    });

    this.worker.addEventListener("message", (event: MessageEvent<WorkerToMainMessage>) => {
      this.handleMessage(event.data);
    });

    this.worker.addEventListener("error", (event: ErrorEvent) => {
      this.readyReject(new Error(event.message));
    });
  }

  /** Wait for WASM to be loaded and the worker to be ready. */
  ready(): Promise<void> {
    return this.readyPromise;
  }

  /** Terminate the worker and reject all pending calls. */
  terminate(): void {
    this.worker.terminate();
    for (const [, pending] of this.pending) {
      pending.reject(new Error("Worker terminated"));
    }
    this.pending.clear();
  }

  // ── Internal message handling ────────────────────────────────────────────

  private handleMessage(msg: WorkerToMainMessage): void {
    switch (msg.type) {
      case "ready":
        this.readyResolve();
        break;

      case "init_error":
        this.readyReject(new Error(msg.error));
        break;

      case "result": {
        const pending = this.pending.get(msg.id);
        if (!pending) return;
        this.pending.delete(msg.id);
        if ("error" in msg) {
          pending.reject(new Error(msg.error));
        } else {
          pending.resolve(msg.ok);
        }
        break;
      }

      case "progress": {
        const pending = this.pending.get(msg.id);
        if (pending?.onProgress) {
          pending.onProgress({
            current: msg.current,
            total: msg.total,
            message: msg.message,
          });
        }
        break;
      }
    }
  }

  // ── Internal call dispatch ───────────────────────────────────────────────

  private call<T>(method: string, args: unknown[], options?: CallOptions): Promise<T> {
    const id = crypto.randomUUID();

    const promise = new Promise<T>((resolve, reject) => {
      this.pending.set(id, {
        resolve: resolve as (value: unknown) => void,
        reject,
        onProgress: options?.onProgress,
      });

      const msg: MainToWorkerMessage = { type: "call", id, method, args };
      this.worker.postMessage(msg);
    });

    // Wire up AbortSignal for cancellation
    if (options?.signal) {
      const signal = options.signal;
      if (signal.aborted) {
        this.pending.delete(id);
        const cancelMsg: MainToWorkerMessage = { type: "cancel", id };
        this.worker.postMessage(cancelMsg);
        return Promise.reject(new Error("Aborted"));
      }
      const onAbort = () => {
        const cancelMsg: MainToWorkerMessage = { type: "cancel", id };
        this.worker.postMessage(cancelMsg);
        const pending = this.pending.get(id);
        if (pending) {
          this.pending.delete(id);
          pending.reject(new Error("Aborted"));
        }
      };
      signal.addEventListener("abort", onAbort, { once: true });
      // Clean up listener when promise settles
      void promise.finally(() => signal.removeEventListener("abort", onAbort));
    }

    return promise;
  }

  // ── Namespace proxies ────────────────────────────────────────────────────

  /** Sequence parsing and manipulation. */
  readonly Seq = {
    parseFasta: (data: string) =>
      this.call<FastaStats>("Seq.parseFasta", [data]),

    parseFastq: (data: string) =>
      this.call<FastqRecord[]>("Seq.parseFastq", [data]),

    gcContent: (seq: string) =>
      this.call<number>("Seq.gcContent", [seq]),

    reverseComplement: (seq: string) =>
      this.call<string>("Seq.reverseComplement", [seq]),

    transcribe: (seq: string) =>
      this.call<string>("Seq.transcribe", [seq]),

    translate: (seq: string) =>
      this.call<string>("Seq.translate", [seq]),

    validate: (seq: string, alphabet: Alphabet) =>
      this.call<boolean>("Seq.validate", [seq, alphabet]),
  };

  /** Sequence alignment. */
  readonly Align = {
    alignDna: (query: string, target: string, mode: AlignmentMode) =>
      this.call<AlignmentResult>("Align.alignDna", [query, target, mode]),

    alignDnaCustom: (
      query: string, target: string, mode: AlignmentMode,
      matchScore: number, mismatchScore: number, gapOpen: number, gapExtend: number,
    ) =>
      this.call<AlignmentResult>("Align.alignDnaCustom", [
        query, target, mode, matchScore, mismatchScore, gapOpen, gapExtend,
      ]),

    alignProtein: (query: string, target: string, mode: AlignmentMode, matrix: SubstitutionMatrix) =>
      this.call<AlignmentResult>("Align.alignProtein", [query, target, mode, matrix]),

    alignBatch: (
      pairs: SeqPair[], mode: AlignmentMode,
      matchScore: number, mismatchScore: number, gapOpen: number, gapExtend: number,
      options?: CallOptions,
    ) =>
      this.call<AlignmentResult[]>(
        "Align.alignBatch",
        [pairs, mode, matchScore, mismatchScore, gapOpen, gapExtend],
        options,
      ),
  };

  /** Statistical methods. */
  readonly Stats = {
    describe: (data: number[]) =>
      this.call<DescriptiveStats>("Stats.describe", [data]),

    pearson: (x: number[], y: number[]) =>
      this.call<number>("Stats.pearson", [x, y]),

    spearman: (x: number[], y: number[]) =>
      this.call<number>("Stats.spearman", [x, y]),

    tTest: (data: number[], mu: number) =>
      this.call<TestResult>("Stats.tTest", [data, mu]),

    tTestTwoSample: (x: number[], y: number[], equalVar: boolean) =>
      this.call<TestResult>("Stats.tTestTwoSample", [x, y, equalVar]),

    mannWhitneyU: (x: number[], y: number[]) =>
      this.call<TestResult>("Stats.mannWhitneyU", [x, y]),

    bonferroni: (pValues: number[]) =>
      this.call<number[]>("Stats.bonferroni", [pValues]),

    benjaminiHochberg: (pValues: number[]) =>
      this.call<number[]>("Stats.benjaminiHochberg", [pValues]),
  };

  /** Machine learning primitives. */
  readonly ML = {
    kmerCount: (seq: string, k: number) =>
      this.call<KmerCounts>("ML.kmerCount", [seq, k]),

    euclideanDistance: (a: number[], b: number[]) =>
      this.call<number>("ML.euclideanDistance", [a, b]),

    manhattanDistance: (a: number[], b: number[]) =>
      this.call<number>("ML.manhattanDistance", [a, b]),

    hammingDistance: (a: string, b: string) =>
      this.call<number>("ML.hammingDistance", [a, b]),

    cosineSimilarity: (a: number[], b: number[]) =>
      this.call<number>("ML.cosineSimilarity", [a, b]),

    umap: (
      data: number[], nFeatures: number, nComponents: number,
      nNeighbors: number, minDist: number, nEpochs: number,
      metric: DistanceMetric, options?: CallOptions,
    ) =>
      this.call<UmapResult>(
        "ML.umap",
        [data, nFeatures, nComponents, nNeighbors, minDist, nEpochs, metric],
        options,
      ),
  };

  /** Chemistry / small molecules. */
  readonly Chem = {
    properties: (smiles: string) =>
      this.call<MolecularProperties>("Chem.properties", [smiles]),

    canonical: (smiles: string) =>
      this.call<string>("Chem.canonical", [smiles]),

    fingerprint: (smiles: string, radius: number, nBits: number) =>
      this.call<Fingerprint>("Chem.fingerprint", [smiles, radius, nBits]),

    tanimoto: (smiles1: string, smiles2: string) =>
      this.call<number>("Chem.tanimoto", [smiles1, smiles2]),

    substructure: (molecule: string, pattern: string) =>
      this.call<SubstructureResult>("Chem.substructure", [molecule, pattern]),
  };

  /** Structural biology. */
  readonly StructBio = {
    pdbInfo: (pdbText: string) =>
      this.call<StructureInfo>("StructBio.pdbInfo", [pdbText]),

    secondaryStructure: (pdbText: string) =>
      this.call<SecondaryStructure>("StructBio.secondaryStructure", [pdbText]),

    rmsd: (coords1: [number, number, number][], coords2: [number, number, number][]) =>
      this.call<number>("StructBio.rmsd", [coords1, coords2]),
  };

  /** Phylogenetics. */
  readonly Phylo = {
    newickInfo: (newick: string) =>
      this.call<TreeInfo>("Phylo.newickInfo", [newick]),

    evolutionaryDistance: (seq1: string, seq2: string, model: DistanceModel) =>
      this.call<number>("Phylo.evolutionaryDistance", [seq1, seq2, model]),

    buildUpgma: (labels: string[], matrix: number[][]) =>
      this.call<string>("Phylo.buildUpgma", [labels, matrix]),

    buildNj: (labels: string[], matrix: number[][]) =>
      this.call<string>("Phylo.buildNj", [labels, matrix]),

    rfDistance: (newick1: string, newick2: string) =>
      this.call<RFDistance>("Phylo.rfDistance", [newick1, newick2]),
  };

  /** Core utilities (hashing, compression). */
  readonly Core = {
    sha256: (data: string) =>
      this.call<string>("Core.sha256", [data]),

    zstdCompress: (data: string, level: number) =>
      this.call<number[]>("Core.zstdCompress", [data, level]),

    zstdDecompress: (data: number[]) =>
      this.call<string>("Core.zstdDecompress", [data]),
  };
}
