// ---------------------------------------------------------------------------
// @cyanea/bio -- Web Worker script
// ---------------------------------------------------------------------------
// Runs in a Worker context (browser or Node.js worker_threads).
// Loads WASM on startup, then dispatches calls via the registry.
//
// Usage (browser):
//   const worker = new Worker(new URL("@cyanea/bio/worker", import.meta.url),
//                             { type: "module" });
//
// Usage (Node.js):
//   const { Worker } = require("worker_threads");
//   const worker = new Worker(require.resolve("@cyanea/bio/worker"));
// ---------------------------------------------------------------------------

import type {
  MainToWorkerMessage,
  WorkerToMainMessage,
  WorkerProgress,
  SeqPair,
  AlignmentMode,
  AlignmentResult,
} from "./types.js";
import { registry } from "./registry.js";
import { Align } from "./index.js";

// ── Cross-environment messaging ────────────────────────────────────────────

/** Whether we're running in a Node.js worker_threads context. */
const isNode =
  typeof (globalThis as Record<string, unknown>).process !== "undefined" &&
  typeof ((globalThis as Record<string, unknown>).process as Record<string, unknown>)?.versions === "object";

/**
 * Post a message to the main thread, abstracting over browser Worker
 * (self.postMessage) and Node.js worker_threads (parentPort.postMessage).
 */
/** Cached parentPort for Node.js worker_threads. Set during initialize(). */
let nodeParentPort: import("worker_threads").MessagePort | null | undefined;

function post(msg: WorkerToMainMessage): void {
  if (isNode && nodeParentPort) {
    nodeParentPort.postMessage(msg);
  } else {
    (self as unknown as { postMessage(msg: unknown): void }).postMessage(msg);
  }
}

// ── Cancellation tracking ──────────────────────────────────────────────────

/** Set of request IDs that have been cancelled. */
const cancelled = new Set<string>();

function isCancelled(id: string): boolean {
  return cancelled.has(id);
}

// ── Progressive dispatch helpers ───────────────────────────────────────────

const PROGRESS_INTERVAL = 50;

/**
 * Decompose Align.alignBatch into individual calls with progress + cancellation.
 */
function alignBatchProgressive(
  id: string,
  pairs: SeqPair[],
  mode: AlignmentMode,
  matchScore: number,
  mismatchScore: number,
  gapOpen: number,
  gapExtend: number,
): AlignmentResult[] {
  const total = pairs.length;
  const results: AlignmentResult[] = [];

  for (let i = 0; i < total; i++) {
    if (isCancelled(id)) {
      throw new Error("Cancelled");
    }

    const pair = pairs[i];
    const result = Align.alignDnaCustom(
      pair.query,
      pair.target,
      mode,
      matchScore,
      mismatchScore,
      gapOpen,
      gapExtend,
    );
    results.push(result);

    if ((i + 1) % PROGRESS_INTERVAL === 0 || i === total - 1) {
      const progress: WorkerProgress = {
        type: "progress",
        id,
        current: i + 1,
        total,
      };
      post(progress);
    }
  }

  return results;
}

// ── Message handler ────────────────────────────────────────────────────────

async function handleMessage(msg: MainToWorkerMessage): Promise<void> {
  switch (msg.type) {
    case "init":
      // WASM is already initialized at module load — just confirm ready.
      post({ type: "ready" });
      break;

    case "cancel":
      cancelled.add(msg.id);
      break;

    case "call": {
      const { id, method, args } = msg;

      try {
        let result: unknown;

        // Handle progressive methods specially
        if (method === "Align.alignBatch") {
          result = alignBatchProgressive(
            id,
            args[0] as SeqPair[],
            args[1] as AlignmentMode,
            args[2] as number,
            args[3] as number,
            args[4] as number,
            args[5] as number,
          );
        } else {
          const entry = registry.get(method);
          if (!entry) {
            post({ type: "result", id, error: `Unknown method: ${method}` });
            return;
          }
          result = entry.fn(...args);
        }

        // Don't send result if cancelled
        if (!isCancelled(id)) {
          post({ type: "result", id, ok: result });
        }
      } catch (e) {
        if (!isCancelled(id)) {
          const error = e instanceof Error ? e.message : String(e);
          post({ type: "result", id, error });
        }
      } finally {
        cancelled.delete(id);
      }
      break;
    }
  }
}

// ── WASM initialization + listener setup ───────────────────────────────────

async function initialize(): Promise<void> {
  try {
    // Dynamic import of the WASM init function
    const { default: init } = await import("../pkg/cyanea_wasm.js");
    await init();

    // Register message handler
    if (isNode) {
      const wt = await import("worker_threads");
      nodeParentPort = wt.parentPort;
      nodeParentPort?.on("message", (msg: MainToWorkerMessage) => {
        void handleMessage(msg);
      });
    } else {
      self.addEventListener("message", (event: MessageEvent<MainToWorkerMessage>) => {
        void handleMessage(event.data);
      });
    }

    post({ type: "ready" });
  } catch (e) {
    const error = e instanceof Error ? e.message : String(e);
    post({ type: "init_error", error });
  }
}

void initialize();
