//! CUDA GPU alignment backend.

use std::sync::Arc;

use cudarc::driver::{
    CudaContext, CudaFunction, CudaSlice, CudaStream, LaunchConfig, PushKernelArg,
};
use cudarc::nvrtc::compile_ptx;

use cyanea_core::{CyaneaError, Result};

use crate::scoring::ScoringScheme;
use crate::types::{AlignmentMode, AlignmentResult};

use super::common::{
    encode_pairs, extract_scoring_params, partition_pairs, reconstruct_alignment, GpuAlignConfig,
};

const ALIGN_CUDA: &str = include_str!("kernels/align.cu");

/// CUDA-based batch aligner.
pub struct CudaAligner {
    _ctx: Arc<CudaContext>,
    stream: Arc<CudaStream>,
    align_fn: CudaFunction,
}

impl CudaAligner {
    pub fn new() -> Result<Self> {
        let ctx = CudaContext::new(0)
            .map_err(|e| CyaneaError::Other(format!("CUDA context init: {e}")))?;
        let stream = ctx.default_stream();

        let ptx = compile_ptx(ALIGN_CUDA)
            .map_err(|e| CyaneaError::Other(format!("CUDA align kernel compile: {e}")))?;
        let module = ctx
            .load_module(ptx)
            .map_err(|e| CyaneaError::Other(format!("CUDA module load: {e}")))?;
        let align_fn = module
            .load_function("align_pairs")
            .map_err(|e| CyaneaError::Other(format!("CUDA function load: {e}")))?;

        Ok(Self {
            _ctx: ctx,
            stream,
            align_fn,
        })
    }

    pub fn align_batch(
        &self,
        pairs: &[(&[u8], &[u8])],
        mode: AlignmentMode,
        scoring: &ScoringScheme,
        config: &GpuAlignConfig,
    ) -> Result<Vec<AlignmentResult>> {
        let n = pairs.len();
        if n == 0 {
            return Ok(Vec::new());
        }

        let (gpu_pairs, cpu_pairs) = partition_pairs(pairs, config.max_bandwidth);

        let mut results = vec![None; n];
        for (orig_idx, q, t) in &cpu_pairs {
            let r = crate::align(q, t, mode, scoring)?;
            results[*orig_idx] = Some(r);
        }

        if gpu_pairs.is_empty() {
            return Ok(results.into_iter().map(|r| r.unwrap()).collect());
        }

        let gpu_pair_refs: Vec<(&[u8], &[u8])> =
            gpu_pairs.iter().map(|(_, q, t)| (*q, *t)).collect();
        let (packed, index) = encode_pairs(&gpu_pair_refs);
        let num_gpu_pairs = gpu_pairs.len();
        let bandwidth = config.max_bandwidth;

        let max_traceback: usize = index
            .iter()
            .map(|idx| (idx.query_len as usize + 1) * bandwidth)
            .sum();

        let (match_s, mismatch_s, gap_o, gap_e) = extract_scoring_params(scoring);
        let mode_val: u32 = match mode {
            AlignmentMode::Local => 0,
            AlignmentMode::Global => 1,
            AlignmentMode::SemiGlobal => 2,
        };

        // Upload to GPU
        let seqs_dev: CudaSlice<u8> = self
            .stream
            .clone_htod(&packed)
            .map_err(|e| CyaneaError::Other(format!("CUDA htod seqs: {e}")))?;

        // Pack index as flat u32 array
        let mut index_flat: Vec<u32> = Vec::with_capacity(index.len() * 4);
        for idx in &index {
            index_flat.push(idx.query_offset);
            index_flat.push(idx.query_len);
            index_flat.push(idx.target_offset);
            index_flat.push(idx.target_len);
        }
        let index_dev: CudaSlice<u32> = self
            .stream
            .clone_htod(&index_flat)
            .map_err(|e| CyaneaError::Other(format!("CUDA htod index: {e}")))?;

        let scores_dev: CudaSlice<i32> = self
            .stream
            .alloc_zeros(num_gpu_pairs)
            .map_err(|e| CyaneaError::Other(format!("CUDA alloc scores: {e}")))?;
        let endpos_dev: CudaSlice<u32> = self
            .stream
            .alloc_zeros(num_gpu_pairs * 2)
            .map_err(|e| CyaneaError::Other(format!("CUDA alloc endpos: {e}")))?;
        let tb_dev: CudaSlice<u8> = self
            .stream
            .alloc_zeros(max_traceback.max(1))
            .map_err(|e| CyaneaError::Other(format!("CUDA alloc traceback: {e}")))?;

        let cfg = LaunchConfig {
            grid_dim: (((num_gpu_pairs as u32) + 255) / 256, 1, 1),
            block_dim: (256.min(num_gpu_pairs as u32).max(1), 1, 1),
            shared_mem_bytes: 0,
        };

        let bw_u32 = bandwidth as u32;
        let npairs_u32 = num_gpu_pairs as u32;

        unsafe {
            self.stream
                .launch_builder(&self.align_fn)
                .arg(&seqs_dev)
                .arg(&index_dev)
                .arg(&scores_dev)
                .arg(&endpos_dev)
                .arg(&tb_dev)
                .arg(&match_s)
                .arg(&mismatch_s)
                .arg(&gap_o)
                .arg(&gap_e)
                .arg(&bw_u32)
                .arg(&mode_val)
                .arg(&npairs_u32)
                .launch(cfg)
        }
        .map_err(|e| CyaneaError::Other(format!("CUDA launch align: {e}")))?;

        // Download results
        let scores = self
            .stream
            .clone_dtoh(&scores_dev)
            .map_err(|e| CyaneaError::Other(format!("CUDA dtoh scores: {e}")))?;
        let end_positions = self
            .stream
            .clone_dtoh(&endpos_dev)
            .map_err(|e| CyaneaError::Other(format!("CUDA dtoh endpos: {e}")))?;
        let traceback_data = self
            .stream
            .clone_dtoh(&tb_dev)
            .map_err(|e| CyaneaError::Other(format!("CUDA dtoh traceback: {e}")))?;

        // Reconstruct
        let mut tb_offset = 0;
        for (pair_idx, (orig_idx, q, t)) in gpu_pairs.iter().enumerate() {
            let score = scores[pair_idx];
            let end_i = end_positions[pair_idx * 2] as usize;
            let end_j = end_positions[pair_idx * 2 + 1] as usize;
            let tb_size = (index[pair_idx].query_len as usize + 1) * bandwidth;
            let tb_slice = &traceback_data
                [tb_offset..tb_offset + tb_size.min(traceback_data.len() - tb_offset)];

            let result = reconstruct_alignment(
                q, t, score, end_i, end_j, tb_slice, bandwidth, mode, scoring,
            );
            results[*orig_idx] = Some(result);
            tb_offset += tb_size;
        }

        Ok(results.into_iter().map(|r| r.unwrap()).collect())
    }
}
