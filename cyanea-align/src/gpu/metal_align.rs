//! Metal GPU alignment backend.

use std::ffi::c_void;

use metal_rs::{
    Buffer as MtlBuffer, CommandQueue, CompileOptions, ComputePipelineState, Device, MTLSize,
    MTLResourceOptions,
};

use cyanea_core::{CyaneaError, Result};

use crate::scoring::ScoringScheme;
use crate::types::{AlignmentMode, AlignmentResult};

use super::common::{
    encode_pairs, extract_scoring_params, partition_pairs, reconstruct_alignment, GpuAlignConfig,
    SeqIndex,
};

const ALIGN_MSL: &str = include_str!("kernels/align.metal");

/// Metal-based batch aligner.
pub struct MetalAligner {
    device: Device,
    queue: CommandQueue,
    pipeline: ComputePipelineState,
}

impl MetalAligner {
    pub fn new() -> Result<Self> {
        let device = Device::system_default()
            .ok_or_else(|| CyaneaError::Other("no Metal device".into()))?;
        let queue = device.new_command_queue();

        let lib = device
            .new_library_with_source(ALIGN_MSL, &CompileOptions::new())
            .map_err(|e| CyaneaError::Other(format!("Metal align shader compile: {e}")))?;
        let func = lib
            .get_function("align_pairs", None)
            .map_err(|e| CyaneaError::Other(format!("Metal align function: {e}")))?;
        let pipeline = device
            .new_compute_pipeline_state_with_function(&func)
            .map_err(|e| CyaneaError::Other(format!("Metal align pipeline: {e}")))?;

        Ok(Self {
            device,
            queue,
            pipeline,
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

        // Partition pairs
        let (gpu_pairs, cpu_pairs) = partition_pairs(pairs, config.max_bandwidth);

        // CPU fallback results
        let mut results = vec![None; n];
        for (orig_idx, q, t) in &cpu_pairs {
            let r = crate::align(q, t, mode, scoring)?;
            results[*orig_idx] = Some(r);
        }

        if gpu_pairs.is_empty() {
            return Ok(results.into_iter().map(|r| r.unwrap()).collect());
        }

        // Encode GPU pairs
        let gpu_pair_refs: Vec<(&[u8], &[u8])> =
            gpu_pairs.iter().map(|(_, q, t)| (*q, *t)).collect();
        let (packed, index) = encode_pairs(&gpu_pair_refs);
        let num_gpu_pairs = gpu_pairs.len();
        let bandwidth = config.max_bandwidth;

        // Calculate traceback size
        let max_traceback: usize = index
            .iter()
            .map(|idx| (idx.query_len as usize + 1) * bandwidth)
            .sum();

        // Scoring params
        let (match_s, mismatch_s, gap_o, gap_e) = extract_scoring_params(scoring);

        // Create Metal buffers
        let seqs_buf = self.make_buffer(&packed);
        let index_buf = self.make_index_buffer(&index);
        let scores_buf = self.alloc_buffer::<i32>(num_gpu_pairs);
        let endpos_buf = self.alloc_buffer::<u32>(num_gpu_pairs * 2);
        let tb_buf = self.alloc_buffer::<u8>(max_traceback.max(1));

        let match_buf = self.scalar_i32(match_s);
        let mismatch_buf = self.scalar_i32(mismatch_s);
        let gap_open_buf = self.scalar_i32(gap_o);
        let gap_extend_buf = self.scalar_i32(gap_e);
        let bw_buf = self.scalar_u32(bandwidth as u32);
        let mode_val: u32 = match mode {
            AlignmentMode::Local => 0,
            AlignmentMode::Global => 1,
            AlignmentMode::SemiGlobal => 2,
        };
        let mode_buf = self.scalar_u32(mode_val);
        let npairs_buf = self.scalar_u32(num_gpu_pairs as u32);

        // Dispatch
        let cmd = self.queue.new_command_buffer();
        let enc = cmd.new_compute_command_encoder();
        enc.set_compute_pipeline_state(&self.pipeline);
        enc.set_buffer(0, Some(&seqs_buf), 0);
        enc.set_buffer(1, Some(&index_buf), 0);
        enc.set_buffer(2, Some(&scores_buf), 0);
        enc.set_buffer(3, Some(&endpos_buf), 0);
        enc.set_buffer(4, Some(&tb_buf), 0);
        enc.set_buffer(5, Some(&match_buf), 0);
        enc.set_buffer(6, Some(&mismatch_buf), 0);
        enc.set_buffer(7, Some(&gap_open_buf), 0);
        enc.set_buffer(8, Some(&gap_extend_buf), 0);
        enc.set_buffer(9, Some(&bw_buf), 0);
        enc.set_buffer(10, Some(&mode_buf), 0);
        enc.set_buffer(11, Some(&npairs_buf), 0);

        let grid = MTLSize::new(num_gpu_pairs as u64, 1, 1);
        let tg = MTLSize::new(1.min(num_gpu_pairs as u64).max(1), 1, 1);
        enc.dispatch_threads(grid, tg);
        enc.end_encoding();
        cmd.commit();
        cmd.wait_until_completed();

        // Read results
        let scores = self.read_i32(&scores_buf, num_gpu_pairs);
        let end_positions = self.read_u32(&endpos_buf, num_gpu_pairs * 2);
        let traceback_data = self.read_u8(&tb_buf, max_traceback.max(1));

        // Reconstruct alignments
        let mut tb_offset = 0;
        for (pair_idx, (orig_idx, q, t)) in gpu_pairs.iter().enumerate() {
            let score = scores[pair_idx];
            let end_i = end_positions[pair_idx * 2] as usize;
            let end_j = end_positions[pair_idx * 2 + 1] as usize;
            let tb_size = (index[pair_idx].query_len as usize + 1) * bandwidth;
            let tb_slice = &traceback_data[tb_offset..tb_offset + tb_size.min(traceback_data.len() - tb_offset)];

            let result =
                reconstruct_alignment(q, t, score, end_i, end_j, tb_slice, bandwidth, mode, scoring);
            results[*orig_idx] = Some(result);
            tb_offset += tb_size;
        }

        Ok(results.into_iter().map(|r| r.unwrap()).collect())
    }

    fn make_buffer(&self, data: &[u8]) -> MtlBuffer {
        if data.is_empty() {
            return self
                .device
                .new_buffer(4, MTLResourceOptions::StorageModeShared);
        }
        self.device.new_buffer_with_data(
            data.as_ptr() as *const c_void,
            data.len() as u64,
            MTLResourceOptions::StorageModeShared,
        )
    }

    fn make_index_buffer(&self, index: &[SeqIndex]) -> MtlBuffer {
        // Pack as uint4 (4 x u32 per entry)
        let mut data = Vec::with_capacity(index.len() * 4);
        for idx in index {
            data.push(idx.query_offset);
            data.push(idx.query_len);
            data.push(idx.target_offset);
            data.push(idx.target_len);
        }
        self.device.new_buffer_with_data(
            data.as_ptr() as *const c_void,
            (data.len() * std::mem::size_of::<u32>()) as u64,
            MTLResourceOptions::StorageModeShared,
        )
    }

    fn alloc_buffer<T>(&self, count: usize) -> MtlBuffer {
        let byte_len = ((count * std::mem::size_of::<T>()) as u64).max(4);
        self.device
            .new_buffer(byte_len, MTLResourceOptions::StorageModeShared)
    }

    fn scalar_i32(&self, val: i32) -> MtlBuffer {
        self.device.new_buffer_with_data(
            &val as *const i32 as *const c_void,
            std::mem::size_of::<i32>() as u64,
            MTLResourceOptions::StorageModeShared,
        )
    }

    fn scalar_u32(&self, val: u32) -> MtlBuffer {
        self.device.new_buffer_with_data(
            &val as *const u32 as *const c_void,
            std::mem::size_of::<u32>() as u64,
            MTLResourceOptions::StorageModeShared,
        )
    }

    fn read_i32(&self, buf: &MtlBuffer, count: usize) -> Vec<i32> {
        let ptr = buf.contents() as *const i32;
        unsafe { std::slice::from_raw_parts(ptr, count) }.to_vec()
    }

    fn read_u32(&self, buf: &MtlBuffer, count: usize) -> Vec<u32> {
        let ptr = buf.contents() as *const u32;
        unsafe { std::slice::from_raw_parts(ptr, count) }.to_vec()
    }

    fn read_u8(&self, buf: &MtlBuffer, count: usize) -> Vec<u8> {
        let ptr = buf.contents() as *const u8;
        unsafe { std::slice::from_raw_parts(ptr, count) }.to_vec()
    }
}

// Safety: Metal objects are thread-safe Obj-C wrappers.
unsafe impl Send for MetalAligner {}
unsafe impl Sync for MetalAligner {}
