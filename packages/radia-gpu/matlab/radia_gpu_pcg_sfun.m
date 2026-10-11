function radia_gpu_pcg_sfun(block)
%RADIA_GPU_PCG_SFUN Readable entry point for the CUDA-only CSR PCG block.
radia.simulink.cudaPCGSFunction(block);
end
