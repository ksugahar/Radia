function [solution, info] = csrPCG(rowPointers, columns, values, rhs, inverseDiagonal, tolerance, maxIterations)
%CSRPCG Solve real SPD CSR data on the NVIDIA GPU. No CPU solver fallback.
% CSR indices are zero-based int32 gpuArray. All numeric data are FP64 gpuArray.
% The caller owns assembly, uploads and the supplied arrays.
persistent configured
if isempty(configured)
    radia.gpu.requireRuntime();
    configured = true;
end
[solution, info] = radia_cuda_csr_pcg(rowPointers, columns, values, rhs, inverseDiagonal, tolerance, maxIterations);
end
