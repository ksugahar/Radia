function [solution, info] = csrPCG(rowPointers, columns, values, rhs, inverseDiagonal, tolerance, maxIterations, backend)
%CSRPCG Explicit CPU/GPU SPD CSR solve through the same mexcuda module.
% CPU: full host int32 indices and double arrays. GPU: corresponding gpuArray.
% Zero-based CSR numbering; no array transfer or fallback is performed here.
% CPU host implementation is serial; MATLAB computational-thread settings do
% not change that fact. CPU no-CUDA-context acceptance remains unverified.
arguments
    rowPointers
    columns
    values
    rhs
    inverseDiagonal
    tolerance (1,1) double
    maxIterations (1,1) double
    backend (1,1) string {mustBeMember(backend,["cpu","gpu"])}
end
if backend == "gpu"
    radia.gpu.requireRuntime();
elseif exist('radia_cuda_csr_pcg','file') ~= 3
    error('radia:csr:missingNative','The mexcuda CSR module is required.');
end
[solution,info] = radia_cuda_csr_pcg(rowPointers,columns,values,rhs,inverseDiagonal,tolerance,maxIterations,char(backend));
end
