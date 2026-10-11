function result = benchmark_gpu(inputPath, outputPath, binaryDir)
% NGSolve owns assembly; only the reduced SPD solve runs on the GPU.
addpath(binaryDir);
s = load(inputPath);
d = gpuDevice();
freeBefore = d.AvailableMemory;
t = tic;
rp = gpuArray(s.rowptr); ci = gpuArray(s.columns); av = gpuArray(s.values);
b = gpuArray(s.rhs); di = gpuArray(s.diagonal_inverse);
wait(d); transfer = toc(t);
tol = 1e-10; maxiter = 10000;
t = tic;
[x, info] = radia_cuda_csr_pcg(rp,ci,av,b,di,tol,maxiter);
wait(d); cold = toc(t);
times = zeros(1,5);
for k = 1:5
    t = tic;
    [x, info] = radia_cuda_csr_pcg(rp,ci,av,b,di,tol,maxiter);
    wait(d); times(k) = toc(t);
end
t = tic; solution = gather(x); download = toc(t);
result = struct('passed',true,'schema','radia.cuda-poisson-gpu.v1', ...
    'backend','custom CUDA CSR PCG; NGSolve CPU assembly; not upstream ngscuda', ...
    'device',d.Name,'compute_capability',d.ComputeCapability, ...
    'total_memory_bytes',d.TotalMemory,'available_memory_before',freeBefore, ...
    'available_memory_after',d.AvailableMemory,'precision','float64', ...
    'upload_seconds',transfer,'cold_solve_seconds',cold,'warm_seconds',times, ...
    'warm_median_seconds',median(times),'download_seconds',download, ...
    'warm_with_transfers_seconds',transfer+median(times)+download, ...
    'info',info);
save(outputPath,'solution','result','-v7');
end
