function requireRuntime()
%REQUIRERUNTIME Fail explicitly when the supported CUDA MEX runtime is absent.
if ~strcmp(version('-release'), '2026a')
    error('radia:gpu:MatlabRelease', 'This alpha was validated with MATLAB R2026a.');
end
if ~license('test', 'Distrib_Computing_Toolbox')
    error('radia:gpu:Toolbox', 'Parallel Computing Toolbox is required.');
end
if exist('radia_cuda_csr_pcg', 'file') ~= 3
    error('radia:gpu:MissingMex', 'The installed CUDA MEX binary is missing.');
end
gpuDevice(); % Validate availability without resetting caller-owned GPU data.
end
