function native_hcurl(volPath,binaryPath,matlabPath,outputPath,gamma)
% End-to-end .vol -> native FE -> CUDA assembly -> native AMS -> GPU solve.
if nargin<5;gamma=0;end
current=1;if gamma>0;current=20;end
addpath(binaryPath);addpath(matlabPath);device=gpuDevice();wait(device);
before=double(device.AvailableMemory);initial=radia_cuda_ams('resident_status');total=tic;
[model,fem,assembly]=radia.gpu.assembleHCurl(volPath, ...
    ["back","front","left","right","top","bottom"],Gamma=gamma,CurrentDensity=[current 0 0]);
[solution,solverReport]=radia.gpu.solveHCurl(volPath, ...
    ["back","front","left","right","top","bottom"],Gamma=gamma,CurrentDensity=[current 0 0],SetupThreads=4,ReturnOnGPU=false);
rhs=gather(fem.rhs);wait(device);elapsed=toc(total);
mass_values=gather(fem.mass.values);
maps=struct('dofs',gather(fem.dofs),'curls',gather(fem.curls),'volume',gather(fem.volume));
assert(isequal(radia_cuda_ams('resident_status'),initial),'radia:cuda:leak','AMS contexts leaked');
after=double(device.AvailableMemory);
report=struct('total_seconds',elapsed,'assembly',assembly,'solver',solverReport, ...
    'gpu_name',device.Name,'available_before',before,'available_after',after, ...
    'retained_device_memory_sample_bytes',max(0,before-after),'matlab_version',version, ...
    'timing_classification','100 build-host diagnostic; not formal performance acceptance');
% Bounded cases only. A larger validator must select checked HDF5 beforehand.
save(outputPath,'solution','rhs','model','mass_values','maps','report');
end
