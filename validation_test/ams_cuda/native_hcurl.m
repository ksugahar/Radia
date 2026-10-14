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
warmed=zeros(3,1);
for k=1:3
    wait(device);timer=tic;
    [repeat,~]=radia.gpu.solveHCurl(volPath, ...
        ["back","front","left","right","top","bottom"],Gamma=gamma,CurrentDensity=[current 0 0],SetupThreads=4,ReturnOnGPU=false);
    wait(device);warmed(k)=toc(timer);
    assert(norm(repeat-solution)<=1e-10*max(norm(solution),realmin),'radia:cuda:repeat','Repeated solve changed solution');
end
assert(isequal(radia_cuda_ams('resident_status'),initial),'radia:cuda:leak','Repeated AMS contexts leaked');
report=struct('total_seconds',elapsed,'assembly',assembly,'solver',solverReport, ...
    'gpu_name',device.Name,'available_before',before,'available_after',after, ...
    'retained_device_memory_sample_bytes',max(0,before-after),'matlab_version',version, ...
    'timing_classification','100 build-host diagnostic; not formal performance acceptance', ...
    'warmed_transfer_inclusive_seconds',warmed,'warmed_transfer_inclusive_median',median(warmed), ...
    'cold_production_seconds',solverReport.total_seconds, ...
    'total_seconds_definition','Validation includes a separate assembly export in addition to the production solve');
% Bounded cases only. A larger validator must select checked HDF5 beforehand.
save(outputPath,'solution','rhs','model','mass_values','maps','report');
end
