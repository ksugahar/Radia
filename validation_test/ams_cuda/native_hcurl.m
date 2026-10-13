function native_hcurl(volPath,binaryPath,matlabPath,outputPath)
% End-to-end .vol -> native FE -> CUDA assembly -> native AMS -> GPU solve.
addpath(binaryPath);addpath(matlabPath);device=gpuDevice();wait(device);
before=double(device.AvailableMemory);initial=radia_cuda_ams('resident_status');total=tic;
[model,fem,assembly]=radia.gpu.assembleHCurl(volPath, ...
    ["back","front","left","right","top","bottom"]);
[hierarchy,setup]=radia.gpu.prepareAMS(model,SetupThreads=4);
context=radia.gpu.AMSContext(hierarchy);cleanup=onCleanup(@()delete(context));
[x,solve]=context.solve(hierarchy.A,fem.rhs,1e-10,4000,true);
solution=gather(x);rhs=gather(fem.rhs);wait(device);elapsed=toc(total);
mass_values=gather(fem.mass.values);
maps=struct('dofs',gather(fem.dofs),'curls',gather(fem.curls),'volume',gather(fem.volume));
clear cleanup context;
assert(isequal(radia_cuda_ams('resident_status'),initial),'radia:cuda:leak','AMS contexts leaked');
after=double(device.AvailableMemory);
report=struct('total_seconds',elapsed,'assembly',assembly,'setup',setup,'solve',solve, ...
    'gpu_name',device.Name,'available_before',before,'available_after',after, ...
    'retained_device_memory_sample_bytes',max(0,before-after),'matlab_version',version, ...
    'timing_classification','100 build-host diagnostic; not formal performance acceptance');
% Bounded cases only. A larger validator must select checked HDF5 beforehand.
save(outputPath,'solution','rhs','model','mass_values','maps','report');
end
