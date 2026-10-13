function fmm_numeric(inputPath,binaryPath,matlabPath,outputPath)
addpath(binaryPath);addpath(matlabPath);s=load(inputPath);device=gpuDevice();wait(device);
initial=radia_cuda_fmm('live_count');before=device.AvailableMemory;total=tic;
fmm=radia.gpu.LaplaceFMM(s.sources,s.targets,Order=s.order,LeafCapacity=16,Separation=1.5);
cleanup=onCleanup(@()delete(fmm));setupSeconds=toc(total);
charge=gpuArray(s.charge);dipole=gpuArray(s.dipole);
[p,e,info]=fmm.apply(charge,dipole);potential=gather(p);field=gather(e);wait(device);elapsed=toc(total);
[p2,e2,~]=fmm.apply(charge,dipole);
assert(isequal(potential,gather(p2))&&isequal(field,gather(e2)),'radia:fmm:repeat','FMM repeat changed output');
clear cleanup fmm;
assert(radia_cuda_fmm('live_count')==initial,'radia:fmm:leak','FMM handle leaked');
report=struct('setup_seconds',setupSeconds,'total_seconds',elapsed,'info',info, ...
    'device',device.Name,'available_before',before,'available_after',device.AvailableMemory,'matlab',version);
save(outputPath,'potential','field','report');
end
