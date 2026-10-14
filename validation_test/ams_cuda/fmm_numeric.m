function fmm_numeric(inputPath,binaryPath,matlabPath,outputPath)
addpath(binaryPath);addpath(matlabPath);s=load(inputPath);device=gpuDevice();wait(device);
initial=radia_cuda_fmm('live_count');before=device.AvailableMemory;total=tic;
fmm=radia.gpu.LaplaceFMM(s.sources,s.targets,Order=s.order,LeafCapacity=s.leaf,Separation=s.separation);
cleanup=onCleanup(@()delete(fmm));setupSeconds=toc(total);
charge=gpuArray(s.charge);dipole=gpuArray(s.dipole);
[p,e,info]=fmm.apply(charge,dipole);potential=gather(p);field=gather(e);wait(device);elapsed=toc(total);
[p2,e2,~]=fmm.apply(charge,dipole);
assert(isequal(potential,gather(p2))&&isequal(field,gather(e2)),'radia:fmm:repeat','FMM repeat changed output');
warmed=zeros(3,1);
for k=1:3
    wait(device);timer=tic;
    q=gpuArray(s.charge);d=gpuArray(s.dipole);
    [a,b]=fmm.apply(q,d);a=gather(a);b=gather(b);wait(device);
    warmed(k)=toc(timer);
end
clear cleanup fmm;
assert(radia_cuda_fmm('live_count')==initial,'radia:fmm:leak','FMM handle leaked');
report=struct('setup_seconds',setupSeconds,'total_seconds',elapsed,'info',info, ...
    'device',device.Name,'available_before',before,'available_after',device.AvailableMemory,'matlab',version, ...
    'warmed_transfer_inclusive_seconds',warmed,'warmed_transfer_inclusive_median',median(warmed), ...
    'cold_includes','native context setup, upload, first apply and download; excludes MATLAB startup', ...
    'memory_definition','Native owned allocations plus sampled device free memory; not a process peak');
save(outputPath,'potential','field','report');
end
