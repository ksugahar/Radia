function ams_tuning_sample(job,binaryPath,output)
% Numerical timing excludes fixture loading and output serialization.
addpath(binaryPath,'-begin');addpath(fullfile(job,'matlab'));g=gpuDevice;wait(g);
initial=radia_cuda_ams('resident_status');labels=["back","front","left","right","top","bottom"];
before=g.AvailableMemory;timer=tic;
[nonlinear,nonlinearInfo]=radia.gpu.solveHCurl(fullfile(job,'model.vol'),labels, ...
    Gamma=4,CurrentDensity=[20 0 0],SetupThreads=1,LinearTolerance=1e-12,ReturnOnGPU=false);
wait(g);nonlinearSeconds=toc(timer);timer=tic;
[linear,linearInfo]=radia.gpu.solveHCurl(fullfile(job,'model.vol'),labels, ...
    Gamma=0,CurrentDensity=[1 0 0],SetupThreads=1,LinearTolerance=1e-12,ReturnOnGPU=false);
wait(g);linearSeconds=toc(timer);timer=tic;
system=radia.gpu.HCurlSystem(fullfile(job,'model.vol'),labels,SetupThreads=1);
clean=onCleanup(@()delete(system));wait(g);preparedSeconds=toc(timer);
retainedBefore=g.AvailableMemory;scales=[1 2 .5 1 0 1];reuse=zeros(numel(linear),numel(scales));
reuseSeconds=zeros(numel(scales),1);reuseInfo=cell(numel(scales),1);
for k=1:numel(scales)
    wait(g);timer=tic;
    [reuse(:,k),reuseInfo{k}]=system.solve(CurrentScale=scales(k),Tolerance=1e-12,ReturnOnGPU=false);
    wait(g);reuseSeconds(k)=toc(timer);
    assert(norm(reuse(:,k)-scales(k)*linear)<=1e-9*max(1,norm(scales(k)*linear)));
end
retainedAfter=g.AvailableMemory;setupInfo=system.SetupInfo;clear clean;
try;system.solve();error('radia:test:missingError','Released system accepted solve');
catch ex;assert(strcmp(ex.identifier,'radia:cuda:released'));end
assert(isequal(radia_cuda_ams('resident_status'),initial),'radia:test:leak','AMS state leaked');
report=struct('nonlinear_seconds',nonlinearSeconds,'linear_fresh_seconds',linearSeconds, ...
    'prepared_setup_seconds',preparedSeconds,'reused_seconds',reuseSeconds,'nonlinear',nonlinearInfo, ...
    'linear',linearInfo,'reuse', {reuseInfo},'setup',setupInfo,'scales',scales, ...
    'device_free_before',before,'device_free_retained_before',retainedBefore, ...
    'device_free_retained_after',retainedAfter,'device_free_after',g.AvailableMemory, ...
    'loaded_binary',which('radia_cuda_ams'),'matlab',version,'device',g.Name,'registry_restored',true);
save(output,'nonlinear','linear','reuse','report');
clear system;clear radia_cuda_ams;rmpath(binaryPath);
end
