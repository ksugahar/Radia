function measure_scale_installed(volPath,payload,output,backend,icPath)
% Fresh linear setup plus zero-initial-guess scaled-RHS solves; FP64.
if nargin<4;backend='ams';end
if nargin<5;icPath='';end
addpath(payload);g=gpuDevice;wait(g);initial=radia_cuda_ams('resident_status');
if strcmp(backend,'ic');addpath(icPath,'-begin');initialIC=radia_cuda_ic('live_count');end
freeBefore=g.AvailableMemory;timer=tic;
labels=["back","front","left","right","top","bottom"];
if strcmp(backend,'ic')
    system=radia.gpu.HCurlICSystem(volPath,labels,Shift=1.05);
else
    assert(strcmp(backend,'ams'),'radia:test:backend','Unknown backend');
    system=radia.gpu.HCurlSystem(volPath,labels,SetupThreads=1);
end
cleanup=onCleanup(@()delete(system));wait(g);setupSeconds=toc(timer);
freeRetained=g.AvailableMemory;setupInfo=system.SetupInfo;
scales=[1 2 .5 1];seconds=zeros(4,1);info=cell(4,1);solution=[];
for k=1:4
    wait(g);timer=tic;
    [x,info{k}]=system.solve(CurrentScale=scales(k),Tolerance=1e-12,ReturnOnGPU=false);
    wait(g);seconds(k)=toc(timer);solution(:,k)=x; %#ok<AGROW>
end
freeAfterSolves=g.AvailableMemory;clear cleanup system;wait(g);
assert(isequal(radia_cuda_ams('resident_status'),initial),'radia:test:leak','AMS state leaked');
if strcmp(backend,'ic');assert(radia_cuda_ic('live_count')==initialIC,'radia:test:leak','IC state leaked');end
report=struct('setup_seconds',setupSeconds,'fresh_total_seconds',setupSeconds+seconds(1), ...
    'solve_seconds',seconds,'scales',scales,'info',{info},'setup',setupInfo, ...
    'device_free_before',freeBefore,'device_free_retained',freeRetained, ...
    'device_free_after_solves',freeAfterSolves,'device_free_after_release',g.AvailableMemory, ...
    'registry_restored',true,'matlab',version,'device',g.Name,'loaded_binary',which('radia_cuda_ams'));
save(output,'solution','report','-v7');
end
