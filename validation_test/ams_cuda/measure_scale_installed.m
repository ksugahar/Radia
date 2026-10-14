function measure_scale_installed(volPath,payload,output)
% Fresh linear setup plus zero-initial-guess scaled-RHS solves; FP64.
addpath(payload);g=gpuDevice;wait(g);initial=radia_cuda_ams('resident_status');
freeBefore=g.AvailableMemory;timer=tic;
system=radia.gpu.HCurlSystem(volPath,["back","front","left","right","top","bottom"],SetupThreads=1);
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
report=struct('setup_seconds',setupSeconds,'fresh_total_seconds',setupSeconds+seconds(1), ...
    'solve_seconds',seconds,'scales',scales,'info',{info},'setup',setupInfo, ...
    'device_free_before',freeBefore,'device_free_retained',freeRetained, ...
    'device_free_after_solves',freeAfterSolves,'device_free_after_release',g.AvailableMemory, ...
    'registry_restored',true,'matlab',version,'device',g.Name,'loaded_binary',which('radia_cuda_ams'));
save(output,'solution','report','-v7');
end
