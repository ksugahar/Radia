function measure_installed_cpu_gpu(job,payload,output)
% Whole nonlinear API, including native setup, upload and solution download.
addpath(payload,'-begin');g=gpuDevice;wait(g);initial=radia_cuda_ams('resident_status');
labels=["back","front","left","right","top","bottom"];
before=g.AvailableMemory;timer=tic;
[x,info]=radia.gpu.solveHCurl(fullfile(job,'model.vol'),labels,Gamma=4, ...
    CurrentDensity=[20 0 0],SetupThreads=1,LinearTolerance=1e-12,ReturnOnGPU=false);
wait(g);seconds=toc(timer);
assert(all(isfinite(x)) && isequal(radia_cuda_ams('resident_status'),initial));
report=struct('whole_api_seconds',seconds,'info',info,'device_free_before',before, ...
    'device_free_after',g.AvailableMemory,'device',g.Name,'matlab',version, ...
    'host_setup_threads',1,'registry_restored',true);
save(output,'x','report','-v7');
end
