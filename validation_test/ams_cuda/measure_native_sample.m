function measure_native_sample(job,name,output)
% Timed public API construction, all transfers, execution, recovery and release.
switch name
    case 'fmm-small';s=load(fullfile(job,'numeric-input.mat'));
    case 'fmm-large';s=load(fullfile(job,'numeric-large.mat'));
    case 'physical';s=load(fullfile(job,'physics-input.mat'));
end
d=gpuDevice();wait(d);before=d.AvailableMemory;timer=tic;info=struct;
switch name
    case 'ams'
        [solution,info]=radia.gpu.solveHCurl(fullfile(job,'model.vol'), ...
            ["back","front","left","right","top","bottom"],Gamma=4,CurrentDensity=[20 0 0], ...
            SetupThreads=1,LinearTolerance=1e-12,ReturnOnGPU=false);
        values={solution};
    case {'fmm-small','fmm-large'}
        f=radia.gpu.LaplaceFMM(s.sources,s.targets,Order=s.order,LeafCapacity=s.leaf,Separation=s.separation);
        clean=onCleanup(@()delete(f));
        [p,e,info]=f.apply(gpuArray(s.charge),gpuArray(s.dipole));values={gather(p),gather(e)};
        clear clean f;
    case 'physical'
        map=gpuArray(spdiags(s.weights,0,numel(s.weights),numel(s.weights)));
        test=gpuArray(speye(size(s.targets,1)));correction=gpuArray(sparse(size(s.targets,1),numel(s.weights)));
        f=radia.gpu.LaplaceBoundaryFMM(s.sphere,s.targets,map,test,correction,Order=8,LeafCapacity=32,Separation=2);
        clean=onCleanup(@()delete(f));[v,info.sl]=f.apply(gpuArray(ones(numel(s.weights),1)));sl=gather(v);clear clean f;
        f=radia.gpu.LaplaceBoundaryFMM(s.sphere,s.targets,map,test,correction,Kind="double-layer",Normals=s.sphere,Order=8,LeafCapacity=32,Separation=2);
        clean=onCleanup(@()delete(f));[v,info.dl]=f.apply(gpuArray(ones(numel(s.weights),1)));dl=gather(v);clear clean f;
        f=radia.gpu.LaplaceFMM(s.coil,s.coilTargets,Order=8,LeafCapacity=32,Separation=2);
        clean=onCleanup(@()delete(f));[b,a,info.bs]=radia.gpu.biotSavartFMM(f,gpuArray(s.current));b=gather(b);a=gather(a);clear clean f;
        f=radia.gpu.LaplaceFMM(s.sphere,s.targets,Order=8,LeafCapacity=32,Separation=2);
        clean=onCleanup(@()delete(f));[v,info.cloud]=radia.gpu.fieldCloudFMM(f,gpuArray(s.weights));cloud=gather(v);clear clean f;
        values={sl,dl,b,a,cloud};
    case 'boundary'
        labels=["back","front","left","right","top","bottom"];
        [f,info.sl_setup]=radia.gpu.prepareBoundaryFMM(fullfile(job,'boundary.vol'),labels);
        clean=onCleanup(@()delete(f));n=info.sl_setup.boundary_dofs;x=sin((1:n)'*.37)+.1*cos((1:n)'*.13);
        [v,info.sl]=f.apply(gpuArray(x));sl=gather(v);clear clean f;
        [f,info.dl_setup]=radia.gpu.prepareBoundaryFMM(fullfile(job,'boundary.vol'),labels,Kind="double-layer");
        clean=onCleanup(@()delete(f));[v,info.dl]=f.apply(gpuArray(x));dl=gather(v);clear clean f;
        values={sl,dl};
    otherwise;error('Unknown benchmark case');
end
wait(d);whole_seconds=toc(timer);after=d.AvailableMemory;
assert(radia_cuda_fmm('live_count')==0);
if strcmp(name,'ams')
    assert(all(structfun(@(x)x==0,radia_cuda_ams('resident_status'))));
end
report=struct('whole_seconds',whole_seconds,'native',info,'device_free_before',before, ...
    'device_free_after',after,'available_memory_delta',before-after);
% Output serialization is outside the numerical API timing, equally for CPU.
save(output,'values','report');
end
