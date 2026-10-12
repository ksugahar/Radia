function nonlinear_resident(inputPath,outputPath,binaryPath,matlabPath,setupThreads)
% Sparse nonlinear assembly and persistent PCG; no dense global operators.
if nargin<5;setupThreads=1;end
addpath(binaryPath);addpath(matlabPath);s=load(inputPath);device=gpuDevice();
initialStatus=radia_cuda_ams('resident_status');
inputBytes=8*(numel(s.fem.curls)+numel(s.fem.volume)+numel(s.fem.mass.values)+numel(s.rhs))+ ...
    4*(numel(s.fem.dofs)+numel(s.fem.slots)+numel(s.fem.mass.rowptr)+numel(s.fem.mass.columns)+numel(s.fem.free));
assert(12*inputBytes<device.AvailableMemory-max(2*1024^3,.2*device.TotalMemory), ...
    'radia:nonlinear:memory','Sparse problem exceeds conservative memory budget');
wait(device);timer=tic;fem=upload(s.fem,'');rhs=gpuArray(s.rhs);
op=s.model.A;op.rowptr=gpuArray(op.rowptr);op.columns=gpuArray(op.columns);
wait(device);uploadSeconds=toc(timer);
routes=s.routes;reports=cell(1,numel(routes));solutions=cell(1,numel(routes));
for sample=1:numel(routes)
    route=routes(sample);wait(device);total=tic;
    nonlinear=radia_cuda_ams('nonlinear_create',fem);
    nonlinearCleanup=onCleanup(@()radia_cuda_ams('nonlinear_release',nonlinear));
    context=[];contextCleanup=[]; %#ok<NASGU>
    x=zeros(size(rhs),'like',rhs);history=[];
    builds=0;setupSeconds=0;solveSeconds=0;assemblySeconds=0;
    firstIterations=0;rebuildNext=false;peakSampleBytes=0;before=device.AvailableMemory;
    for outer=1:30
        wait(device);timer=tic;
        [residual,op.values,oldEnergy]=radia_cuda_ams('nonlinear_evaluate',nonlinear,x);
        wait(device);assemblySeconds=assemblySeconds+toc(timer);
        relative=gather(norm(residual)/norm(rhs));
        fprintf('sample=%d route=%d outer=%d residual=%.3e\n',sample,route,outer,relative);
        if relative<=1e-9;break;end
        rebuild=isempty(context)||route==1||(route==3&&rebuildNext);
        if rebuild
            wait(device);timer=tic;model=s.model;
            model.A.values=radia_cuda_ams('nonlinear_raw',nonlinear);
            [hierarchy,~]=radia.gpu.prepareAMS(model,SetupThreads=setupThreads);
            if isempty(context)
                context=radia.gpu.AMSContext(hierarchy);
                contextCleanup=onCleanup(@()delete(context));
            else
                context.update(hierarchy);
            end
            clear hierarchy;
            wait(device);setupSeconds=setupSeconds+toc(timer);builds=builds+1;firstIterations=0;
        end
        wait(device);timer=tic;
        [delta,info]=context.solve(op,-residual,1e-10,4000,true);
        wait(device);solveSeconds=solveSeconds+toc(timer);
        assert(info.workspace_allocations==1&&info.native_buffer_allocations_during_solve==0,'radia:nonlinear:workspace','PCG workspace was reallocated');
        if firstIterations==0;firstIterations=info.iterations;end
        rebuildNext=info.iterations>max(30,2*firstIterations);
        slope=gather(residual'*delta);alpha=1;
        while radia_cuda_ams('nonlinear_energy',nonlinear,x+alpha*delta)>oldEnergy+1e-4*alpha*slope
            alpha=alpha/2;
            assert(alpha>=2^-20,'radia:nonlinear:lineSearch','Line search failed');
        end
        x=x+alpha*delta;
        peakSampleBytes=max(peakSampleBytes,double(before-device.AvailableMemory));
        history=[history;outer relative info.iterations info.true_relative_residual alpha]; %#ok<AGROW>
    end
    assert(relative<=1e-9,'radia:nonlinear:convergence','Nonlinear residual failed');
    solutions{sample}=gather(x);wait(device);computeSeconds=toc(total);
    reports{sample}=struct('route',route,'sample',sample,'total_seconds',computeSeconds+uploadSeconds, ...
        'compute_seconds',computeSeconds,'upload_seconds',uploadSeconds,'input_owned_bytes',inputBytes, ...
        'device_memory_sample_delta_bytes',peakSampleBytes,'setup_seconds',setupSeconds,'assembly_seconds',assemblySeconds, ...
        'solve_seconds',solveSeconds,'builds',builds,'history',history,'final_residual',relative, ...
        'workspace_bytes',info.workspace_bytes,'workspace_allocations',info.workspace_allocations, ...
        'native_buffer_allocations_during_solve',info.native_buffer_allocations_during_solve,'resident_solve_count',info.solve_count,'resident_hierarchy_builds',info.hierarchy_builds,'graph_cached',info.graph_cached);
    clear contextCleanup context nonlinearCleanup;
    status=radia_cuda_ams('resident_status');
    assert(isequal(status,initialStatus),'radia:nonlinear:lifecycle','Resident contexts leaked');
end
gpuInfo=struct('name',device.Name,'compute_capability',device.ComputeCapability,'total_memory',device.TotalMemory);
save(outputPath,'solutions','reports','gpuInfo');
end
function out=upload(value,name)
if isstruct(value)
    out=value;names=fieldnames(value);
    for k=1:numel(names);out.(names{k})=upload(value.(names{k}),names{k});end
elseif isnumeric(value)&&~ismember(name,{'rows','cols','elements','nu0','gamma'})
    out=gpuArray(value);
else;out=value;
end
end
