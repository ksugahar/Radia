function geometry_reuse(volPath,payload,output,options)
% Paired path/retained-geometry material sweep; independent CPU gate is external.
% Reset tests require a fresh owned MATLAB process; never run in a shared session.
arguments
    volPath (1,1) string
    payload (1,1) string
    output (1,1) string
    options.AllowDeviceReset (1,1) logical = false
    options.LifecycleTestBinary (1,1) logical = false
end
addpath(payload);g=gpuDevice;wait(g);initial=radia_cuda_ams('resident_status');
labels=["back","front","left","right","top","bottom"];
copyPath=string(output)+".geometry.vol";
assert(~isfile(copyPath),'radia:test:output','Owned snapshot path already exists');
copyfile(volPath,copyPath);fileCleanup=onCleanup(@()removeOwned(copyPath));
timer=tic;geometry=radia.gpu.HCurlGeometry(copyPath,labels);wait(g);
geometrySeconds=toc(timer);geometryCleanup=onCleanup(@()delete(geometry));
geometryInfo=geometry.Info;delete(copyPath); % The live immutable snapshot needs no file reads.
dirichlet=char(geometry.Dirichlet);
parameters=[1 1;2 1;.75 2];threads=[1 2 16];records=cell(3,1);solutions=cell(3,1);
for k=1:3
    nu=parameters(k,1);mass=parameters(k,2);current=[k -.25*k .5];
    [coldModel,coldFem]=radia.gpu.assembleHCurl(volPath,labels, ...
        Reluctivity=nu,MassCoefficient=mass,CurrentDensity=current);
    [warmModel,warmFem,assembly]=radia.gpu.assembleHCurl(geometry,labels, ...
        Reluctivity=nu,MassCoefficient=mass,CurrentDensity=current);
    assert(assembly.geometry_reused && assembly.upload_bytes==0);
    % Ten mapped inputs plus mass CSR rowptr/columns; rowptr is checked twice.
    assert(assembly.geometry_liveness_check_calls==12);
    assert(all(cell2mat(struct2cell(assembly.host_staging_phases))==0));
    exact={'rowptr','columns'};
    for field=exact
        assert(isequal(coldModel.A.(field{1}),warmModel.A.(field{1})));
        assert(isequal(coldModel.G.(field{1}),warmModel.G.(field{1})));
    end
    assert(isequal(coldModel.free,warmModel.free));
    errors=zeros(1,7);
    errors(1:3)=[difference(coldModel.A.values,warmModel.A.values), ...
        difference(coldModel.G.values,warmModel.G.values), ...
        difference(coldModel.coordinates,warmModel.coordinates)];
    for field={'dofs','slots','free'}
        assert(isequal(gather(coldFem.(field{1})),gather(warmFem.(field{1}))));
    end
    numericFields={'curls','volume','rhs'};
    for index=1:numel(numericFields)
        name=numericFields{index};errors(3+index)=difference(coldFem.(name),warmFem.(name));
    end
    errors(7)=difference(coldFem.mass.values,warmFem.mass.values);
    assert(max(errors)<=1e-14,'radia:test:geometryParity','Retained geometry changed element/operator data');
    % Exercise MATLAB copy-on-write on arrays cloned from persistent geometry.
    warmFem.dofs(1)=int32(-1);warmFem.slots(1)=int32(-1);warmFem.free(1)=int32(1)-warmFem.free(1);
    warmFem.curls(1)=warmFem.curls(1)+1;warmFem.volume(1)=2*warmFem.volume(1);
    warmFem.mass.rowptr(1)=int32(1);warmFem.mass.columns(1)=int32(-1);
    warmModel.A.values(1)=warmModel.A.values(1)+1;
    warmModel.G.values(1)=warmModel.G.values(1)+1;
    [afterModel,afterFem]=geometry.assemble(Reluctivity=nu,MassCoefficient=mass,CurrentDensity=current);
    assert(bitwise(coldModel.A.values,afterModel.A.values));
    assert(bitwise(coldModel.G.values,afterModel.G.values));
    assert(bitwise(coldFem.mass.rowptr,afterFem.mass.rowptr));
    assert(bitwise(coldFem.mass.columns,afterFem.mass.columns));
    for field={'dofs','slots','free','curls','volume','rhs'}
        assert(bitwise(coldFem.(field{1}),afterFem.(field{1})), ...
            'radia:test:geometryMutation','Returned array mutation corrupted persistent geometry');
    end
    clear coldModel coldFem warmModel warmFem afterModel afterFem;
    timer=tic;
    baseline=radia.gpu.HCurlSystem(volPath,labels,Reluctivity=nu, ...
        MassCoefficient=mass,CurrentDensity=current,SetupThreads=threads(k));
    baselineCleanup=onCleanup(@()delete(baseline));
    [reference,baselineSolve]=baseline.solve(Tolerance=1e-12,ReturnOnGPU=false);
    wait(g);baselineTotal=toc(timer);baselineSetup=baseline.SetupInfo;
    clear baselineCleanup baseline;
    timer=tic;
    retained=radia.gpu.HCurlSystem(geometry,labels,Reluctivity=nu, ...
        MassCoefficient=mass,CurrentDensity=current,SetupThreads=threads(k));
    retainedCleanup=onCleanup(@()delete(retained));
    [value,first]=retained.solve(Tolerance=1e-12,ReturnOnGPU=false);
    wait(g);materialTotal=toc(timer);setup=retained.SetupInfo;
    [reused,reuse]=retained.solve(CurrentScale=2,Tolerance=1e-12,ReturnOnGPU=false);
    assert(norm(value-reference)/norm(reference)<=1e-10);
    assert(norm(reused-2*reference)/norm(2*reference)<=1e-10);
    solutions{k}=struct('current',current,'reluctivity',nu,'mass',mass,'value',value,'scaled_value',reused);
    records{k}=struct('threads',threads(k),'element_operator_max_absolute_difference',max(errors), ...
        'baseline_fresh_seconds',baselineTotal,'retained_material_fresh_seconds',materialTotal, ...
        'baseline_setup',baselineSetup,'setup',setup,'first',first,'reuse',reuse,'baseline_solve',baselineSolve, ...
        'retained_assembly',assembly);
    clear retainedCleanup retained;
end
clear geometryCleanup geometry fileCleanup;
lockedBefore=mislocked('radia_cuda_ams');assert(~lockedBefore);
freeBeforeLifecycle=g.AvailableMemory;
for repeat=1:5
    [owned,~]=radia_cuda_ams('geometry_create',char(volPath),cellstr(labels),'.*');
    assert(mislocked('radia_cuda_ams'));
    radia_cuda_ams('geometry_release',owned);
    assert(~mislocked('radia_cuda_ams'));
end
% Raw gateway ownership/type/generation checks use the same binary.
[h,~]=radia_cuda_ams('geometry_create',char(volPath),cellstr(labels),'.*');
reject(@()badAssembly(double(h)));reject(@()badAssembly(uint64(0)));
radia_cuda_ams('geometry_release',h);reject(@()badAssembly(h));
reject(@()radia_cuda_ams('geometry_release',h));
deviceBefore=g.Index;mismatchChecked=false;allocationRollbackChecked=false;
if options.LifecycleTestBinary
    [h,~]=radia_cuda_ams('geometry_create',char(volPath),cellstr(labels),'.*');
    radia_cuda_ams('geometry_test_device_mismatch',h);
    reject(@()badAssembly(h));radia_cuda_ams('geometry_release',h);
    assert(gpuDevice().Index==deviceBefore && ~mislocked('radia_cuda_ams'));
    mismatchChecked=true;
    radia_cuda_ams('geometry_test_fail_create');reject(@()badCreate(volPath,labels));
    assert(~mislocked('radia_cuda_ams'));allocationRollbackChecked=true;
end
wait(g);freeAfterLifecycle=g.AvailableMemory;
assert(freeAfterLifecycle>=freeBeforeLifecycle-16*2^20, ...
    'radia:test:geometryLeak','Repeated geometry create/release retained unexpected GPU memory');
resetChecked=false;
if options.AllowDeviceReset
    [h,~]=radia_cuda_ams('geometry_create',char(volPath),cellstr(labels),'.*');
    reset(gpuDevice);
    rejectReset(@()badAssembly(h));rejectReset(@()badAssembly(h));
    radia_cuda_ams('geometry_release',h);resetChecked=true;
end
wait(g);assert(isequal(initial,radia_cuda_ams('resident_status')));
report=struct('schema','radia.hcurl-geometry-reuse-validation.v1','internal_parity_accepted',true, ...
    'independent_cpu_accepted',false,'timing_accepted',false,'policy_closed',false, ...
    'geometry_creation_seconds',geometrySeconds,'geometry',geometryInfo,'records',{records}, ...
    'parameters',parameters,'current_sources',[(1:3).' -.25*(1:3).' .5*ones(3,1)], ...
    'boundary_labels',{cellstr(labels)},'dirichlet',dirichlet, ...
    'returned_array_mutation_bitwise_verified',true,'persistent_arrays_across_calls_verified',true, ...
    'lock_before',lockedBefore,'lock_after',mislocked('radia_cuda_ams'), ...
    'lifecycle_free_before',freeBeforeLifecycle,'lifecycle_free_after',freeAfterLifecycle, ...
    'device_mismatch_test_executed',mismatchChecked,'allocation_rollback_test_executed',allocationRollbackChecked, ...
    'device_reset_test_executed',resetChecked, ...
    'registry_restored',true,'matlab',version,'device',g.Name,'loaded_binary',which('radia_cuda_ams'), ...
    'scope','Paired GPU parity and lifecycle only; independent CPU/source/module hash gates remain external');
save(output,'solutions','report','-v7');
file=fopen(string(output)+".json",'w');assert(file>=0);jsonCleanup=onCleanup(@()fclose(file));
fprintf(file,'%s',jsonencode(report,PrettyPrint=true));
end
function equal=bitwise(a,b)
a=gather(a);b=gather(b);
equal=isequal(size(a),size(b)) && strcmp(class(a),class(b)) && ...
    isequal(typecast(a(:),'uint8'),typecast(b(:),'uint8'));
end
function badCreate(path,labels)
[~,~]=radia_cuda_ams('geometry_create',char(path),cellstr(labels),'.*');
end
function error=difference(a,b)
a=gather(a);b=gather(b);assert(isequal(size(a),size(b)));
error=max(abs(a(:)-b(:)),[],'all');assert(isfinite(error));
end
function rejectReset(operation)
failed=false;
try
    operation();
catch failure
    assert(contains(failure.message,'geometry invalidated by device reset'));
    failed=true;
end
assert(failed,'radia:test:geometryReset','Reset geometry was accepted');
end
function badAssembly(h)
[~,~,~]=radia_cuda_ams('assemble_hcurl_geometry',h,1,1,[1 0 0],0);
end
function reject(operation)
failed=false;
try
    operation();
catch
    failed=true;
end
assert(failed,'radia:test:geometryReject','Unsupported geometry handle was accepted');
end
function removeOwned(path)
if isfile(path);delete(path);end
end
