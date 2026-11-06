function deviceRefresh(out,binaryPath,matlabPath)
%DEVICEREFRESH Device AMS refresh against the CPU HypreBasedAMS reuse Update.
% The CPU oracle builds the same hierarchy with reuse_hierarchy and updates it
% with the new values on the unchanged pattern. The device refreshes the resident
% hierarchy in place. Level operators and l1 norms are compared value by value;
% the coarsest inverse is checked through PCG on the refreshed operator.
addpath(binaryPath);addpath(matlabPath);
g=gpuDevice();wait(g);initial=radia_cuda_ams('resident_status');
s=radia.gpu.readAMSExchange(fullfile(out,'input.h5'));model=s.model;
[h1,setup]=radia.gpu.prepareAMS(model,SetupThreads=1,SetupOnGPU=false,ReuseHierarchy=true);
assert(h1.reuse_hierarchy,'Setup did not keep the reuse_hierarchy flag');
for name=["rows","cols","rowptr","columns"]
    assert(isequal(gather(h1.A.(name)(:)),model.A.(name)(:)),'BC operator pattern differs from the input pattern');
end
refreshed=radia_cuda_ams('create',h1);frozen=radia_cuda_ams('create',h1);
cleanup=onCleanup(@()release([refreshed frozen]));
n=double(h1.A.rows);free=gather(h1.free(:))~=0;b=gpuArray(sin((1:n)').*free);
updates=s.updates;results=cell(1,numel(updates));
for k=1:numel(updates)
    op=h1.A;op.values=gpuArray(updates{k}(:));
    wait(g);timer=tic;info=radia_cuda_ams('refresh',refreshed,op);wait(g);wall=toc(timer);
    device=radia_cuda_ams('resident_values',refreshed);
    timer=tic;[reference,referenceInfo]=radia_cuda_ams('setup',model,1,false,true,updates{k}(:));cpu=toc(timer);
    levels=compareLevels(device,reference);
    % Every route solves the same BC-modified operator that the CPU oracle exports.
    bc=h1.A;bc.values=gpuArray(reference.A.values(:));
    oracle=radia_cuda_ams('create',upload(reference,''));
    oracleCleanup=onCleanup(@()radia_cuda_ams('release',oracle));
    [xr,ir]=radia_cuda_ams('solve',refreshed,bc,b,1e-10,4000,false);
    [xo,io]=radia_cuda_ams('solve',oracle,bc,b,1e-10,4000,false);
    [~,iFrozen]=radia_cuda_ams('solve',frozen,bc,b,1e-10,4000,false);
    assert(ir.true_relative_residual<=1e-9 && io.true_relative_residual<=1e-9,'Refreshed or oracle PCG missed the true residual gate');
    difference=gather(norm(xr-xo))/gather(norm(xo));
    results{k}=struct('update',k,'refresh',info,'refresh_wall_seconds',wall, ...
        'cpu_setup_and_update_seconds',cpu,'cpu_info',referenceInfo,'levels',levels, ...
        'iterations_refreshed',ir.iterations,'iterations_cpu_oracle',io.iterations,'iterations_frozen',iFrozen.iterations, ...
        'residual_refreshed',ir.true_relative_residual,'residual_cpu_oracle',io.true_relative_residual, ...
        'residual_frozen',iFrozen.true_relative_residual,'solution_relative_difference',difference);
    clear oracleCleanup xr xo;
end
clear cleanup;wait(g);
assert(isequal(initial,radia_cuda_ams('resident_status')),'Native owners not restored');
result=struct('setup',setup,'updates',{results},'owners_restored',true);
fid=fopen(fullfile(out,'gpu-result.json'),'w');assert(fid~=-1);fwrite(fid,jsonencode(result),'char');fclose(fid);
end
function release(handles)
for h=handles;radia_cuda_ams('release',h);end
end
function rows=compareLevels(device,reference)
rows=struct('space',{},'level',{},'quantity',{},'count',{},'mismatches',{},'max_abs',{},'max_relative',{},'coarse_kind',{});
rows(end+1)=entry(0,0,'A',device.A,reference.A.values,'');
rows(end+1)=entry(0,0,'l1',device.l1,reference.l1,'');
for k=1:numel(device.spaces)
    levels=device.spaces{k};expected=reference.spaces{k}.levels;
    assert(numel(levels)==numel(expected),'Level count differs from the CPU hierarchy');
    for l=1:numel(levels)
        kind='';if ~isempty(levels{l}.coarse_kind);kind=levels{l}.coarse_kind;end
        rows(end+1)=entry(k,l,'A',levels{l}.A,expected{l}.A.values,kind); %#ok<AGROW>
        rows(end+1)=entry(k,l,'l1',levels{l}.l1,expected{l}.l1,kind); %#ok<AGROW>
    end
end
end
function row=entry(space,level,quantity,device,expected,kind)
d=gather(device(:));e=expected(:);assert(numel(d)==numel(e),'Hierarchy value count differs');
assert(all(isfinite(d)),'Nonfinite refreshed hierarchy value');
delta=abs(d-e);scale=max(abs(e));if scale==0;scale=1;end
row=struct('space',space,'level',level,'quantity',quantity,'count',numel(e),'mismatches',nnz(d~=e), ...
    'max_abs',max([delta;0]),'max_relative',max([delta;0])/scale,'coarse_kind',kind);
end
function out=upload(value,name)
if isstruct(value)
    out=value;names=fieldnames(value);
    for k=1:numel(names);out.(names{k})=upload(value.(names{k}),names{k});end
elseif iscell(value)
    out=cell(size(value));for k=1:numel(value);out{k}=upload(value{k},'');end
elseif isnumeric(value) && ~ismember(name,{'rows','cols','num_smooth','cycle_type','correction_weight'})
    out=gpuArray(value);
else
    out=value;
end
end
