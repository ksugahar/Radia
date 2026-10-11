function benchmark(inputPath, outputPath, binaryDir)
addpath(binaryDir);
d=load(inputPath);
g=gpuDevice();
before=g.AvailableMemory;
memoryPlan=memory_plan(d.snapshot,g);
check_memory_budget(memoryPlan);
wait(g); transfer=tic;
s=upload(d.snapshot, '');
rhs=gpuArray(d.rhs);
wait(g);
uploadSeconds=toc(transfer);
wait(g); timer=tic;
[first,info]=radia_cuda_ams(s,rhs,1e-9,2000,false);
wait(g); coldSeconds=toc(timer);
probeErrors=zeros(size(d.probes,1),1);
for k=1:size(d.probes,1)
    [actual,~]=radia_cuda_ams(s,gpuArray(d.probes(k,:)'),1e-9,2000,true);
    actual=gather(actual);
    expected=d.expected(k,:)';
    probeErrors(k)=norm(actual-expected)/max(norm(expected),realmin);
end
times=zeros(5,1);
for k=1:5
    timer=tic; [x,info]=radia_cuda_ams(s,rhs,1e-9,2000,false);wait(g);times(k)=toc(timer);
end
timer=tic;
s2=upload(d.snapshot,'');b2=gpuArray(d.rhs);
[x2,~]=radia_cuda_ams(s2,b2,1e-9,2000,false);
inclusiveSolution=gather(x2);wait(g);inclusiveSeconds=toc(timer);
solution=gather(x);
memoryBytes=before-g.AvailableMemory;
% Zero RHS keeps gateway setup but skips PCG iterations. Apply-only includes
% one AMS cycle and its own gateway setup. These are diagnostic wall times,
% not exact allocations/kernel-only counters; retain their first samples.
zeroRhs=zeros(size(rhs),'like',rhs);
zeroRhsTimes=zeros(6,1);applyOnlyTimes=zeros(6,1);
for k=1:6
    timer=tic;[z,zeroInfo]=radia_cuda_ams(s,zeroRhs,1e-9,2000,false);wait(g);zeroRhsTimes(k)=toc(timer);
    assert(zeroInfo.iterations==0,'Zero RHS unexpectedly ran PCG iterations');
    timer=tic;[applied,~]=radia_cuda_ams(s,rhs,1e-9,2000,true);wait(g);applyOnlyTimes(k)=toc(timer);
end
zeroRhsWarmMedianSeconds=median(zeroRhsTimes(2:end));
applyOnlyWarmMedianSeconds=median(applyOnlyTimes(2:end));
checks=check_errors(s,rhs);
% Exercise rejection with no allocation and no deliberately oversized GPU job.
limited=memory_plan(d.snapshot,struct('TotalMemory',g.TotalMemory,'AvailableMemory',0));
checks.memory_budget_rejection=false;
try;check_memory_budget(limited);catch e;checks.memory_budget_rejection=strcmp(e.identifier,'radia:cuda:memoryBudget');end
assert(checks.memory_budget_rejection,'Memory budget guard failed');
device=struct('name',g.Name,'compute_capability',g.ComputeCapability, ...
    'total_memory',g.TotalMemory,'driver',g.GraphicsDriverVersion,'toolkit',g.ToolkitVersion);
matlabVersion=version;
save(outputPath,'solution','info','probeErrors','coldSeconds','times', ...
    'inclusiveSeconds','inclusiveSolution','memoryBytes','uploadSeconds','checks','device','matlabVersion', ...
    'zeroRhsTimes','applyOnlyTimes','zeroRhsWarmMedianSeconds','applyOnlyWarmMedianSeconds','memoryPlan');
end

function check_memory_budget(plan)
assert(plan.estimatedAdditionalBytes<=plan.budgetBytes, ...
    'radia:cuda:memoryBudget','AMS benchmark exceeds the conservative GPU memory budget');
end

function plan=memory_plan(s,g)
% Allow duplicate uploads, diagnostic copy-on-write arrays, solver work and
% library headroom. This is an admission estimate, not a measured peak.
inputBytes=numeric_bytes(s,'')+8*double(s.A.rows);
workspaceBytes=8*10*double(s.A.rows);
for k=1:numel(s.spaces)
    sp=s.spaces{k};
    workspaceBytes=workspaceBytes+8*(double(sp.P.rows)+double(sp.P.cols));
    for j=1:numel(sp.levels)
        l=sp.levels{j};workspaceBytes=workspaceBytes+8*4*double(l.A.rows);
        if isfield(l,'coarse') && strcmp(l.coarse.kind,'ngsolve_sparse_ldlt')
            workspaceBytes=workspaceBytes+12*numel(l.coarse.diagonal_inverse);
        end
    end
end
requested=3*inputBytes+workspaceBytes+8*6*double(s.A.rows);
headroom=max(256*1024^2,0.25*requested);
reserve=max(2*1024^3,0.2*double(g.TotalMemory));
plan=struct('inputBytes',inputBytes,'workspaceEstimateBytes',workspaceBytes, ...
    'estimatedAdditionalBytes',requested+headroom,'headroomBytes',headroom, ...
    'reservedFreeBytes',reserve,'availableBeforeBytes',double(g.AvailableMemory), ...
    'budgetBytes',max(0,double(g.AvailableMemory)-reserve), ...
    'method','Conservative admission estimate; includes duplicate benchmark uploads; not allocator peak');
end
function bytes=numeric_bytes(a,name)
bytes=0;
if isstruct(a)
    names=fieldnames(a);
    for k=1:numel(names);bytes=bytes+numeric_bytes(a.(names{k}),names{k});end
elseif iscell(a)
    for k=1:numel(a);bytes=bytes+numeric_bytes(a{k},'');end
elseif isnumeric(a) && ~ismember(name,{'rows','cols','num_smooth','cycle_type','correction_weight'})
    item=whos('a');bytes=double(item.bytes);
end
end

function out=upload(a,name)
if isstruct(a)
    out=a;names=fieldnames(a);
    for k=1:numel(names);out.(names{k})=upload(a.(names{k}),names{k});end
elseif iscell(a)
    out=cell(size(a));for k=1:numel(a);out{k}=upload(a{k},'');end
elseif isnumeric(a) && ~ismember(name,{'rows','cols','num_smooth','cycle_type','correction_weight'})
    out=gpuArray(a);
else;out=a;
end
end

function checks=check_errors(s,b)
checks=struct();
bad=s;bad.A.columns(1)=int32(s.A.cols);
checks.invalid_csr=reject(@()invoke(bad,b,1e-9,2000));
checks.cpu_rhs=reject(@()invoke(s,gather(b),1e-9,2000));
bad=b;bad(1)=NaN;checks.nan_rhs=reject(@()invoke(s,bad,1e-9,2000));
checks.nonconvergence=reject(@()invoke(s,b,1e-15,1));
bad=s;bad.A.rowptr(2)=int32(-1);checks.invalid_rowptr=reject(@()invoke(bad,b,1e-9,2000));
bad=s;bad.l1(1)=0;checks.invalid_smoother=reject(@()invoke(bad,b,1e-9,2000));
checks.wrong_rhs_shape=reject(@()invoke(s,b(1:end-1),1e-9,2000));
checks.invalid_tolerance=reject(@()invoke(s,b,0,2000));
bad=s;bad.spaces{1}.levels{end}.coarse.order(1)=int32(s.A.rows);
checks.invalid_coarse_order=reject(@()invoke(bad,b,1e-9,2000));
[z,info]=radia_cuda_ams(s,zeros(size(b),'like',b),1e-9,2000,false);
checks.zero_rhs=gather(norm(z))==0 && info.iterations==0;
clear radia_cuda_ams;
[x,info]=radia_cuda_ams(s,b,1e-9,2000,false);
checks.clear_reload=all(isfinite(gather(x))) && info.true_relative_residual<=1e-9;
assert(all(structfun(@(x)logical(x),checks)),'AMS gateway checks failed');
end
function invoke(s,b,tol,iterations)
[~,~]=radia_cuda_ams(s,b,tol,iterations,false);
end
function yes=reject(f)
yes=false;
try;f();catch e;yes=strcmp(e.identifier,'radia:cuda:ams');end
end
