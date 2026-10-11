function benchmark(inputPath, outputPath, binaryDir)
addpath(binaryDir);
d=load(inputPath);
g=gpuDevice();
before=g.AvailableMemory;
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
device=struct('name',g.Name,'compute_capability',g.ComputeCapability, ...
    'total_memory',g.TotalMemory,'driver',g.GraphicsDriverVersion,'toolkit',g.ToolkitVersion);
matlabVersion=version;
save(outputPath,'solution','info','probeErrors','coldSeconds','times', ...
    'inclusiveSeconds','inclusiveSolution','memoryBytes','uploadSeconds','checks','device','matlabVersion', ...
    'zeroRhsTimes','applyOnlyTimes','zeroRhsWarmMedianSeconds','applyOnlyWarmMedianSeconds');
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
