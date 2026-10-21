function result = checkCsrBackends(inputPath, outputPath, mode)
%CHECKCSRBACKENDS Independent CPU oracle, true residual, errors and repeated use.
% mode='cpu' calls no GPU initialization API; it does not substitute for API tracing.
assert(~isfile(outputPath) && ~isfile(string(outputPath)+".samples.jsonl"));
s = radia.gpu.readAMSExchange(inputPath);
result = struct('schema','radia.csr-backend-validation.v1','accepted',false, ...
    'mode',mode,'cuda_context_gate','unverified','cases',[]);
if strcmp(mode,'gpu')
    initialization=tic;device=gpuDevice; wait(device);
    result.gpu_initialization_seconds=toc(initialization);
elseif ~strcmp(mode,'cpu')
    error('radia:validation:mode','Expected cpu or gpu');
end
for k=1:numel(s.cases)
    f=s.cases{k};started=tic;args={int32(f.row_pointers),int32(f.columns),double(f.values),double(f.rhs),double(f.inverse_diagonal)};
    if strcmp(mode,'gpu'),args=cellfun(@gpuArray,args,'UniformOutput',false);end
    [x,info]=radia.csrPCG(args{:},1e-12,1000,string(mode));
    if strcmp(mode,'gpu'),x=gather(x);wait(device);end
    elapsed=toc(started);
    relative=norm(x(:)-f.expected(:))/max(norm(f.expected),realmin);
    assert(relative<=1e-11 && info.true_relative_residual<=1e-11);
    assert(strcmp(info.backend,mode));
    result.cases=[result.cases;struct('relative_l2',relative,'true_relative_residual',info.true_relative_residual, ...
        'iterations',info.iterations,'transfer_inclusive_seconds',elapsed)]; %#ok<AGROW>
end
f=s.cases{end};args={int32(f.row_pointers),int32(f.columns),double(f.values),double(f.rhs),double(f.inverse_diagonal)};
if strcmp(mode,'gpu'),args=cellfun(@gpuArray,args,'UniformOutput',false);end
process=System.Diagnostics.Process.GetCurrentProcess();
cleanup=onCleanup(@()process.Dispose());
for k=1:5,[x,info]=radia.csrPCG(args{:},1e-12,1000,string(mode));end
samples=zeros(100,1);
for k=1:100
    [x,info]=radia.csrPCG(args{:},1e-12,1000,string(mode));
    if strcmp(mode,'gpu'),wait(device);end
    process.Refresh();samples(k)=double(process.PrivateMemorySize64);
    fid=fopen(string(outputPath)+".samples.jsonl",'a');assert(fid~=-1);
    fprintf(fid,'%s\n',jsonencode(struct('iteration',k,'private_bytes',samples(k),'mode',mode)));fclose(fid);
end
fit=polyfit((1:25)',samples(end-24:end),1);
result.memory=struct('samples',samples,'final25_slope',fit(1),'final25_delta',samples(end)-samples(end-24));
radia.gpu.writeHDivH5(string(outputPath)+".memory.h5",struct('private_samples',samples), ...
    struct('units','SI','labels',{{}},'mode',mode,'warmup',5,'measured',100));
assert(fit(1)<=65536 && result.memory.final25_delta<=1048576);
% Zero RHS has no recursive-residual special case that can bypass validation.
zero=args;zero{4}=zeros(size(f.rhs));if strcmp(mode,'gpu'),zero{4}=gpuArray(zero{4});end
[x,info]=radia.csrPCG(zero{:},1e-12,1000,string(mode));if strcmp(mode,'gpu'),x=gather(x);end
assert(norm(x)==0 && info.true_relative_residual==0);
bad=args;bad{5}=-ones(size(f.inverse_diagonal));if strcmp(mode,'gpu'),bad{5}=gpuArray(bad{5});end
caught=false;
try,[~,~]=radia.csrPCG(bad{:},1e-12,1000,string(mode));catch err,caught=strcmp(err.identifier,'radia:cuda:solve');end
assert(caught);
clear radia_cuda_csr_pcg
result.accepted=true;
fid=fopen(outputPath,'w');assert(fid~=-1);guard=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(result),'char');
end
