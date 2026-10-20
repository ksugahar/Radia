function checkSharedAMS(out,binaryPath,matlabPath,repeats)
% Common initialization; nontrivial graph/direct parity and final-25 memory.
assert(repeats>=100 && repeats==fix(repeats));addpath(binaryPath);addpath(matlabPath);
g=gpuDevice();wait(g);assert(g.AvailableMemory>=32*1024^3);
process=System.Diagnostics.Process.GetCurrentProcess();processCleanup=onCleanup(@()process.Dispose());
assert(~isfile(fullfile(out,'memory-samples.jsonl')));
initial=radia_cuda_ams('resident_status');s=radia.gpu.readAMSExchange(fullfile(out,'input.h5'));
total=tic;setupTimer=tic;[hierarchy,setup]=radia.gpu.prepareAMS(s.model,SetupThreads=1);
levels=cellfun(@(space)numel(space.levels),hierarchy.spaces);assert(max(levels)>=2,'Multilevel AMG not exercised');
b=gpuArray(s.rhs(:));ctx=radia.gpu.AMSContext(hierarchy);cleanup=onCleanup(@()delete(ctx));wait(g);setupSeconds=toc(setupTimer);
solutions=zeros(numel(s.rhs),2);infos=cell(1,2);seconds=zeros(1,2);
for j=1:2
    timer=tic;[x,info]=ctx.solve(hierarchy.A,b,1e-11,4000,j==2);wait(g);seconds(j)=toc(timer);
    solutions(:,j)=gather(x);infos{j}=info;clear x;
    assert(info.true_relative_residual<=1e-9 && info.native_buffer_allocations_during_solve==0);
end
for j=1:5;[x,~]=ctx.solve(hierarchy.A,b,1e-11,4000,true);clear x;end
host=zeros(1,repeats);device=zeros(1,repeats);
for j=1:repeats
    [x,info]=ctx.solve(hierarchy.A,b,1e-11,4000,true);
    assert(info.true_relative_residual<=1e-9 && info.native_buffer_allocations_during_solve==0);clear x;wait(g);
    process.Refresh();host(j)=double(process.PrivateMemorySize64);device(j)=g.AvailableMemory;
    fid=fopen(fullfile(out,'memory-samples.jsonl'),'a');assert(fid~=-1);
    fprintf(fid,'%s\n',jsonencode(struct('iteration',j,'private_bytes',host(j),'device_free_bytes',device(j))));fclose(fid);
end
last=host(end-24:end);fit=polyfit(1:25,last,1);delta=last(end)-last(1);
% Preserve numerical outputs and every sample before the memory verdict.
radia.gpu.writeHDivH5(fullfile(out,'solutions.h5'), ...
    struct('solutions',solutions,'host_samples',host,'device_free_samples',device), ...
    struct('units','SI','labels',{{'unit-cube'}},'warmup',5,'measured',repeats));
assert(fit(1)<=64*1024 && delta<=1024^2,'Host memory plateau failed');
assert(max(device(end-24:end))-min(device(end-24:end))<=1024^2,'Device plateau failed');
clear cleanup;clear hierarchy b;wait(g);assert(isequal(initial,radia_cuda_ams('resident_status')),'Native owners not restored');
result=struct('accepted',true,'setup',setup,'setup_and_upload_seconds',setupSeconds,'solve_seconds',seconds,'total_seconds',toc(total), ...
 'infos',{infos},'levels',levels,'warmup_solves',5,'repeat_solves',repeats,'host_samples',host,'device_free_samples',device, ...
 'private_counter','Windows Process.PrivateMemorySize64','final25_host_slope',fit(1),'final25_host_delta',delta,'owners_restored',true);
fid=fopen(fullfile(out,'gpu-result.json'),'w');assert(fid~=-1);fwrite(fid,jsonencode(result),'char');fclose(fid);
end
