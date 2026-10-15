function highorder_ic_fixture(fixture,binaryDir,output)
% Solver-only validation input: Python assembly is not a GPU product adapter.
addpath(binaryDir,'-begin');data=load(fixture);dev=gpuDevice;wait(dev);
initial=radia_cuda_ic('live_count');timer=tic;
op=struct('rows',double(data.rows),'cols',double(data.rows), ...
    'rowptr',gpuArray(int32(data.rowptr)),'columns',gpuArray(int32(data.columns)), ...
    'values',gpuArray(data.values),'free',gpuArray(int32(data.free)));
rhs=gpuArray(data.rhs);wait(dev);uploadSeconds=toc(timer);
timer=tic;[h,setup]=radia_cuda_ic('create',op,1.05);wait(dev);factorSeconds=toc(timer);
cleanup=onCleanup(@()radia_cuda_ic('release',h));
scales=[1 2 .5 1];seconds=zeros(4,1);info=cell(4,1);solution=[];
for k=1:4
    wait(dev);timer=tic;[x,info{k}]=radia_cuda_ic('solve',h,scales(k)*rhs,1e-12,20000);
    value=gather(x);wait(dev);seconds(k)=toc(timer);solution(:,k)=value; %#ok<AGROW>
end
clear cleanup;assert(radia_cuda_ic('live_count')==initial);
report=struct('scope','Generic CSR GPU IC solver-only; CPU-assembled validation fixture', ...
    'upload_seconds',uploadSeconds,'preconditioner_setup_seconds',factorSeconds, ...
    'solve_seconds',seconds,'scales',scales,'info',{info},'setup',setup, ...
    'matlab',version,'device',dev.Name,'cuda_toolkit',dev.ToolkitVersion, ...
    'graphics_driver',dev.GraphicsDriverVersion,'registry_restored',true);
save(output,'solution','report','-v7');
end
