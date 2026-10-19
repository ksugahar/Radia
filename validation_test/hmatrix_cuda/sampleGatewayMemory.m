function result=sampleGatewayMemory(action,device,warmCount,measureCount,progressPath,observerOnly)
% Snapshot inside one Engine call, after operation-local owners have exited.
% Windows Process.PrivateMemorySize64 is an OS counter, not MATLAB memory().
if nargin<6,observerOnly=false;end
assert(ispc && warmCount==5 && measureCount==100);
process=System.Diagnostics.Process.GetCurrentProcess();
owner=onCleanup(@()process.Dispose());
warm=zeros(warmCount,3);samples=zeros(measureCount,3);
% Initialize the common sampler and its durable progress IO before warmup.
for k=1:5,snapshot();progress('sampler',k);end
for k=1:warmCount
    action();warm(k,:)=snapshot();progress('warmup',k);
end
baseline=snapshot();
for k=1:measureCount
    if ~observerOnly,action();end
    samples(k,:)=snapshot();
    if mod(k,10)==0,progress('measurement',k);end
end
result=struct('schema','radia.gateway-inprocess-samples.v1', ...
    'columns',{{'private_committed_bytes','working_set_bytes','device_available_bytes'}}, ...
    'sampler_warmup_count',5,'baseline',baseline,'warmup',warm,'samples',samples, ...
    'observer_only',observerOnly,'engine_calls_during_measurement',0, ...
    'private_counter','System.Diagnostics.Process.PrivateMemorySize64', ...
    'scope','Process snapshots, including runtime caches; not exact owned or global peak');
    function value=snapshot()
        wait(device);process.Refresh();
        value=[double(process.PrivateMemorySize64),double(process.WorkingSet64),double(device.AvailableMemory)];
    end
    function progress(stage,step)
        fid=fopen(progressPath,'a');assert(fid>=0);
        closer=onCleanup(@()fclose(fid));
        fprintf(fid,'{"stage":"%s","step":%d}\n',stage,step);
    end
end
