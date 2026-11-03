function result=checkHDivDeviceMemory(inputsDir,outputDir)
% 100-run process and device memory slope of the HDiv device solve, under the
% limits of gateway_memory_protocol.json. inputsDir holds hdiv-device-setup.h5
% from `python -m radia_simulink.hdiv`. Writes result.json in outputDir.
% Passing the limits is not proof of zero leaks or an exact global peak.
protocol=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'gateway_memory_protocol.json')));
device=gpuDevice;wait(device);
[action,description]=makeGatewayMemoryWorkload('hdiv-device',inputsDir,outputDir);
samples=sampleGatewayMemory(action,device,protocol.warmup_iterations,protocol.measurement_iterations, ...
    fullfile(outputDir,'progress.jsonl'));
private=samples.samples(:,1);available=samples.samples(:,3);runs=(1:numel(private))';
window=runs>numel(runs)-protocol.plateau_window_iterations;
slope=@(values,keep)polyfit(runs(keep),values(keep),1)*[1;0];
limits=struct('slope',protocol.positive_slope_limit_bytes_per_run,'delta',protocol.positive_total_delta_limit_bytes, ...
    'plateau_slope',protocol.plateau_positive_slope_limit_bytes_per_run,'plateau_delta',protocol.plateau_positive_delta_limit_bytes);
all=true(size(runs));
measured=struct('private_slope_bytes_per_run',slope(private,all),'private_delta_bytes',private(end)-samples.baseline(1), ...
    'private_plateau_slope_bytes_per_run',slope(private,window),'private_plateau_delta_bytes',private(end)-private(find(window,1)), ...
    'device_available_delta_bytes',available(end)-samples.baseline(3));
result=struct('schema','radia.hdiv-device-memory.v1','description',description,'limits',limits,'measured',measured, ...
    'samples',samples,'final_live_count',radia_hdiv_device('live_count'),'accepted', ...
    measured.private_slope_bytes_per_run<=limits.slope && measured.private_delta_bytes<=limits.delta ...
    && measured.private_plateau_slope_bytes_per_run<=limits.plateau_slope ...
    && measured.private_plateau_delta_bytes<=limits.plateau_delta && radia_hdiv_device('live_count')==0);
stream=fopen(fullfile(outputDir,'result.json'),'w');assert(stream>0,'Cannot write the memory result');
closer=onCleanup(@()fclose(stream));fwrite(stream,jsonencode(result,PrettyPrint=true),'char');
end
