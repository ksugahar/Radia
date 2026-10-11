function receipt = validate_cuda(checkout, binaryDirectory, outputJson)
%VALIDATE_CUDA Reproduce focused upstream parity and transfer-inclusive diagnostics.
arguments
    checkout (1,1) string
    binaryDirectory (1,1) string
    outputJson (1,1) string
end
addpath(fullfile(checkout,'packages','radia-optuna','matlab'));
addpath(binaryDirectory,'-begin');
binary=string(which('optuna_mex'));
assert(startsWith(binary,binaryDirectory+filesep), 'Unexpected MEX origin.');
info=optuna_mex('api.info');
assert(string(info.build_route)=="mexcuda" && ~info.taskmanager);
suite=testsuite(fullfile(checkout,'tests','matlab','test_optuna_upstream_oracle.m'));
selected=suite(~cellfun(@isempty,regexp({suite.Name}, ...
    'TPE|Random|Sobol|Mixed|MultiObjective|SingleDistribution')));
extra=testsuite(fullfile(checkout,'tests','matlab','test_optuna_cuda.m'));
results=run([selected(:);extra(:)]);
device=gpuDevice;
% Diagnostic timings on 100 are not formal performance acceptance.
samples=linspace(-2,3,16384)';
weights=ones(64,1)/64; mu=linspace(-1,2,64)'; sigma=ones(64,1);
call=@()optuna_mex('optuna.parzen.log_pdf_numerical',samples,weights,mu,sigma,-2,3,false,NaN);
t=tic; values=call(); first=toc(t);
times=zeros(5,1);
for index=1:numel(times), t=tic; unused=call(); times(index)=toc(t); end
assert(all(isfinite(values)));
[status,host]=system('hostname'); assert(status==0);
receipt=struct('schema','radia.validation.optuna-cuda.v1', ...
    'host',strtrim(host),'matlab',version,'gpu',device.Name, ...
    'compute_capability',device.ComputeCapability,'toolkit',device.ToolkitVersion, ...
    'driver',device.GraphicsDriverVersion,'api',info, ...
    'passed',sum([results.Passed]),'failed',sum([results.Failed]), ...
    'incomplete',sum([results.Incomplete]),'test_names',string({results.Name}), ...
    'timing_kind','100-host diagnostics; not performance acceptance', ...
    'density_samples',numel(samples),'mixture_components',numel(weights), ...
    'first_density_call_s',first,'warmed_median_s',median(times), ...
    'warmed_samples_s',times,'includes_allocation_upload_kernel_download',true, ...
    'isolated_transfer_time_measured',false,'speedup_claimed',false, ...
    'optuna_oracle','5.0.0','taskmanager_used',false);
folder=fileparts(outputJson); if ~isfolder(folder), mkdir(folder); end
fid=fopen(outputJson,'w'); assert(fid~=-1); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(receipt,PrettyPrint=true));
assert(receipt.failed==0 && receipt.incomplete==0, 'Optuna CUDA acceptance failed.');
end
