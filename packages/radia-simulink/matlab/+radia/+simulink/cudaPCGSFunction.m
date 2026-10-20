function cudaPCGSFunction(block)
%CUDAPCGSFUNCTION Upload fixed CSR data at Start and solve each sampled RHS on GPU.
block.NumDialogPrms = 5;
block.DialogPrmsTunable = repmat({'Nontunable'},1,5);
system = block.DialogPrm(1).Data;
required = {'row_pointers','columns','values','inverse_diagonal'};
if ~isstruct(system) || ~isscalar(system) || ~all(isfield(system,required))
    error('radia:gpu:CSRSystem','Expected a scalar CSR system struct.');
end
validateattributes(system.row_pointers,{'int32'},{'vector','nonempty'});
validateattributes(system.columns,{'int32'},{'vector','nonempty'});
validateattributes(system.values,{'double'},{'vector','real','finite','nonempty'});
validateattributes(system.inverse_diagonal,{'double'},{'vector','real','finite','positive','nonempty'});
validateattributes(block.DialogPrm(2).Data,{'double'},{'scalar','real','finite','positive','<=',1e-6});
validateattributes(block.DialogPrm(3).Data,{'double'},{'scalar','integer','positive','<=',100000});
validateattributes(block.DialogPrm(4).Data,{'double'},{'scalar','real','finite','positive'});
n = numel(system.inverse_diagonal);
block.NumInputPorts = 1; block.NumOutputPorts = 1;
block.InputPort(1).Dimensions = n; block.InputPort(1).DatatypeID = 0;
block.InputPort(1).Complexity = 'Real'; block.InputPort(1).DirectFeedthrough = true;
block.OutputPort(1).Dimensions = n; block.OutputPort(1).DatatypeID = 0;
block.OutputPort(1).Complexity = 'Real';
block.SampleTimes = [block.DialogPrm(4).Data,0];
block.SimStateCompliance = 'DisallowSimState';
block.RegBlockMethod('Start',@start);
block.RegBlockMethod('Outputs',@outputs);
block.RegBlockMethod('Terminate',@terminate);
end

function start(block)
folder = string(block.DialogPrm(5).Data);
if ~isfolder(folder), mkdir(folder); end
entry = struct('folder',folder,'calls',0,'max_true_residual',0,'passed',false);
created = false;
try
    radia.gpu.requireRuntime();
    host = block.DialogPrm(1).Data;
    entry.system = struct('row_pointers',gpuArray(host.row_pointers), ...
        'columns',gpuArray(host.columns),'values',gpuArray(host.values), ...
        'inverse_diagonal',gpuArray(host.inverse_diagonal));
    entry.passed = true;
    radia.gpu.internal.systemStore('create',getfullname(block.BlockHandle),entry);
    created = true;
    writeLog(folder,'CUDA MEX PCG started; matrix uploaded once; no CPU solver fallback.');
catch exception
    if created, radia.gpu.internal.systemStore('remove',getfullname(block.BlockHandle)); end
    writeResult(entry,false,exception.message);
    rethrow(exception);
end
end

function outputs(block)
key = getfullname(block.BlockHandle);
entry = radia.gpu.internal.systemStore('get',key);
try
    system = entry.system;
    [solution,info] = radia.gpu.csrPCG(system.row_pointers,system.columns, ...
        system.values,gpuArray(block.InputPort(1).Data),system.inverse_diagonal, ...
        block.DialogPrm(2).Data,block.DialogPrm(3).Data);
    block.OutputPort(1).Data = gather(solution);
    entry.calls = entry.calls+1;
    entry.max_true_residual = max(entry.max_true_residual,info.true_relative_residual);
    radia.gpu.internal.systemStore('update',key,entry);
catch exception
    entry.passed = false;
    radia.gpu.internal.systemStore('update',key,entry);
    writeResult(entry,false,exception.message);
    rethrow(exception);
end
end

function terminate(block)
entry = radia.gpu.internal.systemStore('remove',getfullname(block.BlockHandle));
if isempty(entry), return; end
if entry.passed, writeResult(entry,true,''); end
end

function writeResult(entry,passed,message)
record = struct('schema','radia.gpu.simulink-pcg.v1','passed',passed, ...
    'backend','MATLAB mexcuda','calls',entry.calls, ...
    'max_true_residual',entry.max_true_residual,'error',message);
file = fopen(fullfile(entry.folder,'result.json'),'w');
if file < 0, error('radia:gpu:Evidence','Cannot write CUDA block result.'); end
cleanup = onCleanup(@() fclose(file));
fprintf(file,'%s\n',jsonencode(record));
clear cleanup
writeLog(entry.folder,sprintf('passed=%d calls=%d max_true_residual=%.17g %s', ...
    passed,entry.calls,entry.max_true_residual,message));
end

function writeLog(folder,message)
file = fopen(fullfile(folder,'run.log'),'a');
if file < 0, error('radia:gpu:Evidence','Cannot write CUDA block log.'); end
cleanup = onCleanup(@() fclose(file));
fprintf(file,'%s\n',message);
end
