function [hierarchy,info] = prepareAMS(model,options)
%PREPAREAMS Construct AMS within mexcuda; no Python adapter or snapshot input.
% model: radia.ams-native-input.v1, A/G int32 CSR, coordinates (vertices x 3),
% free int32 mask, and cycle (1 or 7). NGSolve owns FE orientation/assembly.
% By default returns caller-owned GPU arrays ready for radia.gpu.amsPCG.
arguments
    model (1,1) struct
    options.SetupThreads (1,1) double {mustBeInteger,mustBeInRange(options.SetupThreads,1,64)} = 1
    options.Upload (1,1) logical = true
    options.SetupOnGPU (1,1) logical = true
end
[hierarchy,info]=radia_cuda_ams('setup',model,options.SetupThreads,options.SetupOnGPU);
info.upload_seconds=0;
if options.Upload
    device=gpuDevice();
    inputBytes=numeric_bytes(hierarchy,'');
    % Conservative admission estimate, including input duplication and work.
    estimated=1.25*(3*inputBytes+8*16*double(hierarchy.A.rows));
    reserve=max(2*1024^3,0.2*double(device.TotalMemory));
    assert(estimated<=max(0,double(device.AvailableMemory)-reserve), ...
        'radia:cuda:memoryBudget','AMS hierarchy exceeds the conservative GPU memory budget');
    wait(device);timer=tic;hierarchy=upload(hierarchy,'');wait(device);
    info.upload_seconds=toc(timer);
end
end
function bytes=numeric_bytes(value,name)
bytes=0;
if isstruct(value)
    names=fieldnames(value);
    for k=1:numel(names);bytes=bytes+numeric_bytes(value.(names{k}),names{k});end
elseif iscell(value)
    for k=1:numel(value);bytes=bytes+numeric_bytes(value{k},'');end
elseif isnumeric(value) && ~ismember(name,{'rows','cols','num_smooth','cycle_type','correction_weight'})
    item=whos('value');bytes=double(item.bytes);
end
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
