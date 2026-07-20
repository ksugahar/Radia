function stateSpaceMexSFunction(block)
%STATESPACEMEXSFUNCTION Native MEX-backed discrete state-space block.
%   The matrices and initial state cross the MATLAB/MEX boundary once at
%   Start. Each sample then advances the native state and returns only y.

setup(block);
end

function setup(block)
block.NumDialogPrms = 1;
model = block.DialogPrm(1).Data;
validateModel(model);

block.NumInputPorts = 1;
block.NumOutputPorts = 1;
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;
block.InputPort(1).Dimensions = size(model.B, 2);
block.InputPort(1).DatatypeID = 0;
block.InputPort(1).Complexity = "Real";
block.InputPort(1).DirectFeedthrough = true;
block.OutputPort(1).Dimensions = size(model.C, 1);
block.OutputPort(1).DatatypeID = 0;
block.OutputPort(1).Complexity = "Real";
block.SampleTimes = [model.sample_time_s, 0];
block.SimStateCompliance = "DefaultSimState";

block.RegBlockMethod("PostPropagationSetup", @postPropagationSetup);
block.RegBlockMethod("Start", @start);
block.RegBlockMethod("InitializeConditions", @initializeConditions);
block.RegBlockMethod("Outputs", @outputs);
block.RegBlockMethod("Terminate", @terminate);
end

function postPropagationSetup(block)
block.NumDworks = 1;
block.Dwork(1).Name = "native_handle";
block.Dwork(1).Dimensions = 1;
block.Dwork(1).DatatypeID = 0;
block.Dwork(1).Complexity = "Real";
block.Dwork(1).UsedAsDiscState = false;
end

function start(block)
model = block.DialogPrm(1).Data;
radia.setup();
nativeHandle = radia_mex('simulink.state_space.create', ...
    model.A, model.B, model.C, model.D, model.x0);
block.Dwork(1).Data = double(nativeHandle);
end

function initializeConditions(block)
nativeHandle = getNativeHandle(block);
if nativeHandle ~= 0
    radia_mex('simulink.state_space.reset', nativeHandle);
end
end

function outputs(block)
nativeHandle = getNativeHandle(block);
if nativeHandle == 0
    error("radia:simulink:StateSpaceHandle", ...
        "The native state-space handle is not initialized.");
end
u = double(block.InputPort(1).Data(:));
block.OutputPort(1).Data = radia_mex( ...
    'simulink.state_space.step', nativeHandle, u);
end

function terminate(block)
nativeHandle = getNativeHandle(block);
if nativeHandle ~= 0
    try
        radia_mex('simulink.state_space.destroy', nativeHandle);
    catch
    end
    block.Dwork(1).Data = 0;
end
end

function nativeHandle = getNativeHandle(block)
nativeHandle = uint64(block.Dwork(1).Data);
end

function validateModel(model)
if ~isstruct(model) || ~isfield(model, "A") || ~isfield(model, "B") || ...
        ~isfield(model, "C") || ~isfield(model, "D") || ...
        ~isfield(model, "x0") || ~isfield(model, "sample_time_s")
    error("radia:simulink:StateSpaceModel", ...
        "The S-function parameter must contain A, B, C, D, x0, and sample_time_s.");
end
n = size(model.A, 1);
if ndims(model.A) ~= 2 || size(model.A, 2) ~= n || n < 1 || ...
        size(model.B, 1) ~= n || size(model.C, 2) ~= n || ...
        size(model.D, 1) ~= size(model.C, 1) || ...
        size(model.D, 2) ~= size(model.B, 2) || numel(model.x0) ~= n || ...
        ~isscalar(model.sample_time_s) || ~isfinite(model.sample_time_s) || ...
        model.sample_time_s <= 0 || any(~isfinite(model.A), "all") || ...
        any(~isfinite(model.B), "all") || any(~isfinite(model.C), "all") || ...
        any(~isfinite(model.D), "all") || any(~isfinite(model.x0), "all")
    error("radia:simulink:StateSpaceModel", ...
        "The native state-space dimensions or values are invalid.");
end
end
