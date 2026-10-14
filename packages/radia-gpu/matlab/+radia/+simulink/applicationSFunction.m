function applicationSFunction(block)
%APPLICATIONSFUNCTION Explicit-trigger adapter for Radia batch applications.

setup(block);
end

function setup(block)
block.NumDialogPrms = 5;
block.NumInputPorts = 1;
block.NumOutputPorts = 3;
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions = 1;
block.InputPort(1).DatatypeID = 8; % boolean
block.InputPort(1).Complexity = "Real";
block.InputPort(1).DirectFeedthrough = false;

block.OutputPort(1).Dimensions = 1;
block.OutputPort(1).DatatypeID = 6; % int32 status: 0 idle, 2 passed, -1 failed
block.OutputPort(1).Complexity = "Real";
block.OutputPort(2).Dimensions = 1;
block.OutputPort(2).DatatypeID = 0; % primary result value
block.OutputPort(2).Complexity = "Real";
block.OutputPort(3).Dimensions = 1;
block.OutputPort(3).DatatypeID = 0; % wall time [s]
block.OutputPort(3).Complexity = "Real";

block.SampleTimes = [-1, 0];
block.SimStateCompliance = "DefaultSimState";
block.RegBlockMethod("PostPropagationSetup", @postPropagationSetup);
block.RegBlockMethod("InitializeConditions", @initializeConditions);
block.RegBlockMethod("Outputs", @outputs);
block.RegBlockMethod("Update", @update);
end

function postPropagationSetup(block)
block.NumDworks = 4;
configureDwork(block.Dwork(1), "previous_trigger", 8);
configureDwork(block.Dwork(2), "status", 6);
configureDwork(block.Dwork(3), "primary", 0);
configureDwork(block.Dwork(4), "elapsed_s", 0);
end

function configureDwork(dwork, name, datatypeID)
dwork.Name = name;
dwork.Dimensions = 1;
dwork.DatatypeID = datatypeID;
dwork.Complexity = "Real";
dwork.UsedAsDiscState = true;
end

function initializeConditions(block)
block.Dwork(1).Data = false;
block.Dwork(2).Data = int32(0);
block.Dwork(3).Data = NaN;
block.Dwork(4).Data = 0.0;
end

function outputs(block)
block.OutputPort(1).Data = int32(block.Dwork(2).Data);
block.OutputPort(2).Data = double(block.Dwork(3).Data);
block.OutputPort(3).Data = double(block.Dwork(4).Data);
end

function update(block)
trigger = logical(block.InputPort(1).Data);
previous = logical(block.Dwork(1).Data);
block.Dwork(1).Data = trigger;
if ~trigger || previous
    return
end

application = string(block.DialogPrm(1).Data);
configFile = string(block.DialogPrm(2).Data);
runRoot = string(block.DialogPrm(3).Data);
timeout_s = double(block.DialogPrm(4).Data);
pythonExecutable = string(block.DialogPrm(5).Data);

try
    result = radia.simulink.runApplication(application, configFile, ...
        RunRoot=runRoot, Timeout_s=timeout_s, ...
        PythonExecutable=pythonExecutable, ThrowOnFailure=false);
    if string(result.status) == "passed"
        block.Dwork(2).Data = int32(2);
    else
        block.Dwork(2).Data = int32(-1);
    end
    block.Dwork(3).Data = primaryValue(result);
    block.Dwork(4).Data = double(result.elapsed_s);
    set_param(block.BlockHandle, "UserData", result, ...
        "UserDataPersistent", "off");
catch exception
    block.Dwork(2).Data = int32(-1);
    block.Dwork(3).Data = NaN;
    block.Dwork(4).Data = 0.0;
    diagnostic = struct("status", "failed", ...
        "error", string(exception.message), ...
        "identifier", string(exception.identifier));
    set_param(block.BlockHandle, "UserData", diagnostic, ...
        "UserDataPersistent", "off");
    warning("radia:simulink:ApplicationBlockFailed", "%s", exception.message);
end
end

function value = primaryValue(result)
value = NaN;
if ~isfield(result, "primary") || ~isstruct(result.primary) || ...
        ~isfield(result.primary, "value") || isempty(result.primary.value)
    return
end
candidate = double(result.primary.value);
if isscalar(candidate) && isfinite(candidate)
    value = candidate;
end
end
