function motorAngleFamilyMexSFunction(block)
%MOTORANGLEFAMILYMEXSFUNCTION Native periodic motor ROM block.
%   Input layout: [mechanical_angle_rad; model_inputs].
%   Output layout: [linear_outputs; electromagnetic_torque_Nm].

setup(block);
end

function setup(block)
block.NumDialogPrms = 1;
family = block.DialogPrm(1).Data;
validateFamily(family);

block.NumInputPorts = 1;
block.NumOutputPorts = 1;
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;
block.InputPort(1).Dimensions = family.input_count + 1;
block.InputPort(1).DatatypeID = 0;
block.InputPort(1).Complexity = "Real";
block.InputPort(1).DirectFeedthrough = true;
block.OutputPort(1).Dimensions = family.output_count;
block.OutputPort(1).DatatypeID = 0;
block.OutputPort(1).Complexity = "Real";
block.SampleTimes = [family.sample_time_s, 0];
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
family = block.DialogPrm(1).Data;
radia.setup();
nativeHandle = radia_mex('simulink.state_space.create', ...
    family.angle_grid_rad, family.period_rad, ...
    family.A, family.B, family.C, family.D, ...
    family.Q, family.R, family.S, family.x0);
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
    error("radia:simulink:MotorAngleFamilyHandle", ...
        "The native motor angle-family handle is not initialized.");
end
input = double(block.InputPort(1).Data(:));
block.OutputPort(1).Data = radia_mex( ...
    'simulink.state_space.step', nativeHandle, input(1), input(2:end));
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

function validateFamily(family)
required = ["schema", "angle_grid_rad", "period_rad", "sample_time_s", ...
    "snapshot_count", "state_order", "input_count", ...
    "linear_output_count", "output_count", "A", "B", "C", "D", ...
    "Q", "R", "S", "x0"];
if ~isstruct(family) || ~all(isfield(family, required)) || ...
        string(family.schema) ~= "radia.motor.periodic-angle-family.v1"
    error("radia:simulink:MotorAngleFamily", ...
        "The parameter must come from makeMotorAngleFamily.");
end
end
