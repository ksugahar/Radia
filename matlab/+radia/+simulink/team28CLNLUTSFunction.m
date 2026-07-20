function team28CLNLUTSFunction(block)
%TEAM28CLNLUTSFUNCTION Level-2 S-function for the TEAM 28 CLN LUT.

setup(block);
end

function setup(block)
block.NumDialogPrms = 1;
lut = block.DialogPrm(1).Data;
validateLUT(lut);

block.NumInputPorts = 1;
block.NumOutputPorts = 1;
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;
block.InputPort(1).Dimensions = 2;
block.InputPort(1).DatatypeID = 0;
block.InputPort(1).Complexity = "Real";
block.InputPort(1).DirectFeedthrough = true;
block.OutputPort(1).Dimensions = 3;
block.OutputPort(1).DatatypeID = 0;
block.OutputPort(1).Complexity = "Real";
block.SampleTimes = [lut.sample_time_s, 0];
block.SimStateCompliance = "DefaultSimState";
block.RegBlockMethod("Outputs", @outputs);
end

function outputs(block)
lut = block.DialogPrm(1).Data;
u = double(block.InputPort(1).Data(:));
[force_N, lift_N, slope_N_per_m] = radia.simulink.evaluateTeam28CLNForce( ...
    lut, u(1), u(2));
block.OutputPort(1).Data = [force_N; lift_N; slope_N_per_m];
end

function validateLUT(lut)
if ~isstruct(lut) || ~isfield(lut, "schema") || ...
        lut.schema ~= "radia.team28.cln_lut.v1" || ...
        ~isfield(lut, "sample_time_s")
    error("radia:simulink:Team28LUT", ...
        "The S-function parameter must be a TEAM 28 CLN LUT struct.");
end
end
