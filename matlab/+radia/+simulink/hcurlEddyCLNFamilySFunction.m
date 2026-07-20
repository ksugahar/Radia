function hcurlEddyCLNFamilySFunction(block)
%HCURLEDDYCLNFAMILYSFUNCTION Moving common-basis HCurl CLN block.
%   Input layout: [minus_dI_dt(port_count); height_m; coil_current(port_count)].
%   Output layout: [port_response(port_count); force_x; force_y; force_z].

setup(block);
end

function setup(block)
block.NumDialogPrms = 1;
family = block.DialogPrm(1).Data;
validateFamily(family);
nPort = family.port_count;
nState = family.state_order;

block.NumInputPorts = 1;
block.NumOutputPorts = 1;
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;
block.InputPort(1).Dimensions = 2 * nPort + 1;
block.InputPort(1).DatatypeID = 0;
block.InputPort(1).Complexity = "Real";
block.InputPort(1).DirectFeedthrough = true;
block.OutputPort(1).Dimensions = nPort + 3;
block.OutputPort(1).DatatypeID = 0;
block.OutputPort(1).Complexity = "Real";
block.SampleTimes = [family.sample_time_s, 0];
block.SimStateCompliance = "DefaultSimState";

block.RegBlockMethod("PostPropagationSetup", @postPropagationSetup);
block.RegBlockMethod("InitializeConditions", @initializeConditions);
block.RegBlockMethod("Outputs", @outputs);
block.RegBlockMethod("Update", @update);
end

function postPropagationSetup(block)
family = block.DialogPrm(1).Data;
block.NumDworks = 1;
block.Dwork(1).Name = "state";
block.Dwork(1).Dimensions = family.state_order;
block.Dwork(1).DatatypeID = 0;
block.Dwork(1).Complexity = "Real";
block.Dwork(1).UsedAsDiscState = true;
end

function initializeConditions(block)
family = block.DialogPrm(1).Data;
block.Dwork(1).Data = family.models{1}.x0(:);
end

function outputs(block)
family = block.DialogPrm(1).Data;
u = double(block.InputPort(1).Data(:));
nPort = family.port_count;
height = u(nPort + 1);
coilCurrent = u(nPort + 2:end);
x = double(block.Dwork(1).Data(:));
model = radia.simulink.interpolateHCurlEddyCLNFamily( ...
    family, height, BuildStateSpace=false);
response = model.Cd * x + model.Dd * u(1:nPort);
force = radia.simulink.evaluateHCurlEddyCLNForce(model, x, coilCurrent);
block.OutputPort(1).Data = [response(:); force(:)];
end

function update(block)
family = block.DialogPrm(1).Data;
u = double(block.InputPort(1).Data(:));
nPort = family.port_count;
height = u(nPort + 1);
model = radia.simulink.interpolateHCurlEddyCLNFamily( ...
    family, height, BuildStateSpace=false);
x = double(block.Dwork(1).Data(:));
block.Dwork(1).Data = model.Ad * x + model.Bd * u(1:nPort);
end

function validateFamily(family)
if ~isstruct(family) || ~isfield(family, "schema") || ...
        family.schema ~= "radia.hcurl.eddy_cln.family.v1" || ...
        ~isfield(family, "models") || family.snapshot_count < 1
    error("radia:simulink:HCurlCLNFamily", ...
        "The S-function parameter must be a loaded HCurl CLN family.");
end
if numel(family.models) ~= family.snapshot_count || ...
        family.state_order < 1 || family.port_count < 1
    error("radia:simulink:HCurlCLNFamily", ...
        "The HCurl CLN family dimensions are inconsistent.");
end
end
