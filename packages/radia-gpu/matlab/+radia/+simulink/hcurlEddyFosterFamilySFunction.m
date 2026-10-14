function hcurlEddyFosterFamilySFunction(block)
%HCURLEDDYFOSTERFAMILYSFUNCTION Moving shared-mode Foster HCurl block.
%   Input layout: [minus_dI_dt(port_count); height_m; coil_current(port_count)].
%   Output layout: [port_response(port_count); force_x; force_y; force_z].
%   The discrete state matrix is shared by every height; Outputs reads the
%   modal state and Update advances it with the height-interpolated input.
%   force is the instantaneous product K*z(t)*i(t); its cycle average is the
%   phasor time-average force.

setup(block);
end

function setup(block)
block.NumDialogPrms = 1;
family = block.DialogPrm(1).Data;
validateFamily(family);
nPort = family.port_count;

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
block.Dwork(1).Name = "modal_state";
block.Dwork(1).Dimensions = family.state_order;
block.Dwork(1).DatatypeID = 0;
block.Dwork(1).Complexity = "Real";
block.Dwork(1).UsedAsDiscState = true;
end

function initializeConditions(block)
family = block.DialogPrm(1).Data;
block.Dwork(1).Data = family.x0(:);
end

function outputs(block)
family = block.DialogPrm(1).Data;
u = double(block.InputPort(1).Data(:));
nPort = family.port_count;
model = radia.simulink.interpolateHCurlEddyFosterFamily(family, u(nPort + 1));
z = double(block.Dwork(1).Data(:));
response = model.Cd * z + model.Dd * u(1:nPort);
force = radia.internal.hcurlEddyFosterInstantaneousForce( ...
    model.force_operator, z, u(nPort + 2:end));
block.OutputPort(1).Data = [response(:); force(:)];
end

function update(block)
family = block.DialogPrm(1).Data;
u = double(block.InputPort(1).Data(:));
nPort = family.port_count;
model = radia.simulink.interpolateHCurlEddyFosterFamily(family, u(nPort + 1));
z = double(block.Dwork(1).Data(:));
block.Dwork(1).Data = model.Ad * z + model.Bd * u(1:nPort);
end

function validateFamily(family)
if ~isstruct(family) || ~isfield(family, "schema") || ...
        family.schema ~= "radia.hcurl.eddy_foster.family.v1" || ...
        ~isfield(family, "shared_modes") || ~family.shared_modes || ...
        family.snapshot_count < 1
    error("radia:simulink:HCurlFosterFamily", ...
        "The S-function parameter must be a loaded Foster HCurl family.");
end
if family.state_order < 1 || family.port_count < 1 || ...
        numel(family.positions_m) ~= family.snapshot_count
    error("radia:simulink:HCurlFosterFamily", ...
        "The Foster HCurl family dimensions are inconsistent.");
end
end
