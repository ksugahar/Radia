function ihPlantSFunction(block)
%IHPLANTSFUNCTION Simulink block for one discrete Radia IH plant step.
%   The dialog parameter is the struct returned by makeIHPlant.  This
%   Level-2 MATLAB S-function keeps the Simulink boundary thin; the plant
%   matrices remain owned by the Radia MATLAB API and can later be replaced
%   by a VIM/FEM/ESIM state update without changing the model ports.

setup(block);
end

function setup(block)
block.NumDialogPrms = 1;
plant = block.DialogPrm(1).Data;
validatePlant(plant);

block.NumInputPorts = 1;
block.NumOutputPorts = 1;
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions = 2;
block.InputPort(1).DatatypeID = 0;
block.InputPort(1).Complexity = "Real";
block.InputPort(1).DirectFeedthrough = true;

block.OutputPort(1).Dimensions = size(plant.C, 1);
block.OutputPort(1).DatatypeID = 0;
block.OutputPort(1).Complexity = "Real";

block.SampleTimes = [plant.sample_time_s, 0];
block.SimStateCompliance = "DefaultSimState";

block.RegBlockMethod("PostPropagationSetup", @postPropagationSetup);
block.RegBlockMethod("InitializeConditions", @initializeConditions);
block.RegBlockMethod("Outputs", @outputs);
block.RegBlockMethod("Update", @update);
end

function postPropagationSetup(block)
plant = block.DialogPrm(1).Data;
block.NumDworks = 1;
block.Dwork(1).Name = "state";
block.Dwork(1).Dimensions = numel(plant.x0);
block.Dwork(1).DatatypeID = 0;
block.Dwork(1).Complexity = "Real";
block.Dwork(1).UsedAsDiscState = true;
end

function initializeConditions(block)
plant = block.DialogPrm(1).Data;
block.Dwork(1).Data = plant.x0(:);
end

function outputs(block)
plant = block.DialogPrm(1).Data;
u = double(block.InputPort(1).Data(:));
x = double(block.Dwork(1).Data(:));
block.OutputPort(1).Data = plant.C * x + plant.D * u;
end

function update(block)
plant = block.DialogPrm(1).Data;
u = double(block.InputPort(1).Data(:));
x = double(block.Dwork(1).Data(:));
block.Dwork(1).Data = plant.A * x + plant.B * u;
end

function validatePlant(plant)
if ~isstruct(plant) || ~isfield(plant, "schema") || ...
        plant.schema ~= "radia.ih.simulink.plant.v1"
    error("radia:simulink:InvalidPlant", ...
        "The S-function parameter must be a Radia IH plant struct.");
end
if ~isequal(size(plant.A), [2, 2]) || ~isequal(size(plant.B), [2, 2]) || ...
        size(plant.C, 2) ~= 2 || size(plant.D, 2) ~= 2 || ...
        numel(plant.x0) ~= 2
    error("radia:simulink:InvalidPlant", ...
        "The Radia IH plant matrices have incompatible dimensions.");
end
end
