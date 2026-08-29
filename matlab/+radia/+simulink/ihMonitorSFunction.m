function ihMonitorSFunction(block)
%IHMONITORSFUNCTION Reduce distributed IH fields to scalar telemetry.
setup(block);
end

function setup(block)
block.NumDialogPrms = 1;
block.DialogPrmsTunable = {'Nontunable'};
config = radia.simulink.validateIHNativeConfig(block.DialogPrm(1).Data);
block.NumInputPorts = 6;
block.NumOutputPorts = 15;
for index = 1:6
    block.InputPort(index).DatatypeID = 0;
    block.InputPort(index).Complexity = "Real";
    block.InputPort(index).DirectFeedthrough = true;
end
block.InputPort(1).Dimensions = config.n_heat;
block.InputPort(2).Dimensions = config.n_temperature;
for index = 3:6
    block.InputPort(index).Dimensions = 1;
end
for index = 1:15
    block.OutputPort(index).Dimensions = 1;
    block.OutputPort(index).DatatypeID = 0;
    block.OutputPort(index).Complexity = "Real";
end
block.SampleTimes = [config.sample_time_s, 0];
block.SimStateCompliance = 'DisallowSimState';
block.RegBlockMethod("PostPropagationSetup", @postPropagationSetup);
block.RegBlockMethod("Start", @start);
block.RegBlockMethod("Outputs", @outputs);
block.RegBlockMethod("Update", @update);
end

function postPropagationSetup(block)
block.NumDworks = 1;
block.Dwork(1).Name = "sample_index";
block.Dwork(1).Dimensions = 1;
block.Dwork(1).DatatypeID = 0;
block.Dwork(1).Complexity = "Real";
block.Dwork(1).UsedAsDiscState = true;
end

function start(block)
block.Dwork(1).Data = 0;
end

function outputs(block)
config = radia.simulink.validateIHNativeConfig(block.DialogPrm(1).Data);
heat = double(block.InputPort(1).Data(:));
temperature = double(block.InputPort(2).Data(:));
heatWeights = double(config.heat_cell_weights(:));
temperatureWeights = double(config.temperature_cell_weights(:));
revision = double(block.InputPort(6).Data);
if ~isfinite(revision)
    revision = 0;
end

block.OutputPort(1).Data = 2; % running / healthy
block.OutputPort(2).Data = double(block.CurrentTime);
block.OutputPort(3).Data = revision;
block.OutputPort(4).Data = 0; % last successfully emitted sample
block.OutputPort(5).Data = block.Dwork(1).Data;
block.OutputPort(6).Data = double(block.InputPort(3).Data);
block.OutputPort(7).Data = double(block.InputPort(4).Data);
block.OutputPort(8).Data = double(block.InputPort(5).Data);
block.OutputPort(9).Data = min(heat);
block.OutputPort(10).Data = weightedMean(heat, heatWeights);
block.OutputPort(11).Data = max(heat);
block.OutputPort(12).Data = dot(heat, heatWeights);
block.OutputPort(13).Data = min(temperature);
block.OutputPort(14).Data = weightedMean(temperature, temperatureWeights);
block.OutputPort(15).Data = max(temperature);
end

function update(block)
block.Dwork(1).Data = block.Dwork(1).Data + 1;
end

function value = weightedMean(field, weights)
denominator = sum(weights);
if ~(isfinite(denominator) && denominator > 0)
    error("radia:simulink:IHMonitorWeights", ...
        "IH monitor weights must have a finite positive sum.");
end
value = dot(field, weights) / denominator;
end
