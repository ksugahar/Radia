function ihThermalSFunction(block)
%IHTHERMALSFUNCTION Advance accepted temperature from distributed heating.
setup(block);
end

function setup(block)
block.NumDialogPrms = 1;
block.DialogPrmsTunable = {'Nontunable'};
config = radia.simulink.validateIHNativeConfig(block.DialogPrm(1).Data);
block.NumInputPorts = 3;
block.NumOutputPorts = 1;
for k = 1:3
    block.InputPort(k).DatatypeID = 0;
    block.InputPort(k).Complexity = "Real";
    block.InputPort(k).DirectFeedthrough = false;
end
block.InputPort(1).Dimensions = config.n_heat; % heat_density_W_per_m3
block.InputPort(2).Dimensions = 1; % ambient_temperature_K
block.InputPort(3).Dimensions = 1; % angle_rad
block.OutputPort(1).Dimensions = config.n_temperature; % temperature_K
block.OutputPort(1).DatatypeID = 0;
block.OutputPort(1).Complexity = "Real";
block.SampleTimes = [config.sample_time_s, 0];
block.SimStateCompliance = 'DisallowSimState';
block.RegBlockMethod("PostPropagationSetup", @postPropagationSetup);
block.RegBlockMethod("Start", @start);
block.RegBlockMethod("InitializeConditions", @initializeConditions);
block.RegBlockMethod("Outputs", @outputs);
block.RegBlockMethod("Update", @update);
block.RegBlockMethod("Terminate", @terminate);
end

function postPropagationSetup(block)
block.NumDworks = 2;
for k = 1:2
    block.Dwork(k).Dimensions = 1;
    block.Dwork(k).DatatypeID = 7; % uint32
    block.Dwork(k).Complexity = "Real";
    block.Dwork(k).UsedAsDiscState = false;
end
block.Dwork(1).Name = "native_handle_low";
block.Dwork(2).Name = "native_handle_high";
end

function start(block)
radia.setup(RequireMex=true);
radia.simulink.requireIHNativeRuntime();
destroyHandle(block);
config = radia.simulink.validateIHNativeConfig(block.DialogPrm(1).Data);
setHandle(block, radia_mex('ih.thermal.create', config));
end

function initializeConditions(block)
h = getHandle(block);
if h ~= 0, radia_mex('ih.thermal.reset', h); end
end

function outputs(block)
block.OutputPort(1).Data = radia_mex('ih.thermal.output', requireHandle(block));
end

function update(block)
radia_mex('ih.thermal.update', requireHandle(block), ...
    double(block.InputPort(1).Data(:)), double(block.InputPort(2).Data), ...
    double(block.InputPort(3).Data));
end

function terminate(block)
destroyHandle(block);
end

function destroyHandle(block)
h = getHandle(block);
if h ~= 0
    try
        radia_mex('ih.thermal.destroy', h);
    catch
    end
    setHandle(block, uint64(0));
end
end

function h = requireHandle(block)
h = getHandle(block);
if h == 0, error("radia:simulink:IHThermalHandle", "IH Thermal handle is not initialized."); end
end

function h = getHandle(block)
h = bitor(uint64(block.Dwork(1).Data), ...
    bitshift(uint64(block.Dwork(2).Data), 32));
end

function setHandle(block, h)
block.Dwork(1).Data = uint32(bitand(h, uint64(4294967295)));
block.Dwork(2).Data = uint32(bitshift(h, -32));
end
