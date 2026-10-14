function nonlinearReactorSFunction(block)
%NONLINEARREACTORSFUNCTION Level-2 wrapper for the native HDiv-MMM reactor.
setup(block);
end

function setup(block)
block.NumDialogPrms = 1;
block.DialogPrmsTunable = {'Nontunable'};
config = radia.simulink.validateNonlinearReactorConfig(block.DialogPrm(1).Data);
block.NumInputPorts = 1;
block.NumOutputPorts = 8;
block.InputPort(1).Dimensions = 1;
block.InputPort(1).DatatypeID = 0;
block.InputPort(1).Complexity = "Real";
block.InputPort(1).DirectFeedthrough = true;
for port = 1:7
    block.OutputPort(port).Dimensions = 1;
    block.OutputPort(port).DatatypeID = 0;
    block.OutputPort(port).Complexity = "Real";
end
block.OutputPort(8).Dimensions = double(config.n_samples);
block.OutputPort(8).DatatypeID = 0;
block.OutputPort(8).Complexity = "Real";
block.SampleTimes = [double(config.sample_time_s),0];
block.SimStateCompliance = "CustomSimState";
block.RegBlockMethod("PostPropagationSetup",@postPropagationSetup);
block.RegBlockMethod("Start",@start);
block.RegBlockMethod("InitializeConditions",@initializeConditions);
block.RegBlockMethod("Outputs",@outputs);
block.RegBlockMethod("Update",@update);
block.RegBlockMethod("Terminate",@terminate);
block.RegBlockMethod("GetSimState",@getSimState);
block.RegBlockMethod("SetSimState",@setSimState);
end

function postPropagationSetup(block)
block.NumDworks = 2;
for index = 1:2
    block.Dwork(index).Dimensions = 1;
    block.Dwork(index).DatatypeID = 7; % uint32
    block.Dwork(index).Complexity = "Real";
    block.Dwork(index).UsedAsDiscState = false;
end
block.Dwork(1).Name = "native_handle_low";
block.Dwork(2).Name = "native_handle_high";
end

function start(block)
radia.setup(RequireMex=true);
destroyHandle(block);
config = radia.simulink.validateNonlinearReactorConfig(block.DialogPrm(1).Data);
setHandle(block,radia_mex('reactor.create',config));
end

function initializeConditions(block)
handle = getHandle(block);
if handle ~= 0
    radia_mex('reactor.reset',handle);
end
end

function outputs(block)
result = radia_mex('reactor.output',requireHandle(block), ...
    double(block.InputPort(1).Data));
block.OutputPort(1).Data = result.voltage_V;
block.OutputPort(2).Data = result.flux_linkage_Wb_turn;
block.OutputPort(3).Data = result.differential_inductance_H;
block.OutputPort(4).Data = result.peak_flux_density_T;
block.OutputPort(5).Data = result.magnetic_energy_J;
block.OutputPort(6).Data = result.nonlinear_iterations;
block.OutputPort(7).Data = result.residual_relative_norm;
block.OutputPort(8).Data = result.flux_density_T;
end

function update(block)
radia_mex('reactor.update',requireHandle(block), ...
    double(block.InputPort(1).Data));
end

function terminate(block)
destroyHandle(block);
end

function state = getSimState(block)
state = radia_mex('reactor.snapshot',requireHandle(block));
end

function setSimState(block,state)
radia_mex('reactor.restore',requireHandle(block),state);
end

function handle = requireHandle(block)
handle = getHandle(block);
if handle == 0
    error("radia:simulink:NonlinearReactorHandle", ...
        "The nonlinear reactor native handle is not initialized.");
end
end

function handle = getHandle(block)
handle = bitor(uint64(block.Dwork(1).Data), ...
    bitshift(uint64(block.Dwork(2).Data),32));
end

function setHandle(block,handle)
block.Dwork(1).Data = uint32(bitand(handle,uint64(4294967295)));
block.Dwork(2).Data = uint32(bitshift(handle,-32));
end

function destroyHandle(block)
handle = getHandle(block);
if handle ~= 0
    try
        radia_mex('reactor.destroy',handle);
    catch
    end
    setHandle(block,uint64(0));
end
end
