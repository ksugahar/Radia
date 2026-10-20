function adjointOptimizationSFunction(block)
%ADJOINTOPTIMIZATIONSFUNCTION Execute one checked MMA/SQP study on a trigger.
setup(block);
end

function setup(block)
block.NumDialogPrms = 2;
block.DialogPrmsTunable = {'Nontunable','Nontunable'};
block.NumInputPorts = 1;
block.NumOutputPorts = 4;
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;
block.InputPort(1).Dimensions = 1;
block.InputPort(1).DirectFeedthrough = true;

runner = block.DialogPrm(1).Data;
if ~isa(runner,"radia.topopt.AdjointRunner")
    error("radia:simulink:AdjointRunner", ...
        "The first dialog parameter must be a radia.topopt.AdjointRunner.");
end
n = numel(runner.InitialDesign);
block.OutputPort(1).Dimensions = 1;
block.OutputPort(2).Dimensions = 1;
block.OutputPort(3).Dimensions = 1;
block.OutputPort(4).Dimensions = n;
for index = 1:block.NumOutputPorts
    block.OutputPort(index).DatatypeID = 0;
end
block.SampleTimes = [double(block.DialogPrm(2).Data), 0];
block.RegBlockMethod('PostPropagationSetup', @postSetup);
block.RegBlockMethod('Start', @start);
block.RegBlockMethod('Outputs', @outputs);
block.RegBlockMethod('Terminate', @terminate);
end

function postSetup(block)
runner = block.DialogPrm(1).Data;
n = numel(runner.InitialDesign);
names = {'previous_trigger','objective','status','iterations','design'};
block.NumDworks = numel(names);
for index = 1:numel(names)
    block.Dwork(index).Name = names{index};
    if index == 5
        dimensions = n;
    else
        dimensions = 1;
    end
    block.Dwork(index).Dimensions = dimensions;
    block.Dwork(index).DatatypeID = 0;
    block.Dwork(index).Complexity = 'Real';
    block.Dwork(index).UsedAsDiscState = true;
end
end

function start(block)
for index = 1:block.NumDworks
    block.Dwork(index).Data = zeros(block.Dwork(index).Dimensions,1);
end
block.Dwork(2).Data = NaN;
block.Dwork(3).Data = 0;
block.Dwork(4).Data = 0;
block.Dwork(5).Data = block.DialogPrm(1).Data.InitialDesign;
end

function outputs(block)
trigger = double(block.InputPort(1).Data);
if trigger > 0 && block.Dwork(1).Data <= 0
    runner = block.DialogPrm(1).Data;
    block.Dwork(3).Data = 1;
    try
        result = runner.run();
        block.Dwork(2).Data = result.objective;
        block.Dwork(4).Data = height(result.history) - 1;
        block.Dwork(5).Data = result.design;
        block.Dwork(3).Data = 2;
    catch exception
        block.Dwork(3).Data = -1;
        warning("radia:simulink:AdjointFailed", ...
            "Adjoint optimization failed: %s", exception.message);
    end
end
block.Dwork(1).Data = trigger;
block.OutputPort(1).Data = block.Dwork(2).Data;
block.OutputPort(2).Data = block.Dwork(3).Data;
block.OutputPort(3).Data = block.Dwork(4).Data;
block.OutputPort(4).Data = block.Dwork(5).Data;
end

function terminate(block)
block.Dwork(3).Data = 0;
end
