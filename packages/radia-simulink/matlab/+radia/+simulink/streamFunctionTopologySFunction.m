function streamFunctionTopologySFunction(block)
%STREAMFUNCTIONTOPOLOGYSFUNCTION Run and expose one analytic SF study.
setup(block);
end

function setup(block)
block.NumDialogPrms = 3;
block.DialogPrmsTunable = {'Nontunable','Nontunable','Nontunable'};
block.NumInputPorts = 1;
block.NumOutputPorts = 8;
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;
block.InputPort(1).Dimensions = 1;
block.InputPort(1).DirectFeedthrough = true;

runner = block.DialogPrm(1).Data;
validateRunner(runner);
designCount = numel(runner.InitialDesign);
historyCapacity = validateHistoryCapacity(block.DialogPrm(3).Data);
dimensions = [1,1,1,designCount,1,historyCapacity, ...
    historyCapacity,historyCapacity];
for index = 1:block.NumOutputPorts
    block.OutputPort(index).Dimensions = dimensions(index);
    block.OutputPort(index).DatatypeID = 0;
end
block.SampleTimes = [double(block.DialogPrm(2).Data),0];
block.RegBlockMethod('PostPropagationSetup',@postSetup);
block.RegBlockMethod('Start',@start);
block.RegBlockMethod('Outputs',@outputs);
end

function postSetup(block)
runner = block.DialogPrm(1).Data;
designCount = numel(runner.InitialDesign);
historyCapacity = validateHistoryCapacity(block.DialogPrm(3).Data);
names = {'previous_trigger','objective','status','iterations','design', ...
    'history_count','history_iteration','history_objective', ...
    'history_constraint'};
dimensions = [1,1,1,1,designCount,1,historyCapacity, ...
    historyCapacity,historyCapacity];
block.NumDworks = numel(names);
for index = 1:numel(names)
    block.Dwork(index).Name = names{index};
    block.Dwork(index).Dimensions = dimensions(index);
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
block.Dwork(5).Data = block.DialogPrm(1).Data.InitialDesign;
for index = 7:9
    block.Dwork(index).Data = NaN(block.Dwork(index).Dimensions,1);
end
end

function outputs(block)
trigger = double(block.InputPort(1).Data);
if trigger > 0 && block.Dwork(1).Data <= 0
    runner = block.DialogPrm(1).Data;
    block.Dwork(3).Data = 1;
    try
        result = runner.run();
        history = result.history;
        capacity = block.Dwork(7).Dimensions;
        if height(history) > capacity
            error("radia:simulink:StreamFunctionTopologyHistoryCapacity", ...
                "Optimization history has %d rows but HistoryCapacity is %d.", ...
                height(history),capacity);
        end
        count = height(history);
        iteration = NaN(block.Dwork(7).Dimensions,1);
        objective = NaN(block.Dwork(8).Dimensions,1);
        constraint = NaN(block.Dwork(9).Dimensions,1);
        iteration(1:count) = history.Iteration(1:count);
        objective(1:count) = history.Objective(1:count);
        constraint(1:count) = history.MaxConstraint(1:count);
        block.Dwork(2).Data = result.objective;
        block.Dwork(3).Data = 2*double(result.converged) - ...
            2*double(~result.converged);
        block.Dwork(4).Data = height(history)-1;
        block.Dwork(5).Data = result.design;
        block.Dwork(6).Data = count;
        block.Dwork(7).Data = iteration;
        block.Dwork(8).Data = objective;
        block.Dwork(9).Data = constraint;
    catch exception
        block.Dwork(3).Data = -1;
        warning("radia:simulink:StreamFunctionTopologyFailed", ...
            "Stream-function topology optimization failed: %s", ...
            exception.message);
    end
end
block.Dwork(1).Data = trigger;
for index = 1:block.NumOutputPorts
    block.OutputPort(index).Data = block.Dwork(index+1).Data;
end
end

function validateRunner(runner)
if ~isa(runner,"radia.topopt.AdjointRunner") || ...
        ~isfield(runner.Metadata,"domain") || ...
        string(runner.Metadata.domain) ~= "stream-function"
    error("radia:simulink:StreamFunctionTopologyRunner", ...
        "Runner must be created by makeStreamFunctionAdjointRunner.");
end
end

function value = validateHistoryCapacity(value)
value = double(value);
if ~isscalar(value) || ~isfinite(value) || value < 2 || value ~= floor(value)
    error("radia:simulink:StreamFunctionTopologyHistory", ...
        "HistoryCapacity must be an integer greater than one.");
end
end
