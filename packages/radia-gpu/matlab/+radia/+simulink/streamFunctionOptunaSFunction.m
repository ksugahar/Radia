function streamFunctionOptunaSFunction(block)
%STREAMFUNCTIONOPTUNASFUNCTION Run one explicit Stream Function Optuna study.
setup(block);
end

function setup(block)
block.NumDialogPrms = 4;
block.DialogPrmsTunable = repmat({'Nontunable'}, 1, 4);
block.NumInputPorts = 1;
block.NumOutputPorts = 7;
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;
block.InputPort(1).Dimensions = 1;
block.InputPort(1).DatatypeID = 8;
block.InputPort(1).Complexity = "Real";
block.InputPort(1).DirectFeedthrough = false;
block.OutputPort(1).DatatypeID = 6;
for index = 1:block.NumOutputPorts
    block.OutputPort(index).Dimensions = 1;
    block.OutputPort(index).Complexity = "Real";
    if index > 1
        block.OutputPort(index).DatatypeID = 0;
    end
end
block.SampleTimes = [validateSampleTime(block.DialogPrm(4).Data), 0];
block.SimStateCompliance = "DefaultSimState";
block.RegBlockMethod("PostPropagationSetup", @postSetup);
block.RegBlockMethod("InitializeConditions", @initializeConditions);
block.RegBlockMethod("Outputs", @outputs);
block.RegBlockMethod("Update", @update);
end

function postSetup(block)
names = {'previous_trigger', 'status', 'best_value', 'best_trial', ...
    'completed', 'failed', 'pareto_count', 'elapsed_s'};
types = [8, 6, 0, 0, 0, 0, 0, 0];
block.NumDworks = numel(names);
for index = 1:numel(names)
    block.Dwork(index).Name = names{index};
    block.Dwork(index).Dimensions = 1;
    block.Dwork(index).DatatypeID = types(index);
    block.Dwork(index).Complexity = "Real";
    block.Dwork(index).UsedAsDiscState = true;
end
end

function initializeConditions(block)
block.Dwork(1).Data = false;
block.Dwork(2).Data = int32(0);
block.Dwork(3).Data = NaN;
block.Dwork(4).Data = NaN;
for index = 5:8
    block.Dwork(index).Data = 0;
end
end

function outputs(block)
for index = 1:block.NumOutputPorts
    block.OutputPort(index).Data = block.Dwork(index + 1).Data;
end
end

function update(block)
trigger = logical(block.InputPort(1).Data);
previous = logical(block.Dwork(1).Data);
block.Dwork(1).Data = trigger;
if ~trigger || previous
    return
end

block.Dwork(2).Data = int32(1);
timer = tic;
try
    runnerVariable = string(block.DialogPrm(1).Data);
    studyVariable = string(block.DialogPrm(2).Data);
    nTrials = validateTrialCount(block.DialogPrm(3).Data);
    modelName = string(get_param(bdroot(block.BlockHandle), "Name"));
    [runner, study] = radia.simulink.resolveStreamFunctionOptunaObjects( ...
        modelName, runnerVariable, studyVariable);
    runner.optimize(study, nTrials, ContinueOnError=true);

    trials = study.TrialTable;
    complete = trials.State == "COMPLETE";
    failed = trials.State == "FAIL";
    block.Dwork(5).Data = sum(complete);
    block.Dwork(6).Data = sum(failed);
    front = study.paretoFront();
    block.Dwork(7).Data = height(front);
    if isscalar(study.Directions) && any(complete)
        best = study.bestSolution();
        block.Dwork(3).Data = best.value;
        block.Dwork(4).Data = best.trial_number;
    else
        block.Dwork(3).Data = NaN;
        block.Dwork(4).Data = NaN;
    end
    block.Dwork(8).Data = toc(timer);
    block.Dwork(2).Data = int32(2 * any(complete) - ~any(complete));
    summary = struct( ...
        "schema", "radia.simulink.streamfunction-optuna.v1", ...
        "status", double(block.Dwork(2).Data), ...
        "completed_trials", double(block.Dwork(5).Data), ...
        "failed_trials", double(block.Dwork(6).Data), ...
        "pareto_count", double(block.Dwork(7).Data), ...
        "best_value", double(block.Dwork(3).Data), ...
        "best_trial", double(block.Dwork(4).Data), ...
        "elapsed_s", double(block.Dwork(8).Data), ...
        "study_name", study.Name, ...
        "storage_path", study.StoragePath);
    set_param(block.BlockHandle, "UserData", summary, ...
        "UserDataPersistent", "off");
catch exception
    block.Dwork(2).Data = int32(-1);
    block.Dwork(8).Data = toc(timer);
    diagnostic = struct( ...
        "schema", "radia.simulink.streamfunction-optuna.v1", ...
        "status", -1, ...
        "error", string(exception.message), ...
        "identifier", string(exception.identifier), ...
        "elapsed_s", double(block.Dwork(8).Data));
    set_param(block.BlockHandle, "UserData", diagnostic, ...
        "UserDataPersistent", "off");
    warning("radia:simulink:StreamFunctionOptunaFailed", ...
        "Stream Function Optuna study failed: %s", exception.message);
end
end

function value = validateTrialCount(value)
value = double(value);
if ~isscalar(value) || ~isfinite(value) || value < 1 || value ~= floor(value)
    error("radia:simulink:StreamFunctionOptunaTrials", ...
        "Number of trials must be a positive integer.");
end
end

function value = validateSampleTime(value)
value = double(value);
if ~isscalar(value) || ~isfinite(value) || value <= 0
    error("radia:simulink:StreamFunctionOptunaSampleTime", ...
        "Sample time must be finite and positive.");
end
end
