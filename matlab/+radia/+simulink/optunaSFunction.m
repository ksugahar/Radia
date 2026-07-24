function optunaSFunction(block)
%OPTUNASFUNCTION Incremental Optuna study with Simulink-native telemetry.
setup(block);
end

function setup(block)
block.NumDialogPrms = 7;
block.DialogPrmsTunable = repmat({'Nontunable'}, 1, 7);
block.NumInputPorts = 1;
block.NumOutputPorts = 11;
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;
block.InputPort(1).Dimensions = 1;
block.InputPort(1).DirectFeedthrough = true;

nTrials = max(1, round(double(block.DialogPrm(2).Data)));
for index = 1:8
    block.OutputPort(index).Dimensions = 1;
    block.OutputPort(index).DatatypeID = 0;
end
for index = 9:10
    block.OutputPort(index).Dimensions = nTrials;
    block.OutputPort(index).DatatypeID = 0;
end
block.OutputPort(11).Dimensions = 1;
block.OutputPort(11).DatatypeID = 0;

block.SampleTimes = [block.DialogPrm(5).Data 0];
block.RegBlockMethod('PostPropagationSetup', @postSetup);
block.RegBlockMethod('Start', @start);
block.RegBlockMethod('Outputs', @outputs);
block.RegBlockMethod('Terminate', @terminate);
end

function postSetup(block)
nTrials = max(1, round(double(block.DialogPrm(2).Data)));
names = {'previous_trigger','best_value','best_trial','status', ...
    'completed_trials','last_value','elapsed_s','best_updated', ...
    'pareto_count','pareto_x','pareto_y','pareto_revision'};
block.NumDworks = numel(names);
for index = 1:numel(names)
    block.Dwork(index).Name = names{index};
    dimensions = 1;
    if index == 10 || index == 11
        dimensions = nTrials;
    end
    block.Dwork(index).Dimensions = dimensions;
    block.Dwork(index).DatatypeID = 0;
    block.Dwork(index).Complexity = 'Real';
    block.Dwork(index).UsedAsDiscState = true;
end
end

function start(block)
for index = 1:block.NumDworks
    block.Dwork(index).Data = zeros(block.Dwork(index).Dimensions, 1);
end
block.Dwork(2).Data = NaN;
block.Dwork(3).Data = NaN;
block.Dwork(6).Data = NaN;
block.Dwork(10).Data = NaN(block.Dwork(10).Dimensions, 1);
block.Dwork(11).Data = NaN(block.Dwork(11).Dimensions, 1);
radia.simulink.optunaRuntimeStore("remove", runtimeKey(block));
end

function outputs(block)
trigger = double(block.InputPort(1).Data);
block.Dwork(8).Data = 0;
key = runtimeKey(block);

if trigger > 0 && block.Dwork(1).Data <= 0
    initializeStudy(block, key);
end

runtime = radia.simulink.optunaRuntimeStore("get", key);
if ~isempty(runtime) && runtime.active
    try
        runtime = runOneTrial(block, runtime);
        radia.simulink.optunaRuntimeStore("set", key, runtime);
    catch exception
        runtime.active = false;
        radia.simulink.optunaRuntimeStore("set", key, runtime);
        block.Dwork(4).Data = -1;
        warning('radia:simulink:OptunaFailed', ...
            'Optuna block failed: %s', exception.message);
    end
end

block.Dwork(1).Data = trigger;
for index = 1:8
    block.OutputPort(index).Data = block.Dwork(index + 1).Data;
end
block.OutputPort(9).Data = block.Dwork(10).Data;
block.OutputPort(10).Data = block.Dwork(11).Data;
block.OutputPort(11).Data = block.Dwork(12).Data;
end

function initializeStudy(block, key)
objectiveParameter = block.DialogPrm(1).Data;
if isa(objectiveParameter, 'function_handle')
    objective = objectiveParameter;
else
    objectiveName = string(objectiveParameter);
    if ~isscalar(objectiveName) || strlength(objectiveName) == 0
        error('radia:simulink:OptunaObjective', ...
            'Objective must be a function handle or a nonempty function name.');
    end
    objective = str2func(char(objectiveName));
end
nTrials = max(1, round(double(block.DialogPrm(2).Data)));
directions = reshape(string(block.DialogPrm(3).Data), 1, []);
storage = string(block.DialogPrm(4).Data);
showExternalMonitor = logical(block.DialogPrm(6).Data);
samplerName = normalizeSamplerName(block.DialogPrm(7).Data);
sampler = makeSampler(samplerName, directions);

previous = radia.simulink.optunaRuntimeStore("get", key);
if ~isempty(previous) && ~isempty(previous.monitor) && isvalid(previous.monitor)
    delete(previous.monitor);
end

monitor = [];
progress = [];
if showExternalMonitor
    monitor = radia.optuna.LiveMonitor();
    progress = @monitor.update;
end
study = radia.optuna.createStudy(directions=directions, ...
    sampler=sampler, ...
    StoragePath=storage, AutoSave=strlength(storage) > 0, ...
    ProgressFcn=progress);
runtime = struct( ...
    "study", study, ...
    "objective", objective, ...
    "monitor", monitor, ...
    "started", tic, ...
    "n_trials", nTrials, ...
    "sampler", samplerName, ...
    "active", true, ...
    "best_trial", NaN, ...
    "pareto_signature", "");
radia.simulink.optunaRuntimeStore("set", key, runtime);

block.Dwork(2).Data = NaN;
block.Dwork(3).Data = NaN;
block.Dwork(4).Data = 2;
block.Dwork(5).Data = 0;
block.Dwork(6).Data = NaN;
block.Dwork(7).Data = 0;
block.Dwork(8).Data = 0;
block.Dwork(9).Data = 0;
block.Dwork(10).Data = NaN(block.Dwork(10).Dimensions, 1);
block.Dwork(11).Data = NaN(block.Dwork(11).Dimensions, 1);
block.Dwork(12).Data = 0;
end

function name = normalizeSamplerName(value)
choices = ["auto","random","tpe","cmaes","motpe","nsgaii"];
if isnumeric(value) && isscalar(value) && isfinite(value) && ...
        value == floor(value) && value >= 1 && value <= numel(choices)
    name = choices(value);
else
    name = lower(string(value));
end
end

function sampler = makeSampler(name, directions)
choices = ["auto","random","tpe","cmaes","motpe","nsgaii"];
if ~isscalar(name) || ~ismember(name, choices)
    error('radia:simulink:OptunaSampler', ...
        'Sampler must be auto, random, tpe, cmaes, motpe, or nsgaii.');
end
isMultiObjective = numel(directions) > 1;
if name == "auto"
    if isMultiObjective
        name = "nsgaii";
    else
        name = "cmaes";
    end
end
if isMultiObjective && ismember(name, ["tpe","cmaes"])
    error('radia:simulink:OptunaSampler', ...
        'Use random, motpe, or nsgaii for multiple objectives.');
end
if ~isMultiObjective && ismember(name, ["motpe","nsgaii"])
    error('radia:simulink:OptunaSampler', ...
        'Use random, tpe, or cmaes for a single objective.');
end
switch name
    case "random"
        sampler = radia.optuna.RandomSampler(0);
    case "tpe"
        sampler = radia.optuna.TPESampler(Seed=0, NStartupTrials=10);
    case "cmaes"
        sampler = radia.optuna.CmaEsSampler(Seed=0, NStartupTrials=1);
    case "motpe"
        sampler = radia.optuna.MOTPESampler(Seed=0, NStartupTrials=20);
    case "nsgaii"
        sampler = radia.optuna.NSGAIISampler(Seed=0, PopulationSize=24);
end
end

function runtime = runOneTrial(block, runtime)
trial = runtime.study.ask();
try
    value = runtime.objective(trial);
    if trial.State == "RUNNING"
        runtime.study.tell(trial, value);
    end
catch exception
    if trial.State == "RUNNING"
        runtime.study.fail(trial, exception.message);
    end
    rethrow(exception);
end

values = reshape(double(value), 1, []);
block.Dwork(6).Data = values(1);
complete = runtime.study.TrialTable.State == "COMPLETE";
completedCount = sum(complete);
block.Dwork(5).Data = completedCount;
block.Dwork(7).Data = toc(runtime.started);

[bestValue, bestTrial] = primaryBest(runtime.study);
block.Dwork(2).Data = bestValue;
block.Dwork(3).Data = bestTrial;
if ~isequaln(bestTrial, runtime.best_trial)
    block.Dwork(8).Data = 1;
    runtime.best_trial = bestTrial;
end

pareto = runtime.study.paretoFront();
nPoints = min(height(pareto), block.Dwork(10).Dimensions);
x = NaN(block.Dwork(10).Dimensions, 1);
y = NaN(block.Dwork(11).Dimensions, 1);
if nPoints > 0
    matrix = vertcat(pareto.Values{1:nPoints});
    x(1:nPoints) = matrix(:, 1);
    if size(matrix, 2) >= 2
        y(1:nPoints) = matrix(:, 2);
    end
end
block.Dwork(9).Data = nPoints;
block.Dwork(10).Data = x;
block.Dwork(11).Data = y;
signature = string(jsonencode(struct( ...
    "trial", pareto.TrialNumber, "x", x(1:nPoints), "y", y(1:nPoints))));
if signature ~= runtime.pareto_signature
    block.Dwork(12).Data = block.Dwork(12).Data + 1;
    runtime.pareto_signature = signature;
end

if completedCount >= runtime.n_trials
    runtime.active = false;
    block.Dwork(4).Data = 1;
else
    block.Dwork(4).Data = 2;
end
end

function [value, number] = primaryBest(study)
complete = study.TrialTable.State == "COMPLETE" & ...
    isfinite(study.TrialTable.Value);
rows = find(complete);
if isempty(rows)
    value = NaN;
    number = NaN;
    return;
end
values = study.TrialTable.Value(rows);
if study.Directions(1) == "minimize"
    [value, index] = min(values);
else
    [value, index] = max(values);
end
number = study.TrialTable.TrialNumber(rows(index));
end

function terminate(block)
key = runtimeKey(block);
runtime = radia.simulink.optunaRuntimeStore("get", key);
if ~isempty(runtime) && ~isempty(runtime.monitor) && isvalid(runtime.monitor)
    delete(runtime.monitor);
end
radia.simulink.optunaRuntimeStore("remove", key);
end

function key = runtimeKey(block)
key = string(sprintf('%.0f', block.BlockHandle));
end
