function optunaSFunction(block)
%OPTUNASFUNCTION Incremental Optuna study with Simulink-native telemetry.
setup(block);
end

function setup(block)
block.NumDialogPrms = 7;
block.DialogPrmsTunable = repmat({'Nontunable'}, 1, 7);
block.NumInputPorts = 2;
block.NumOutputPorts = 14;
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;
for index = 1:2
    block.InputPort(index).Dimensions = 1;
    block.InputPort(index).DirectFeedthrough = true;
end

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
for index = 12:14
    block.OutputPort(index).Dimensions = 1;
    block.OutputPort(index).DatatypeID = 0;
end

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
    'pareto_count','pareto_x','pareto_y','pareto_revision', ...
    'previous_cancel','failed_trials','last_failure_code', ...
    'attempted_trials'};
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
cancel = double(block.InputPort(2).Data);
block.Dwork(8).Data = 0;
key = runtimeKey(block);

if trigger > 0 && block.Dwork(1).Data <= 0
    initializeStudy(block, key);
end

runtime = radia.simulink.optunaRuntimeStore("get", key);
if ~isempty(runtime) && runtime.active && ...
        cancel > 0 && block.Dwork(13).Data <= 0
    runtime.active = false;
    radia.simulink.optunaRuntimeStore("set", key, runtime);
    block.Dwork(4).Data = 3;
elseif ~isempty(runtime) && runtime.active
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
block.Dwork(13).Data = cancel;
for index = 1:8
    block.OutputPort(index).Data = block.Dwork(index + 1).Data;
end
block.OutputPort(9).Data = block.Dwork(10).Data;
block.OutputPort(10).Data = block.Dwork(11).Data;
block.OutputPort(11).Data = block.Dwork(12).Data;
block.OutputPort(12).Data = block.Dwork(14).Data;
block.OutputPort(13).Data = block.Dwork(16).Data;
block.OutputPort(14).Data = block.Dwork(15).Data;
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
samplerSpec = normalizeSamplerSpec(block.DialogPrm(7).Data);
[sampler, samplerName, samplerDecision] = ...
    makeSampler(samplerSpec, directions, nTrials);

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
study.setUserAttr("auto_sampler_decision", samplerDecision);
runtime = struct( ...
    "study", study, ...
    "objective", objective, ...
    "monitor", monitor, ...
    "started", tic, ...
    "n_trials", nTrials, ...
    "sampler", samplerName, ...
    "sampler_decision", samplerDecision, ...
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
block.Dwork(14).Data = 0;
block.Dwork(15).Data = 0;
block.Dwork(16).Data = 0;
end

function spec = normalizeSamplerSpec(value)
choices = ["auto","random","tpe","cmaes","motpe","nsgaii", ...
    "gp","nsgaiii","bruteforce","qmc"];
spec = struct( ...
    "name", "auto", ...
    "fixed_numeric", false, ...
    "dimensions", NaN, ...
    "has_constraints", false, ...
    "constraints_declared", false, ...
    "has_categorical", false, ...
    "is_conditional", false);
if isstruct(value)
    if ~isscalar(value)
        error('radia:simulink:OptunaSamplerMetadata', ...
            'Auto sampler metadata must be one scalar struct.');
    end
    spec.name = lower(string(metadataField(value, ...
        {"Name","Sampler"}, "auto")));
    spec.fixed_numeric = metadataLogical(value, {"FixedNumeric"}, false);
    spec.has_categorical = metadataLogical(value, ...
        {"HasCategorical"}, false);
    spec.is_conditional = metadataLogical(value, {"IsConditional"}, false);
    [constraintValue, constraintsDeclared] = metadataField(value, ...
        {"HasConstraints"}, false);
    spec.constraints_declared = constraintsDeclared;
    spec.has_constraints = logicalScalar(constraintValue, "HasConstraints");
    [dimensionValue, dimensionDeclared] = metadataField(value, ...
        {"Dimensions","Dimension"}, NaN);
    if dimensionDeclared
        if ~(isnumeric(dimensionValue) && isscalar(dimensionValue) && ...
                isfinite(dimensionValue) && dimensionValue >= 1 && ...
                dimensionValue == floor(dimensionValue))
            error('radia:simulink:OptunaSamplerMetadata', ...
                'Dimensions must be a positive integer.');
        end
        spec.dimensions = double(dimensionValue);
    end
    if spec.fixed_numeric && ...
            (spec.has_categorical || spec.is_conditional)
        error('radia:simulink:OptunaSamplerMetadata', ...
            ['FixedNumeric cannot be combined with HasCategorical or ' ...
            'IsConditional.']);
    end
elseif isnumeric(value) && isscalar(value) && isfinite(value) && ...
        value == floor(value) && value >= 1 && value <= numel(choices)
    spec.name = choices(value);
else
    spec.name = lower(string(value));
end
end

function [sampler, selectedName, decision] = makeSampler(spec, directions, nTrials)
choices = ["auto","random","tpe","cmaes","motpe","nsgaii", ...
    "gp","nsgaiii","bruteforce","qmc"];
if ~isscalar(spec.name) || ~ismember(spec.name, choices)
    error('radia:simulink:OptunaSampler', ...
        ['Sampler must be auto, random, tpe, cmaes, gp, motpe, ' ...
        'nsgaii, nsgaiii, bruteforce, or qmc.']);
end
isMultiObjective = numel(directions) > 1;
selectedName = spec.name;
reason = "explicit_sampler";
if selectedName == "auto"
    [selectedName, reason] = chooseAutoSampler( ...
        spec, numel(directions), nTrials);
end
if spec.has_constraints && selectedName == "cmaes"
    error('radia:simulink:OptunaSampler', ...
        ['This MATLAB sampler does not implement constrained ranking; ' ...
        'use GP, TPE, MOTPE, NSGA-II, or NSGA-III.']);
end
if isMultiObjective && ismember(selectedName, ["tpe","cmaes"])
    error('radia:simulink:OptunaSampler', ...
        ['Use random, gp, motpe, nsgaii, nsgaiii, bruteforce, ' ...
        'or qmc for multiple objectives.']);
end
if ~isMultiObjective && ismember(selectedName, ...
        ["motpe","nsgaii","nsgaiii"])
    error('radia:simulink:OptunaSampler', ...
        ['Use random, tpe, cmaes, gp, bruteforce, or qmc for ' ...
        'a single objective.']);
end
switch selectedName
    case "random"
        sampler = radia.optuna.RandomSampler(0);
    case "tpe"
        sampler = radia.optuna.TPESampler(Seed=0, NStartupTrials=10);
    case "cmaes"
        sampler = radia.optuna.CmaEsSampler(Seed=0, NStartupTrials=1);
    case "gp"
        sampler = radia.optuna.GPSampler(Seed=0,NStartupTrials=10, ...
            DeterministicObjective=true);
    case "motpe"
        sampler = radia.optuna.MOTPESampler(Seed=0, NStartupTrials=20);
    case "nsgaii"
        sampler = radia.optuna.NSGAIISampler(Seed=0, PopulationSize=24);
    case "nsgaiii"
        sampler = radia.optuna.NSGAIIISampler(Seed=0,PopulationSize=24);
    case "bruteforce"
        sampler = radia.optuna.BruteForceSampler(Seed=0);
    case "qmc"
        sampler = radia.optuna.QMCSampler( ...
            QMCType="sobol",Scramble=true,Seed=0);
end
decision = struct( ...
    "schema", "radia.optuna.auto-sampler-lite.v2", ...
    "requested", spec.name, ...
    "selected", selectedName, ...
    "reason", reason, ...
    "n_trials", double(nTrials), ...
    "n_objectives", double(numel(directions)), ...
    "fixed_numeric", spec.fixed_numeric, ...
    "dimensions", spec.dimensions, ...
    "has_constraints", spec.has_constraints, ...
    "constraints_declared", spec.constraints_declared, ...
    "has_categorical", spec.has_categorical, ...
    "is_conditional", spec.is_conditional);
end

function [name, reason] = chooseAutoSampler(spec, nObjectives, nTrials)
[name,reason]=radia.optuna.internal.AutoSamplerPolicy. ...
    choose(spec,nObjectives,nTrials);
end

function value = metadataLogical(metadata, names, fallback)
[candidate, found] = metadataField(metadata, names, fallback);
if ~found
    value = fallback;
else
    value = logicalScalar(candidate, names{1});
end
end

function value = logicalScalar(candidate, name)
if ~((islogical(candidate) || isnumeric(candidate)) && ...
        isscalar(candidate) && isfinite(double(candidate)) && ...
        ismember(double(candidate), [0, 1]))
    error('radia:simulink:OptunaSamplerMetadata', ...
        '%s must be one logical scalar.', name);
end
value = logical(candidate);
end

function [value, found] = metadataField(metadata, names, fallback)
fields = fieldnames(metadata);
found = false;
value = fallback;
for index = 1:numel(fields)
    if any(strcmpi(fields{index}, names))
        value = metadata.(fields{index});
        found = true;
        return
    end
end
end

function runtime = runOneTrial(block, runtime)
trial = runtime.study.ask();
failed = false;
try
    value = runtime.objective(trial);
    if trial.State == "RUNNING"
        runtime.study.tell(trial, value);
    end
catch exception
    if trial.State == "RUNNING"
        runtime.study.fail(trial, exception.message);
    end
    failed = true;
    value = NaN;
    block.Dwork(14).Data = block.Dwork(14).Data + 1;
    block.Dwork(15).Data = failureCode(exception);
end

values = reshape(double(value), 1, []);
block.Dwork(6).Data = values(1);
complete = runtime.study.TrialTable.State == "COMPLETE";
completedCount = sum(complete);
attemptedCount = height(runtime.study.TrialTable);
block.Dwork(5).Data = completedCount;
block.Dwork(16).Data = attemptedCount;
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

if attemptedCount >= runtime.n_trials
    runtime.active = false;
    block.Dwork(4).Data = 1;
else
    block.Dwork(4).Data = 2;
end
if failed && attemptedCount < runtime.n_trials
    block.Dwork(4).Data = 2;
end
end

function code = failureCode(exception)
text = lower(string(exception.identifier) + " " + string(exception.message));
if contains(text,"timeout") || contains(text,"timed out")
    code = 1;
elseif contains(text,"license") || contains(text,"resource") || ...
        contains(text,"checkout")
    code = 2;
elseif contains(text,"mesh") || contains(text,"element quality") || ...
        contains(text,"jacobian")
    code = 3;
elseif contains(text,"converge") || contains(text,"solver") || ...
        contains(text,"singular")
    code = 4;
elseif contains(text,"observable") || contains(text,"nan") || ...
        contains(text,"inf")
    code = 5;
elseif contains(text,"config") || contains(text,"parameter") || ...
        contains(text,"model")
    code = 6;
else
    code = 99;
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
