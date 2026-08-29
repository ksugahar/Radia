function optunaSFunction(block)
%OPTUNASFUNCTION Incremental Optuna study with Simulink-native telemetry.
setup(block);
end

function setup(block)
block.NumDialogPrms = 7;
block.DialogPrmsTunable = repmat({'Nontunable'}, 1, 7);
block.NumInputPorts = 6;
block.NumOutputPorts = 18;
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;
for index = 1:6
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
for index = 11:18
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
    'attempted_trials','previous_pause','previous_resume', ...
    'selected_trial','pruned_trials','current_trial', ...
    'checkpoint_revision','previous_apply'};
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
block.Dwork(19).Data = NaN;
block.Dwork(21).Data = NaN;
radia.simulink.optunaRuntimeStore("remove", runtimeKey(block));
end

function outputs(block)
trigger = double(block.InputPort(1).Data);
cancel = double(block.InputPort(2).Data);
pauseSignal = double(block.InputPort(3).Data);
resumeSignal = double(block.InputPort(4).Data);
selectedTrial = double(block.InputPort(5).Data);
applySignal = double(block.InputPort(6).Data);
block.Dwork(8).Data = 0;
key = runtimeKey(block);

if trigger > 0 && block.Dwork(1).Data <= 0
    initializeStudy(block, key);
end

runtime = radia.simulink.optunaRuntimeStore("get", key);
if ~isempty(runtime) && runtime.active && ...
        cancel > 0 && block.Dwork(13).Data <= 0
    runtime = cancelRuntime(runtime);
    radia.simulink.optunaRuntimeStore("set", key, runtime);
    block.Dwork(4).Data = 3;
elseif ~isempty(runtime) && runtime.active && ...
        pauseSignal > 0 && block.Dwork(17).Data <= 0
    runtime = pauseRuntime(runtime);
    radia.simulink.optunaRuntimeStore("set", key, runtime);
    block.Dwork(4).Data = 4;
elseif ~isempty(runtime) && runtime.paused && ...
        resumeSignal > 0 && block.Dwork(18).Data <= 0
    runtime = resumeRuntime(runtime);
    radia.simulink.optunaRuntimeStore("set", key, runtime);
    block.Dwork(4).Data = 2;
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

if ~isempty(runtime) && applySignal > 0 && block.Dwork(23).Data <= 0
    try
        runtime = applySelectedTrial(block, runtime, selectedTrial);
        radia.simulink.optunaRuntimeStore("set", key, runtime);
    catch exception
        block.Dwork(15).Data = failureCode(exception);
        warning('radia:simulink:OptunaApplyFailed', ...
            'Could not apply selected trial: %s', exception.message);
    end
end

block.Dwork(1).Data = trigger;
block.Dwork(13).Data = cancel;
block.Dwork(17).Data = pauseSignal;
block.Dwork(18).Data = resumeSignal;
block.Dwork(23).Data = applySignal;
for index = 1:8
    block.OutputPort(index).Data = block.Dwork(index + 1).Data;
end
block.OutputPort(9).Data = block.Dwork(10).Data;
block.OutputPort(10).Data = block.Dwork(11).Data;
block.OutputPort(11).Data = block.Dwork(12).Data;
block.OutputPort(12).Data = block.Dwork(14).Data;
block.OutputPort(13).Data = block.Dwork(16).Data;
block.OutputPort(14).Data = block.Dwork(15).Data;
block.OutputPort(15).Data = block.Dwork(19).Data;
block.OutputPort(16).Data = block.Dwork(20).Data;
block.OutputPort(17).Data = block.Dwork(21).Data;
block.OutputPort(18).Data = block.Dwork(22).Data;
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
rawSamplerSpec = block.DialogPrm(7).Data;
samplerSpec = normalizeSamplerSpec(rawSamplerSpec);
requestedSampler = samplerSpec.name;
[samplerSpec, resumedConfiguration] = resumeSamplerConfiguration( ...
    samplerSpec, storage, directions, objective);
[sampler, samplerName, samplerDecision] = ...
    makeSampler(samplerSpec, directions, nTrials, samplerSpec.seed);
if resumedConfiguration && requestedSampler == "auto"
    samplerDecision.requested = "auto";
    samplerDecision.reason = "resumed_saved_sampler";
end

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
session = [];
mode = "legacy-trial";
if ~isempty(samplerSpec.parameters)
    options = radia.optuna.optimoptions( ...
        MaxTrials=nTrials, Sampler=samplerName, ...
        Pruner=samplerSpec.pruner, Seed=samplerSpec.seed, ...
        StoragePath=storage, Directions=directions, Display="off");
    session = radia.optuna.OptimizationSession( ...
        objective, samplerSpec.parameters, options, ...
        ModelName=samplerSpec.model_name, ProgressFcn=progress);
    session.start();
    study = session.Study;
    mode = "parameter-session";
else
    settings = {"directions", directions, "sampler", sampler, ...
        "StoragePath", storage, "AutoSave", strlength(storage) > 0, ...
        "ProgressFcn", progress, ...
        "load_if_exists", strlength(storage) > 0};
    pruner = radia.optuna.internal.prunerFromName(samplerSpec.pruner);
    if ~isempty(pruner)
        settings = [settings, {"pruner", pruner}];
    end
    study = radia.optuna.createStudy(settings{:});
end
study.setUserAttr("auto_sampler_decision", samplerDecision);
runtime = struct( ...
    "study", study, ...
    "session", session, ...
    "mode", mode, ...
    "objective", objective, ...
    "monitor", monitor, ...
    "started", tic, ...
    "n_trials", nTrials, ...
    "sampler", samplerName, ...
    "sampler_decision", samplerDecision, ...
    "active", true, ...
    "paused", false, ...
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
block.Dwork(17).Data = 0;
block.Dwork(18).Data = 0;
block.Dwork(19).Data = NaN;
block.Dwork(20).Data = 0;
block.Dwork(21).Data = NaN;
block.Dwork(22).Data = 0;
block.Dwork(23).Data = 0;
end

function [spec, resumed] = resumeSamplerConfiguration( ...
        spec, storage, directions, objective)
resumed = false;
if strlength(storage) == 0 || ...
        (~isfile(storage) && ~isfile(storage + ".bak"))
    return
end

saved = radia.optuna.loadStudy(storage=storage);
validateSessionExperiment(storage, spec.parameters, objective);
if ~isequal(reshape(string(saved.Directions), 1, []), ...
        reshape(string(directions), 1, []))
    error('radia:simulink:OptunaResumeDirections', ...
        ['The saved study uses different objective directions. ' ...
        'Select a new Study MAT file for this experiment.']);
end
if ~isfield(saved.UserAttrs, "auto_sampler_decision")
    return
end

decision = saved.UserAttrs.auto_sampler_decision;
required = ["selected", "pruner"];
if ~isstruct(decision) || ~isscalar(decision) || ...
        any(~isfield(decision, required))
    error('radia:simulink:OptunaResumeMetadata', ...
        'The saved study has invalid sampler-decision metadata.');
end
savedSampler = string(decision.selected);
if spec.name == "auto"
    spec.name = savedSampler;
elseif spec.name ~= savedSampler
    error('radia:simulink:OptunaResumeSampler', ...
        ['The saved study uses sampler ''%s'', but the block requests ' ...
        '''%s''. Select a new Study MAT file to compare samplers.'], ...
        savedSampler, spec.name);
end
if string(decision.pruner) ~= spec.pruner
    error('radia:simulink:OptunaResumePruner', ...
        ['The saved study uses pruner ''%s'', but the block requests ' ...
        '''%s''. Select a new Study MAT file to compare pruners.'], ...
        string(decision.pruner), spec.pruner);
end

savedSeed = samplerStateSeed(saved, savedSampler);
if isempty(spec.seed)
    spec.seed = savedSeed;
elseif ~isempty(savedSeed) && double(spec.seed) ~= double(savedSeed)
    error('radia:simulink:OptunaResumeSeed', ...
        ['The saved study uses seed %.0f, but the block requests %.0f. ' ...
        'Select a new Study MAT file to compare seeds.'], ...
        double(savedSeed), double(spec.seed));
end
resumed = true;
end

function validateSessionExperiment(storage, parameters, objective)
sessionPath = storage + ".session.mat";
if ~isfile(sessionPath)
    return
end
loaded = builtin("load", sessionPath, "SessionData", "-mat");
if ~isfield(loaded, "SessionData")
    error('radia:simulink:OptunaResumeMetadata', ...
        "Session checkpoint '%s' has no SessionData.", sessionPath);
end
data = loaded.SessionData;
if ~isfield(data, "Parameters") || ...
        ~isequaln(data.Parameters, parameters)
    error('radia:simulink:OptunaResumeSearchSpace', ...
        ['The saved study uses a different OptimizationParameter ' ...
        'definition. Select a new Study MAT file for this experiment.']);
end
if isfield(data, "ObjectiveFcn") && ...
        string(func2str(data.ObjectiveFcn)) ~= string(func2str(objective))
    error('radia:simulink:OptunaResumeObjective', ...
        ['The saved study uses a different objective function. ' ...
        'Select a new Study MAT file for this experiment.']);
end
end

function seed = samplerStateSeed(study, samplerName)
seed = [];
states = study.SamplerStateTable;
rows = states.Sampler == samplerName;
if ~any(rows)
    % Composite samplers may still be in their independent-random startup
    % phase, so only the child sampler has emitted state. All child streams
    % are initialized from the resolved parent seed.
    rows = true(height(states), 1);
end
indices = find(rows);
candidates = NaN(numel(indices), 1);
for position = 1:numel(indices)
    state = states.State{indices(position)};
    if isstruct(state) && isscalar(state) && isfield(state, "seed") && ...
            isnumeric(state.seed) && isscalar(state.seed) && ...
            isfinite(double(state.seed))
        candidates(position) = double(state.seed);
    end
end
candidates = unique(candidates(isfinite(candidates)));
if isscalar(candidates)
    seed = candidates;
end
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
    "is_conditional", false, ...
    "seed", [], ...
    "pruner", "none", ...
    "parameters", [], ...
    "model_name", "");
if isstruct(value)
    if ~isscalar(value)
        error('radia:simulink:OptunaSamplerMetadata', ...
            'Auto sampler metadata must be one scalar struct.');
    end
    samplerValue = metadataField(value, {"Name","Sampler"}, "auto");
    if isnumeric(samplerValue) && isscalar(samplerValue) && ...
            isfinite(samplerValue) && samplerValue == floor(samplerValue) && ...
            samplerValue >= 1 && samplerValue <= numel(choices)
        spec.name = choices(samplerValue);
    else
        spec.name = lower(string(samplerValue));
    end
    spec.seed = metadataField(value, {"Seed"}, []);
    if ~(isempty(spec.seed) || (isnumeric(spec.seed) && ...
            isscalar(spec.seed) && isfinite(spec.seed) && ...
            spec.seed == floor(spec.seed) && spec.seed >= 0))
        error('radia:simulink:OptunaSamplerMetadata', ...
            'Seed must be empty or one nonnegative integer.');
    end
    prunerValue = metadataField(value, {"Pruner"}, "none");
    prunerChoices = ["none","median","hyperband","percentile", ...
        "patient","successivehalving","threshold"];
    if isnumeric(prunerValue) && isscalar(prunerValue) && ...
            isfinite(prunerValue) && prunerValue == floor(prunerValue) && ...
            prunerValue >= 1 && prunerValue <= numel(prunerChoices)
        spec.pruner = prunerChoices(prunerValue);
    else
        spec.pruner = lower(string(prunerValue));
    end
    if ~ismember(spec.pruner, ["none","median","hyperband", ...
            "percentile","patient","successivehalving","threshold"])
        error('radia:simulink:OptunaSamplerMetadata', ...
            'Unknown pruner ''%s''.', spec.pruner);
    end
    spec.parameters = metadataField(value, {"Parameters"}, []);
    spec.model_name = string(metadataField(value, {"ModelName"}, ""));
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
    if ~isempty(spec.parameters)
        if ~isa(spec.parameters, "radia.optuna.OptimizationParameter")
            error('radia:simulink:OptunaSamplerMetadata', ...
                'Parameters must contain OptimizationParameter objects.');
        end
        free = spec.parameters([spec.parameters.Free]);
        spec.dimensions = numel(free);
        spec.has_categorical = any([free.Type] == "categorical");
        spec.fixed_numeric = all([free.Type] ~= "categorical");
    end
elseif isnumeric(value) && isscalar(value) && isfinite(value) && ...
        value == floor(value) && value >= 1 && value <= numel(choices)
    spec.name = choices(value);
else
    spec.name = lower(string(value));
end
end

function [sampler, selectedName, decision] = ...
        makeSampler(spec, directions, nTrials, seed)
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
% The short-name mapping lives in one place, shared with
% radia.optuna.optimize. A new seed-less sampler draws private entropy; a
% resumed study supplies the exact seed recorded with its sampler state.
sampler = radia.optuna.internal.samplerFromName(selectedName, seed);
decision = struct( ...
    "schema", "radia.optuna.auto-sampler-lite.v2", ...
    "requested", spec.name, ...
    "selected", selectedName, ...
    "reason", reason, ...
    "n_trials", double(nTrials), ...
    "n_objectives", double(numel(directions)), ...
    "seed", seed, ...
    "pruner", spec.pruner, ...
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
    if any(strcmpi(fields{index}, cellstr(string(names))))
        value = metadata.(fields{index});
        found = true;
        return
    end
end
end

function runtime = runOneTrial(block, runtime)
failed = false;
if runtime.mode == "parameter-session"
    frozen = runtime.session.runNext();
    if isempty(frozen)
        value = NaN;
    else
        value = frozen.Values;
        failed = frozen.State == "FAIL";
        if failed
            block.Dwork(15).Data = failureCode(MException( ...
                "radia:simulink:OptunaObjective", ...
                "%s", frozen.ErrorMessage));
        end
    end
else
    trial = runtime.study.ask();
    try
        value = runtime.objective(trial);
        if trial.State == "RUNNING"
            runtime.study.tell(trial, value);
        end
    catch exception
        if trial.State == "RUNNING"
            if exception.identifier == "radia:optuna:TrialPruned"
                runtime.study.tell(trial, State="PRUNED");
            else
                runtime.study.fail(trial, exception.message);
            end
        end
        failed = trial.State == "FAIL";
        value = NaN;
        if failed
            block.Dwork(15).Data = failureCode(exception);
        end
    end
end
if failed
    block.Dwork(14).Data = block.Dwork(14).Data + 1;
end

values = reshape(double(value), 1, []);
if isempty(values)
    block.Dwork(6).Data = NaN;
else
    block.Dwork(6).Data = values(1);
end
complete = runtime.study.TrialTable.State == "COMPLETE";
completedCount = sum(complete);
attemptedCount = height(runtime.study.TrialTable);
prunedCount = sum(runtime.study.TrialTable.State == "PRUNED");
block.Dwork(5).Data = completedCount;
block.Dwork(16).Data = attemptedCount;
block.Dwork(20).Data = prunedCount;
if attemptedCount > 0
    block.Dwork(21).Data = ...
        runtime.study.TrialTable.TrialNumber(end);
end
block.Dwork(7).Data = toc(runtime.started);
if runtime.mode == "parameter-session"
    block.Dwork(22).Data = runtime.session.CheckpointRevision;
end

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

if runtime.mode == "parameter-session" && ...
        runtime.session.State == "completed"
    runtime.active = false;
    block.Dwork(4).Data = 1;
elseif attemptedCount >= runtime.n_trials
    runtime.active = false;
    block.Dwork(4).Data = 1;
else
    block.Dwork(4).Data = 2;
end
if failed && attemptedCount < runtime.n_trials
    block.Dwork(4).Data = 2;
end
end

function runtime = cancelRuntime(runtime)
if runtime.mode == "parameter-session" && ...
        ismember(runtime.session.State, ["configured","running","paused"])
    runtime.session.cancel();
end
runtime.active = false;
runtime.paused = false;
end

function runtime = pauseRuntime(runtime)
if runtime.mode == "parameter-session"
    runtime.session.pause();
end
runtime.active = false;
runtime.paused = true;
end

function runtime = resumeRuntime(runtime)
if runtime.mode == "parameter-session"
    runtime.session.resume();
end
runtime.active = true;
runtime.paused = false;
end

function runtime = applySelectedTrial(block, runtime, selectedTrial)
if runtime.mode ~= "parameter-session"
    error("radia:simulink:OptunaSelection", ...
        "Trial selection requires a nonempty Parameters specification.");
end
if ~isfinite(selectedTrial) || selectedTrial < 0 || ...
        selectedTrial ~= floor(selectedTrial)
    error("radia:simulink:OptunaSelection", ...
        "Selected trial must be one nonnegative integer.");
end
runtime.session.selectTrial(selectedTrial);
runtime.session.applySelectedToModel();
block.Dwork(19).Data = selectedTrial;
block.Dwork(22).Data = runtime.session.CheckpointRevision;
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
if ~isempty(runtime) && isfield(runtime, "session") && ...
        ~isempty(runtime.session) && runtime.session.State == "running"
    runtime.session.pause();
end
if ~isempty(runtime) && ~isempty(runtime.monitor) && isvalid(runtime.monitor)
    delete(runtime.monitor);
end
radia.simulink.optunaRuntimeStore("remove", key);
end

function key = runtimeKey(block)
key = string(sprintf('%.0f', block.BlockHandle));
end
