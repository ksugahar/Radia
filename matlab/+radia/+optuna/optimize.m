function [x, fval, exitflag, output] = optimize(fun, parameters, options)
%OPTIMIZE Run an Optuna study with a Global Optimization Toolbox signature.
%   [X, FVAL, EXITFLAG, OUTPUT] = RADIA.OPTUNA.OPTIMIZE(FUN, PARAMETERS)
%   searches PARAMETERS for the parameter set minimizing FUN, and returns the
%   four values ga, particleswarm and simulannealbnd return.
%
%   FUN is called as FVAL = FUN(X), where X is a struct whose fields are the
%   parameter names.  This is the toolbox's contract, not Optuna's: the
%   objective receives values, never a Trial.  That matters for UseParallel --
%   suggesting is a serial operation on the study, so only an objective that
%   takes values can be evaluated on workers without making the sampler's
%   sequence depend on scheduling.
%
%   PARAMETERS is an array of radia.optuna.OptimizationParameter, or a struct
%   of bounds accepted by OptimizationParameter.fromStruct.
%
%   OPTIONS comes from radia.optuna.optimoptions.
%
%   EXITFLAG uses ga's codes so the meaning carries over:
%      1  FunctionTolerance satisfied over MaxStallTrials
%      0  MaxTrials reached
%     -1  stopped by an OutputFcn or PlotFcn
%     -2  no feasible point found (no trial completed)
%     -5  MaxTime exceeded
%
%   OUTPUT reports funccount, trials, message and rngstate as the toolbox
%   solvers do, and adds what only a study can report: the study handle, the
%   best trial number, and the pruned and failed counts.
%
%   OutputFcn and PlotFcn are called as STOP = FCN(X, OPTIMVALUES, STATE) with
%   STATE "init", "iter" and "done", and a true STOP ends the run -- the same
%   protocol the toolbox uses.
%
%   Example:
%       p = radia.optuna.OptimizationParameter.fromStruct( ...
%           struct("Kp", [0 10], "Ki", [1e-3 1]));
%       opts = radia.optuna.optimoptions(MaxTrials=50, Display="iter");
%       [x, fval, flag] = radia.optuna.optimize(@(v) cost(v), p, opts);
%
%   See also radia.optuna.optimoptions, radia.optuna.OptimizationParameter,
%   radia.optuna.getParameterFromModel.

arguments
    fun (1,1) function_handle
    parameters
    options (1,1) radia.optuna.OptimizeOptions = radia.optuna.OptimizeOptions()
end

if isstruct(parameters)
    parameters = radia.optuna.OptimizationParameter.fromStruct(parameters);
end
if ~isa(parameters, "radia.optuna.OptimizationParameter")
    error("radia:optuna:Optimize", ...
        "PARAMETERS must be radia.optuna.OptimizationParameter objects " + ...
        "or a struct of bounds.");
end
parameters = reshape(parameters, 1, []);
searchable = parameters([parameters.Free]);
if isempty(searchable)
    error("radia:optuna:Optimize", ...
        "Every parameter has Free=false, so there is nothing to search.");
end
arrayfun(@(p) p.mustBeSearchable(), searchable);

fixed = parameters(~[parameters.Free]);
for index = 1:numel(fixed)
    if isempty(fixed(index).Value)
        error("radia:optuna:Optimize", ...
            "Parameter '%s' is fixed (Free=false) but has no Value.", ...
            fixed(index).Name);
    end
end

study = buildStudy(options, searchable);
distributions = arrayfun(@(p) p.distribution(), searchable, ...
    UniformOutput=false);
names = [searchable.Name];

started = tic;
best = struct("x", struct(), "fval", [], "trialNumber", []);
stalled = 0;
funccount = 0;
prunedCount = 0;
failedCount = 0;
exitflag = 0;

stopRequested = callHooks(options, best, emptyOptimValues(options, started), "init");
if stopRequested
    exitflag = -1;
end
printHeader(options);

trial = 0;
while ~stopRequested && trial < options.MaxTrials
    if toc(started) >= options.MaxTime
        exitflag = -5;
        break
    end
    batch = batchSize(options, options.MaxTrials - trial);
    [trials, values] = askBatch(study, names, distributions, fixed, batch);
    scores = evaluateBatch(fun, values, options);

    for index = 1:numel(trials)
        trial = trial + 1;
        funccount = funccount + 1;
        outcome = scores(index);
        if ~isempty(outcome.constraints)
            study.recordConstraints(trials(index), outcome.constraints);
        end
        if outcome.failed
            failedCount = failedCount + 1;
            study.fail(trials(index), outcome.errorMessage);
            fvalThisTrial = NaN;
        elseif outcome.pruned
            prunedCount = prunedCount + 1;
            study.tell(trials(index), State="PRUNED");
            fvalThisTrial = NaN;
        else
            study.tell(trials(index), outcome.value);
            fvalThisTrial = outcome.value;
        end

        improved = false;
        if ~isnan(fvalThisTrial)
            if isempty(best.fval) || isImprovement(fvalThisTrial, best.fval, ...
                    options.Directions(1))
                improvement = Inf;
                if ~isempty(best.fval)
                    improvement = abs(best.fval - fvalThisTrial);
                end
                best.fval = fvalThisTrial;
                best.x = values(index);
                best.trialNumber = trials(index).Number;
                improved = improvement > options.FunctionTolerance;
            end
        end
        if improved
            stalled = 0;
        else
            stalled = stalled + 1;
        end

        optimValues = struct( ...
            "trial", trial, ...
            "fval", fvalThisTrial, ...
            "bestfval", pick(best.fval, NaN), ...
            "bestx", best.x, ...
            "funccount", funccount, ...
            "elapsedtime", toc(started), ...
            "stalled", stalled, ...
            "state", outcome.state, ...
            "sampler", options.Sampler, ...
            "pruned", prunedCount, ...
            "failed", failedCount);
        printIteration(options, optimValues);
        if callHooks(options, best, optimValues, "iter")
            stopRequested = true;
            exitflag = -1;
            break
        end
        if isfinite(options.MaxStallTrials) && stalled >= options.MaxStallTrials
            exitflag = 1;
            stopRequested = true;
            break
        end
    end
end

if isempty(best.fval)
    exitflag = -2;
end

x = best.x;
fval = pick(best.fval, NaN);
for index = 1:numel(fixed)
    x.(matlab.lang.makeValidName(fixed(index).Name)) = fixed(index).Value;
end

output = struct( ...
    "funccount", funccount, ...
    "trials", trial, ...
    "message", exitMessage(exitflag, trial, options), ...
    "rngstate", struct("seed", options.Seed, "sampler", options.Sampler), ...
    "study", study, ...
    "bestTrialNumber", pick(best.trialNumber, NaN), ...
    "prunedcount", prunedCount, ...
    "failedcount", failedCount, ...
    "elapsedtime", toc(started));

finalValues = struct( ...
    "trial", trial, "fval", fval, "bestfval", fval, "bestx", x, ...
    "funccount", funccount, "elapsedtime", toc(started), ...
    "stalled", stalled, "state", "DONE", "sampler", options.Sampler, ...
    "pruned", prunedCount, "failed", failedCount);
callHooks(options, best, finalValues, "done");
printFooter(options, output, exitflag);
end

function study = buildStudy(options, parameters)
name = options.Sampler;
if name == "auto"
    % Resolving "auto" needs the study's shape, which only the caller knows.
    spec = struct( ...
        "fixed_numeric", all([parameters.Type] ~= "categorical"), ...
        "dimensions", numel(parameters), ...
        "has_constraints", ~isempty(options.ConstraintFcn), ...
        "constraints_declared", ~isempty(options.ConstraintFcn), ...
        "has_categorical", any([parameters.Type] == "categorical"), ...
        "is_conditional", false);
    name = radia.optuna.internal.AutoSamplerPolicy.choose( ...
        spec, numel(options.Directions), options.MaxTrials);
end
sampler = radia.optuna.internal.samplerFromName(name, options.Seed);
settings = {"Sampler", sampler, "Directions", options.Directions};
if strlength(options.StudyName) > 0
    settings = [settings, {"Name", options.StudyName}];
end
if options.StoragePath == ""
    settings = [settings, {"AutoSave", false}];
else
    if ~options.Resume && (isfile(options.StoragePath) || ...
            isfile(options.StoragePath + ".bak"))
        error("radia:optuna:Resume", ...
            "Storage '%s' already exists and Resume=false.", ...
            options.StoragePath);
    end
    settings = [settings, {"StoragePath", options.StoragePath}];
end
if options.Pruner ~= "none"
    settings = [settings, {"Pruner", ...
        radia.optuna.internal.prunerFromName(options.Pruner)}];
end
study = radia.optuna.Study(settings{:});
end

function [trials, values] = askBatch(study, names, distributions, fixed, count)
trials = radia.optuna.Trial.empty(1, 0);
values = struct([]);
for index = 1:count
    trial = study.ask();
    entry = struct();
    for k = 1:numel(names)
        entry.(matlab.lang.makeValidName(names(k))) = ...
            radia.optuna.internal.suggestFromDistribution( ...
                trial, names(k), distributions{k});
    end
    for k = 1:numel(fixed)
        entry.(matlab.lang.makeValidName(fixed(k).Name)) = fixed(k).Value;
    end
    trials(index) = trial;
    if isempty(values)
        values = entry;
    else
        values(index) = entry;
    end
end
end

function scores = evaluateBatch(fun, values, options)
count = numel(values);
scores = repmat(struct("value", NaN, "failed", false, "pruned", false, ...
    "state", "COMPLETE", "constraints", double.empty(1,0), ...
    "errorMessage", ""), 1, count);
if options.UseParallel && count > 1
    parfor index = 1:count
        scores(index) = evaluateOne(fun, values(index), options);
    end
else
    for index = 1:count
        scores(index) = evaluateOne(fun, values(index), options);
    end
end
end

function score = evaluateOne(fun, value, options)
score = struct("value", NaN, "failed", false, "pruned", false, ...
    "state", "COMPLETE", "constraints", double.empty(1,0), ...
    "errorMessage", "");
try
    if ~isempty(options.ConstraintFcn)
        score.constraints = reshape(double(options.ConstraintFcn(value)), 1, []);
        if any(~isfinite(score.constraints))
            score.failed = true;
            score.state = "FAIL";
            score.errorMessage = "ConstraintFcn returned a non-finite value.";
            return
        end
    end
    score.value = double(fun(value));
    if ~isscalar(score.value) || ~isfinite(score.value)
        score.failed = true;
        score.state = "FAIL";
        score.errorMessage = ...
            "Objective function must return one finite scalar.";
    end
catch problem
    if problem.identifier == "radia:optuna:TrialPruned"
        score.pruned = true;
        score.state = "PRUNED";
    else
        if ~options.CatchObjectiveErrors
            rethrow(problem)
        end
        score.failed = true;
        score.state = "FAIL";
        score.errorMessage = string(problem.message);
    end
end
end

function count = batchSize(options, remaining)
if options.UseParallel
    pool = gcp("nocreate");
    if isempty(pool)
        error("radia:optuna:UseParallel", ...
            "UseParallel=true but no parallel pool is open. " + ...
            "Open one with parpool, or set UseParallel=false.");
    end
    count = pool.NumWorkers;
    if options.BatchSize > 0
        count = min(count, options.BatchSize);
    end
    count = min(count, remaining);
else
    count = 1;
end
count = max(1, count);
end

function tf = isImprovement(candidate, incumbent, direction)
if direction == "maximize"
    tf = candidate > incumbent;
else
    tf = candidate < incumbent;
end
end

function stop = callHooks(options, best, optimValues, state)
stop = false;
hooks = [options.OutputFcn, options.PlotFcn];
for index = 1:numel(hooks)
    result = hooks{index}(best.x, optimValues, state);
    if ~isempty(result) && islogical(result) && any(result)
        stop = true;
    end
end
end

function values = emptyOptimValues(options, started)
values = struct("trial", 0, "fval", NaN, "bestfval", NaN, ...
    "bestx", struct(), "funccount", 0, "elapsedtime", toc(started), ...
    "stalled", 0, "state", "INIT", "sampler", options.Sampler, ...
    "pruned", 0, "failed", 0);
end

function value = pick(candidate, fallback)
if isempty(candidate)
    value = fallback;
else
    value = candidate;
end
end

function printHeader(options)
if options.Display ~= "iter"
    return
end
fprintf("\n%8s %14s %14s %10s %10s\n", ...
    "Trial", "f(x)", "Best f(x)", "State", "Time (s)");
end

function printIteration(options, optimValues)
if options.Display ~= "iter"
    return
end
fprintf("%8d %14.6g %14.6g %10s %10.2f\n", optimValues.trial, ...
    optimValues.fval, optimValues.bestfval, optimValues.state, ...
    optimValues.elapsedtime);
end

function printFooter(options, output, exitflag)
if options.isQuiet()
    return
end
fprintf("\n%s\n", output.message);
fprintf("Trials: %d (pruned %d, failed %d), exitflag %d, %.2f s\n\n", ...
    output.trials, output.prunedcount, output.failedcount, exitflag, ...
    output.elapsedtime);
end

function message = exitMessage(exitflag, trial, options)
switch exitflag
    case 1
        message = sprintf("Optimization stopped: the best objective value " + ...
            "improved by less than FunctionTolerance (%g) over the last " + ...
            "%d trials.", options.FunctionTolerance, options.MaxStallTrials);
    case 0
        message = sprintf("Optimization stopped: MaxTrials (%d) reached.", ...
            options.MaxTrials);
    case -1
        message = "Optimization stopped by an output or plot function.";
    case -2
        message = sprintf("No feasible point found: none of the %d trials " + ...
            "completed with a finite objective value.", trial);
    case -5
        message = sprintf("Optimization stopped: MaxTime (%g s) exceeded.", ...
            options.MaxTime);
    otherwise
        message = "Optimization stopped.";
end
end
