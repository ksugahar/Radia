function tests = test_optuna_got_interface
%TEST_OPTUNA_GOT_INTERFACE The Global Optimization Toolbox-shaped surface.
%   These are MATLAB-integration tests. They lock the toolbox idiom that
%   radia.optuna.optimize adds on top of the study -- the options object, the
%   [x, fval, exitflag, output] contract, the "init"/"iter"/"done" output
%   protocol, and the sdo-style parameter array. Upstream Optuna has no
%   optimoptions and no exitflag, so none of this is evidence of upstream
%   parity, and none of it may be cited as such.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabDirectory = fullfile(root, "matlab");
entries = string(strsplit(path, pathsep));
testCase.TestData.RemovePath = ~any(strcmpi(entries, string(matlabDirectory)));
if testCase.TestData.RemovePath
    addpath(matlabDirectory);
end
testCase.TestData.MatlabDirectory = matlabDirectory;
end

function teardownOnce(testCase)
if testCase.TestData.RemovePath
    rmpath(testCase.TestData.MatlabDirectory);
end
end

function testOptimoptionsMatchesToolboxUsage(testCase)
% optimoptions accepts the three call shapes the toolbox accepts: bare
% Name=Value, a leading solver name, and an existing options object to copy.
defaults = radia.optuna.optimoptions();
verifyEqual(testCase, defaults.MaxTrials, 100);
verifyEqual(testCase, defaults.Display, "final");
verifyEqual(testCase, defaults.MaxTime, Inf);
verifyFalse(testCase, defaults.UseParallel);

named = radia.optuna.optimoptions(MaxTrials=7, Sampler="tpe", Display="off");
verifyEqual(testCase, named.MaxTrials, 7);
verifyEqual(testCase, named.Sampler, "tpe");
verifyTrue(testCase, named.isQuiet());

withSolver = radia.optuna.optimoptions("optuna", MaxTrials=9);
verifyEqual(testCase, withSolver.MaxTrials, 9);

copied = radia.optuna.optimoptions(named, MaxTrials=11);
verifyEqual(testCase, copied.MaxTrials, 11);
verifyEqual(testCase, copied.Sampler, "tpe", ...
    "copying an options object must preserve the options not overridden");

% A misspelled option is an error naming the available options, not a
% silently ignored argument.
verifyError(testCase, @() radia.optuna.optimoptions(MaxTrails=5), ...
    "radia:optuna:OptimOptions");
verifyError(testCase, @() radia.optuna.optimoptions(Display="verbose"), ...
    ?MException);
end

function testOptimizeReturnsToolboxOutputContract(testCase)
parameters = radia.optuna.OptimizationParameter.fromStruct( ...
    struct("a", [-2 2], "b", [-2 2]));
options = radia.optuna.optimoptions(MaxTrials=25, Sampler="random", ...
    Seed=17, Display="off");
[x, fval, exitflag, output] = radia.optuna.optimize( ...
    @(v) (v.a - 0.5)^2 + (v.b + 0.25)^2, parameters, options);

verifyEqual(testCase, sort(string(fieldnames(x)))', ["a" "b"]);
verifyTrue(testCase, isfinite(fval));
verifyEqual(testCase, fval, (x.a - 0.5)^2 + (x.b + 0.25)^2, ...
    "fval must be the objective at the returned x", AbsTol=1e-12);
verifyEqual(testCase, exitflag, 0, "MaxTrials reached is exitflag 0");
verifyEqual(testCase, output.funccount, 25);
verifyEqual(testCase, output.trials, 25);
verifyEqual(testCase, output.prunedcount, 0);
verifyEqual(testCase, output.failedcount, 0);
verifyTrue(testCase, contains(output.message, "MaxTrials"));
verifyClass(testCase, output.study, "radia.optuna.Study");
verifyEqual(testCase, numel(output.study.get_trials()), 25, ...
    "the study behind the toolbox surface holds every trial");

% The reported best really is the best completed trial in the study.
values = arrayfun(@(t) t.Values, output.study.get_trials());
verifyEqual(testCase, fval, min(values), AbsTol=0);
end

function testOutputFcnSeesInitIterDoneAndCanStop(testCase)
states = strings(0, 1);
trials = [];
    function stop = record(~, optimValues, state)
        states(end + 1, 1) = string(state);
        trials(end + 1, 1) = optimValues.trial;
        stop = false;
    end
parameters = radia.optuna.OptimizationParameter.fromStruct(struct("a", [0 1]));
options = radia.optuna.optimoptions(MaxTrials=4, Sampler="random", ...
    Seed=3, Display="off", OutputFcn=@record);
radia.optuna.optimize(@(v) v.a, parameters, options);

verifyEqual(testCase, states(1), "init");
verifyEqual(testCase, states(end), "done");
verifyEqual(testCase, sum(states == "iter"), 4);
verifyEqual(testCase, trials(states == "iter"), (1:4)');

% Returning true stops the run, and the toolbox reports that with -1.
    function stop = stopAtThree(~, optimValues, ~)
        stop = optimValues.trial >= 3;
    end
stopped = radia.optuna.optimoptions(MaxTrials=50, Sampler="random", ...
    Seed=3, Display="off", OutputFcn=@stopAtThree);
[~, ~, exitflag, output] = radia.optuna.optimize(@(v) v.a, parameters, stopped);
verifyEqual(testCase, exitflag, -1);
verifyEqual(testCase, output.trials, 3);
verifyTrue(testCase, contains(output.message, "output or plot function"));
end

function testStallAndTimeAndInfeasibleExitFlags(testCase)
parameters = radia.optuna.OptimizationParameter.fromStruct(struct("a", [0 1]));

% A constant objective never improves, so the stall counter fires.
stall = radia.optuna.optimoptions(MaxTrials=50, MaxStallTrials=5, ...
    FunctionTolerance=1e-6, Sampler="random", Seed=5, Display="off");
[~, ~, exitflag, output] = radia.optuna.optimize(@(~) 1.0, parameters, stall);
verifyEqual(testCase, exitflag, 1);
verifyEqual(testCase, output.trials, 6, ...
    "the first trial sets the incumbent and is not a stall, so five " + ...
    "stalled trials means six in total -- ga counts MaxStallGenerations " + ...
    "the same way");
verifyTrue(testCase, contains(output.message, "FunctionTolerance"));

% An objective that always throws leaves no feasible point.
infeasible = radia.optuna.optimoptions(MaxTrials=3, Sampler="random", ...
    Seed=5, Display="off");
[x, fval, exitflag, output] = radia.optuna.optimize( ...
    @(~) error("cost:blew", "no"), parameters, infeasible);
verifyEqual(testCase, exitflag, -2);
verifyTrue(testCase, isnan(fval));
verifyEmpty(testCase, fieldnames(x));
verifyEqual(testCase, output.failedcount, 3);
verifyTrue(testCase, contains(output.message, "No feasible point"));
verifyTrue(testCase, all(arrayfun(@(t) t.State == "FAIL", ...
    output.study.get_trials())), ...
    "a failed objective must be recorded as FAIL, not dropped");
end

function testFixedParametersAreNotSearchedButAreReturned(testCase)
searched = radia.optuna.OptimizationParameter("a", Minimum=0, Maximum=1);
pinned = radia.optuna.OptimizationParameter("b", Value=0.75, Free=false);
options = radia.optuna.optimoptions(MaxTrials=5, Sampler="random", ...
    Seed=11, Display="off");
[x, ~, ~, output] = radia.optuna.optimize(@(v) v.a + v.b, ...
    [searched pinned], options);
verifyEqual(testCase, x.b, 0.75, ...
    "a fixed parameter must come back with the value the caller pinned");
frozen = output.study.get_trials();
verifyEqual(testCase, sort(string(fieldnames(frozen(1).Params)))', "a", ...
    "a fixed parameter must not enter the study's search space");

% A fixed parameter with no value, and a searched one with no bounds, are
% both refused rather than defaulted to something plausible.
verifyError(testCase, @() radia.optuna.optimize(@(v) v.a, ...
    [searched radia.optuna.OptimizationParameter("c", Free=false)], options), ...
    "radia:optuna:Optimize");
verifyError(testCase, @() radia.optuna.optimize(@(v) v.a, ...
    radia.optuna.OptimizationParameter("unbounded"), options), ...
    "radia:optuna:ParameterBounds");
end

function testParameterTypesMapOntoDistributions(testCase)
continuous = radia.optuna.OptimizationParameter("f", Minimum=0, Maximum=2);
logarithmic = radia.optuna.OptimizationParameter("g", Minimum=1, Maximum=100, ...
    Transform="log");
integral = radia.optuna.OptimizationParameter("k", Minimum=0, Maximum=9, ...
    Type="integer");
categorical = radia.optuna.OptimizationParameter("c", Type="categorical", ...
    Choices=["low" "high"]);

verifyEqual(testCase, continuous.distribution().kind, "float");
verifyFalse(testCase, logical(logarithmic.distribution().log) == false);
verifyEqual(testCase, integral.distribution().kind, "integer");
verifyEqual(testCase, categorical.distribution().kind, "categorical");

% A log-scaled parameter must not reach zero, and a categorical one needs
% choices; both are refused up front rather than at sample time.
bad = radia.optuna.OptimizationParameter("h", Minimum=0, Maximum=1, ...
    Transform="log");
verifyError(testCase, @() bad.distribution(), "radia:optuna:ParameterBounds");
verifyError(testCase, ...
    @() radia.optuna.OptimizationParameter("e", Type="categorical").distribution(), ...
    "radia:optuna:ParameterChoices");

options = radia.optuna.optimoptions(MaxTrials=6, Sampler="random", ...
    Seed=23, Display="off");
[x, ~, ~, output] = radia.optuna.optimize(@(v) double(v.k), ...
    [integral categorical], options);
verifyEqual(testCase, x.k, round(x.k), "an integer parameter must stay integral");
verifyTrue(testCase, ismember(string(x.c), ["low" "high"]));
verifyEqual(testCase, output.trials, 6);
end

function testSamplerFactoryIsSharedAndRefusesAuto(testCase)
% The short names map onto samplers in exactly one place, and that place
% refuses "auto" instead of guessing without the study's shape.
verifyClass(testCase, radia.optuna.internal.samplerFromName("random", 5), ...
    "radia.optuna.RandomSampler");
verifyClass(testCase, radia.optuna.internal.samplerFromName("tpe", 5), ...
    "radia.optuna.TPESampler");
verifyClass(testCase, radia.optuna.internal.samplerFromName("nsgaii", 5), ...
    "radia.optuna.NSGAIISampler");
verifyError(testCase, @() radia.optuna.internal.samplerFromName("auto", 1), ...
    "radia:optuna:SamplerName");
verifyError(testCase, @() radia.optuna.internal.samplerFromName("bogus", 1), ...
    "radia:optuna:SamplerName");

% An equal seed through the factory reproduces the study's own sequence.
first = radia.optuna.Study(Sampler=radia.optuna.internal.samplerFromName( ...
    "random", 31), AutoSave=false);
second = radia.optuna.Study(Sampler=radia.optuna.RandomSampler(31), ...
    AutoSave=false);
verifyEqual(testCase, first.ask().suggest_float("x", 0, 1), ...
    second.ask().suggest_float("x", 0, 1), AbsTol=0);
end

function testUseParallelRefusesWithoutAPool(testCase)
% UseParallel evaluates a batch of already-suggested parameter sets, so it
% needs a pool. Without one it says so rather than quietly running serially,
% which would make a timing comparison meaningless.
if ~isempty(gcp("nocreate"))
    return
end
parameters = radia.optuna.OptimizationParameter.fromStruct(struct("a", [0 1]));
options = radia.optuna.optimoptions(MaxTrials=4, Sampler="random", ...
    Seed=7, Display="off", UseParallel=true);
verifyError(testCase, @() radia.optuna.optimize(@(v) v.a, parameters, options), ...
    "radia:optuna:UseParallel");
end
