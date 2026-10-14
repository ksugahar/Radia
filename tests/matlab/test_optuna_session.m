function tests = test_optuna_session
%TEST_OPTUNA_SESSION MATLAB-only persistent student-workflow integration.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabDirectory = fullfile(root, "matlab");
entries = string(strsplit(path, pathsep));
testCase.TestData.RemovePath = ...
    ~any(strcmpi(entries, string(matlabDirectory)));
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

function testLifecycleSelectAndApply(testCase)
[storage, cleanup] = temporaryStorage(); %#ok<ASGLU>
parameter = radia.optuna.OptimizationParameter("session_x", ...
    Minimum=-1, Maximum=1);
options = radia.optuna.optimoptions(MaxTrials=6, Sampler="random", ...
    Seed=17, Display="off", StoragePath=storage);
session = radia.optuna.OptimizationSession( ...
    @(v) (v.session_x - 0.2)^2, parameter, options);
verifyEqual(testCase, session.State, "configured");

session.start();
first = session.runNext();
verifyEqual(testCase, first.State, "COMPLETE");
session.pause();
verifyEqual(testCase, session.State, "paused");
snapshot = session.snapshot();
verifyEqual(testCase, snapshot.complete, 1);
verifyEqual(testCase, snapshot.attempted, 1);

session.resume();
session.run();
verifyEqual(testCase, session.State, "completed");
snapshot = session.snapshot();
verifyEqual(testCase, snapshot.complete, 6);
verifyEqual(testCase, snapshot.failed, 0);

best = session.Study.bestSolution();
session.selectTrial(best.trial_number);
selected = session.selectedParameters();
verifyEqual(testCase, selected.session_x, best.params.session_x, AbsTol=0);
baseCleanup = onCleanup(@() evalin("base", ...
    "clear('session_x')"));
session.applySelectedToModel(Target="base");
verifyEqual(testCase, evalin("base", "session_x"), ...
    selected.session_x, AbsTol=0);
end

function testCheckpointRestoreAndResume(testCase)
[storage, cleanup] = temporaryStorage(); %#ok<ASGLU>
parameter = radia.optuna.OptimizationParameter("x", ...
    Minimum=0, Maximum=1);
options = radia.optuna.optimoptions(MaxTrials=4, Sampler="random", ...
    Seed=23, Display="off", StoragePath=storage);
session = radia.optuna.OptimizationSession( ...
    @(v) (v.x - 0.4)^2, parameter, options);
session.start();
session.runNext();
session.pause();
sessionPath = session.save();
verifyTrue(testCase, isfile(sessionPath));

restored = radia.optuna.OptimizationSession.load(sessionPath, ...
    ObjectiveFcn=@(v) (v.x - 0.4)^2);
verifyEqual(testCase, restored.State, "paused");
verifyEqual(testCase, restored.snapshot().complete, 1);
restored.resume();
restored.run();
verifyEqual(testCase, restored.State, "completed");
verifyEqual(testCase, restored.snapshot().complete, 4);
end

function testRestoreFailsOrphanedRunningTrial(testCase)
[storage, cleanup] = temporaryStorage(); %#ok<ASGLU>
parameter = radia.optuna.OptimizationParameter("x", ...
    Minimum=0, Maximum=1);
options = radia.optuna.optimoptions(MaxTrials=3, Sampler="random", ...
    Seed=29, Display="off", StoragePath=storage);
session = radia.optuna.OptimizationSession(@(v) v.x, parameter, options);
session.start();
session.Study.ask();

restored = radia.optuna.OptimizationSession.load( ...
    storage + ".session.mat", ObjectiveFcn=@(v) v.x);
snapshot = restored.snapshot();
verifyEqual(testCase, restored.State, "paused");
verifyEqual(testCase, restored.StopReason, "recovered-stale-running");
verifyEqual(testCase, snapshot.failed, 1);
verifyEqual(testCase, snapshot.running, 0);
end

function testStructuredObjectiveConstraintsAndIntermediateValues(testCase)
parameter = radia.optuna.OptimizationParameter("x", ...
    Minimum=0, Maximum=1);
options = radia.optuna.optimoptions(MaxTrials=5, Sampler="random", ...
    Seed=31, Display="off");
session = radia.optuna.OptimizationSession( ...
    @structuredObjective, parameter, options);
session.run();

verifyEqual(testCase, height(session.Study.IntermediateTable), 10);
verifyEqual(testCase, height(session.Study.ConstraintTable), 5);
verifyTrue(testCase, all(session.Study.TrialTable.State == "COMPLETE"), ...
    "Constraint violations remain completed trials; feasibility is " + ...
    "recorded separately, as in Optuna.");
best = session.Study.bestSolution();
verifyLessThanOrEqual(testCase, best.params.x, 0.5);
end

function testExplicitPruningIsRecorded(testCase)
parameter = radia.optuna.OptimizationParameter("x", ...
    Minimum=0, Maximum=1);
options = radia.optuna.optimoptions(MaxTrials=3, Sampler="random", ...
    Seed=37, Display="off");
session = radia.optuna.OptimizationSession( ...
    @alwaysPruned, parameter, options);
session.run();
verifyEqual(testCase, session.snapshot().pruned, 3);
verifyTrue(testCase, all(session.Study.TrialTable.State == "PRUNED"));
verifyEqual(testCase, height(session.Study.IntermediateTable), 3);
end

function testMultiobjectiveParetoAndSelection(testCase)
parameter = radia.optuna.OptimizationParameter("x", ...
    Minimum=0, Maximum=1);
options = radia.optuna.optimoptions(MaxTrials=8, Sampler="nsgaii", ...
    Seed=41, Display="off", ...
    Directions=["minimize", "minimize"]);
session = radia.optuna.OptimizationSession( ...
    @(v) [v.x^2, (v.x - 1)^2], parameter, options);
session.run();
pareto = session.Study.paretoFront();
verifyGreaterThan(testCase, height(pareto), 0);
session.selectTrial(pareto.TrialNumber(1));
verifyTrue(testCase, isfinite(session.selectedParameters().x));
verifyTrue(testCase, isnan(session.snapshot().bestValue), ...
    "A multiobjective session exposes Pareto results, not a fake scalar best.");
end

function testEqualSeedProducesEqualSequence(testCase)
parameter = radia.optuna.OptimizationParameter("x", ...
    Minimum=-2, Maximum=2);
options = radia.optuna.optimoptions(MaxTrials=7, Sampler="random", ...
    Seed=43, Display="off");
first = radia.optuna.OptimizationSession(@(v) v.x^2, parameter, options);
second = radia.optuna.OptimizationSession(@(v) v.x^2, parameter, options);
first.run();
second.run();
verifyEqual(testCase, first.Study.ParamTable.ValueNumeric, ...
    second.Study.ParamTable.ValueNumeric, AbsTol=0);
end

function testCancelIsIdempotent(testCase)
parameter = radia.optuna.OptimizationParameter("x", ...
    Minimum=0, Maximum=1);
options = radia.optuna.optimoptions(MaxTrials=2, Display="off");
session = radia.optuna.OptimizationSession(@(v) v.x, ...
    parameter, options);
session.cancel();
session.cancel();
verifyEqual(testCase, session.State, "cancelled");
verifyError(testCase, @() session.runNext(), ...
    "radia:optuna:SessionState");
end

function result = structuredObjective(values, ~)
objective = (values.x - 0.25)^2;
result = struct( ...
    "Objectives", objective, ...
    "InequalityConstraints", values.x - 0.5, ...
    "IntermediateValues", [0, objective + 1; 1, objective]);
end

function result = alwaysPruned(values, context)
context.report(values.x, 0);
context.prune();
result = values.x;
end

function [pathValue, cleanup] = temporaryStorage()
directory = string(tempname("C:\temp"));
mkdir(directory);
pathValue = directory + "\study.mat";
cleanup = onCleanup(@() removeTemporaryDirectory(directory));
end

function removeTemporaryDirectory(directory)
if isfolder(directory)
    rmdir(directory, "s");
end
end
