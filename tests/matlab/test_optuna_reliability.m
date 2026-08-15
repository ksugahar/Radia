function tests = test_optuna_reliability
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabDirectory = fullfile(repositoryRoot, "matlab");
entries = string(strsplit(path, pathsep));
testCase.TestData.RemoveMatlabDirectory = ...
    ~any(strcmpi(entries, string(matlabDirectory)));
if testCase.TestData.RemoveMatlabDirectory
    addpath(matlabDirectory);
end
testCase.TestData.MatlabDirectory = matlabDirectory;
end

function teardownOnce(testCase)
if testCase.TestData.RemoveMatlabDirectory
    rmpath(testCase.TestData.MatlabDirectory);
end
end

function testConstrainedBestTrialRequiresFeasibility(testCase)
study = radia.optuna.Study( ...
    Sampler=radia.optuna.RandomSampler(1), AutoSave=false);

infeasible = study.ask();
study.tell(infeasible, 0);
study.recordConstraints(infeasible, [1, -1]);

feasible = study.ask();
study.tell(feasible, 2);
study.recordConstraints(feasible, [-1, 0]);

nonfinite = study.ask();
study.tell(nonfinite, -100);
study.recordConstraints(nonfinite, [Inf, -Inf]);

verifyEqual(testCase, study.bestTrial().TrialNumber, feasible.Number);
verifyEqual(testCase, study.bestValue(), 2);

allInfeasible = radia.optuna.Study(AutoSave=false);
trial = allInfeasible.ask();
allInfeasible.tell(trial, -10);
allInfeasible.recordConstraints(trial, 0.1);
verifyError(testCase, @()allInfeasible.bestTrial(), ...
    "radia:optuna:NoFeasibleTrial");
end

function testConstrainedParetoFrontExcludesUncheckedAndInfeasible(testCase)
study = radia.optuna.Study(Directions=["minimize", "minimize"], ...
    Sampler=radia.optuna.RandomSampler(2), AutoSave=false);
values = [0, 0; 1, 2; 2, 1; -1, -1];
constraints = {1, -1, 0, []};
for index = 1:size(values, 1)
    trial = study.ask();
    study.tell(trial, values(index, :));
    if ~isempty(constraints{index})
        study.recordConstraints(trial, constraints{index});
    end
end

front = study.paretoFront();
verifyEqual(testCase, sort(front.TrialNumber), [1; 2]);
end

function testConstrainedTPERanksFeasibleTrialsBeforeLowerObjectives(testCase)
sampler = radia.optuna.TPESampler( ...
    Seed=9, NStartupTrials=8, Gamma=0.25, ...
    ConstraintsFcn=@(trial)0.8 - trial.Params.x);
study = radia.optuna.Study(Sampler=sampler, AutoSave=false);
study.optimize(@(trial)trial.suggestFloat("x", 0, 1)^2, 40);

postStartup = study.ParamTable.ValueNumeric(21:end);
verifyGreaterThan(testCase, mean(postStartup), 0.75);
verifyGreaterThanOrEqual(testCase, sum(postStartup >= 0.8), 10);
verifyLessThanOrEqual(testCase, ...
    study.constraintsForTrial(study.bestTrial().TrialNumber), 0);
end

function testTellIsStrictByDefaultAndExplicitlyIdempotent(testCase)
study = radia.optuna.Study(AutoSave=false);
trial = study.ask();
study.tell(trial, 1);

verifyError(testCase, @()study.tell(trial, 1), ...
    "radia:optuna:TrialState");
study.tell(trial, 999, SkipIfFinished=true);
study.fail(trial, "late failure", SkipIfFinished=true);
verifyEqual(testCase, study.TrialTable.Value, 1);
verifyEqual(testCase, study.TrialTable.State, "COMPLETE");

other = radia.optuna.Study(AutoSave=false);
foreign = other.ask();
verifyError(testCase, @()study.tell(foreign, 0), ...
    "radia:optuna:StudyMismatch");
end

function testConstraintCallbackFailureStillCommitsCompleteTrial(testCase)
storagePath = string(tempname("C:\temp")) + ".mat";
cleanup = onCleanup(@()deleteStudyStorage(storagePath));
sampler = radia.optuna.TPESampler( ...
    ConstraintsFcn=@(~)constraintCallbackFailure());
study = radia.optuna.Study( ...
    StoragePath=storagePath, Sampler=sampler, AutoSave=true);

verifyError(testCase, ...
    @()study.optimize(@(trial)trial.suggestFloat("x", 0, 1), 1), ...
    "radia:test:ConstraintCallback");
verifyEqual(testCase, study.TrialTable.State, "COMPLETE");
verifyTrue(testCase, isfinite(study.TrialTable.Value));
verifyEqual(testCase, height(study.ObjectiveTable), 1);
verifyEqual(testCase, height(study.ConstraintTable), 0);
reloaded = radia.optuna.Study(StoragePath=storagePath, ...
    Sampler=radia.optuna.TPESampler(), AutoSave=false);
verifyEqual(testCase, reloaded.TrialTable.State, "COMPLETE");
verifyEqual(testCase, height(reloaded.ObjectiveTable), 1);
clear cleanup
deleteStudyStorage(storagePath);
end

function testUnsupportedConstrainedCMAFailsClosed(testCase)
study = radia.optuna.Study(Sampler=radia.optuna.CmaEsSampler(Seed=3), ...
    AutoSave=false);
trial = study.ask();
verifyError(testCase, @()study.recordConstraints(trial, 1), ...
    "radia:optuna:SamplerConstraints");
study.fail(trial, "unsupported constrained sampler");
end

function testStaleRunningTrialRecoveryLeavesNoRunningRow(testCase)
study = radia.optuna.Study(AutoSave=false);
trial = study.ask();
recovered = study.recoverStaleRunning(0, Message="worker lease expired");

verifyEqual(testCase, recovered.TrialNumber, trial.Number);
verifyEqual(testCase, recovered.State, "FAIL");
verifyEqual(testCase, recovered.ErrorMessage, "worker lease expired");
verifyEqual(testCase, sum(study.TrialTable.State == "RUNNING"), 0);
study.tell(trial, 1, SkipIfFinished=true);
verifyError(testCase, @()study.tell(trial, 1), ...
    "radia:optuna:TrialState");
end

function testAtomicStorageSchemaAndBackupRecovery(testCase)
storagePath = string(tempname("C:\temp")) + ".mat";
cleanup = onCleanup(@()deleteStudyStorage(storagePath));
study = radia.optuna.Study(Name="recoverable", ...
    StoragePath=storagePath, Sampler=radia.optuna.RandomSampler(3));
trial = study.ask();
trial.suggestFloat("x", -1, 1);
study.tell(trial, 0.25);

verifyTrue(testCase, isfile(storagePath));
verifyTrue(testCase, isfile(storagePath + ".bak"));
loaded = load(storagePath, "StudyData", "-mat");
verifyEqual(testCase, string(loaded.StudyData.Schema), ...
    "radia.optuna.study");
verifyEqual(testCase, loaded.StudyData.Version, 4);

overwriteFile(storagePath, "damaged primary");
lastwarn("");
warning("on", "radia:optuna:RecoveredStorage");
reloaded = radia.optuna.Study(StoragePath=storagePath, ...
    Sampler=radia.optuna.RandomSampler(3), AutoSave=false);
[~, warningId] = lastwarn;
verifyEqual(testCase, string(warningId), ...
    "radia:optuna:RecoveredStorage");
verifyEqual(testCase, reloaded.Name, "recoverable");
verifyEqual(testCase, reloaded.TrialTable.State, "COMPLETE");

clear cleanup
deleteStudyStorage(storagePath);
end

function testEmptyConstraintRecordSurvivesStorageRoundTrip(testCase)
storagePath = string(tempname("C:\temp")) + ".mat";
cleanup = onCleanup(@()deleteStudyStorage(storagePath));
study = radia.optuna.Study(StoragePath=storagePath, ...
    Sampler=radia.optuna.NSGAIISampler(PopulationSize=2));
trial = study.ask();
study.tell(trial,0);
study.recordConstraints(trial,zeros(1,0));
reloaded = radia.optuna.Study(StoragePath=storagePath, ...
    Sampler=radia.optuna.NSGAIISampler(PopulationSize=2), ...
    AutoSave=false);
[present,values] = reloaded.constraintRecord(trial.Number);
verifyTrue(testCase,present);
verifyEmpty(testCase,values);
verifyEqual(testCase,reloaded.ConstraintCountTable.Count,0);
clear cleanup
deleteStudyStorage(storagePath);
end

function testSamplerRandomStatesResumeDeterministically(testCase)
verifySamplerResume(testCase, "random");
verifySamplerResume(testCase, "tpe");
verifySamplerResume(testCase, "motpe");
verifySamplerResume(testCase, "nsgaii");
end

function testReusedSamplerStartsEachNewStudyFromConfiguredSeed(testCase)
for samplerName = ["random", "tpe", "motpe", "nsgaii"]
    [sharedSampler, directions] = makeSampler(samplerName, 23);
    firstStudy = radia.optuna.Study( ...
        Directions=directions, Sampler=sharedSampler, AutoSave=false);
    firstTrial = firstStudy.ask();
    firstTrial.suggestFloat("x", -1, 1);

    secondStudy = radia.optuna.Study( ...
        Directions=directions, Sampler=sharedSampler, AutoSave=false);
    secondTrial = secondStudy.ask();
    reusedValue = secondTrial.suggestFloat("x", -1, 1);

    [freshSampler, ~] = makeSampler(samplerName, 23);
    referenceStudy = radia.optuna.Study( ...
        Directions=directions, Sampler=freshSampler, AutoSave=false);
    referenceTrial = referenceStudy.ask();
    referenceValue = referenceTrial.suggestFloat("x", -1, 1);
    verifyEqual(testCase, reusedValue, referenceValue, AbsTol=0);
end
end

function testOptunaLikeDefaultSamplerSelection(testCase)
singleObjective = radia.optuna.Study(AutoSave=false);
multiObjective = radia.optuna.Study( ...
    Directions=["minimize", "maximize"], AutoSave=false);

verifyClass(testCase, singleObjective.Sampler, ...
    "radia.optuna.TPESampler");
verifyClass(testCase, multiObjective.Sampler, ...
    "radia.optuna.MOTPESampler");
end

function testIntegerStepDistributionCompatibilityAndNameCollision(testCase)
study = radia.optuna.Study( ...
    Sampler=radia.optuna.RandomSampler(4), AutoSave=false);
trial = study.ask();
value = trial.suggest_int("mesh order", 2, 10, Step=3);
verifyTrue(testCase, ismember(value, [2, 5, 8]));
verifyEqual(testCase, trial.suggest_int("mesh order", 2, 10, Step=3), ...
    value);
verifyError(testCase, ...
    @()trial.suggest_int("mesh order", 2, 10, Step=2), ...
    "radia:optuna:IncompatibleDistribution");
verifyError(testCase, ...
    @()trial.suggestFloat("mesh order", 2, 10), ...
    "radia:optuna:IncompatibleDistribution");

collision = study.ask();
collision.suggestFloat("a-b", 0, 1);
verifyError(testCase, @()collision.suggestFloat("a_b", 0, 1), ...
    "radia:optuna:ParameterNameCollision");
end

function testMedianPrunerUsesCurrentTrialBestIntermediate(testCase)
study = radia.optuna.Study( ...
    Sampler=radia.optuna.RandomSampler(5), ...
    Pruner=radia.optuna.MedianPruner( ...
        NStartupTrials=0, MinCompletedTrials=1), ...
    AutoSave=false);
reference = study.ask();
reference.report(1, 1);
study.tell(reference, 1);

candidate = study.ask();
candidate.report(0, 0);
candidate.report(2, 1);
verifyFalse(testCase, candidate.shouldPrune());
end

function verifySamplerResume(testCase, samplerName)
storagePath = string(tempname("C:\temp")) + ".mat";
cleanup = onCleanup(@()deleteStudyStorage(storagePath));
[continuousSampler, directions] = makeSampler(samplerName, 17);
[persistedSampler, ~] = makeSampler(samplerName, 17);
continuous = radia.optuna.Study(Directions=directions, ...
    Sampler=continuousSampler, AutoSave=false);
persisted = radia.optuna.Study(Directions=directions, ...
    StoragePath=storagePath, Sampler=persistedSampler, AutoSave=true);

for index = 1:5
    continuousTrial = continuous.ask();
    persistedTrial = persisted.ask();
    continuousX = continuousTrial.suggestFloat("x", -2, 2);
    persistedX = persistedTrial.suggestFloat("x", -2, 2);
    verifyEqual(testCase, persistedX, continuousX, AbsTol=0);
    if isscalar(directions)
        continuousValue = (continuousX - 0.3)^2;
        persistedValue = (persistedX - 0.3)^2;
    else
        continuousValue = [continuousX^2, (continuousX - 1)^2];
        persistedValue = [persistedX^2, (persistedX - 1)^2];
    end
    continuous.tell(continuousTrial, continuousValue);
    persisted.tell(persistedTrial, persistedValue);
end

[reloadedSampler, ~] = makeSampler(samplerName, 17);
reloaded = radia.optuna.Study(Directions=directions, ...
    StoragePath=storagePath, Sampler=reloadedSampler, AutoSave=false);
continuousTrial = continuous.ask();
reloadedTrial = reloaded.ask();
continuousX = continuousTrial.suggestFloat("x", -2, 2);
reloadedX = reloadedTrial.suggestFloat("x", -2, 2);
verifyEqual(testCase, reloadedX, continuousX, AbsTol=0);

clear cleanup
deleteStudyStorage(storagePath);
end

function [sampler, directions] = makeSampler(name, seed)
switch name
    case "random"
        sampler = radia.optuna.RandomSampler(seed);
        directions = "minimize";
    case "tpe"
        sampler = radia.optuna.TPESampler( ...
            Seed=seed, NStartupTrials=2, Multivariate=true);
        directions = "minimize";
    case "motpe"
        sampler = radia.optuna.MOTPESampler( ...
            Seed=seed, NStartupTrials=2);
        directions = ["minimize", "minimize"];
    case "nsgaii"
        sampler = radia.optuna.NSGAIISampler( ...
            Seed=seed, PopulationSize=2);
        directions = ["minimize", "minimize"];
    otherwise
        error("radia:test:Sampler", "Unknown sampler '%s'.", name);
end
end

function overwriteFile(path, text)
stream = fopen(path, "w");
if stream < 0
    error("radia:test:Storage", "Cannot overwrite '%s'.", path);
end
cleanup = onCleanup(@()fclose(stream));
fwrite(stream, char(text), "char");
clear cleanup
end

function value = constraintCallbackFailure()
value = NaN; %#ok<NASGU>
error("radia:test:ConstraintCallback", ...
    "Constraint evaluation failed after the objective returned.");
end

function deleteStudyStorage(path)
paths = [string(path), string(path) + ".bak"];
for candidate = paths
    if isfile(candidate)
        delete(candidate);
    end
end
end
