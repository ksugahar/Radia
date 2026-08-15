function tests=test_optuna_pruner_parity
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabDirectory=fullfile(root,"matlab");
entries=string(strsplit(path,pathsep));
testCase.TestData.RemovePath= ...
    ~any(strcmpi(entries,string(matlabDirectory)));
if testCase.TestData.RemovePath, addpath(matlabDirectory); end
testCase.TestData.MatlabDirectory=matlabDirectory;
end

function teardownOnce(testCase)
if testCase.TestData.RemovePath
    rmpath(testCase.TestData.MatlabDirectory);
end
end

function testReportMatchesOptunaSemantics(testCase)
study=radia.optuna.Study(Pruner=radia.optuna.NopPruner(),AutoSave=false);
trial=study.ask();
trial.report(NaN,0);
verifyWarning(testCase,@()trial.report(7,0), ...
    "radia:optuna:DuplicateReport");
verifyTrue(testCase,isnan(trial.IntermediateValues.Value(1)));
verifyError(testCase,@()trial.report(1,-1),"radia:optuna:Report");
verifyFalse(testCase,trial.shouldPrune());
end

function testPercentileUsesBestAndSparseIntervals(testCase)
pruner=radia.optuna.PercentilePruner(50,NStartupTrials=0, ...
    NWarmupSteps=0,IntervalSteps=2,NMinTrials=1);
study=radia.optuna.Study(Pruner=pruner,AutoSave=false);
addCompleted(study,[1 3],[1 1]);
addCompleted(study,[1 3],[3 3]);
trial=study.ask();
trial.report(5,1);
trial.report(4,3);
verifyTrue(testCase,trial.shouldPrune());

maximize=radia.optuna.Study(Directions="maximize",Pruner= ...
    radia.optuna.PercentilePruner(50,NStartupTrials=0),AutoSave=false);
addCompleted(maximize,0,1);
addCompleted(maximize,0,3);
candidate=maximize.ask();
candidate.report(1.5,0);
verifyTrue(testCase,candidate.shouldPrune());
end

function testThresholdUsesLatestValueAndNaN(testCase)
study=radia.optuna.Study(Pruner=radia.optuna.ThresholdPruner( ...
    Lower=0,Upper=10),AutoSave=false);
trial=study.ask();
trial.report(5,0);
verifyFalse(testCase,trial.shouldPrune());
trial.report(11,1);
verifyTrue(testCase,trial.shouldPrune());
nanTrial=study.ask();
nanTrial.report(NaN,0);
verifyTrue(testCase,nanTrial.shouldPrune());
end

function testPatientRequiresAFullNoImprovementWindow(testCase)
study=radia.optuna.Study(Pruner=radia.optuna.PatientPruner([], ...
    Patience=1,MinDelta=0),AutoSave=false);
trial=study.ask();
values=[10 4 5 6];
for index=1:4
    trial.report(values(index),index-1);
end
verifyTrue(testCase,trial.shouldPrune());

wrapped=radia.optuna.PatientPruner(radia.optuna.NopPruner(),Patience=1);
study2=radia.optuna.Study(Pruner=wrapped,AutoSave=false);
trial2=study2.ask();
values=[10 4 5 6];
for index=1:4
    trial2.report(values(index),index-1);
end
verifyFalse(testCase,trial2.shouldPrune());
end

function testSuccessiveHalvingPromotionAndRungState(testCase)
pruner=radia.optuna.SuccessiveHalvingPruner( ...
    MinResource=1,ReductionFactor=2);
study=radia.optuna.Study(Pruner=pruner,AutoSave=false);
first=study.ask();
first.report(1,1);
verifyFalse(testCase,first.shouldPrune());
verifyTrue(testCase,isfield(first.SystemAttrs,"completed_rung_0"));
study.tell(first,1);

second=study.ask();
second.report(2,1);
verifyTrue(testCase,second.shouldPrune());
verifyEqual(testCase,second.SystemAttrs.completed_rung_0,2);
end

function testSuccessiveHalvingBootstrapPrunesFirstRung(testCase)
study=radia.optuna.Study(Pruner= ...
    radia.optuna.SuccessiveHalvingPruner(MinResource=1, ...
    ReductionFactor=2,BootstrapCount=1),AutoSave=false);
trial=study.ask();
trial.report(1,1);
verifyTrue(testCase,trial.shouldPrune());
end

function testHyperbandUsesOptunaCrc32BracketAssignment(testCase)
pruner=radia.optuna.HyperbandPruner(MinResource=1, ...
    MaxResource=9,ReductionFactor=3);
study=radia.optuna.Study(Name="hb",Pruner=pruner,AutoSave=false);
trial=study.ask();
trial.report(1,0);
verifyFalse(testCase,trial.shouldPrune());
actual=zeros(1,10);
for number=0:9
    actual(number+1)=pruner.bracketId(study,number);
end
verifyEqual(testCase,actual,[0 0 2 2 1 1 1 2 1 0]);
end

function testWilcoxonPrunesPairedConsistentlyWorseTrial(testCase)
study=radia.optuna.Study(Pruner=radia.optuna.WilcoxonPruner( ...
    PThreshold=0.1,NStartupSteps=2),AutoSave=false);
addCompleted(study,0:5,zeros(1,6));
trial=study.ask();
for step=0:5
    trial.report(10,step);
end
verifyTrue(testCase,trial.shouldPrune());

nonfinite=study.ask();
nonfinite.report(Inf,0);
verifyWarning(testCase,@()nonfinite.shouldPrune(), ...
    "radia:optuna:WilcoxonNonfinite");
verifyFalse(testCase,nonfinite.shouldPrune());
end

function addCompleted(study,steps,values)
trial=study.ask();
for index=1:numel(steps)
    trial.report(values(index),steps(index));
end
study.tell(trial,values(end));
end
