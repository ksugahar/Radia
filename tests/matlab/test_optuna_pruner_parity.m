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
fixture=fullfile(root,"tests","matlab","fixtures","optuna49_oracle.json");
testCase.TestData.Oracle=jsondecode(fileread(fixture));
end

function teardownOnce(testCase)
if testCase.TestData.RemovePath
    rmpath(testCase.TestData.MatlabDirectory);
end
end

function testPrunerDecisionsMatchUpstream(testCase)
expected=testCase.TestData.Oracle.pruners;
percentile=radia.optuna.Study(Pruner=radia.optuna.PercentilePruner( ...
    50,NStartupTrials=0,NWarmupSteps=0,IntervalSteps=2,NMinTrials=1), ...
    AutoSave=false);
addCompletedTrial(percentile,[1,3],[1,1]);
addCompletedTrial(percentile,[1,3],[3,3]);
trial=percentile.ask();
trial.report(5,1); trial.report(4,3);
verifyEqual(testCase,trial.shouldPrune(), ...
    logical(expected.percentile_minimize));

maximize=radia.optuna.Study(Directions="maximize",Pruner= ...
    radia.optuna.PercentilePruner(50,NStartupTrials=0),AutoSave=false);
addCompletedTrial(maximize,0,1); addCompletedTrial(maximize,0,3);
trial=maximize.ask(); trial.report(1.5,0);
verifyEqual(testCase,trial.shouldPrune(), ...
    logical(expected.percentile_maximize));

threshold=radia.optuna.Study(Pruner=radia.optuna.ThresholdPruner( ...
    Lower=0,Upper=10),AutoSave=false);
trial=threshold.ask(); trial.report(5,0);
actual(1)=trial.shouldPrune();
trial.report(11,1); actual(2)=trial.shouldPrune();
nanTrial=threshold.ask(); nanTrial.report(NaN,0);
actual(3)=nanTrial.shouldPrune();
verifyEqual(testCase,actual,reshape(logical(expected.threshold),1,[]));

patient=radia.optuna.Study(Pruner=radia.optuna.PatientPruner([], ...
    Patience=1,MinDelta=0),AutoSave=false);
trial=patient.ask();
patientValues=[10,4,5,6];
for index=1:4, trial.report(patientValues(index),index-1); end
wrapped=radia.optuna.Study(Pruner=radia.optuna.PatientPruner( ...
    radia.optuna.NopPruner(),Patience=1),AutoSave=false);
wrappedTrial=wrapped.ask();
for index=1:4, wrappedTrial.report(patientValues(index),index-1); end
verifyEqual(testCase,[trial.shouldPrune(),wrappedTrial.shouldPrune()], ...
    reshape(logical(expected.patient),1,[]));

halvingPruner=radia.optuna.SuccessiveHalvingPruner( ...
    MinResource=1,ReductionFactor=2);
halving=radia.optuna.Study(Pruner=halvingPruner,AutoSave=false);
first=halving.ask(); first.report(1,1);
firstDecision=first.shouldPrune(); halving.tell(first,1);
second=halving.ask(); second.report(2,1);
secondDecision=second.shouldPrune();
verifyEqual(testCase,[firstDecision,secondDecision], ...
    reshape(logical(expected.successive_halving.decisions),1,[]));
verifyEqual(testCase,[first.SystemAttrs.completed_rung_0, ...
    second.SystemAttrs.completed_rung_0], ...
    reshape(double(expected.successive_halving.rung_values),1,[]));
bootstrap=radia.optuna.Study(Pruner= ...
    radia.optuna.SuccessiveHalvingPruner(MinResource=1, ...
    ReductionFactor=2,BootstrapCount=1),AutoSave=false);
bootstrapTrial=bootstrap.ask(); bootstrapTrial.report(1,1);
verifyEqual(testCase,bootstrapTrial.shouldPrune(), ...
    logical(expected.successive_halving.bootstrap));

hyperbandPruner=radia.optuna.HyperbandPruner(MinResource=1, ...
    MaxResource=9,ReductionFactor=3);
hyperband=radia.optuna.Study(Name="hb",Pruner=hyperbandPruner,AutoSave=false);
hyperbandTrial=hyperband.ask(); hyperbandTrial.report(1,0);
verifyEqual(testCase,hyperbandTrial.shouldPrune(), ...
    logical(expected.hyperband.first_decision));
brackets=zeros(1,10);
for number=0:9, brackets(number+1)=hyperbandPruner.bracketId(hyperband,number); end
verifyEqual(testCase,brackets, ...
    reshape(double(expected.hyperband.bracket_ids),1,[]));

wilcoxon=radia.optuna.Study(Pruner=radia.optuna.WilcoxonPruner( ...
    PThreshold=0.1,NStartupSteps=2),AutoSave=false);
addCompletedTrial(wilcoxon,0:5,zeros(1,6));
wilcoxonTrial=wilcoxon.ask();
for step=0:5, wilcoxonTrial.report(10,step); end
nonfinite=wilcoxon.ask(); nonfinite.report(Inf,0);
warning("off","radia:optuna:WilcoxonNonfinite");
cleanup=onCleanup(@()warning("on","radia:optuna:WilcoxonNonfinite"));
verifyEqual(testCase,[wilcoxonTrial.shouldPrune(),nonfinite.shouldPrune()], ...
    reshape(logical(expected.wilcoxon),1,[]));
clear cleanup
end

function addCompletedTrial(study,steps,values)
trial=study.ask();
for index=1:numel(steps), trial.report(values(index),steps(index)); end
study.tell(trial,values(end));
end
