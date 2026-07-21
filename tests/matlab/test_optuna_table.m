function tests = test_optuna_table
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabDirectory = fullfile(repositoryRoot, "matlab");
addpath(matlabDirectory);
testCase.TestData.MatlabDirectory = matlabDirectory;
end

function teardownOnce(testCase)
rmpath(testCase.TestData.MatlabDirectory);
end

function testDefineByRunAndPersistence(testCase)
path = string(tempname("C:\temp")) + ".mat";
cleanup = onCleanup(@() deleteIfPresent(path));
study = radia.optuna.Study(Name="table-study", StoragePath=path, ...
    Sampler=radia.optuna.RandomSampler(7));
trial = study.ask();
x = trial.suggestFloat("frequency", 1e3, 1e5, Log=true);
p = trial.suggestInteger("p_order", 2, 6);
mode = trial.suggestCategorical("mode", ["FEM", "VIM"]);
trial.report(x / 1e5, 0);
study.tell(trial, x / 1e5);

verifyTrue(testCase, isfinite(x));
verifyTrue(testCase, ismember(p, 2:6));
verifyTrue(testCase, ismember(mode, ["FEM", "VIM"]));
verifyEqual(testCase, height(study.TrialTable), 1);
verifyEqual(testCase, height(study.ParamTable), 3);
verifyEqual(testCase, height(study.IntermediateTable), 1);
verifyEqual(testCase, study.TrialTable.State(1), "COMPLETE");

reloaded = radia.optuna.Study(Name="ignored", StoragePath=path);
verifyEqual(testCase, reloaded.Name, "table-study");
verifyEqual(testCase, height(reloaded.TrialTable), 1);
best = reloaded.bestTrial();
verifyEqual(testCase, best.TrialNumber, 0);
verifyEqual(testCase, reloaded.bestValue(), x / 1e5, "AbsTol", 1e-12);
bestParams = reloaded.bestParams();
verifyEqual(testCase, bestParams.frequency, x, "AbsTol", 1e-12);
verifyEqual(testCase, bestParams.p_order, p);
verifyEqual(testCase, bestParams.mode, mode);
snapshot = reloaded.bestSolution();
verifyTrue(testCase, snapshot.available);
verifyEqual(testCase, snapshot.trial_number, 0);
verifyEqual(testCase, snapshot.value, x / 1e5, "AbsTol", 1e-12);
verifyEqual(testCase, snapshot.params, bestParams);
verifyEqual(testCase, reloaded.best_value(), reloaded.bestValue(), ...
    "AbsTol", 1e-12);
verifyEqual(testCase, reloaded.best_params(), bestParams);
verifyEqual(testCase, reloaded.best_solution().trial_number, 0);
end

function testMedianPruningAndOptimize(testCase)
study = radia.optuna.Study( ...
    Sampler=radia.optuna.RandomSampler(11), ...
    Pruner=radia.optuna.MedianPruner(NStartupTrials=0, MinCompletedTrials=1));
first = study.ask();
first.report(1, 1);
study.tell(first, 1);

second = study.ask();
second.report(2, 1);
verifyTrue(testCase, second.shouldPrune());
second.prune();
verifyEqual(testCase, second.State, "PRUNED");

study2 = radia.optuna.Study(Sampler=radia.optuna.RandomSampler(3));
results = study2.optimize(@localObjective, 3);
verifyEqual(testCase, height(results), 3);
verifyEqual(testCase, sum(results.State == "COMPLETE"), 3);
end

function testSimulinkRunnerContract(testCase)
verifyError(testCase, @() radia.optuna.SimulinkRunner("model", ScoreFcn=[]), ...
    "radia:optuna:SimulinkScore");
runner = radia.optuna.SimulinkRunner("model", ...
    ScoreFcn=@(simOut, trial) 0, StopTime="1.0");
verifyTrue(testCase, isa(runner, "radia.optuna.SimulinkRunner"));
end

function testSimulinkRunnerOptimizesStudy(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end

modelName = "radia_optuna_runner_test";
cleanup = onCleanup(@() closeIfLoaded(modelName));
plant = radia.simulink.makeIHPlant( ...
    HeatCapacity_J_per_K=10, ThermalConductance_W_per_K=2, SampleTime_s=0.1);
radia.simulink.buildIHControlModel(modelName, plant, ...
    StopTime_s=1, Save=false, Open=false);
time_s = reshape(0:0.1:1, [], 1);
inputData = [time_s, 100 * ones(size(time_s)), 293.15 * ones(size(time_s))];
runner = radia.optuna.SimulinkRunner(modelName, ...
    ConfigureFcn=@(simInput, trial) simInput.setExternalInput(inputData), ...
    ScoreFcn=@(simOut, trial) simOut.get("yout").getElement(1).Values.Data(end));
study = radia.optuna.create_study( ...
    study_name="simulink-runner", direction="minimize", ...
    sampler=radia.optuna.TPESampler(Seed=3, NStartupTrials=1));
results = runner.optimize(study, 2);
verifyEqual(testCase, height(results), 2);
verifyTrue(testCase, all(results.State == "COMPLETE"));
verifyTrue(testCase, all(isfinite(results.Value)));
clear cleanup
closeIfLoaded(modelName);
end

function testTPEAndCmaEsSamplers(testCase)
tpe = radia.optuna.create_study( ...
    study_name="tpe-study", direction="minimize", ...
    sampler=radia.optuna.TPESampler(Seed=4, NStartupTrials=2));
tpe.optimize(@(trial) tpeObjective(trial), 8);
verifyEqual(testCase, height(tpe.get_trials()), 8);
verifyEqual(testCase, tpe.best_trial().State, "COMPLETE");

cma = radia.optuna.create_study( ...
    study_name="cma-study", direction="minimize", ...
    sampler=radia.optuna.CmaEsSampler(Seed=5, NStartupTrials=1));
cma.optimize(@(trial) cmaObjective(trial), 8);
verifyEqual(testCase, height(cma.get_trials()), 8);
verifyTrue(testCase, all(cma.get_trials().State == "COMPLETE"));
end

function testMultiObjectiveParetoFrontAndPersistence(testCase)
path=string(tempname("C:\temp"))+".mat"; cleanup=onCleanup(@()deleteIfPresent(path));
study=radia.optuna.create_study(directions=["minimize","minimize"], ...
    storage=path, sampler=radia.optuna.RandomSampler(2));
values=[0 3;1 2;2 1;3 3];
for k=1:height(values)
    trial=study.ask(); trial.suggestFloat("x",0,1); study.tell(trial,values(k,:));
end
front=study.paretoFront();
verifyEqual(testCase,sort(front.TrialNumber),[0;1;2]);
verifyEqual(testCase,height(study.ObjectiveTable),8);
verifyError(testCase,@()study.bestTrial(),"radia:optuna:MultiObjectiveBest");
reloaded=radia.optuna.Study(StoragePath=path);
verifyEqual(testCase,height(reloaded.paretoFront()),3);
clear cleanup; deleteIfPresent(path);
end

function testSimulinkRunnerParallelTrials(testCase)
if isempty(ver("parallel")), testCase.assumeFail("Parallel Computing Toolbox is unavailable."); end
name="radia_optuna_parallel_gate"; file="C:\temp\radia_optuna_parallel_gate.slx";
cleanup=onCleanup(@()closeAndDelete(name,file));
plant=radia.simulink.makeIHPlant(HeatCapacity_J_per_K=10, ...
    ThermalConductance_W_per_K=2,SampleTime_s=0.1);
radia.simulink.buildIHControlModel(name,plant,StopTime_s=0.2,Save=false,Open=false);
save_system(name,file); close_system(name,0);
time=(0:0.1:0.2)'; input=[time,100*ones(size(time)),293.15*ones(size(time))];
runner=radia.optuna.SimulinkRunner(file, ...
    ConfigureFcn=@(simInput,trial)simInput.setExternalInput(input), ...
    ScoreFcn=@(simOut,trial)simOut.get("yout").getElement(1).Values.Data(end));
study=radia.optuna.createStudy(AutoSave=false);
result=runner.optimizeParallel(study,2,ShowProgress=false);
verifyEqual(testCase,result.State,["COMPLETE";"COMPLETE"]);
clear cleanup; closeAndDelete(name,file);
end

function testMOTPEAndNSGAIIMultiObjectiveSamplers(testCase)
motpe=radia.optuna.createStudy(directions=["minimize","minimize"], ...
    sampler=radia.optuna.MOTPESampler(Seed=8,NStartupTrials=4),AutoSave=false);
motpe.optimize(@multiObjective,20);
verifyEqual(testCase,sum(motpe.TrialTable.State=="COMPLETE"),20);
verifyGreaterThan(testCase,height(motpe.paretoFront()),1);

nsga=radia.optuna.createStudy(directions=["minimize","minimize"], ...
    sampler=radia.optuna.NSGAIISampler(Seed=9,PopulationSize=4),AutoSave=false);
nsga.optimize(@multiObjective,20);
verifyEqual(testCase,sum(nsga.TrialTable.State=="COMPLETE"),20);
verifyGreaterThan(testCase,height(nsga.paretoFront()),1);
verifyTrue(testCase,all(nsga.ParamTable.ValueNumeric>=-1 & nsga.ParamTable.ValueNumeric<=3));
end

function testParetoRankAndCrowdingContract(testCase)
values=[0 3;1 2;2 1;3 3;2 4];
[rank,crowding]=radia.optuna.internal.ParetoSupport.rankAndCrowding( ...
    values,["minimize","minimize"]);
verifyEqual(testCase,rank(1:3),ones(3,1));
verifyGreaterThan(testCase,rank(4),1);
verifyTrue(testCase,isinf(crowding(1)));
end

function value = localObjective(trial)
x = trial.suggestFloat("x", -1, 1);
trial.report(x^2, 0);
value = x^2;
end

function value = tpeObjective(trial)
x = trial.suggest_float("x", -2, 2);
mode = trial.suggest_categorical("mode", {"a", "b"});
value = x^2 + double(string(mode) == "b");
end

function value = cmaObjective(trial)
x = trial.suggest_float("x", 0.05, 2, Log=true);
y = trial.suggest_float("y", -1, 3);
value = (log(x) - log(0.4))^2 + (y - 1.2)^2;
end

function values=multiObjective(trial)
x=trial.suggestFloat("x",-1,3);
values=[x^2,(x-2)^2];
end

function closeIfLoaded(modelName)
if bdIsLoaded(modelName)
    close_system(modelName, 0);
end
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end

function closeAndDelete(name,file)
if bdIsLoaded(name), close_system(name,0); end
deleteIfPresent(file);
end
