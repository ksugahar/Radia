function tests = test_optuna_table
tests = functiontests(localfunctions);
end

function testCollidingParameterNamesKeepMatchingFrozenKeys(testCase)
study=radia.optuna.Study(Sampler=radia.optuna.RandomSampler(3), ...
    AutoSave=false);
trial=study.ask();
first=trial.suggest_float("x-1",0,1);
second=trial.suggest_float("x.1",2,3);
study.tell(trial,first+second);
frozen=study.get_trials();
parameterKeys=sort(string(fieldnames(frozen(1).Params)));
distributionKeys=sort(string(fieldnames(frozen(1).Distributions)));
verifyEqual(testCase,numel(parameterKeys),2);
verifyEqual(testCase,distributionKeys,parameterKeys);
target=radia.optuna.Study(Sampler=radia.optuna.RandomSampler(3), ...
    AutoSave=false);
target.add_trial(frozen(1));
verifyEqual(testCase,height(target.TrialTable),1);
verifyEqual(testCase,target.TrialTable.State(1),"COMPLETE");
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabDirectory = fullfile(repositoryRoot, "matlab");
entries = string(strsplit(path,pathsep));
testCase.TestData.RemoveMatlabDirectory = ...
    ~any(strcmpi(entries,string(matlabDirectory)));
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

function testTrialCompatibilityMetadata(testCase)
trial = radia.optuna.Study(AutoSave=false).ask();
x = trial.suggest_float("positive", 1, 100, Log=true);
trial.set_user_attr("role", "compatibility-test");
trial.set_system_attr("source", "matlab");

verifyTrue(testCase, isfinite(x));
verifyTrue(testCase, isfield(trial.Distributions, "positive"));
verifyEqual(testCase, trial.Distributions.positive.name, "FloatDistribution");
verifyEqual(testCase, trial.UserAttrs.role, "compatibility-test");
verifyEqual(testCase, trial.SystemAttrs.source, "matlab");
end

function testStudyUserAttributesPersist(testCase)
path = string(tempname("C:\temp")) + ".mat";
cleanup = onCleanup(@() deleteIfPresent(path));
study = radia.optuna.Study(StoragePath=path, AutoSave=true);
study.set_user_attr("owner", "radia");
reloaded = radia.optuna.Study(StoragePath=path);
verifyEqual(testCase, reloaded.UserAttrs.owner, "radia");
clear cleanup;
deleteIfPresent(path);
end

function testDistributionCodecReadsVersionedAndLegacyMetadata(testCase)
spec = radia.optuna.internal.DistributionCodec.float(0.1, 10, true, NaN);
roundTrip = radia.optuna.internal.DistributionCodec.decode( ...
    "float", radia.optuna.internal.DistributionCodec.encode(spec));
legacyFloat = radia.optuna.internal.DistributionCodec.decode( ...
    "float", "[-2,3]");
legacyInteger = radia.optuna.internal.DistributionCodec.decode( ...
    "integer", "[1,9] log");
legacyCategorical = radia.optuna.internal.DistributionCodec.decode( ...
    "categorical", string(jsonencode({"a","b"})));
singleCategorical = ...
    radia.optuna.internal.DistributionCodec.categorical("steel");
singleRoundTrip = radia.optuna.internal.DistributionCodec.decode( ...
    "categorical", ...
    radia.optuna.internal.DistributionCodec.encode(singleCategorical));

verifyTrue(testCase, ...
    radia.optuna.internal.DistributionCodec.equivalent(spec, roundTrip));
verifyEqual(testCase, [legacyFloat.low,legacyFloat.high], [-2,3]);
verifyFalse(testCase, legacyFloat.log);
verifyEqual(testCase, [legacyInteger.low,legacyInteger.high], [1,9]);
verifyTrue(testCase, legacyInteger.log);
verifyEqual(testCase, ...
    radia.optuna.internal.DistributionCodec.choiceTokens( ...
    legacyCategorical.choices), ...
    [string(jsonencode("a"));string(jsonencode("b"))]);
verifyTrue(testCase, ...
    radia.optuna.internal.DistributionCodec.isSingle(singleRoundTrip));
verifyTrue(testCase, ...
    radia.optuna.internal.DistributionCodec.equivalent( ...
    singleCategorical, singleRoundTrip));
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
buildRunnerFixture(modelName, 1.0);
time_s = reshape(0:0.1:1, [], 1);
inputData = [time_s, 10 * ones(size(time_s))];
runner = radia.optuna.SimulinkRunner(modelName, ...
    ConfigureFcn=@(simInput, trial) configureRunnerInput( ...
        simInput, trial, inputData), ...
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

function testSimulinkRunnerRecordsCAETrialContract(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end
modelName = "radia_optuna_cae_contract";
modelFile = "C:\temp\radia_optuna_cae_contract.slx";
cleanup = onCleanup(@() closeAndDelete(modelName, modelFile));
buildRunnerFixture(modelName, 0.2);
save_system(modelName, modelFile);
close_system(modelName, 0);
time_s = reshape(0:0.1:0.2, [], 1);
inputData = [time_s, 10 * ones(size(time_s))];
runner = radia.optuna.SimulinkRunner(modelFile, ...
    ConfigureFcn=@(simInput, trial) configureRunnerInput( ...
        simInput, trial, inputData), ...
    ScoreFcn=@runnerScore, ...
    ConstraintFcn=@runnerConstraint, ...
    ValidationFcn=@runnerValidation, ...
    ResultFcn=@runnerArtifacts, ...
    Context=struct( ...
        "geometry_id", "fixture-v1", ...
        "mesh_id", "analytic-no-mesh", ...
        "material_id", "unit-gain", ...
        "excitation_id", "constant-10"));
study = radia.optuna.createStudy(AutoSave=false, ...
    sampler=radia.optuna.RandomSampler(71));
result = runner.optimize(study, 2);

verifyEqual(testCase, result.State, ["COMPLETE";"COMPLETE"]);
verifyEqual(testCase, height(study.ConstraintTable), 2);
verifyTrue(testCase, all(isfinite(study.ConstraintTable.Value)));
verifyTrue(testCase, any(study.UserAttrTable.Name == "cae_execution"));
verifyTrue(testCase, any(study.UserAttrTable.Name == "cae_validation"));
verifyTrue(testCase, any(study.UserAttrTable.Name == "cae_artifacts"));
row = find(study.UserAttrTable.Name == "cae_execution", 1);
record = jsondecode(study.UserAttrTable.ValueJSON(row));
verifyEqual(testCase, string(record.schema), "radia.optuna.cae-trial.v1");
verifyEqual(testCase, string(record.status), "complete");
verifyEqual(testCase, string(record.context.geometry_id), "fixture-v1");
verifyEqual(testCase, strlength(string(record.model.sha256)), 64);
verifyGreaterThanOrEqual(testCase, record.timing_s.configuration, 0);
verifyGreaterThanOrEqual(testCase, record.timing_s.simulation, 0);
verifyGreaterThanOrEqual(testCase, record.timing_s.postprocess, 0);
verifyGreaterThanOrEqual(testCase, record.timing_s.total, ...
    record.timing_s.simulation);
clear cleanup
closeAndDelete(modelName, modelFile);
end

function testSimulinkRunnerClassifiesCAEFailureAndContinues(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is not installed on this MATLAB runtime.");
    return
end
modelName = "radia_optuna_cae_failure";
cleanup = onCleanup(@() closeIfLoaded(modelName));
buildRunnerFixture(modelName, 0.1);
time_s = reshape(0:0.1:0.1, [], 1);
inputData = [time_s, 10 * ones(size(time_s))];
runner = radia.optuna.SimulinkRunner(modelName, ...
    ConfigureFcn=@(simInput, trial) configureRunnerInput( ...
        simInput, trial, inputData), ...
    ScoreFcn=@runnerMeshFailure);
study = radia.optuna.createStudy(AutoSave=false);
result = runner.optimize(study, 2, ContinueOnError=true);

verifyEqual(testCase, result.State, ["FAIL";"FAIL"]);
rows = study.UserAttrTable.Name == "cae_failure";
verifyEqual(testCase, sum(rows), 2);
failure = jsondecode(study.UserAttrTable.ValueJSON(find(rows, 1)));
verifyEqual(testCase, string(failure.class), "mesh_invalid");
verifyEqual(testCase, string(failure.stage), "postprocess");
verifyFalse(testCase, failure.retryable);
verifyEqual(testCase, string(failure.identifier), "radia:test:MeshQuality");
clear cleanup
closeIfLoaded(modelName);
end

function testCmaPartialGenerationSurvivesStudyReload(testCase)
path = string(tempname("C:\temp")) + ".mat";
cleanup = onCleanup(@() deleteIfPresent(path));
continuous = radia.optuna.Study(Sampler=radia.optuna.CmaEsSampler( ...
    Seed=53, NStartupTrials=1, PopulationSize=4), AutoSave=false);
persisted = radia.optuna.Study(StoragePath=path, ...
    Sampler=radia.optuna.CmaEsSampler( ...
    Seed=53, NStartupTrials=1, PopulationSize=4), AutoSave=true);

for index = 1:3
    continuous.optimize(@cmaObjective, 1);
    persisted.optimize(@cmaObjective, 1);
end
reloaded = radia.optuna.Study(StoragePath=path, AutoSave=false, ...
    Sampler=radia.optuna.CmaEsSampler( ...
    Seed=53, NStartupTrials=1, PopulationSize=4));

expected = continuous.ask();
actual = reloaded.ask();
expectedX = expected.suggest_float("x", 0.05, 2, Log=true);
expectedY = expected.suggest_float("y", -1, 3);
actualX = actual.suggest_float("x", 0.05, 2, Log=true);
actualY = actual.suggest_float("y", -1, 3);
verifyEqual(testCase, [actualX,actualY], [expectedX,expectedY], AbsTol=0);
row=reloaded.SamplerStateTable.Sampler=="cmaes";
verifyEqual(testCase,sum(row),1);
verifyEqual(testCase, ...
    reloaded.SamplerStateTable.State{row}.population_trial_numbers, [1;2]);
clear cleanup
deleteIfPresent(path);
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
buildRunnerFixture(name,0.2);
save_system(name,file); close_system(name,0);
time=(0:0.1:0.2)'; input=[time,10*ones(size(time))];
runner=radia.optuna.SimulinkRunner(file, ...
    ConfigureFcn=@(simInput,trial)configureRunnerInput(simInput,trial,input), ...
    ScoreFcn=@(simOut,trial)simOut.get("yout").getElement(1).Values.Data(end), ...
    UseFastRestart=true);
study=radia.optuna.createStudy(AutoSave=false);
result=runner.optimizeParallel(study,2,ShowProgress=false);
verifyEqual(testCase,result.State,["COMPLETE";"COMPLETE"]);
clear cleanup; closeAndDelete(name,file);
end

function testSimulinkRunnerSteadyStateFailuresLeaveNoRunningTrials(testCase)
if isempty(ver("parallel"))
    testCase.assumeFail("Parallel Computing Toolbox is unavailable.");
end
pool = gcp();
if isa(pool, "parallel.ThreadPool")
    testCase.assumeFail("Steady-state Simulink execution requires a process pool.");
end
name = "radia_optuna_parallel_failure_gate";
file = "C:\temp\radia_optuna_parallel_failure_gate.slx";
cleanup = onCleanup(@()closeAndDelete(name, file));
buildRunnerFixture(name, 0.1);
save_system(name, file);
close_system(name, 0);
time = (0:0.1:0.1)';
input = [time, 10 * ones(size(time))];
runner = radia.optuna.SimulinkRunner(file, ...
    ConfigureFcn=@(simInput, trial)configureRunnerInput( ...
        simInput, trial, input), ...
    ScoreFcn=@runnerMeshFailure);
study = radia.optuna.createStudy(AutoSave=false);
result = runner.optimizeParallel(study, 3, BatchSize=2, ...
    ShowProgress=false, ExecutionMode="steady_state");

verifyEqual(testCase, result.State, repmat("FAIL", 3, 1));
verifyEqual(testCase, sum(result.State == "RUNNING"), 0);
contract = study.UserAttrs.cae_execution_contract;
verifyEqual(testCase, contract.execution_mode, ...
    "parfeval-steady-state");
verifyEqual(testCase, contract.parallel_decision.selected, ...
    "steady_state");
clear cleanup
closeAndDelete(name, file);
end

function testConstraintTablePersists(testCase)
path=string(tempname("C:\temp"))+".mat";
cleanup=onCleanup(@()deleteIfPresent(path));
sampler=radia.optuna.TPESampler(Seed=31,NStartupTrials=1, ...
    ConstraintsFcn=@(trial)trial.Params.x-0.25);
study=radia.optuna.Study(StoragePath=path,Sampler=sampler,AutoSave=true);
study.optimize(@constraintObjective,3);
verifyEqual(testCase,height(study.ConstraintTable),3);
reloaded=radia.optuna.Study(StoragePath=path);
verifyEqual(testCase,reloaded.ConstraintTable,study.ConstraintTable);
clear cleanup; deleteIfPresent(path);
end

function value = localObjective(trial)
x = trial.suggestFloat("x", -1, 1);
trial.report(x^2, 0);
value = x^2;
end

function value=constraintObjective(trial)
x=trial.suggestFloat("x",-1,1);
value=x^2;
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

function buildRunnerFixture(modelName, stopTime)
new_system(modelName);
add_block("simulink/Ports & Subsystems/In1",modelName+"/Input", ...
    Position=[30 45 60 65]);
add_block("simulink/Math Operations/Gain",modelName+"/Trial Gain", ...
    Gain="runner_gain",Position=[110 35 190 75]);
add_block("simulink/Ports & Subsystems/Out1",modelName+"/Response", ...
    Position=[250 45 280 65]);
add_line(modelName,"Input/1","Trial Gain/1");
add_line(modelName,"Trial Gain/1","Response/1");
set_param(modelName,"SolverType","Fixed-step","Solver","FixedStepDiscrete", ...
    "FixedStep","0.1","StopTime",num2str(stopTime), ...
    "SaveOutput","on","OutputSaveName","yout","SaveFormat","Dataset");
assignin("base","runner_gain",1);
end

function simInput=configureRunnerInput(simInput,trial,inputData)
gain=trial.suggestFloat("gain",0.5,1.5);
simInput=simInput.setVariable("runner_gain",gain);
simInput=simInput.setExternalInput(inputData);
end

function value=runnerScore(simOut,trial) %#ok<INUSD>
value=simOut.get("yout").getElement(1).Values.Data(end);
end

function values=runnerConstraint(simOut,trial) %#ok<INUSD>
values=trial.Params.gain-1.1;
end

function result=runnerValidation(simOut,trial) %#ok<INUSD>
result=struct("reference","analytic-gain", ...
    "relative_error",0,"passed",true);
end

function result=runnerArtifacts(simOut,trial) %#ok<INUSD>
result=struct("role","simulation-output", ...
    "trial_number",trial.Number,"embedded",true);
end

function value=runnerMeshFailure(simOut,trial) %#ok<INUSD>
value=NaN; %#ok<NASGU>
error("radia:test:MeshQuality", ...
    "Mesh element quality is invalid for this CAE trial.");
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
