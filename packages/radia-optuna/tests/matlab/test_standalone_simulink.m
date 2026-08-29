function result = test_standalone_simulink(installedMatlabPath)
%TEST_STANDALONE_SIMULINK Exercise Simulink using only an installed wheel.
%
% The caller must run restoredefaultpath before adding installedMatlabPath and
% this test folder.  This test deliberately builds a scratch model rather
% than relying on a Radia model or on the monorepo MATLAB tree.

arguments
    installedMatlabPath (1,1) string
end

studyPath = string(which("radia.optuna.Study"));
runnerPath = string(which("radia.optuna.SimulinkRunner"));
sessionPath = string(which("radia.optuna.OptimizationSession"));
teachingBuilderPath = string(which( ...
    "radia.simulink.buildOptunaTeachingModel"));
mexPath = string(which("optuna_mex"));
assertPathComesFromWheel(studyPath, installedMatlabPath, "Study");
assertPathComesFromWheel(runnerPath, installedMatlabPath, "SimulinkRunner");
assertPathComesFromWheel(sessionPath, installedMatlabPath, ...
    "OptimizationSession");
assertPathComesFromWheel(teachingBuilderPath, installedMatlabPath, ...
    "buildOptunaTeachingModel");
assertPathComesFromWheel(mexPath, installedMatlabPath, "optuna_mex");
assert(~isempty(ver("simulink")), ...
    "radia:optuna:StandaloneSimulinkUnavailable", ...
    "Simulink is required for the installed-wheel acceptance test.");

scratchDirectory = string(tempname("C:\temp"));
mkdir(scratchDirectory);
[~, token] = fileparts(scratchDirectory);
modelName = "radia_optuna_wheel_" + string(token);
modelFile = fullfile(scratchDirectory, modelName + ".slx");
blockModelName = modelName + "_block";
blockModelFile = fullfile(scratchDirectory, blockModelName + ".slx");
teachingModelName = modelName + "_teaching";
teachingModelFile = fullfile(scratchDirectory, ...
    teachingModelName + ".slx");
cleanup = onCleanup(@()cleanupFixture( ...
    [modelName, blockModelName, teachingModelName], scratchDirectory));

buildFixture(modelName);
save_system(modelName, modelFile);
close_system(modelName, 0);

time = (0:0.1:0.2)';
inputData = [time, 10 * ones(size(time))];
runner = radia.optuna.SimulinkRunner(modelFile, ...
    ConfigureFcn=@(simInput, trial)configureInput( ...
        simInput, trial, inputData), ...
    ScoreFcn=@scoreOutput, ...
    ConstraintFcn=@constraintValue, ...
    ValidationFcn=@validationRecord, ...
    ResultFcn=@artifactRecord, ...
    UseFastRestart=false, ...
    Context=struct("distribution", "radia-optuna", ...
        "installation", "wheel", "simulink", "standalone"));
study = radia.optuna.createStudy(AutoSave=false, ...
    Sampler=radia.optuna.RandomSampler(41));
trials = runner.optimize(study, 4, ContinueOnError=false);

assert(height(trials) == 4);
assert(all(trials.State == "COMPLETE"));
assert(all(isfinite(trials.Value)));
gainRows = study.ParamTable.Name == "gain";
assert(sum(gainRows) == 4);
assert(numel(unique(study.ParamTable.ValueNumeric(gainRows))) > 1);
assert(height(study.ConstraintTable) == 4);
assert(any(study.UserAttrTable.Name == "cae_execution"));
assert(any(study.UserAttrTable.Name == "cae_validation"));
assert(any(study.UserAttrTable.Name == "cae_artifacts"));
executionRow = find(study.UserAttrTable.Name == "cae_execution", 1);
execution = jsondecode(study.UserAttrTable.ValueJSON(executionRow));
assert(string(execution.schema) == "radia.optuna.cae-trial.v1");
assert(string(execution.status) == "complete");
assert(string(execution.context.distribution) == "radia-optuna");
assert(strlength(string(execution.model.sha256)) == 64);

failureRunner = radia.optuna.SimulinkRunner(modelFile, ...
    ConfigureFcn=@(simInput, trial)configureInput( ...
        simInput, trial, inputData), ...
    ScoreFcn=@standaloneFailureScore, ...
    UseFastRestart=false);
failureStudy = radia.optuna.createStudy(AutoSave=false, ...
    Sampler=radia.optuna.RandomSampler(43));
failureTrials = failureRunner.optimize(failureStudy, 1, ContinueOnError=true);
assert(height(failureTrials) == 1);
assert(failureTrials.State == "FAIL");
failureRow = find(failureStudy.UserAttrTable.Name == "cae_failure", 1);
assert(~isempty(failureRow));
failure = jsondecode(failureStudy.UserAttrTable.ValueJSON(failureRow));
assert(string(failure.identifier) == "radia:optuna:StandaloneFailure");
assert(string(failure.stage) == "postprocess");

blockStudyFile = fullfile(scratchDirectory, "block-study.mat");
new_system(blockModelName);
blockWorkspace = get_param(blockModelName, "ModelWorkspace");
assignin(blockWorkspace, "x", NaN);
optimizationBlock = radia.simulink.buildOptunaStudyBlock(blockModelName, ...
    BlockName="Optuna Study", ...
    ObjectiveFcn="radia_optuna_standalone_objective", ...
    NumTrials=4, StoragePath=blockStudyFile, ...
    SampleTime_s=0.1, Sampler="tpe", Save=false);
monitorBlock = radia.simulink.addOptunaMonitor(blockModelName, ...
    Position=[480 300 710 440]);
assert(string(get_param(monitorBlock + "/Optimization History", ...
    "BlockType")) == "Scope");
monitorContract = get_param(monitorBlock, "UserData");
assert(~monitorContract.browser_required);
assert(string(monitorContract.visualization) == "simulink-scope-xy");

add_block("simulink/Sources/Constant", blockModelName + "/Start", ...
    Value="1", Position=[30 100 80 130]);
add_block("simulink/Sources/Constant", blockModelName + "/Cancel", ...
    Value="0", Position=[30 155 80 185]);
add_block("simulink/Sinks/To Workspace", blockModelName + "/Best", ...
    VariableName="wheel_best", SaveFormat="Array", ...
    Position=[430 75 530 105]);
add_block("simulink/Sinks/To Workspace", blockModelName + "/Status", ...
    VariableName="wheel_status", SaveFormat="Array", ...
    Position=[430 125 530 155]);
add_block("simulink/Sinks/To Workspace", blockModelName + "/Attempted", ...
    VariableName="wheel_attempted", SaveFormat="Array", ...
    Position=[430 175 530 205]);
add_line(blockModelName, "Start/1", "Optuna Study/1");
add_line(blockModelName, "Cancel/1", "Optuna Study/2");
add_line(blockModelName, "Optuna Study/1", "Best/1");
add_line(blockModelName, "Optuna Study/2", "Status/1");
add_line(blockModelName, "Optuna Study/3", "Attempted/1");
set_param(blockModelName, "SolverType", "Fixed-step", ...
    "Solver", "FixedStepDiscrete", "FixedStep", "0.1", ...
    "StopTime", "0.4");
save_system(blockModelName, blockModelFile);
blockOutput = sim(blockModelName, ReturnWorkspaceOutputs="on");
blockBest = blockOutput.get("wheel_best");
blockStatus = blockOutput.get("wheel_status");
blockAttempted = blockOutput.get("wheel_attempted");
assert(isfinite(blockBest(end)));
assert(blockStatus(end) == 1);
assert(blockAttempted(end) == 4);
assert(isfile(blockStudyFile));
blockReview = radia.simulink.reviewOptunaStudy(optimizationBlock);
blockStudy = radia.optuna.loadStudy(storage=blockStudyFile);
assert(height(blockStudy.TrialTable) == 4);
assert(height(blockReview.trials) == 4);
assert(all(blockStudy.TrialTable.State == "COMPLETE"));
tableNames = ["TrialTable", "ObjectiveTable", "ParamTable", ...
    "IntermediateTable", "UserAttrTable", "ConstraintTable", ...
    "SamplerStateTable"];
for tableName = tableNames
    assert(isprop(blockStudy, tableName), ...
        "radia:optuna:StandaloneResumeTable", ...
        "Restored Study is missing %s.", tableName);
end
tableResume = struct( ...
    "schema", "radia-optuna.table-resume-test.v1", ...
    "ok", true, ...
    "source_trial_count", 4, ...
    "restored_trial_count", height(blockStudy.TrialTable), ...
    "table_names", tableNames);
blockContract = get_param(optimizationBlock, "UserData");
assert(~blockContract.python_per_trial);
assert(string(blockContract.interface) == "compact-student");
portsBefore = get_param(optimizationBlock, "PortConnectivity");
set_param(optimizationBlock, "num_trials", "6");
set_param(blockModelName, "StopTime", "0.6");
resumedOutput = sim(blockModelName, ReturnWorkspaceOutputs="on");
resumedReview = radia.simulink.reviewOptunaStudy(optimizationBlock);
resumedAttempted = resumedOutput.get("wheel_attempted");
assert(resumedAttempted(end) == 6);
assert(height(resumedReview.trials) == 6);
applied = radia.simulink.applyOptunaTrial(optimizationBlock, -1);
assert(blockWorkspace.getVariable("x") == applied.x);
portsAfter = get_param(optimizationBlock, "PortConnectivity");
assert(isequaln(portsAfter, portsBefore));

sessionStudyFile = fullfile(scratchDirectory, "session-study.mat");
parameter = radia.optuna.OptimizationParameter("x", ...
    Minimum=-1, Maximum=1);
options = radia.optuna.optimoptions(MaxTrials=4, Sampler="random", ...
    Seed=47, Display="off", StoragePath=sessionStudyFile);
session = radia.optuna.OptimizationSession( ...
    @radia.optuna.teachingQuadraticObjective, parameter, options);
session.start();
session.runNext();
session.pause();
checkpointPath = session.save();
assert(isfile(checkpointPath));
restoredSession = radia.optuna.OptimizationSession.load( ...
    checkpointPath, ...
    ObjectiveFcn=@radia.optuna.teachingQuadraticObjective);
restoredSession.resume();
restoredSession.run();
sessionSnapshot = restoredSession.snapshot();
assert(restoredSession.State == "completed");
assert(sessionSnapshot.complete == 4);
best = restoredSession.Study.bestSolution();
restoredSession.selectTrial(best.trial_number);
assert(isfinite(restoredSession.selectedParameters().x));

teachingStudyFile = fullfile(scratchDirectory, "teaching-study.mat");
createdTeachingModel = radia.simulink.buildOptunaTeachingModel( ...
    OutputPath=teachingModelFile, Exercise="quadratic", ...
    StoragePath=teachingStudyFile, Save=true);
assert(string(createdTeachingModel) == string(teachingModelFile));
load_system(teachingModelFile);
teachingOutput = sim(teachingModelName, ReturnWorkspaceOutputs="on");
teachingAttempted = teachingOutput.get("teaching_attempted");
teachingStatus = teachingOutput.get("teaching_status");
assert(teachingAttempted(end) == 12);
assert(teachingStatus(end) == 1);
assert(isfile(teachingStudyFile));

result = struct( ...
    "schema", "radia-optuna.standalone-simulink-test.v1", ...
    "ok", true, ...
    "study_path", studyPath, ...
    "runner_path", runnerPath, ...
    "session_path", sessionPath, ...
    "teaching_builder_path", teachingBuilderPath, ...
    "mex_path", mexPath, ...
    "complete_trials", sum(trials.State == "COMPLETE"), ...
    "failed_trials", sum(failureTrials.State == "FAIL"), ...
    "simulink_block_trials", height(blockStudy.TrialTable), ...
    "simulink_block_resumed_trials", height(resumedReview.trials), ...
    "simulink_block_best", blockBest(end), ...
    "simulink_applied_x", applied.x, ...
    "simulink_topology_stable", isequaln(portsAfter, portsBefore), ...
    "simulink_monitor", string(monitorContract.visualization), ...
    "session_trials", sessionSnapshot.complete, ...
    "teaching_model_trials", teachingAttempted(end), ...
    "table_resume", tableResume, ...
    "distinct_parameter_values", ...
        numel(unique(study.ParamTable.ValueNumeric(gainRows))));
disp(jsonencode(result));

clear cleanup
cleanupFixture([modelName, blockModelName, teachingModelName], ...
    scratchDirectory);
end

function assertPathComesFromWheel(actualPath, installedMatlabPath, label)
assert(strlength(actualPath) > 0, ...
    "radia:optuna:StandaloneMissingEntry", ...
    "%s was not found on the MATLAB path.", label);
actual = lower(replace(actualPath, "/", "\"));
expected = lower(replace(installedMatlabPath, "/", "\"));
if ~endsWith(expected, "\")
    expected = expected + "\";
end
assert(startsWith(actual, expected), ...
    "radia:optuna:StandalonePathLeak", ...
    "%s resolved outside the installed wheel: %s", label, actualPath);
end

function buildFixture(modelName)
new_system(modelName);
add_block("simulink/Ports & Subsystems/In1", modelName + "/Input", ...
    Position=[30 45 60 65]);
add_block("simulink/Math Operations/Gain", modelName + "/Trial Gain", ...
    Gain="standalone_gain", Position=[110 35 190 75]);
add_block("simulink/Ports & Subsystems/Out1", modelName + "/Response", ...
    Position=[250 45 280 65]);
add_line(modelName, "Input/1", "Trial Gain/1");
add_line(modelName, "Trial Gain/1", "Response/1");
set_param(modelName, "SolverType", "Fixed-step", ...
    "Solver", "FixedStepDiscrete", "FixedStep", "0.1", ...
    "StopTime", "0.2", "SaveOutput", "on", ...
    "OutputSaveName", "yout", "SaveFormat", "Dataset");
end

function simInput = configureInput(simInput, trial, inputData)
gain = trial.suggestFloat("gain", 0.5, 1.5);
simInput = simInput.setVariable("standalone_gain", gain);
simInput = simInput.setExternalInput(inputData);
end

function value = scoreOutput(simOut, trial) %#ok<INUSD>
response = simOut.get("yout").getElement(1).Values.Data(end);
value = (response - 8)^2;
end

function values = constraintValue(simOut, trial) %#ok<INUSD>
values = trial.Params.gain - 1.25;
end

function record = validationRecord(simOut, trial) %#ok<INUSD>
record = struct("reference", "analytic-gain", "passed", true);
end

function record = artifactRecord(simOut, trial) %#ok<INUSD>
record = struct("role", "scratch-simulation-output", ...
    "trial_number", trial.Number);
end

function value = standaloneFailureScore(simOut, trial) %#ok<INUSD>
value = NaN; %#ok<NASGU>
error("radia:optuna:StandaloneFailure", ...
    "Intentional installed-wheel Simulink failure-path test.");
end

function cleanupFixture(modelNames, scratchDirectory)
for modelName = reshape(string(modelNames), 1, [])
    if bdIsLoaded(modelName)
        close_system(modelName, 0);
    end
end
scratch = string(java.io.File(char(scratchDirectory)).getCanonicalPath());
tempRoot = string(java.io.File("C:\temp").getCanonicalPath()) + filesep;
assert(startsWith(lower(scratch + filesep), lower(tempRoot)), ...
    "radia:optuna:UnsafeScratchCleanup", ...
    "Refusing to remove a scratch directory outside C:\temp: %s", scratch);
if isfolder(scratch)
    rmdir(scratch, "s");
end
end
