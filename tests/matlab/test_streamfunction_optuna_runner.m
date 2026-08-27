function tests = test_streamfunction_optuna_runner
tests = functiontests(localfunctions);
end

function testRunnerRecordsMultiobjectiveTrialsAndArtifacts(testCase)
root = freshRoot("runner");
cleanup = onCleanup(@() removeTree(root));
[runner, study] = makeRunner(root, 11);

results = runner.optimize(study, 3, ContinueOnError=false);

verifyEqual(testCase, height(results), 3);
verifyTrue(testCase, all(results.State == "COMPLETE"));
verifyEqual(testCase, height(study.ObjectiveTable), 6);
verifyEqual(testCase, height(study.ConstraintCountTable), 3);
verifyEqual(testCase, study.ConstraintCountTable.Count, ones(3, 1));
verifyTrue(testCase, all(isfinite(study.ObjectiveTable.Value)));
verifyGreaterThanOrEqual(testCase, height(study.paretoFront()), 1);

trialDirectories = dir(fullfile(root, "trial_*"));
verifyEqual(testCase, numel(trialDirectories), 3);
for index = 1:numel(trialDirectories)
    trialRoot = fullfile(trialDirectories(index).folder, ...
        trialDirectories(index).name);
    verifyTrue(testCase, isfile(fullfile(trialRoot, "config.json")));
    verifyTrue(testCase, isfile(fullfile(trialRoot, "optuna_trial.json")));
    config = jsondecode(fileread(fullfile(trialRoot, "config.json")));
    verifyEqual(testCase, string(config.application), "streamfunction");
    verifyTrue(testCase, isfield(config.settings, "alpha"));
    verifyTrue(testCase, isfield(config.settings, "order"));
    verifyTrue(testCase, isfield(config.settings, "confine"));
end
clear cleanup
end

function testExplicitTriggerSimulinkBlockRunsWholeStudy(testCase)
root = freshRoot("simulink");
cleanupRoot = onCleanup(@() removeTree(root));
[runner, study] = makeSingleObjectiveRunner(root, 17);
token = char(java.util.UUID.randomUUID());
modelName = "radia_sf_optuna_" + erase(string(token), "-");
new_system(modelName);
cleanupModel = onCleanup(@() closeModel(modelName));
workspace = get_param(modelName, "ModelWorkspace");
workspace.assignin("radia_streamfunction_optuna_runner", runner);
workspace.assignin("radia_streamfunction_optuna_study", study);
add_block("simulink/Sources/Constant", modelName + "/Run", ...
    Value="true", OutDataTypeStr="boolean");
add_block("simulink/User-Defined Functions/Level-2 MATLAB S-Function", ...
    modelName + "/Stream Function Optuna", ...
    FunctionName="radia_streamfunction_optuna_sfun", ...
    Parameters="'radia_streamfunction_optuna_runner'," + ...
        "'radia_streamfunction_optuna_study',2,1");
add_line(modelName, "Run/1", "Stream Function Optuna/1");
set_param(modelName, "Solver", "FixedStepDiscrete", ...
    "FixedStep", "1", "StopTime", "1");

sim(modelName);

summary = get_param(modelName + "/Stream Function Optuna", "UserData");
verifyEqual(testCase, string(summary.schema), ...
    "radia.simulink.streamfunction-optuna.v1");
verifyEqual(testCase, summary.status, 2);
verifyEqual(testCase, summary.completed_trials, 2);
verifyEqual(testCase, summary.failed_trials, 0);
verifyTrue(testCase, isfinite(summary.best_value));
verifyGreaterThanOrEqual(testCase, summary.best_trial, 0);
verifyEqual(testCase, height(study.TrialTable), 2);
clear cleanupModel cleanupRoot
end

function [runner, study] = makeSingleObjectiveRunner(root, seed)
settings = struct( ...
    "method", "Design", ...
    "coil_vol", "coil.vol", ...
    "eval_vol", "eval.vol", ...
    "target_cf", "1", ...
    "order", 2, ...
    "confine", "abe", ...
    "aca_eps", 1.0e-10);
space = struct();
space.alpha = struct( ...
    "Kind", "float", "Low", 1.0e-4, "High", 1.0e-2, "Log", true);
runner = radia.stream.OptunaRunner(settings, ...
    SearchSpace=space, OutputRoot=root, RunFcn=@fakeStreamFunctionRun);
study = radia.optuna.createStudy( ...
    direction="minimize", sampler=radia.optuna.RandomSampler(seed), ...
    AutoSave=false);
end

function [runner, study] = makeRunner(root, seed)
settings = struct( ...
    "method", "Design", ...
    "coil_vol", "coil.vol", ...
    "eval_vol", "eval.vol", ...
    "target_cf", "1", ...
    "aca_eps", 1.0e-10);
space = struct();
space.alpha = struct( ...
    "Kind", "float", "Low", 1.0e-4, "High", 1.0e-2, "Log", true);
space.order = struct( ...
    "Kind", "int", "Low", 1, "High", 3, "Step", 1);
space.confine = struct( ...
    "Kind", "categorical", "Choices", ["abe", "gauge"]);
runner = radia.stream.OptunaRunner(settings, ...
    SearchSpace=space, ...
    ObjectiveKeys=["homogeneity_rms", "peak_J"], ...
    ConstraintFcn=@(solver, application, trial) ...
        solver.fit_residual_rms - 0.1, ...
    OutputRoot=root, RunFcn=@fakeStreamFunctionRun);
study = radia.optuna.createStudy( ...
    directions=["minimize", "minimize"], ...
    sampler=radia.optuna.RandomSampler(seed), AutoSave=false);
end

function result = fakeStreamFunctionRun(configPath, runRoot, ~, ~)
config = jsondecode(fileread(configPath));
settings = config.settings;
if ~isfolder(runRoot)
    mkdir(runRoot);
end
solverPath = fullfile(runRoot, "solver_result.json");
resultPath = fullfile(runRoot, "result.json");
alpha = double(settings.alpha);
order = double(settings.order);
confinementPenalty = double(string(settings.confine) == "gauge") * 0.02;
solver = struct( ...
    "schema", "radia.test.streamfunction-result.v1", ...
    "homogeneity_rms", (log10(alpha) + 3.0)^2 + ...
        0.01 / order + confinementPenalty, ...
    "peak_J", 2.0 + 20.0 * alpha + 0.1 * order, ...
    "fit_residual_rms", 0.01, ...
    "factorization", "aca_plus_qr_tsvd", ...
    "aca_eps", double(settings.aca_eps));
solver.rms = solver.homogeneity_rms;
writeJSON(solverPath, solver);
application = struct( ...
    "schema", "radia.simulink.application_run.v1", ...
    "status", "passed", ...
    "solver_result", string(solverPath), ...
    "result_json", string(resultPath), ...
    "run_dir", string(runRoot), ...
    "elapsed_s", 0.001, ...
    "primary", struct( ...
        "key", "homogeneity_rms", ...
        "value", solver.homogeneity_rms));
writeJSON(resultPath, struct("radia_result", application));
result = application;
end

function writeJSON(path, value)
file = fopen(path, "w", "n", "UTF-8");
assert(file >= 0);
cleanup = onCleanup(@() fclose(file));
fprintf(file, "%s\n", jsonencode(value, PrettyPrint=true));
clear cleanup
end

function root = freshRoot(label)
root = fullfile("C:\temp", "radia_sf_optuna_" + label + "_" + ...
    erase(string(java.util.UUID.randomUUID()), "-"));
mkdir(root);
end

function removeTree(path)
if isfolder(path)
    rmdir(path, "s");
end
end

function closeModel(modelName)
if bdIsLoaded(modelName)
    close_system(modelName, 0);
end
end
