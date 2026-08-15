function tests = test_optuna_table
tests = functiontests(localfunctions);
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

function testTPEKeepsBadSetForGammaOne(testCase)
sampler = radia.optuna.TPESampler(Seed=17, NStartupTrials=1, Gamma=1);
study = radia.optuna.Study(Name="tpe-boundary", Sampler=sampler, ...
    AutoSave=false);
study.optimize(@tpeBoundaryObjective, 12);

verifyEqual(testCase, height(study.TrialTable), 12);
verifyTrue(testCase, all(study.TrialTable.State == "COMPLETE"));
verifyTrue(testCase, all(isfinite(study.TrialTable.Value)));
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

function testJointSearchSpaceContract(testCase)
study = radia.optuna.Study(Sampler=radia.optuna.TPESampler( ...
    Seed=21, NStartupTrials=0), AutoSave=false);
trial = study.ask();
values = trial.suggest_vector(["x","y"], [-1,-2], [1,2]);

verifySize(testCase, values, [1 2]);
verifyTrue(testCase, all(values >= [-1,-2] & values <= [1,2]));
verifyTrue(testCase, isfield(trial.Distributions, "x"));
verifyTrue(testCase, isfield(trial.Distributions, "y"));
verifyEqual(testCase, height(study.ParamTable), 2);
end

function testJointTPESamplesSharedMixture(testCase)
study = radia.optuna.Study(Sampler=radia.optuna.TPESampler( ...
    Seed=22, NStartupTrials=3), AutoSave=false);
study.optimize(@jointTPEObjective, 18);
verifyEqual(testCase, height(study.TrialTable), 18);
verifyTrue(testCase, all(study.TrialTable.State == "COMPLETE"));
verifyLessThan(testCase, study.bestValue(), 0.5);
end

function testAutomaticIntersectionMatchesExplicitJointSampling(testCase)
automatic = radia.optuna.Study(Sampler=radia.optuna.TPESampler( ...
    Seed=22, NStartupTrials=3, Multivariate=true), AutoSave=false);
explicit = radia.optuna.Study(Sampler=radia.optuna.TPESampler( ...
    Seed=22, NStartupTrials=3), AutoSave=false);

automatic.optimize(@automaticJointTPEObjective, 18);
explicit.optimize(@jointTPEObjective, 18);

verifyEqual(testCase, automatic.ParamTable.Name, explicit.ParamTable.Name);
verifyEqual(testCase, automatic.ParamTable.ValueNumeric, ...
    explicit.ParamTable.ValueNumeric, AbsTol=0);
verifyEqual(testCase, automatic.TrialTable.Value, ...
    explicit.TrialTable.Value, AbsTol=0);
next = automatic.ask();
verifyEqual(testCase, next.SystemAttrs.tpe_relative_search_space, ["x","y"]);
end

function testAutomaticIntersectionShrinksDynamicSearchSpace(testCase)
study = radia.optuna.Study(Sampler=radia.optuna.TPESampler( ...
    Seed=23, NStartupTrials=0, Multivariate=true), AutoSave=false);
first = study.ask();
first.suggestFloat("x", -1, 1);
first.suggestFloat("conditional", -2, 2);
study.tell(first, 1);

second = study.ask();
second.suggestFloat("x", -1, 1);
study.tell(second, 0.5);

third = study.ask();
verifyEqual(testCase, ...
    third.SystemAttrs.tpe_relative_search_space, "x");
verifyFalse(testCase, isfield(third.Params, "x"));
end

function testAutomaticIntersectionSupportsMixedDistributions(testCase)
sampler = radia.optuna.TPESampler( ...
    Seed=24, NStartupTrials=2, Multivariate=true);
study = radia.optuna.Study(Sampler=sampler, AutoSave=false);
choices = ["steel","copper"];
for index = 1:2
    trial = study.ask();
    x = trial.suggestFloat("x", 0.1, 2, Log=true);
    n = trial.suggest_int("n", 1, 5);
    material = trial.suggestCategorical("material", choices);
    study.tell(trial, (x-0.7)^2 + (n-3)^2 + double(material=="copper"));
end

trial = study.ask();
verifyEqual(testCase, trial.SystemAttrs.tpe_relative_search_space, ...
    ["material","n","x"]);
x = trial.suggestFloat("x", 0.1, 2, Log=true);
n = trial.suggest_int("n", 1, 5);
material = trial.suggestCategorical("material", choices);
verifyTrue(testCase, x >= 0.1 && x <= 2);
verifyTrue(testCase, n >= 1 && n <= 5 && n == floor(n));
verifyTrue(testCase, isstring(material) && isscalar(material));
verifyTrue(testCase, any(material == choices));
end

function testAutomaticIntersectionSurvivesStudyReload(testCase)
path = string(tempname("C:\temp")) + ".mat";
cleanup = onCleanup(@() deleteIfPresent(path));
study = radia.optuna.Study(StoragePath=path, AutoSave=true, ...
    Sampler=radia.optuna.TPESampler( ...
    Seed=25, NStartupTrials=1, Multivariate=true));
trial = study.ask();
trial.suggestFloat("x", -1, 1);
trial.suggestInteger("n", 1, 4);
study.tell(trial, 1);
verifyTrue(testCase, all(startsWith(study.ParamTable.Distribution, "{")));

reloaded = radia.optuna.Study(StoragePath=path, AutoSave=false, ...
    Sampler=radia.optuna.TPESampler( ...
    Seed=25, NStartupTrials=1, Multivariate=true));
next = reloaded.ask();
verifyEqual(testCase, next.SystemAttrs.tpe_relative_search_space, ["n","x"]);
clear cleanup;
deleteIfPresent(path);
end

function value = jointTPEObjective(trial)
xy = trial.suggestVector(["x","y"], [-2,-2], [2,2]);
value = (xy(2) - 0.8 * xy(1))^2 + 0.02 * sum(xy.^2);
end

function value = automaticJointTPEObjective(trial)
x = trial.suggestFloat("x", -2, 2);
y = trial.suggestFloat("y", -2, 2);
value = (y - 0.8 * x)^2 + 0.02 * (x^2 + y^2);
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

function testTPEParzenMixtureHasOptunaPrior(testCase)
estimator = radia.optuna.internal.ParzenEstimator.numerical( ...
    [0.2; 0.4; 0.8], 0, 1, PriorWeight=1);

verifyEqual(testCase, numel(estimator.weights), 4);
verifyEqual(testCase, sum(estimator.weights), 1, AbsTol=1e-14);
verifyEqual(testCase, estimator.mu(end), 0.5, AbsTol=1e-14);
verifyEqual(testCase, estimator.sigma(end), 1, AbsTol=1e-14);
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

function testParzenUsesPredeterminedObservationWeights(testCase)
numerical = radia.optuna.internal.ParzenEstimator.numerical( ...
    [0.2; 0.8], 0, 1, PriorWeight=2, ...
    ObservationWeights=[1;3]);
categorical = radia.optuna.internal.ParzenEstimator.categorical( ...
    [1;2], 3, PriorWeight=2, ObservationWeights=[1;3]);

verifyEqual(testCase, numerical.weights, [1;3;2]/6, AbsTol=1e-14);
verifyEqual(testCase, categorical.weights, [1;3;2]/6, AbsTol=1e-14);
verifyError(testCase, @() ...
    radia.optuna.internal.ParzenEstimator.numerical( ...
    [0.2;0.8], 0, 1, ObservationWeights=1), ...
    "radia:optuna:TPEWeights");
verifyError(testCase, @() ...
    radia.optuna.internal.ParzenEstimator.numerical( ...
    [], 0, 1, ObservationWeights=1), ...
    "radia:optuna:TPEWeights");
end

function testMultivariateParzenUsesOptunaBandwidth(testCase)
observations = linspace(0.05, 0.95, 20).';
dimension = 2;
estimator = radia.optuna.internal.ParzenEstimator.numerical( ...
    observations, 0, 1, MultivariateDimension=dimension);
expected = 0.2 * numel(observations)^(-1/(dimension+4));

verifyEqual(testCase, estimator.sigma(1:end-1), ...
    repmat(expected, numel(observations), 1), AbsTol=1e-14);
verifyEqual(testCase, estimator.sigma(end), 1, AbsTol=1e-14);
end

function testTPEDefaultsMatchOptuna49Lifecycle(testCase)
sampler = radia.optuna.TPESampler();
verifyFalse(testCase, sampler.ConstantLiar);
end

function value = tpeBoundaryObjective(trial)
x = trial.suggestFloat("x", -1, 1);
value = (x - 0.25)^2;
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

function testCanonicalCmaEvolutionStateUpdate(testCase)
engine = radia.optuna.internal.CMAEvolutionStrategy( ...
    [0.5,0.5], 1/6, Bounds=[0,1;0,1], ...
    PopulationSize=4, Seed=51);
points = zeros(engine.PopulationSize, engine.Dimension);
for index = 1:engine.PopulationSize
    points(index,:) = engine.ask();
end
fitness = sum((points - [0.2,0.8]).^2, 2);
engine.tell(points, fitness);

verifyEqual(testCase, engine.Generation, 1);
verifyGreaterThan(testCase, norm(engine.PSigma), 0);
verifyGreaterThan(testCase, norm(engine.PC), 0);
verifyGreaterThan(testCase, min(eig(engine.Covariance)), 0);
verifyGreaterThan(testCase, engine.Sigma, 0);

restored = radia.optuna.internal.CMAEvolutionStrategy.fromSnapshot( ...
    engine.snapshot());
verifyEqual(testCase, restored.ask(), engine.ask(), AbsTol=1e-15);
end

function testCanonicalCmaGenerationGolden(testCase)
engine = radia.optuna.internal.CMAEvolutionStrategy( ...
    [0.5,0.5], 1/6, Bounds=[0,1;0,1], ...
    PopulationSize=4, Seed=61);
points = [0.2,0.8;0.4,0.6;0.7,0.3;0.9,0.1];
fitness = [0.1;0.2;0.8;1.2];
engine.tell(points, fitness);

verifyEqual(testCase, engine.Mean, ...
    [0.2391674280134541,0.7608325719865459], AbsTol=2e-15);
verifyEqual(testCase, engine.Sigma, ...
    0.20539021051529205, AbsTol=2e-15);
verifyEqual(testCase, engine.PSigma, ...
    [-1.5252557671746065,1.5252557671746065], AbsTol=2e-15);
verifyEqual(testCase, engine.PC, ...
    [-1.759696529682537,1.759696529682537], AbsTol=2e-15);
verifyEqual(testCase, engine.Covariance, ...
    [1.367328043231014,-0.5132177541621417; ...
    -0.5132177541621417,1.367328043231014], AbsTol=2e-15);
end

function testCmaUsesAutomaticIntersectionAndFullGenerations(testCase)
study = radia.optuna.Study(Sampler=radia.optuna.CmaEsSampler( ...
    Seed=52, NStartupTrials=1, PopulationSize=4), AutoSave=false);
study.optimize(@cmaObjective, 5);

verifyEqual(testCase, height(study.SamplerStateTable), 1);
verifyEqual(testCase, study.SamplerStateTable.Sampler, "cmaes");
verifyEqual(testCase, study.SamplerStateTable.Generation, 1);
state = study.SamplerStateTable.State{1};
verifyEqual(testCase, state.engine.generation, 1);
verifySize(testCase, state.engine.covariance, [2,2]);
verifyEmpty(testCase, state.population_fitness);

next = study.ask();
verifyEqual(testCase, next.SystemAttrs.cmaes_relative_search_space, ...
    ["x","y"]);
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
verifyEqual(testCase, height(reloaded.SamplerStateTable), 1);
verifyEqual(testCase, ...
    reloaded.SamplerStateTable.State{1}.population_trial_numbers, [1;2]);
clear cleanup
deleteIfPresent(path);
end

function testCmaRejectsChangedDistributionFromGeneration(testCase)
study = radia.optuna.Study(Sampler=radia.optuna.CmaEsSampler( ...
    Seed=54, NStartupTrials=1, PopulationSize=4), AutoSave=false);
first = study.ask();
first.suggestFloat("x", -1, 1);
first.suggestFloat("y", -1, 1);
study.tell(first, 1);

changed = study.ask();
changed.suggestFloat("x", -1, 2);
changed.suggestFloat("y", -1, 1);
study.tell(changed, 0.5);
state = study.SamplerStateTable.State{1};
verifyEmpty(testCase, state.population_fitness);

next = study.ask();
verifyFalse(testCase, ...
    isfield(next.SystemAttrs, "cmaes_relative_search_space"));
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

function testMOTPEQualityBeatsRandomOnZDT1(testCase)
errors = zeros(3, 2);
for seed = 0:2
    errors(seed+1, 1) = zdt1FrontError( ...
        radia.optuna.RandomSampler(seed), 60);
    errors(seed+1, 2) = zdt1FrontError( ...
        radia.optuna.MOTPESampler(Seed=seed, NStartupTrials=10), 60);
end
verifyLessThan(testCase, median(errors(:,2)), median(errors(:,1)));
end

function value = zdt1FrontError(sampler, trialCount)
study = radia.optuna.createStudy(directions=["minimize","minimize"], ...
    Sampler=sampler, AutoSave=false);
for index = 1:trialCount
    trial = study.ask();
    x1 = trial.suggestFloat("x1", 0, 1);
    x2 = trial.suggestFloat("x2", 0, 1);
    g = 1 + 9*x2;
    study.tell(trial, [x1,g*(1-sqrt(x1/g))]);
end
front = study.paretoFront();
values = vertcat(front.Values{:});
ideal = 1 - sqrt(max(0, min(1, values(:,1))));
value = mean(max(0, values(:,2)-ideal));
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

function testConstrainedMOTPESplitPrefersFeasibleTrials(testCase)
sampler=radia.optuna.MOTPESampler(Seed=32,NStartupTrials=2,Gamma=0.5, ...
    ConstraintsFcn=@(trial)trial.Params.x-0.2);
study=radia.optuna.createStudy(directions=["minimize","minimize"], ...
    sampler=sampler,AutoSave=false);
study.optimize(@constrainedMultiObjective,12);
verifyEqual(testCase,height(study.ConstraintTable),12);
constraints=arrayfun(@(n)study.constraintsForTrial(n), ...
    study.TrialTable.TrialNumber,'UniformOutput',false);
verifyTrue(testCase,any(cellfun(@(v)all(v<=0),constraints)));
verifyGreaterThan(testCase,height(study.paretoFront()),0);
end

function testMOTPESplitSelectsFeasibleBeforeInfeasible(testCase)
study=radia.optuna.Study(Directions=["minimize","minimize"],AutoSave=false);
values=[0.1 0.9;0.9 0.1;0.0 0.0;0.2 0.2];
violations=[-1;-1;1;2];
for k=1:4
    trial=study.ask();
    study.tell(trial,values(k,:));
    study.recordConstraints(trial,violations(k));
end
[good,weights]=radia.optuna.internal.ParetoSupport.splitMOTPE( ...
    study,(0:3)',values,2);
verifyEqual(testCase,good,[true;true;false;false]);
verifySize(testCase,weights,[2 1]);
verifyTrue(testCase,all(weights>0));
end

function testParetoRankAndCrowdingContract(testCase)
values=[0 3;1 2;2 1;3 3;2 4];
[rank,crowding]=radia.optuna.internal.ParetoSupport.rankAndCrowding( ...
    values,["minimize","minimize"]);
verifyEqual(testCase,rank(1:3),ones(3,1));
verifyEqual(testCase,rank,[1;1;1;2;2]);
verifyEqual(testCase,crowding,[Inf;2;Inf;Inf;Inf]);
native=radia.optuna.nativeStatus();
verifyEqual(testCase,native.backend,"native-mex");
verifyEmpty(testCase,native.missing_commands);
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

function values=constrainedMultiObjective(trial)
x=trial.suggestFloat("x",-1,1);
values=[x^2,(x-0.5)^2];
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
