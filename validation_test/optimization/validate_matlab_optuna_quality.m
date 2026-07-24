function report = validate_matlab_optuna_quality(options)
%VALIDATE_MATLAB_OPTUNA_QUALITY Measure sampler quality and ask/tell cost.
%   This validation compares seeded samplers at equal trial budgets. It is
%   deliberately separate from unit tests: the output characterizes search
%   quality and overhead, not numerical correctness of a field solver.
arguments
    options.Seeds (1,:) double = 0:4
    options.SingleTrials (1,1) double {mustBeInteger,mustBePositive} = 80
    options.MultiTrials (1,1) double {mustBeInteger,mustBePositive} = 120
    options.OutputFile (1,1) string = ""
end
root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
matlabPath = fullfile(root, "matlab");
matlabPathWasPresent = any( ...
    strcmpi(string(strsplit(path, pathsep)), string(matlabPath)));
if ~matlabPathWasPresent
    addpath(matlabPath);
end
pathCleanup = onCleanup(@() restoreMatlabPath( ...
    matlabPath, matlabPathWasPresent));
if strlength(options.OutputFile) == 0
    options.OutputFile = fullfile(fileparts(mfilename("fullpath")), ...
        "validation_matlab_optuna_quality_summary.json");
end

singleSamplers = ["random", "tpe", "tpe_joint", "cmaes"];
singleProblems = ["branin", "correlated_valley"];
single = struct();
for problem = singleProblems
    for samplerName = singleSamplers
        bestValues = zeros(numel(options.Seeds), 1);
        elapsed = zeros(numel(options.Seeds), 1);
        for index = 1:numel(options.Seeds)
            sampler = makeSingleSampler(samplerName, options.Seeds(index));
            study = radia.optuna.createStudy(direction="minimize", ...
                Sampler=sampler, AutoSave=false);
            timer = tic;
            switch problem
                case "branin"
                    study.optimize(@braninObjective, options.SingleTrials);
                    bestValues(index) = study.bestValue() - ...
                        0.39788735772973816;
                case "correlated_valley"
                    study.optimize(@correlatedValleyObjective, ...
                        options.SingleTrials);
                    bestValues(index) = study.bestValue();
            end
            elapsed(index) = toc(timer);
        end
        key = matlab.lang.makeValidName(problem + "_" + samplerName);
        single.(key) = summarizeSingle(bestValues, elapsed, ...
            options.SingleTrials);
    end
end

multiSamplers = ["random", "motpe", "nsgaii"];
multi = struct();
for samplerName = multiSamplers
    frontError = zeros(numel(options.Seeds), 1);
    coverage = zeros(numel(options.Seeds), 1);
    pointCount = zeros(numel(options.Seeds), 1);
    elapsed = zeros(numel(options.Seeds), 1);
    for index = 1:numel(options.Seeds)
        sampler = makeMultiSampler(samplerName, options.Seeds(index));
        study = radia.optuna.createStudy( ...
            directions=["minimize","minimize"], ...
            Sampler=sampler, AutoSave=false);
        timer = tic;
        study.optimize(@zdt1Objective, options.MultiTrials);
        elapsed(index) = toc(timer);
        front = study.paretoFront();
        values = vertcat(front.Values{:});
        idealF2 = 1 - sqrt(max(0, min(1, values(:,1))));
        frontError(index) = mean(max(0, values(:,2) - idealF2));
        coverage(index) = max(values(:,1)) - min(values(:,1));
        pointCount(index) = size(values, 1);
    end
    key = matlab.lang.makeValidName(samplerName);
    multi.(key) = struct( ...
        "median_front_error", median(frontError), ...
        "front_error", frontError, ...
        "median_f1_coverage", median(coverage), ...
        "f1_coverage", coverage, ...
        "median_pareto_points", median(pointCount), ...
        "pareto_points", pointCount, ...
        "median_elapsed_s", median(elapsed), ...
        "median_ask_tell_ms_per_trial", ...
            1e3 * median(elapsed) / options.MultiTrials);
end

qualityGates = struct( ...
    "branin_tpe_beats_random", ...
        single.branin_tpe.median_regret < ...
        single.branin_random.median_regret, ...
    "branin_joint_tpe_beats_random", ...
        single.branin_tpe_joint.median_regret < ...
        single.branin_random.median_regret, ...
    "branin_cmaes_beats_random", ...
        single.branin_cmaes.median_regret < ...
        single.branin_random.median_regret, ...
    "correlated_tpe_beats_random", ...
        single.correlated_valley_tpe.median_regret < ...
        single.correlated_valley_random.median_regret, ...
    "correlated_joint_tpe_beats_random", ...
        single.correlated_valley_tpe_joint.median_regret < ...
        single.correlated_valley_random.median_regret, ...
    "correlated_cmaes_beats_random", ...
        single.correlated_valley_cmaes.median_regret < ...
        single.correlated_valley_random.median_regret, ...
    "motpe_front_error_beats_random", ...
        multi.motpe.median_front_error < ...
        multi.random.median_front_error, ...
    "motpe_coverage_exceeds_half_domain", ...
        multi.motpe.median_f1_coverage > 0.5);
gateValues = struct2array(qualityGates);
qualityGates.all_passed = all(gateValues);

report = struct( ...
    "schema", "radia.validation.matlab-optuna-quality.v3", ...
    "generated_at_utc", string(datetime("now", "TimeZone", "UTC", ...
        "Format", "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'")), ...
    "hostname", currentHostname(), ...
    "matlab_version", string(version), ...
    "matlab_release", string(version("-release")), ...
    "platform", string(computer), ...
    "seeds", options.Seeds, ...
    "single_trials", options.SingleTrials, ...
    "multi_trials", options.MultiTrials, ...
    "single_objective", single, ...
    "multi_objective_zdt1", multi, ...
    "quality_gates", qualityGates, ...
    "interpretation", struct( ...
        "lower_single_regret_is_better", true, ...
        "lower_front_error_is_better", true, ...
        "higher_f1_coverage_is_better", true, ...
        "tpe_joint_uses_automatic_intersection_search_space", true, ...
        "cmaes_uses_full_covariance_evolution_paths", true, ...
        "sampler_state_storage", "MATLAB SamplerStateTable", ...
        "timing_scope", "MATLAB sampler plus table-backed ask/tell and analytic objective"));
writeJson(options.OutputFile, report);
fprintf("MATLAB Optuna quality report: %s\n", options.OutputFile);
clear pathCleanup
end

function sampler = makeSingleSampler(name, seed)
switch name
    case "random"
        sampler = radia.optuna.RandomSampler(seed);
    case "tpe"
        sampler = radia.optuna.TPESampler(Seed=seed, NStartupTrials=10);
    case "tpe_joint"
        sampler = radia.optuna.TPESampler( ...
            Seed=seed, NStartupTrials=10, Multivariate=true);
    case "cmaes"
        sampler = radia.optuna.CmaEsSampler(Seed=seed, NStartupTrials=1);
end
end

function sampler = makeMultiSampler(name, seed)
switch name
    case "random"
        sampler = radia.optuna.RandomSampler(seed);
    case "motpe"
        sampler = radia.optuna.MOTPESampler(Seed=seed, NStartupTrials=20);
    case "nsgaii"
        sampler = radia.optuna.NSGAIISampler(Seed=seed, PopulationSize=24);
end
end

function value = braninObjective(trial)
x1 = trial.suggestFloat("x1", -5, 10);
x2 = trial.suggestFloat("x2", 0, 15);
value = branin(x1, x2);
end

function value = branin(x1, x2)
a = 1;
b = 5.1/(4*pi^2);
c = 5/pi;
r = 6;
s = 10;
t = 1/(8*pi);
value = a*(x2-b*x1^2+c*x1-r)^2 + s*(1-t)*cos(x1) + s;
end

function value = correlatedValleyObjective(trial)
x = zeros(1, 4);
for index = 1:numel(x)
    x(index) = trial.suggestFloat("x" + index, -2, 2);
end
residual = x(2:end) - 0.85*x(1:end-1);
value = sum(residual.^2) + 0.01*sum(x.^2);
end

function values = zdt1Objective(trial)
x1 = trial.suggestFloat("x1", 0, 1);
x2 = trial.suggestFloat("x2", 0, 1);
g = 1 + 9*x2;
values = [x1, g*(1-sqrt(x1/g))];
end

function summary = summarizeSingle(regret, elapsed, trials)
summary = struct( ...
    "median_regret", median(regret), ...
    "regret", regret, ...
    "median_elapsed_s", median(elapsed), ...
    "median_ask_tell_ms_per_trial", 1e3*median(elapsed)/trials);
end

function writeJson(path, value)
parent = fileparts(path);
if strlength(parent) > 0 && ~isfolder(parent)
    mkdir(parent);
end
file = fopen(path, "w", "n", "UTF-8");
if file < 0
    error("radia:validation:OptunaWrite", ...
        "Cannot write MATLAB Optuna validation report: %s", path);
end
cleanup = onCleanup(@() fclose(file));
count = fprintf(file, "%s\n", jsonencode(value, PrettyPrint=true));
if count <= 0
    error("radia:validation:OptunaWrite", ...
        "Cannot write MATLAB Optuna validation report: %s", path);
end
clear cleanup
end

function restoreMatlabPath(matlabPath, wasPresent)
if wasPresent
    return
end
entries = string(strsplit(path, pathsep));
if any(strcmpi(entries, string(matlabPath)))
    rmpath(matlabPath);
end
end

function name = currentHostname()
name = string(getenv("COMPUTERNAME"));
if strlength(name) > 0
    return
end
[status, output] = system("hostname");
if status == 0
    name = strip(string(output));
end
end
