function result = benchmark_matlab_optuna49(outputPath)
%BENCHMARK_MATLAB_OPTUNA49 Benchmark the MATLAB Optuna 4.9-compatible TPE.
% Run this on the same otherwise-idle host as benchmark_optuna49_python.py.

arguments
    outputPath (1, 1) string = ""
end
here = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(here));
addpath(fullfile(repoRoot, "matlab"));

trials = 100;
repeats = 11;
warmupRepeats = 3;
expectedScalarChecksum = 20.040135043951892;
expectedGroupedChecksum = 104.33176385944043;

[scalarSeconds, scalarChecksums] = measure(@runScalar, trials, repeats);
[groupedSeconds, groupedChecksums] = measure(@runGrouped, trials, repeats);
verifyChecksums(scalarChecksums, expectedScalarChecksum, "scalar");
verifyChecksums(groupedChecksums, expectedGroupedChecksum, "grouped conditional");

commands = string(radia.optuna.internal.NativeKernels.call("api.commands"));
nativeStatus = radia.optuna.nativeStatus();
host = string(getenv("COMPUTERNAME"));
if strlength(host) == 0
    [status, hostText] = system("hostname");
    if status == 0
        host = strtrim(string(hostText));
    end
end
result = struct( ...
    "schema", "radia.validation.optuna49-performance-runtime.v1", ...
    "generated_at", string(datetime("now", "TimeZone", "UTC", ...
        "Format", "yyyy-MM-dd'T'HH:mm:ss.SSSXXX")), ...
    "runtime", "matlab", ...
    "host", host, ...
    "versions", struct( ...
        "matlab", version, ...
        "native_gateway", nativeStatus.gateway, ...
        "optuna_mex_command_count", numel(commands)), ...
    "settings", struct( ...
        "trials", trials, ...
        "total_repeats", repeats, ...
        "warmup_repeats", warmupRepeats, ...
        "reported_repeats", repeats - warmupRepeats), ...
    "scalar", summarize(scalarSeconds, scalarChecksums, ...
        trials, warmupRepeats), ...
    "grouped_conditional", summarize(groupedSeconds, groupedChecksums, ...
        trials, warmupRepeats));

encoded = jsonencode(result, PrettyPrint=true);
if strlength(outputPath) > 0
    fileId = fopen(outputPath, "w");
    if fileId < 0
        error("radia:optuna:BenchmarkOutput", ...
            "Could not open benchmark output: %s", outputPath);
    end
    cleanup = onCleanup(@() fclose(fileId));
    fprintf(fileId, "%s\n", encoded);
    clear cleanup
end
fprintf("%s\n", encoded);
end

function [durations, checksums] = measure(workload, trials, repeats)
durations = zeros(1, repeats);
checksums = zeros(1, repeats);
for repeat = 1:repeats
    started = tic;
    checksums(repeat) = workload(trials);
    durations(repeat) = toc(started);
end
end

function summary = summarize(durations, checksums, trials, warmupRepeats)
medianSeconds = median(durations(warmupRepeats + 1:end));
summary = struct( ...
    "all_seconds", durations, ...
    "median_warmed_seconds", medianSeconds, ...
    "trials_per_second", trials / medianSeconds, ...
    "checksum", checksums(end));
end

function verifyChecksums(checksums, expected, label)
if any(abs(checksums - expected) > 1e-12)
    error("radia:optuna:BenchmarkChecksum", ...
        "The %s checksum changed: %.17g instead of %.17g.", ...
        label, checksums(end), expected);
end
end

function checksum = runScalar(trials)
study = radia.optuna.Study(Sampler=radia.optuna.TPESampler( ...
    Seed=37, NStartupTrials=4), AutoSave=false);
checksum = 0;
for index = 1:trials
    trial = study.ask();
    x = trial.suggest_float("x", -2, 2);
    study.tell(trial, (x - 0.25)^2);
    checksum = checksum + x;
end
end

function checksum = runGrouped(trials)
study = radia.optuna.Study(Sampler=radia.optuna.TPESampler( ...
    Seed=101, NStartupTrials=4, Multivariate=true, Group=true), ...
    AutoSave=false);
checksum = 0;
for index = 1:trials
    trial = study.ask();
    branch = string(trial.suggest_categorical( ...
        "branch", ["left", "right"]));
    x = trial.suggest_float("x", -1, 1);
    if branch == "left"
        y = trial.suggest_float("y", 0, 2);
        value = (x - 0.2)^2 + (y - 0.4)^2;
        checksum = checksum + x + y;
    else
        z = trial.suggest_int("z", 1, 5);
        value = (x + 0.1)^2 + 0.05 * z;
        checksum = checksum + x + z;
    end
    study.tell(trial, value);
end
end
