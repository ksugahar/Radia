function benchmark_optuna_history_store(options)
%BENCHMARK_OPTUNA_HISTORY_STORE Scaling of radia.optuna's history stores.
%   Study keeps the optimization history in column stores behind
%   materialized table views, with one trial-number bucket index per store
%   (radia.optuna.internal.TrialRowIndex). Three separate things are timed
%   because they are limited by different costs:
%
%     write   ask/suggest/report/tell for N trials. Appending to a MATLAB
%             table copies every existing row, so this is what the column
%             stores buy; the exponent of t ~ N^p should stay near 1.
%     freeze  study.get_trials(). Dominated by per-trial FrozenTrial
%             construction and datetime conversion, not by history lookup,
%             so freezeTrials converts each serial timestamp column once
%             and slices it instead of calling datetime() per trial.
%     lookup  a FIXED number of per-trial history probes against a GROWING
%             store. This is the only measurement that isolates the bucket
%             index: indexed time should stay flat while the scan reference
%             grows with the store.
%
%   POLICY: wall-clock numbers are only meaningful on an idle quiet compute
%   host (mdx or hibino). Timings taken on LAB are polluted by concurrent
%   builds and must not be quoted.
%
%   Usage:
%     benchmark_optuna_history_store
%     benchmark_optuna_history_store(TrialCounts=[250 500 1000 2000 4000])

arguments
    options.TrialCounts (1,:) double = [250 500 1000 2000 4000]
    options.ReportedSteps (1,1) double = 4
    options.ProbeCount (1,1) double = 200
    options.Output (1,1) string = ""
end

here = fileparts(mfilename("fullpath"));
repositoryRoot = fileparts(fileparts(here));
addpath(fullfile(repositoryRoot, "matlab"));
if strlength(options.Output) == 0
    options.Output = fullfile(here, "results_optuna_history_store.json");
end

% MATLAB compiles class methods on first use; without a warm-up the
% smallest case absorbs that cost and the fitted exponent is meaningless.
warmup = buildStudy(64, options.ReportedSteps);
warmup.get_trials();
probeHistory(warmup, warmup.TrialTable.TrialNumber(1:16));
scanHistory(warmup, warmup.TrialTable.TrialNumber(1:16));

counts = sort(options.TrialCounts);
results = cell(1, numel(counts));
for index = 1:numel(counts)
    results{index} = measureOne(counts(index), options.ReportedSteps, ...
        options.ProbeCount);
    entry = results{index};
    fprintf("N=%5d  write=%7.3f s  freeze=%7.3f s  " + ...
        "lookup=%8.4f s  scan=%8.4f s  scan/lookup=%5.1fx\n", ...
        entry.iterations, entry.t_setup, entry.t_solve, ...
        entry.t_indexed_lookup, entry.t_scan_lookup, entry.lookup_speedup);
end

payload = struct( ...
    "timestamp", string(datetime("now", TimeZone="local", ...
        Format="uuuu-MM-dd'T'HH:mm:ssXXX")), ...
    "hostname", string(getenv("COMPUTERNAME")), ...
    "benchmark", "radia.optuna history store scaling", ...
    "schema", "radia.optuna.history-store-benchmark.v1", ...
    "matlab_version", string(version), ...
    "problem", struct( ...
        "trial_counts", counts, ...
        "parameters_per_trial", 2, ...
        "reported_steps_per_pruned_trial", options.ReportedSteps, ...
        "pruned_fraction", 1/3, ...
        "lookup_probe_count", options.ProbeCount), ...
    "results", {results});
payload.scaling = struct( ...
    "write_exponent", fitExponent(counts, cellfun(@(e) e.t_setup, results)), ...
    "freeze_exponent", fitExponent(counts, cellfun(@(e) e.t_solve, results)), ...
    "indexed_lookup_exponent", fitExponent(counts, ...
        cellfun(@(e) e.t_indexed_lookup, results)), ...
    "scan_lookup_exponent", fitExponent(counts, ...
        cellfun(@(e) e.t_scan_lookup, results)));
fprintf("fitted exponents: write=%.2f freeze=%.2f " + ...
    "indexed_lookup=%.2f scan_lookup=%.2f\n", ...
    payload.scaling.write_exponent, payload.scaling.freeze_exponent, ...
    payload.scaling.indexed_lookup_exponent, ...
    payload.scaling.scan_lookup_exponent);

fid = fopen(options.Output, "w", "n", "UTF-8");
closer = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(payload, PrettyPrint=true), "char");
clear closer
fprintf("Results saved to %s\n", options.Output);
end

function study = buildStudy(trialCount, reportedSteps)
study = radia.optuna.Study(Name="history-bench", ...
    Sampler=radia.optuna.RandomSampler(17), ...
    Pruner=radia.optuna.NopPruner(), AutoSave=false);
for index = 1:trialCount
    trial = study.ask();
    x = trial.suggest_float("x", -1, 1);
    trial.suggest_int("k", 0, 9);
    if mod(index, 3) == 0
        for step = 0:(reportedSteps - 1)
            trial.report(x * x + step, step);
        end
        study.tell(trial, State="PRUNED");
    else
        study.tell(trial, x * x);
    end
end
end

function entry = measureOne(trialCount, reportedSteps, probeCount)
writeStart = tic;
study = buildStudy(trialCount, reportedSteps);
tWrite = toc(writeStart);

freezeStart = tic;
frozen = study.get_trials();
tFreeze = toc(freezeStart);

% Hold the probe count fixed so only the store size varies.
numbers = study.TrialTable.TrialNumber;
probes = numbers(round(linspace(1, numel(numbers), ...
    min(probeCount, numel(numbers)))));
indexedStart = tic;
indexed = probeHistory(study, probes);
tIndexed = toc(indexedStart);
scanStart = tic;
scanned = scanHistory(study, probes);
tScan = toc(scanStart);

entry = struct( ...
    "iterations", trialCount, ...
    "t_setup", tWrite, ...
    "t_solve", tFreeze, ...
    "t_indexed_lookup", tIndexed, ...
    "t_scan_lookup", tScan, ...
    "lookup_speedup", tScan / max(tIndexed, eps), ...
    "peak_memory_mb", peakMemoryMB(), ...
    "frozen_trials", numel(frozen), ...
    "intermediate_rows", height(study.IntermediateTable), ...
    "parameter_rows", height(study.ParamTable), ...
    "probe_count", numel(probes), ...
    "converged", numel(frozen) == trialCount && isequaln(indexed, scanned));
end

function [steps, values] = probeHistory(study, probes)
% The indexed path: one bucket lookup per probed trial.
[steps, values] = study.lastIntermediateValues(probes);
end

function [steps, values] = scanHistory(study, probes)
% The unindexed reference: scan the public view once per probed trial.
intermediate = study.IntermediateTable;
allNumbers = intermediate.TrialNumber;
allSteps = intermediate.Step;
allValues = intermediate.Value;
probes = reshape(double(probes), [], 1);
steps = NaN(size(probes));
values = NaN(size(probes));
for index = 1:numel(probes)
    rows = find(allNumbers == probes(index));
    if isempty(rows)
        continue
    end
    [steps(index), position] = max(allSteps(rows));
    values(index) = allValues(rows(position));
end
end

function value = peakMemoryMB()
try
    info = memory;
    value = info.MemUsedMATLAB / 1024 / 1024;
catch
    value = NaN;
end
end

function exponent = fitExponent(counts, times)
usable = counts > 0 & times > 0;
if sum(usable) < 2
    exponent = NaN;
    return
end
coefficients = polyfit(log(counts(usable)), log(times(usable)), 1);
exponent = coefficients(1);
end
