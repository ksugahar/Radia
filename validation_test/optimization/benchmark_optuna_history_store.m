function payload=benchmark_optuna_history_store(options)
%BENCHMARK_OPTUNA_HISTORY_STORE Long scaling check for Study history access.
%   Keep this benchmark in validation_test: it deliberately builds several
%   large studies and is not part of the short unit-test lane.

arguments
    options.TrialCounts (1,:) double = [250 500 1000 2000 4000]
    options.ReportedSteps (1,1) double = 4
    options.ProbeCount (1,1) double = 200
    options.Output (1,1) string = ""
end

here=fileparts(mfilename("fullpath"));
root=fileparts(fileparts(here));
addpath(fullfile(root,"matlab"));
warmup=buildStudy(32,options.ReportedSteps);
warmup.get_trials();

counts=sort(options.TrialCounts);
results=cell(1,numel(counts));
for index=1:numel(counts)
    results{index}=measureOne( ...
        counts(index),options.ReportedSteps,options.ProbeCount);
    item=results{index};
    fprintf("N=%5d write=%7.3f freeze=%7.3f " + ...
        "indexed=%8.4f scan=%8.4f speedup=%5.1fx\n", ...
        item.trials,item.write_s,item.freeze_s,item.indexed_lookup_s, ...
        item.scan_lookup_s,item.lookup_speedup);
end
payload=struct( ...
    "schema","radia.optuna.history-store-benchmark.v1", ...
    "timestamp",string(datetime("now",TimeZone="local", ...
        Format="uuuu-MM-dd'T'HH:mm:ssXXX")), ...
    "hostname",string(getenv("COMPUTERNAME")), ...
    "matlab_version",string(version), ...
    "trial_counts",counts,"results",{results}, ...
    "freeze_exponent",fitExponent(counts, ...
        cellfun(@(item)item.freeze_s,results)), ...
    "indexed_lookup_exponent",fitExponent(counts, ...
        cellfun(@(item)item.indexed_lookup_s,results)));
if strlength(options.Output)>0
    fid=fopen(options.Output,"w","n","UTF-8");
    if fid<0
        error("radia:optuna:HistoryBenchmark", ...
            "Cannot open '%s' for writing.",options.Output);
    end
    cleanup=onCleanup(@()fclose(fid));
    fwrite(fid,jsonencode(payload,PrettyPrint=true),"char");
    clear cleanup
end
end

function study=buildStudy(trialCount,reportedSteps)
study=radia.optuna.Study(Name="history-bench", ...
    Sampler=radia.optuna.RandomSampler(17), ...
    Pruner=radia.optuna.NopPruner(),AutoSave=false);
for index=1:trialCount
    trial=study.ask();
    x=trial.suggest_float("x",-1,1);
    trial.suggest_int("k",0,9);
    if mod(index,3)==0
        for step=0:(reportedSteps-1)
            trial.report(x*x+step,step);
        end
        study.tell(trial,State="PRUNED");
    else
        study.tell(trial,x*x);
    end
end
end

function result=measureOne(trialCount,reportedSteps,probeCount)
started=tic;
study=buildStudy(trialCount,reportedSteps);
writeTime=toc(started);
started=tic;
frozen=study.get_trials();
freezeTime=toc(started);
numbers=study.TrialTable.TrialNumber;
probes=numbers(round(linspace(1,numel(numbers), ...
    min(probeCount,numel(numbers)))));
started=tic;
[indexedSteps,indexedValues]=study.lastIntermediateValues(probes);
indexedTime=toc(started);
started=tic;
[scanSteps,scanValues]=scanHistory(study,probes);
scanTime=toc(started);
result=struct( ...
    "trials",trialCount,"write_s",writeTime,"freeze_s",freezeTime, ...
    "indexed_lookup_s",indexedTime,"scan_lookup_s",scanTime, ...
    "lookup_speedup",scanTime/max(indexedTime,eps), ...
    "frozen_trials",numel(frozen), ...
    "intermediate_rows",height(study.IntermediateTable), ...
    "parameter_rows",height(study.ParamTable), ...
    "converged",isequaln(indexedSteps,scanSteps) && ...
        isequaln(indexedValues,scanValues));
end

function [steps,values]=scanHistory(study,probes)
source=study.IntermediateTable;
probes=reshape(double(probes),[],1);
steps=NaN(size(probes));
values=NaN(size(probes));
for index=1:numel(probes)
    rows=find(source.TrialNumber==probes(index));
    if ~isempty(rows)
        [steps(index),position]=max(source.Step(rows));
        values(index)=source.Value(rows(position));
    end
end
end

function exponent=fitExponent(counts,times)
coefficients=polyfit(log(counts),log(times),1);
exponent=coefficients(1);
end
