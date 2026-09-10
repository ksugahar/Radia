function result = benchmark_optuna_storage_scaling(outputPath, withStorage)
%BENCHMARK_OPTUNA_STORAGE_SCALING Attribute scalar TPE and storage costs.
% Diagnostic timings, not compatibility expectations or pure transfer costs.
arguments
    outputPath (1,1) string
    withStorage (1,1) logical = true
end
root=fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(root,"matlab"));
folder=string(tempname(fileparts(outputPath))); mkdir(folder);
cleanup=onCleanup(@() removeScratch(folder));
levels=[100,1000,10000]; repeats=7;
storagePath="";
if withStorage, storagePath=fullfile(folder,"study.mat"); end
study=radia.optuna.Study(Sampler=radia.optuna.TPESampler( ...
    Seed=37,NStartupTrials=4),AutoSave=false, ...
    StoragePath=storagePath);
count=0; entries=cell(1,numel(levels));
for k=1:numel(levels)
    started=tic;
    while count<levels(k)
        trial=study.ask(); x=trial.suggest_float("x",-2,2);
        study.tell(trial,(x-0.25)^2); count=count+1;
    end
    fillSeconds=toc(started);
    timings=zeros(repeats,6);
    for r=1:repeats
        started=tic; trial=study.ask(); timings(r,1)=toc(started);
        started=tic; x=trial.suggest_float("x",-2,2); timings(r,2)=toc(started);
        started=tic; study.tell(trial,(x-0.25)^2); timings(r,3)=toc(started);
        count=count+1;
        started=tic; a=study.TrialTable; b=study.ParamTable;
        c=study.ObjectiveTable; timings(r,4)=toc(started);
        assert(height(a)==count && height(b)==count && height(c)==count);
        clear a b c
        started=tic; a=study.TrialTable; b=study.ParamTable;
        c=study.ObjectiveTable; timings(r,5)=toc(started);
        clear a b c
        if withStorage
            started=tic; study.save(); timings(r,6)=toc(started);
        else
            timings(r,6)=NaN;
        end
    end
    profile clear; profile on;
    for r=1:repeats
        trial=study.ask(); x=trial.suggest_float("x",-2,2);
        study.tell(trial,(x-0.25)^2); count=count+1;
    end
    profile off; info=profile('info');
    trialProfile=profileRows(info.FunctionTable);
    saveProfile=[]; savedBytes=0;
    if withStorage
        profile clear; profile on; study.save(); profile off;
        info=profile('info'); saveProfile=profileRows(info.FunctionTable);
        file=dir(fullfile(folder,"study.mat")); savedBytes=file.bytes;
    end
    entries{k}=struct('history_before_probes',levels(k), ...
        'history_after_probes',count,'fill_seconds',fillSeconds, ...
        'columns',{{'ask','suggest','tell','dirty_tables','cached_tables','save_cached_tables'}}, ...
        'all_seconds',timings,'median_seconds',median(timings(3:end,:),1), ...
        'saved_bytes',savedBytes,'trial_profile',trialProfile,'save_profile',saveProfile);
    result=struct('schema','radia.validation.optuna-storage-scaling.v1', ...
        'host',getenv('COMPUTERNAME'),'matlab',version,'seed',37, ...
        'storage_path_enabled',withStorage,'autosave',false, ...
        'repeats',repeats,'discarded_repeats',2,'entries',{entries(1:k)}, ...
        'notes',{{'Scalar TPE, one float, completed trials, no callbacks.', ...
        'ask includes sampler preparation and native computation; suggest includes recording.', ...
        'MEX profile time includes gateway and native work, not separately measured transfer.', ...
        'save includes complete MAT serialization, read-back validation, and replacement/backup.', ...
        'Profiling runs are separate from uninstrumented wall timings.'}});
    fid=fopen(outputPath,'w'); assert(fid>=0);
    fileCleanup=onCleanup(@() fclose(fid)); fprintf(fid,'%s',jsonencode(result,PrettyPrint=true));
    clear fileCleanup
    fprintf('HISTORY %d MEDIAN_SECONDS %s\n',levels(k),mat2str(entries{k}.median_seconds));
end
clear cleanup
end

function rows=profileRows(functions)
rows=struct('name',{},'calls',{},'inclusive_seconds',{},'self_seconds',{});
for k=1:numel(functions)
    f=functions(k);
    rows(k)=struct('name',f.FunctionName,'calls',f.NumCalls, ...
        'inclusive_seconds',f.TotalTime, ...
        'self_seconds',max(0,f.TotalTime-sum([f.Children.TotalTime])));
end
[~,order]=sort([rows.self_seconds],'descend');
rows=rows(order(1:min(40,numel(order))));
end

function removeScratch(folder)
% Only the fresh, task-owned directory created above is removed.
if isfolder(folder), rmdir(folder,'s'); end
end
