function result=benchmark_optuna_parameter_history(outputPath,options)
%BENCHMARK_OPTUNA_PARAMETER_HISTORY In-memory independent TPE end-to-end cost.
arguments
    outputPath (1,1) string
    options.TrialCounts (1,:) double {mustBeInteger,mustBePositive} = [100 1000 5000]
    options.Dimensions (1,:) double {mustBeInteger,mustBePositive} = [1 10 40]
    options.Repeats (1,1) double {mustBeInteger,mustBeGreaterThan(options.Repeats,4)} = 13
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'packages','radia-optuna','matlab'));
previous=radia.optuna.executionMode("cpu");
restore=onCleanup(@()radia.optuna.executionMode(previous));
entries=cell(0,1);
for dimensions=options.Dimensions
    for count=options.TrialCounts
        sampler=radia.optuna.TPESampler(Seed=37,Multivariate=false);
        study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
        names="x"+string(1:dimensions);
        % Populate identical, deterministic history without timing setup.
        random=radia.optuna.Study(Sampler=radia.optuna.RandomSampler(19), ...
            AutoSave=false);
        for index=1:count
            trial=random.ask(); values=zeros(1,dimensions);
            for column=1:dimensions
                values(column)=trial.suggest_float(names(column),-2,2);
            end
            random.tell(trial,sum((values-0.25).^2));
        end
        study.add_trials(random.get_trials());
        clear random
        memoryBefore=memory;
        elapsed=zeros(options.Repeats,3);
        proposals=zeros(options.Repeats,dimensions);
        for repeat=1:options.Repeats
            started=tic; trial=study.ask(); elapsed(repeat,1)=toc(started);
            started=tic;
            for column=1:dimensions
                proposals(repeat,column)=trial.suggest_float(names(column),-2,2);
            end
            elapsed(repeat,2)=toc(started);
            started=tic;
            study.tell(trial,sum((proposals(repeat,:)-0.25).^2));
            elapsed(repeat,3)=toc(started);
        end
        entries{end+1,1}=struct('history',count,'dimensions',dimensions, ...
            'columns',{{'ask','suggest_all','tell'}},'samples_s',elapsed, ...
            'first_trial_s',sum(elapsed(1,:)), ...
            'warmed_median_s',median(sum(elapsed(5:end,:),2)), ...
            'matlab_memory_before_probes_bytes',memoryBefore.MemUsedMATLAB, ...
            'proposals',proposals,'parameter_rows',height(study.ParamTable));
        fprintf('N=%d D=%d median=%g s\n',count,dimensions, ...
            entries{end}.warmed_median_s);
        result=struct('schema','radia.validation.optuna-parameter-history.v1', ...
            'hostname',string(getenv('COMPUTERNAME')),'matlab',version, ...
            'execution_mode','cpu','storage_enabled',false,'autosave',false, ...
            'sampler','independent TPE','seed',37,'history_seed',19, ...
            'repeats',options.Repeats,'discarded_repeats',4, ...
            'entries',{entries},'api',optuna_mex('api.info'), ...
            'notes',{{'Includes ask, all suggestions, recording and tell.', ...
            'Cold history attachment and first proposal are excluded from warmed medians.', ...
            'History is retained in full; no MAT serialization in the measured region.', ...
            'Compare proposal arrays between source variants before reporting speedup.'}});
        fid=fopen(outputPath,'w'); assert(fid>=0);
        closeFile=onCleanup(@()fclose(fid));
        fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true)); clear closeFile
    end
end
end
