function result=benchmark_optuna_inference(outputPath,options)
%BENCHMARK_OPTUNA_INFERENCE Closed-loop optimization and proposal latency.
arguments
    outputPath (1,1) string
    options.Trials (1,1) double {mustBeInteger,mustBePositive} = 800
    options.Seeds (1,:) double = [11 37]
    options.ExecutionMode (1,1) string {mustBeMember(options.ExecutionMode,["cpu","gpu"])} = "cpu"
    options.BatchSize (1,1) double {mustBeInteger,mustBePositive} = 1
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'packages','radia-optuna','matlab'));
started=tic; previous=radia.optuna.executionMode(options.ExecutionMode);
configureSeconds=toc(started);
restoreMode=onCleanup(@()radia.optuna.executionMode(previous));
entries=cell(0,1);
for problem=["branin","correlated_valley"]
    for seed=options.Seeds
        sampler=radia.optuna.TPESampler(Seed=seed,Multivariate=true);
        study=radia.optuna.Study(Sampler=sampler,AutoSave=false);
        if problem=="branin", dimensions=2; else, dimensions=40; end
        proposals=zeros(options.Trials,dimensions);
        values=zeros(options.Trials,1); best=zeros(options.Trials,1);
        elapsed=zeros(options.Trials,4);
        before=optuna_mex('api.info'); memoryBefore=memory;
        for batchStart=1:options.BatchSize:options.Trials
        batchEnd=min(options.Trials,batchStart+options.BatchSize-1);
        pending=cell(batchEnd-batchStart+1,1);
        for index=batchStart:batchEnd
            started=tic; trial=study.ask(); elapsed(index,1)=toc(started);
            pending{index-batchStart+1}=trial;
            started=tic;
            if problem=="branin"
                proposals(index,1)=trial.suggest_float('x1',-5,10);
                proposals(index,2)=trial.suggest_float('x2',0,15);
            else
                for column=1:dimensions
                    proposals(index,column)=trial.suggest_float("x"+column,-2,2);
                end
            end
            elapsed(index,2)=toc(started);
        end
        for index=batchStart:batchEnd
            trial=pending{index-batchStart+1};
            started=tic; x=proposals(index,:);
            if problem=="branin"
                values(index)=(x(2)-5.1/(4*pi^2)*x(1)^2+5/pi*x(1)-6)^2 + ...
                    10*(1-1/(8*pi))*cos(x(1))+10;
            else
                values(index)=sum((x(2:end)-0.85*x(1:end-1)).^2)+0.01*sum(x.^2);
            end
            elapsed(index,3)=toc(started);
            started=tic; study.tell(trial,values(index)); elapsed(index,4)=toc(started);
            if index==1, best(index)=values(index); else, best(index)=min(best(index-1),values(index)); end
        end
        clear pending trial
        end
        memoryAfter=memory;
        entries{end+1,1}=struct('problem',problem,'seed',seed,'dimensions',dimensions, ...
            'trial_count',options.Trials,'proposals',proposals,'values',values, ...
            'best_curve',best,'best_value',best(end),'initial_ten_best',best(min(10,end)), ...
            'time_columns',{{'ask','suggest_all','objective','tell'}},'samples_s',elapsed, ...
            'inference_seconds',sum(elapsed(:,1:2),'all'), ...
            'objective_seconds',sum(elapsed(:,3)), 'history_update_seconds',sum(elapsed(:,4)), ...
            'first_proposal_seconds',sum(elapsed(1,1:2)), ...
            'native_before',before,'native_after',optuna_mex('api.info'), ...
            'matlab_memory_before_bytes',memoryBefore.MemUsedMATLAB, ...
            'matlab_memory_after_bytes',memoryAfter.MemUsedMATLAB);
        fprintf('%s seed=%d best=%g inference=%g s\n',problem,seed,best(end),entries{end}.inference_seconds);
        result=struct('schema','radia.validation.optuna-closed-loop-inference.v1', ...
            'matlab',version,'execution_mode',options.ExecutionMode, ...
            'multivariate',true,'candidates',24,'batch_size',options.BatchSize, ...
            'module_configure_s',configureSeconds, ...
            'storage_enabled',false,'all_history_retained',true,'entries',{entries}, ...
            'notes',{{'Every trial evaluates its objective and updates the study; no prefilled history.', ...
            'Inference is ask plus all suggestions, including model/candidate preparation and parameter recording.', ...
            'Objective evaluation and tell/history update are separately timed.', ...
            'Batch proposals all precede objective/tell updates, in ordered ask and tell loops; objectives are evaluated serially.', ...
            'Branin and a 40-dimensional correlated quadratic valley are optimization benchmarks, not CAE acceptance models.', ...
            'Preserve proposal, objective and convergence curves when comparing backends.'}});
        fid=fopen(outputPath,'w','n','UTF-8'); assert(fid>=0);
        closeFile=onCleanup(@()fclose(fid));
        fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true)); clear closeFile
        clear study sampler
    end
end
end
