function result=profileResidentHMatrix(inputPath,resultPath)
% Diagnostic resident vs per-application setup/transfer; no FE assembly claim.
arguments
    inputPath (1,:) char
    resultPath (1,:) char
end
inputs=load(inputPath);test=inputs.cases(4);
assert(strcmp(test.name,'mixed-uneven'));
result=struct('accepted',false,'steps',3,'samples',struct([]), ...
    'scope','Numeric leaf construction, transfers, repeated operator apply, recovery and release only');
try
    gpu=gpuDevice();result.sampled_free_before=gpu.AvailableMemory;
    op=radia.gpu.HMatrixOperator(test.config);free=onCleanup(@()delete(op));
    result.memory=op.stats();clear free
    result.cold_rebuild=run('rebuild');result.cold_resident=run('resident');
    run('rebuild');run('resident');
    schedule={'rebuild','resident','resident','rebuild'};
    index=0;
    for round=1:2
        for k=1:4
            sample=run(schedule{k});index=index+1;
            if index==1,result.samples=sample;else,result.samples(index)=sample;end
        end
    end
    result.sampled_free_after=gpu.AvailableMemory;
    result.live_count=radia_cuda_hmatrix('live_count');assert(result.live_count==0);
    result.cold_scope='First payload calls after MATLAB/MEX initialization and a memory probe; not cold process or MEX loading.';
    result.memory_limitations='Owned device buffers exact; CUDA driver/MATLAB free-memory samples are not peaks. Native host/process peak unmeasured.';
    result.accepted=true;
catch failure
    result.error=struct('identifier',failure.identifier,'message',failure.message);
    saveResult();rethrow(failure)
end
saveResult();
    function sample=run(mode)
        value=test.initial;timer=tic;
        if strcmp(mode,'resident')
            object=radia.gpu.HMatrixOperator(test.config);closer=onCleanup(@()delete(object));
            object.upload(value);object.apply(3);value=object.download();clear closer
        else
            for step=1:3
                object=radia.gpu.HMatrixOperator(test.config);closer=onCleanup(@()delete(object));
                object.upload(value);object.apply();value=object.download();clear closer
            end
        end
        elapsed=toc(timer);
        difference=norm(value-test.expected_third)/max(norm(test.expected_third),realmin);
        assert(difference<=1e-12);
        sample=struct('mode',mode,'whole_seconds',elapsed,'amortized_apply_seconds',elapsed/3,'relative_error',difference);
    end
    function saveResult()
        file=fopen(resultPath,'w');assert(file>=0);closer=onCleanup(@()fclose(file)); %#ok<NASGU>
        fprintf(file,'%s\n',jsonencode(result,PrettyPrint=true));
    end
end
