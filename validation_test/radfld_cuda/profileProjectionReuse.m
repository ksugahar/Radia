function [result,vectors] = profileProjectionReuse(volPath,source,cpuValues,jsonPath)
% Continuous three-step comparisons. Source stays unchanged across steps.
% CPU values are validation-only; GPU usage has no CPU-reference prerequisite.
arguments
    volPath (1,:) char
    source (1,1) struct
    cpuValues (:,3) double {mustBeReal,mustBeFinite}
    jsonPath (1,:) char
end
result=struct('accepted',false,'relative_limit',1e-4,'timing_acceptance',false, ...
    'steps',3,'chunk_size',128,'samples',struct([]),'a100_used',false);
vectors=[];
try
    validationTimer=tic;
    h=radia_cuda_projection('create',volPath,'hdiv',1);
    free=onCleanup(@()radia_cuda_projection('destroy',h));
    [reference,~]=radia_cuda_projection('project',h,cpuValues);
    targets=radia_cuda_projection('targets',h);
    clear free
    if strcmp(source.representation,'cloud')
        probeField=radia.gpu.FieldCloud(source.xyz,source.strength,LeafSize=8);
    else
        probeField=radia.gpu.FieldElements(source,LeafSize=8);
    end
    freeProbe=onCleanup(@()delete(probeField));
    baseline=probeField.field(targets,'tree');
    mixed=probeField.field(targets,'tree_mixed32');
    result.field_relative_error=norm(mixed-cpuValues,'fro')/max(norm(cpuValues,'fro'),realmin);
    result.FP64_field_relative_error=norm(baseline-cpuValues,'fro')/max(norm(cpuValues,'fro'),realmin);
    assert(max(result.field_relative_error,result.FP64_field_relative_error)<=1e-4);
    [~,cacheProbe]=evaluateProjectionTargets(targets,@(p)probeField.field(p,'tree'),"source-v1");
    changedSource=source;
    if strcmp(source.representation,'cloud')
        changedSource.strength=2*source.strength;
        changedField=radia.gpu.FieldCloud(changedSource.xyz,changedSource.strength,LeafSize=8);
    else
        vcol=13;scol=10;
        if strcmp(source.representation,'curved_tet'),vcol=31;scol=19;end
        changedSource.volume(:,vcol:end)=2*source.volume(:,vcol:end);
        changedSource.surface(:,scol:end)=2*source.surface(:,scol:end);
        changedField=radia.gpu.FieldElements(changedSource,LeafSize=8);
    end
    freeChanged=onCleanup(@()delete(changedField));
    [updated,~,updatedInfo]=evaluateProjectionTargets(targets,@(p)changedField.field(p,'tree'),"source-v2",cacheProbe);
    assert(~updatedInfo.reused && norm(updated-2*baseline,'fro')<=1e-12*max(norm(2*baseline,'fro'),realmin));
    result.source_revision_invalidation_passed=true;
    clear freeChanged freeProbe
    result.precision_validation_seconds=toc(validationTimer);
    modes={'rebuild64','chunk64','batch64','batch32','cache64','cache32'};
    for k=1:numel(modes)
        [sample,~]=runSequence(modes{k});
        result.warmup_sequences(k)=sample;
    end
    % Baseline/chunk vs batch, precision pair, and exact-value cache pair.
    pairs=[1 3;2 3;3 4;3 5;4 6];
    index=0;
    for pair=1:size(pairs,1)
        schedule=pairs(pair,[1 2 2 1]);
        for k=1:4
            [sample,vector]=runSequence(modes{schedule(k)});
            index=index+1;
            if index==1,result.samples=sample;else,result.samples(index)=sample;end
            vectors(:,index)=vector; %#ok<AGROW>
        end
    end
    result.accepted=true;
    result.cache_scope='Unchanged geometry, strength, options and explicit source revision only; not valid after source updates.';
    result.memory_limitations='Device owned buffer peak is exact; CUDA free-memory samples are not process/global peaks. MATLAB target/cache arrays counted by numel*8; native host peak unmeasured.';
catch failure
    result.error=struct('identifier',failure.identifier,'message',failure.message);
    saveResult();
    rethrow(failure)
end
saveResult();
    function [sample,vector]=runSequence(mode)
        token=uint64(0);field=[];tokenCleanup=[];fieldCleanup=[];
        calls=0;deviceMemory=struct();algorithm='tree';
        if endsWith(mode,'32'),algorithm='tree_mixed32';end
        cache=struct();timer=tic;
        rebuild=strcmp(mode,'rebuild64');
        if ~rebuild,create();end
        for step=1:3
            if rebuild,create();end
            if ~startsWith(mode,'cache'),cache=struct();end
            [vector,cache,info]=projectBatchedField(token,@evaluate,"fixed-source-v1",cache);
            if rebuild,release();end
        end
        if ~rebuild,release();end
        seconds=toc(timer);
        difference=norm(vector-reference)/max(norm(reference),realmin);
        assert(difference<=1e-4,'radia:projection:accuracy','CPU projection difference exceeds limit');
        sample=struct('mode',mode,'whole_seconds',seconds, ...
            'amortized_step_seconds',seconds/3,'relative_coefficient_error',difference, ...
            'field_calls',calls,'targets',info.target_count,'device_memory',deviceMemory, ...
            'matlab_cache_array_bytes',8*(numel(cache.points)+numel(cache.values)));
        if startsWith(mode,'cache'),assert(calls==1);end
        function create()
            token=radia_cuda_projection('create',volPath,'hdiv',1);
            tokenCleanup=onCleanup(@()radia_cuda_projection('destroy',token));
            if strcmp(source.representation,'cloud')
                field=radia.gpu.FieldCloud(source.xyz,source.strength,LeafSize=8);
            else
                field=radia.gpu.FieldElements(source,LeafSize=8);
            end
            fieldCleanup=onCleanup(@()delete(field));
        end
        function release()
            delete(field);
            fieldCleanup=[];tokenCleanup=[];
            field=[];token=uint64(0);
        end
        function values=evaluate(points)
            stride=size(points,1);
            if strcmp(mode,'chunk64'),stride=128;end
            values=zeros(size(points));
            for first=1:stride:size(points,1)
                rows=first:min(first+stride-1,size(points,1));
                [values(rows,:),deviceMemory]=field.field(points(rows,:),algorithm);
                calls=calls+1;
            end
        end
    end
    function saveResult()
        file=fopen(jsonPath,'w');assert(file>=0);close=onCleanup(@()fclose(file));
        fprintf(file,'%s\n',jsonencode(result,PrettyPrint=true));
    end
end
