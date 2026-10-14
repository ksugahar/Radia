function result = compareProjectionPrecision(volPath,xyz,strength,jsonPath,options)
% Diagnostic comparison on an otherwise idle 100; not development/build host acceptance.
% The caller checks CI, computation load and memory before calling this.
arguments
    volPath (1,:) char
    xyz (:,3) double {mustBeReal,mustBeFinite}
    strength (:,1) double {mustBeReal,mustBeFinite}
    jsonPath (1,:) char
    options.ElementSource (1,1) struct = struct()
end
result=struct('accepted',false,'timing_acceptance',false, ...
    'scope','Point-cloud actual NGSolve projection, FP32 far / FP64 near', ...
    'a100_used',false,'relative_limit',1e-4,'samples',struct([]));
if isfield(options.ElementSource,'representation')
    result.scope=['Retained ' options.ElementSource.representation ' actual NGSolve projection, FP32 far / FP64 near'];
end
try
% Warm both algorithms, then balanced ABBA repeats. Each sample includes
% source construction/upload, target capture, evaluation, FE projection,
% result recovery and release. Reuse is measured separately elsewhere.
algorithms={'tree','tree_mixed32'};
for k=1:2
    sample=runOne(algorithms{k});
    result.cold_samples(k)=rmfield(sample,'coefficients');
end
schedule=[1 2 2 1 1 2 2 1];
reference=[];
for k=1:numel(schedule)
    sample=runOne(algorithms{schedule(k)});
    if isempty(reference), reference=sample.coefficients; end
    sample.relative_coefficient_error=norm(sample.coefficients-reference)/max(norm(reference),realmin);
    assert(sample.relative_coefficient_error<=1e-4,'radia:precision:projection','Projection precision limit exceeded');
    sample=rmfield(sample,'coefficients');
    if k==1
        result.samples=sample;
    else
        result.samples(k)=sample;
    end
end
result.accepted=true;
result.memory_limitations='Owned device buffers are exact; cudaMemGetInfo samples do not measure process/device peak or host peak.';
result.precision_limitations='This comparison does not replace independent CPU physical FE L2 acceptance.';
catch failure
    result.error=struct('identifier',failure.identifier,'message',failure.message);
    saveResult();
    rethrow(failure)
end
saveResult();
    function sample=runOne(algorithm)
        timer=tic;
        h=radia_cuda_projection('create',volPath,'hdiv',1);
        release=onCleanup(@()radia_cuda_projection('destroy',h));
        if isfield(options.ElementSource,'representation')
            f=radia.gpu.FieldElements(options.ElementSource,LeafSize=8);
        else
            f=radia.gpu.FieldCloud(xyz,strength,LeafSize=8);
        end
        releaseField=onCleanup(@()delete(f));
        points=radia_cuda_projection('targets',h);
        [values,memory]=f.field(points,algorithm);
        [coefficients,~]=radia_cuda_projection('project',h,values);
        clear releaseField release
        sample=struct('algorithm',algorithm,'whole_seconds',toc(timer), ...
            'coefficients',coefficients,'memory',memory,'target_count',size(points,1));
    end
    function saveResult()
        file=fopen(jsonPath,'w');
        assert(file>=0,'Cannot write result JSON');
        close=onCleanup(@()fclose(file));
        fprintf(file,'%s\n',jsonencode(result,PrettyPrint=true));
    end
end
