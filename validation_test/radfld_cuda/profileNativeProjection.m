function result = profileNativeProjection(gf,coefficient,reference,outputPath)
%PROFILENATIVEPROJECTION Profile actual NGSolve interpolation from MATLAB.
% Caller owns objects and provides an independent direct reference GF vector.
% Run only in a job-owned environment. This does not select an accelerator.
arguments
    gf (1,1) radia.ngsolve.GridFunction
    coefficient (1,1) radia.ngsolve.CoefficientFunction
    reference (:,1) double {mustBeFinite}
    outputPath (1,1) string
end
result=struct('schema','radia.native-projection-profile.v1', ...
    'accepted',false,'scope','NGSolve interpolation through MATLAB native API', ...
    'relative_limit',1e-4,'timing_acceptance',false, ...
    'timing_class','diagnostic only; caller must verify idle host', ...
    'a100_used',false,'matlab',version,'dof_count',gf.DofCount, ...
    'cf_profile_status','native CF call counters not connected', ...
    'memory_limitations','Host snapshots only; device peak not measured');
try
    assert(numel(reference)==gf.DofCount,'radia:projection:ReferenceSize', ...
        'Reference must use the identical NGSolve space and DoF ordering');
    [before,system]=memory;
    timer=tic;
    gf.interpolate(coefficient);
    actual=gf.vector(); % Include result recovery in the continuous interval.
    result.cold_projection_seconds=toc(timer);
    result.host_available_before_bytes=system.PhysicalMemory.Available;
    result.host_used_before_bytes=before.MemUsedMATLAB;
    result.samples_seconds=zeros(1,3);
    for repeat=1:3
        timer=tic; gf.interpolate(coefficient); actual=gf.vector();
        result.samples_seconds(repeat)=toc(timer);
    end
    [after,system]=memory;
    result.host_used_after_bytes=after.MemUsedMATLAB;
    result.host_available_after_bytes=system.PhysicalMemory.Available;
    assert(all(isfinite(actual),'all'),'radia:projection:Nonfinite', ...
        'Native projection returned nonfinite values');
    result.coefficient_vector_absolute_error=norm(actual(:)-reference);
    denominator=norm(reference);
    if denominator==0
        assert(result.coefficient_vector_absolute_error==0, ...
            'radia:projection:ZeroReference','Nonzero result for zero reference');
        result.coefficient_vector_relative_error=0;
    else
        result.coefficient_vector_relative_error= ...
            result.coefficient_vector_absolute_error/denominator;
    end
    assert(result.coefficient_vector_relative_error<=1e-4, ...
        'radia:projection:Parity','Projection coefficient-vector mismatch');
    % A coefficient norm is NOT the physical FE L2/field acceptance gate.
    result.algebra_check_passed=true;
    result.physical_acceptance_pending=true;
    result.pending_reason='Mapped-point batching and FE L2/field validation pending';
catch failure
    result.error_identifier=failure.identifier;result.error_message=failure.message;
    writeResult(outputPath,result);
    rethrow(failure)
end
writeResult(outputPath,result);
end
function writeResult(path,value)
f=fopen(path,'w');
assert(f~=-1,'radia:projection:Output','Cannot open result JSON');
cleanup=onCleanup(@() fclose(f));
fprintf(f,'%s\n',jsonencode(value,PrettyPrint=true));
end
