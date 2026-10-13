function [vector,cache,info] = projectBatchedField(token,evaluate,revision,cache)
%PROJECTBATCHEDFIELD Replay exact NGSolve targets after bulk MATLAB/CUDA field evaluation.
% token belongs to radia_cuda_projection, never another MEX handle registry.
arguments
    token (1,1) uint64
    evaluate (1,1) function_handle
    revision (1,1) string
    cache struct = struct()
end
timer=tic;
points=radia_cuda_projection('targets',token);
[values,cache,evaluation]=evaluateProjectionTargets(points,evaluate,revision,cache);
[vector,nativeSeconds]=radia_cuda_projection('project',token,values);
info=struct('whole_projection_seconds',toc(timer), ...
    'native_interpolation_seconds',nativeSeconds,'evaluation',evaluation, ...
    'source_revision',revision,'target_count',size(points,1), ...
    'device_peak_bytes',[],'memory_limitations','Device peak not instrumented');
end
