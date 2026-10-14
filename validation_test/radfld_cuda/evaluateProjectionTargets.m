function [values,cache,info] = evaluateProjectionTargets(points,evaluate,revision,cache)
%EVALUATEPROJECTIONTARGETS Exact-coordinate batching for projection experiments.
% points must be the actual NGSolve mapped targets, not replacement nodes.
% evaluate(uniquePoints) completes evaluation and all transfers before returning.
% revision must change with source values, geometry, transforms or evaluator options.
arguments
    points (:,3) double {mustBeReal,mustBeFinite}
    evaluate (1,1) function_handle
    revision (1,1) string
    cache struct = struct()
end
assert(strlength(revision)>0,'radia:projection:Revision','Explicit revision required');
timer=tic;
reuse=isfield(cache,'revision') && isequal(cache.revision,revision) && ...
    isfield(cache,'points') && isequal(cache.points,points);
if reuse
    values=cache.values;
    info=struct('reused',true,'requested_points',size(points,1), ...
        'evaluated_points',0,'whole_evaluation_seconds',toc(timer));
    return
end
% Exact equality only: no spatial quantization or tolerance-based merging.
[uniquePoints,~,inverse]=unique(points,'rows');
uniqueValues=evaluate(uniquePoints);
assert(isa(uniqueValues,'double') && isreal(uniqueValues) && ...
    isequal(size(uniqueValues),[size(uniquePoints,1),3]) && ...
    all(isfinite(uniqueValues),'all'),'radia:projection:Values', ...
    'Evaluator must return completed finite real double Nx3 field values');
values=uniqueValues(inverse,:);
cache=struct('revision',revision,'points',points,'values',values);
info=struct('reused',false,'requested_points',size(points,1), ...
    'evaluated_points',size(uniquePoints,1),'whole_evaluation_seconds',toc(timer));
end
