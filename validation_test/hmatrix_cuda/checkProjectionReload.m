function result=checkProjectionReload(meshPath)
%CHECKPROJECTIONRELOAD Reject handles from a previous MEX load.
% Caller supplies a public .vol checked against its label contract.
arguments
    meshPath (1,:) char
end
active=uint64(0);
cleanup=onCleanup(@cleanupOwned);
assert(radia_cuda_projection('live_count')==0,'radia:ProjectionTest','Expected empty registry');
old=radia_cuda_projection('create',meshPath,'hdiv',1);
active=old;
radia_cuda_projection('destroy',old);active=uint64(0);
assert(radia_cuda_projection('live_count')==0,'radia:ProjectionTest','Destroy leaked context');
clear radia_cuda_projection
fresh=radia_cuda_projection('create',meshPath,'hdiv',1);
active=fresh;
assert(fresh~=old,'radia:ProjectionTest','Handle alias across MEX reload');
oldTargetsRejected=false;oldDestroyRejected=false;
try
    unused=radia_cuda_projection('targets',old); %#ok<NASGU>
catch ex
    oldTargetsRejected=strcmp(ex.identifier,'radia:cuda:projection') && contains(ex.message,'Stale or wrong projection handle');
end
try
    radia_cuda_projection('destroy',old);
catch ex
    oldDestroyRejected=strcmp(ex.identifier,'radia:cuda:projection') && contains(ex.message,'Stale or wrong projection handle');
end
assert(oldTargetsRejected && oldDestroyRejected,'radia:ProjectionTest','Stale handle accepted');
points=radia_cuda_projection('targets',fresh);
assert(size(points,2)==3 && ~isempty(points) && all(isfinite(points),'all'), ...
    'radia:ProjectionTest','Fresh context unusable after stale rejection');
clear cleanup
assert(radia_cuda_projection('live_count')==0,'radia:ProjectionTest','Cleanup leaked context');
result=struct('accepted',true,'old_targets_rejected',oldTargetsRejected, ...
    'old_destroy_rejected',oldDestroyRejected,'live_count',0, ...
    'scope','Checked handle lifecycle only; no GPU solver or performance acceptance');
    function cleanupOwned()
        if active~=0
            radia_cuda_projection('destroy',active);
            active=uint64(0);
        end
    end
end
