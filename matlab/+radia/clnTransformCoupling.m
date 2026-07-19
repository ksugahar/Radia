function reduced = clnTransformCoupling(Q, MLS)
%CLNTRANSFORMCOUPLING Transform loop-star coupling into reduced loop space.
reduced = radia.internal.callMex( ...
    'cln.transform_coupling', double(Q), double(MLS));
end
