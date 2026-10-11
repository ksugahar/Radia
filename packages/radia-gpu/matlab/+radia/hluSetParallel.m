function hluSetParallel(on)
%HLUSETPARALLEL Enable or disable HACApK H-LU block parallelism.
radia.internal.callMex('hlu.set_parallel', logical(on));
end
