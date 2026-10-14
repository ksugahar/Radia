function on = hluGetParallel()
%HLUGETPARALLEL Return whether HACApK H-LU block parallelism is enabled.
on = radia.internal.callMex('hlu.get_parallel');
end
