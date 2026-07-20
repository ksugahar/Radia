function n = hluMaxThreads()
%HLUMAXTHREADS Return the available HACApK H-LU thread count.
n = radia.internal.callMex('hlu.max_threads');
end
