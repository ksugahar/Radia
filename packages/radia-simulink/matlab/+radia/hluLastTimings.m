function stats = hluLastTimings()
%HLULASTTIMINGS Return timing counters from the most recent H-LU operation.
stats = radia.internal.callMex('hlu.last_timings');
end
