function result = acceleratorAbeTopopt(functionName, positional, options)
%ACCELERATORABETOPOPT Explicit Python boundary for the callback outer loop.
%   The dense element-fill inverse has a native MEX entry point at
%   radia.topopt.solveAbeElementFillPlan. Use this fallback only for NGSolve
%   mesh callbacks and complete-solve orchestration outside Simulink step time.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.accelerator_abe_topopt",functionName,positional, ...
    Keywords=options.Keywords);
end
