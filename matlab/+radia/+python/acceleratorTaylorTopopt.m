function result = acceleratorTaylorTopopt(functionName, positional, options)
%ACCELERATORTAYLORTOPOPT Call the Python high-order R/T/U topology API.
%   This explicit batch boundary exposes normal/skew multipole observation,
%   forward-AD Taylor-map objectives through octupole, local TSVD
%   reachability certificates, global KAN parameterization of decoupled
%   first-order symplectic targets, rank-6 continuation, and HDiv-MMM
%   optimization. It is not a Simulink step-time backend.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.accelerator_taylor_topopt",functionName,positional, ...
    Keywords=options.Keywords);
end
