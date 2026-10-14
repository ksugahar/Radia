function result = acceleratorLieTopopt(functionName, positional, options)
%ACCELERATORLIETOPOPT Call the canonical high-order Lie topology API.
%   This explicit batch boundary exposes the native fifth-degree body-field
%   Hamiltonian jet, Dragt-Finn f3/f4/f5 factorization and R/T/U/V map,
%   RK-orbit moving-frame bridge, finite-amplitude symplectic map application,
%   forward-AD derivatives, full HCurl transverse A-jet Lie maps, sampled-field
%   objective gradients, reference-orbit constraint Jacobians, reachability
%   certificates, and HDiv-MMM topology optimization. It is not a Simulink
%   step-time backend.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.accelerator_lie_topopt",functionName,positional, ...
    Keywords=options.Keywords);
end
