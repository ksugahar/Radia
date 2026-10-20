function result = p1Newton(functionName, positional, options)
% Explicit batch entry for radia.p1_newton; not a per-step backend.
% First-order (lowest-order Nedelec) nonlinear magnetostatic Newton for
% reduced or total A with a closed-form iron Jacobian and AMS/IC/direct
% inner solves. Batch only: the mesh, the A GridFunction and the field
% CoefficientFunctions stay on the Python side and the native AMS setup owns
% its own parallel regions. Python ecosystem objects remain in result.value;
% native MEX handles cannot be passed to this interface.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.p1_newton", ...
    functionName, positional, Keywords=options.Keywords);
end
