function result = eddyAxisymRing(functionName, positional, options)
% Explicit batch entry for radia.eddy_axisym_ring; not a per-step backend.
% Voltage-driven axisymmetric eddy solves on the axifem operators. Loss is
% reported through radia.axisym_measure's revolved volume measure, so a
% MATLAB caller gets the same measure the Python side names.
% Python ecosystem objects remain in result.value; native MEX handles cannot
% be passed to this interface.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.eddy_axisym_ring", ...
    functionName, positional, Keywords=options.Keywords);
end