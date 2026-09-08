function result = AMS(matrix, space, options)
%AMS Native coordinate-based AMS for a real, order-one HCurl auxiliary matrix.
%   Complex=true returns a complex preconditioner for a matching eddy matrix.
%   No Python objects or matrix copies cross this boundary. NoGrads must be true.
arguments
    matrix (1,1) radia.ngsolve.Matrix
    space (1,1) radia.ngsolve.FESpace
    options.Complex (1,1) logical = false
    options.Cycle (1,1) double {mustBeMember(options.Cycle,[1 7])} = 1
    options.SmoothingSteps (1,1) double {mustBeInteger,mustBePositive} = 1
end
h = radia.internal.callMex('sparsesolv.ams', matrix.nativeHandle(), ...
    space.nativeHandle(), options.Complex, options.Cycle, options.SmoothingSteps);
result = radia.ngsolve.Matrix.fromNativeHandle(h);
end
