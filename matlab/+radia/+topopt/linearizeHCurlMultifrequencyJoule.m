function result=linearizeHCurlMultifrequencyJoule(operator,resistance,cellVertexVelocities,frequenciesHz,rhs,options)
%LINEARIZEHCURLMULTIFREQUENCYJOULE Native complex-adjoint Joule gradient.
arguments
    operator (1,1) radia.topopt.HCurlTopologyOperator
    resistance (1,1) struct
    cellVertexVelocities double {mustBeFinite}
    frequenciesHz (:,1) double {mustBePositive,mustBeFinite}
    rhs double
    options.Weights (:,1) double {mustBeNonnegative,mustBeFinite}=ones(numel(frequenciesHz),1)
    options.RHSJacobian double=double.empty
end
if ~all(isfield(resistance,["matrix","jacobian"]))
    error("radia:topopt:Resistance", ...
        "resistance must contain matrix and jacobian fields.");
end
result=radia.internal.callMex('hcurl.topopt.multifrequency_joule', ...
    operator.nativeHandle(),double(resistance.matrix), ...
    double(resistance.jacobian),double(cellVertexVelocities), ...
    frequenciesHz,double(rhs),options.Weights,double(options.RHSJacobian));
end
