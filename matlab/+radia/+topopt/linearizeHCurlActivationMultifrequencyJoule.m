function result=linearizeHCurlActivationMultifrequencyJoule(operator,cellCurlGrams,activation,frequenciesHz,rhs,conductivity,options)
%LINEARIZEHCURLACTIVATIONMULTIFREQUENCYJOULE Cellwise SIMP conductivity adjoint.
arguments
    operator (1,1) radia.topopt.HCurlTopologyOperator
    cellCurlGrams double {mustBeFinite}
    activation (:,1) double {mustBeBetween(activation,0,1)}
    frequenciesHz (:,1) double {mustBePositive,mustBeFinite}
    rhs double
    conductivity (1,1) struct
    options.Weights (:,1) double {mustBeNonnegative,mustBeFinite}=ones(numel(frequenciesHz),1)
    options.RHSJacobian double=double.empty
    options.InductancePower (1,1) double {mustBeGreaterThanOrEqual(options.InductancePower,1)}=1
end
if ~all(isfield(conductivity,["solid","void"]))
    error("radia:topopt:Conductivity", ...
        "conductivity must contain solid and void fields.");
end
if isfield(conductivity,"power"), power=conductivity.power; else, power=3; end
result=radia.internal.callMex( ...
    'hcurl.topopt.activation_multifrequency_joule',operator.nativeHandle(), ...
    double(cellCurlGrams),activation,frequenciesHz,double(rhs), ...
    options.Weights,double(options.RHSJacobian),double(conductivity.solid), ...
    double(conductivity.void),double(power),options.InductancePower);
end
