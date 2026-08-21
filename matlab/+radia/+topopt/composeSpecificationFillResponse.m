function response = composeSpecificationFillResponse( ...
        specificationFieldJacobian, elementFieldResponse, specificationRows)
%COMPOSESPECIFICATIONFILLRESPONSE Analytic/AD EarlyTimes-to-Abe junction.
arguments
    specificationFieldJacobian (:,:) double
    elementFieldResponse (:,:) double
    specificationRows (:,1) double {mustBeInteger,mustBePositive} = []
end
if size(specificationFieldJacobian,2) ~= size(elementFieldResponse,1)
    error("radia:topopt:SpecificationFillDimension", ...
        "The specification Jacobian and field response do not compose.");
end
if isempty(specificationRows)
    jacobian = specificationFieldJacobian;
else
    jacobian = specificationFieldJacobian(specificationRows,:);
end
response = jacobian * elementFieldResponse;
end
