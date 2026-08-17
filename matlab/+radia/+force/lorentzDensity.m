function forceDensity_N_per_m3 = lorentzDensity( ...
    currentDensity_A_per_m2, magneticFluxDensity_T)
%LORENTZDENSITY Magnetostatic Lorentz force density J cross B in SI units.
%   Inputs are one three-vector or N-by-3 arrays. A one-row input is
%   expanded across the other input. The result is N-by-3 in N/m^3.

[currentDensity, magneticFluxDensity] = localBroadcastVectors( ...
    currentDensity_A_per_m2, "currentDensity_A_per_m2", ...
    magneticFluxDensity_T, "magneticFluxDensity_T");
forceDensity_N_per_m3 = cross(currentDensity, magneticFluxDensity, 2);
end

function [first, second] = localBroadcastVectors( ...
    firstValue, firstName, secondValue, secondName)
first = localVectors(firstValue, firstName);
second = localVectors(secondValue, secondName);
if size(first, 1) == 1 && size(second, 1) > 1
    first = repmat(first, size(second, 1), 1);
elseif size(second, 1) == 1 && size(first, 1) > 1
    second = repmat(second, size(first, 1), 1);
elseif size(first, 1) ~= size(second, 1)
    error("radia:force:SampleCount", ...
        "%s and %s must have the same row count or one row", ...
        firstName, secondName);
end
end

function vectors = localVectors(value, name)
if ~isnumeric(value) || ~isreal(value) || any(~isfinite(value), "all")
    error("radia:force:Vector", "%s must contain finite real values", name);
end
vectors = double(value);
if isvector(vectors) && numel(vectors) == 3
    vectors = reshape(vectors, 1, 3);
elseif ndims(vectors) ~= 2 || size(vectors, 2) ~= 3 %#ok<ISMAT>
    error("radia:force:VectorShape", "%s must be a three-vector or N-by-3", name);
end
end
