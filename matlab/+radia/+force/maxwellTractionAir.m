function traction_Pa = maxwellTractionAir(magneticFluxDensity_T, ...
    outwardNormal, permeability_H_per_m)
%MAXWELLTRACTIONAIR Air-side Maxwell traction T*n in pascals.
%   Inputs are one three-vector or N-by-3 arrays. Normals are normalized
%   sample by sample and must point out of the body whose force is requested.

arguments
    magneticFluxDensity_T {mustBeNumeric,mustBeReal,mustBeFinite}
    outwardNormal {mustBeNumeric,mustBeReal,mustBeFinite}
    permeability_H_per_m (1,1) double {mustBeFinite,mustBePositive} = 4*pi*1e-7
end

[field, normals] = localBroadcastVectors( ...
    magneticFluxDensity_T, outwardNormal);
normalMagnitude = vecnorm(normals, 2, 2);
if any(normalMagnitude == 0)
    error("radia:force:Normal", "outwardNormal must be nonzero");
end
normals = normals ./ normalMagnitude;
normalField = sum(field .* normals, 2);
fieldSquared = sum(field.^2, 2);
traction_Pa = (field .* normalField ...
    - 0.5 * fieldSquared .* normals) / permeability_H_per_m;
end

function [first, second] = localBroadcastVectors(firstValue, secondValue)
first = localVectors(firstValue, "magneticFluxDensity_T");
second = localVectors(secondValue, "outwardNormal");
if size(first, 1) == 1 && size(second, 1) > 1
    first = repmat(first, size(second, 1), 1);
elseif size(second, 1) == 1 && size(first, 1) > 1
    second = repmat(second, size(first, 1), 1);
elseif size(first, 1) ~= size(second, 1)
    error("radia:force:SampleCount", ...
        "field and normal must have the same row count or one row");
end
end

function vectors = localVectors(value, name)
vectors = double(value);
if isvector(vectors) && numel(vectors) == 3
    vectors = reshape(vectors, 1, 3);
elseif ndims(vectors) ~= 2 || size(vectors, 2) ~= 3 %#ok<ISMAT>
    error("radia:force:VectorShape", "%s must be a three-vector or N-by-3", name);
end
end
