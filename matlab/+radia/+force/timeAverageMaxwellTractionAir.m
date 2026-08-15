function traction_Pa = timeAverageMaxwellTractionAir( ...
    magneticFluxDensityPhasor_T, outwardNormal, permeability_H_per_m, amplitude)
%TIMEAVERAGEMAXWELLTRACTIONAIR Cycle-averaged air-side traction <T>*n.

if nargin < 3 || isempty(permeability_H_per_m)
    permeability_H_per_m = 4*pi*1e-7;
end
if nargin < 4 || isempty(amplitude)
    amplitude = "peak";
end
[field, normals] = localBroadcastVectors(magneticFluxDensityPhasor_T, outwardNormal);
normalMagnitude = vecnorm(normals, 2, 2);
if any(normalMagnitude == 0)
    error("radia:force:Normal", "outwardNormal must be nonzero");
end
normals = normals ./ normalMagnitude;
stress = radia.force.timeAverageMaxwellStressAir( ...
    field, permeability_H_per_m, amplitude);
traction_Pa = zeros(size(field, 1), 3);
if size(field, 1) == 1
    traction_Pa(1, :) = (stress * normals(1, :).').';
else
    for index = 1:size(field, 1)
        traction_Pa(index, :) = (stress(:, :, index) * normals(index, :).').';
    end
end
end

function [field, normals] = localBroadcastVectors(fieldValue, normalValue)
if ~isnumeric(fieldValue) || any(~isfinite(real(fieldValue)), "all") || any(~isfinite(imag(fieldValue)), "all")
    error("radia:force:Vector", "magneticFluxDensityPhasor_T must contain finite numeric phasors");
end
field = localShape(double(fieldValue), "magneticFluxDensityPhasor_T");
if ~isnumeric(normalValue) || ~isreal(normalValue) || any(~isfinite(normalValue), "all")
    error("radia:force:Vector", "outwardNormal must contain finite real values");
end
normals = localShape(double(normalValue), "outwardNormal");
if size(field, 1) == 1 && size(normals, 1) > 1
    field = repmat(field, size(normals, 1), 1);
elseif size(normals, 1) == 1 && size(field, 1) > 1
    normals = repmat(normals, size(field, 1), 1);
elseif size(field, 1) ~= size(normals, 1)
    error("radia:force:SampleCount", "field and normal must have the same row count or one row");
end
end

function vectors = localShape(value, name)
vectors = value;
if isvector(vectors) && numel(vectors) == 3
    vectors = reshape(vectors, 1, 3);
elseif ndims(vectors) ~= 2 || size(vectors, 2) ~= 3 %#ok<ISMAT>
    error("radia:force:VectorShape", "%s must be a three-vector or N-by-3", name);
end
end
