function stress_Pa = maxwellStressAir(magneticFluxDensity_T, permeability_H_per_m)
%MAXWELLSTRESSAIR Static magnetic Maxwell stress tensor in air.
%   A single B vector returns a 3-by-3 tensor. N-by-3 samples return a
%   3-by-3-by-N array. The surface used for force integration must lie in air.

arguments
    magneticFluxDensity_T {mustBeNumeric,mustBeReal,mustBeFinite}
    permeability_H_per_m (1,1) double {mustBeFinite,mustBePositive} = 4*pi*1e-7
end

field = localVectors(magneticFluxDensity_T);
count = size(field, 1);
stress_Pa = zeros(3, 3, count);
identity = eye(3);
for index = 1:count
    vector = field(index, :).';
    stress_Pa(:, :, index) = (vector * vector.' ...
        - 0.5 * dot(vector, vector) * identity) / permeability_H_per_m;
end
if count == 1
    stress_Pa = stress_Pa(:, :, 1);
end
end

function vectors = localVectors(value)
vectors = double(value);
if isvector(vectors) && numel(vectors) == 3
    vectors = reshape(vectors, 1, 3);
elseif ndims(vectors) ~= 2 || size(vectors, 2) ~= 3 %#ok<ISMAT>
    error("radia:force:VectorShape", ...
        "magneticFluxDensity_T must be a three-vector or N-by-3");
end
end
