function stress_Pa = timeAverageMaxwellStressAir( ...
    magneticFluxDensityPhasor_T, permeability_H_per_m, amplitude)
%TIMEAVERAGEMAXWELLSTRESSAIR Cycle-averaged magnetic stress in air.
%   Uses Re(B tensor conj(B) - |B|^2 I/2) / mu, multiplied by 1/2 for
%   peak phasors and by 1 for RMS phasors.

if nargin < 2 || isempty(permeability_H_per_m)
    permeability_H_per_m = 4*pi*1e-7;
end
if nargin < 3 || isempty(amplitude)
    amplitude = "peak";
end
if ~isscalar(permeability_H_per_m) || ~isfinite(permeability_H_per_m) || permeability_H_per_m <= 0
    error("radia:force:Permeability", "permeability_H_per_m must be finite and positive");
end
field = localVectors(magneticFluxDensityPhasor_T);
factor = localFactor(amplitude);
count = size(field, 1);
stress_Pa = zeros(3, 3, count);
identity = eye(3);
for index = 1:count
    vector = field(index, :).';
    stress_Pa(:, :, index) = factor * ( ...
        real(vector * vector') ...
        - 0.5 * sum(abs(vector).^2) * identity) / permeability_H_per_m;
end
if count == 1
    stress_Pa = stress_Pa(:, :, 1);
end
end

function factor = localFactor(amplitude)
key = lower(strtrim(string(amplitude)));
if key == "peak"
    factor = 0.5;
elseif key == "rms"
    factor = 1.0;
else
    error("radia:force:PhasorAmplitude", "amplitude must be 'peak' or 'rms'");
end
end

function vectors = localVectors(value)
if ~isnumeric(value) || any(~isfinite(real(value)), "all") || any(~isfinite(imag(value)), "all")
    error("radia:force:Vector", "magneticFluxDensityPhasor_T must contain finite numeric phasors");
end
vectors = double(value);
if isvector(vectors) && numel(vectors) == 3
    vectors = reshape(vectors, 1, 3);
elseif ndims(vectors) ~= 2 || size(vectors, 2) ~= 3 %#ok<ISMAT>
    error("radia:force:VectorShape", "magneticFluxDensityPhasor_T must be a three-vector or N-by-3");
end
end
