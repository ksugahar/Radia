function forceDensity_N_per_m3 = timeAverageLorentzDensity( ...
    currentDensityPhasor_A_per_m2, magneticFluxDensityPhasor_T, amplitude)
%TIMEAVERAGELORENTZDENSITY Cycle-averaged Re(J cross conj(B)).
%   AMPLITUDE is "peak" (factor 1/2) or "rms" (factor 1).

if nargin < 3 || isempty(amplitude)
    amplitude = "peak";
end
[currentDensity, magneticFluxDensity] = localBroadcastVectors( ...
    currentDensityPhasor_A_per_m2, "currentDensityPhasor_A_per_m2", ...
    magneticFluxDensityPhasor_T, "magneticFluxDensityPhasor_T");
factor = localFactor(amplitude);
forceDensity_N_per_m3 = factor * real(cross( ...
    currentDensity, conj(magneticFluxDensity), 2));
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

function [first, second] = localBroadcastVectors(firstValue, firstName, secondValue, secondName)
first = localVectors(firstValue, firstName);
second = localVectors(secondValue, secondName);
if size(first, 1) == 1 && size(second, 1) > 1
    first = repmat(first, size(second, 1), 1);
elseif size(second, 1) == 1 && size(first, 1) > 1
    second = repmat(second, size(first, 1), 1);
elseif size(first, 1) ~= size(second, 1)
    error("radia:force:SampleCount", "%s and %s must have the same row count or one row", firstName, secondName);
end
end

function vectors = localVectors(value, name)
if ~isnumeric(value) || any(~isfinite(real(value)), "all") || any(~isfinite(imag(value)), "all")
    error("radia:force:Vector", "%s must contain finite numeric phasors", name);
end
vectors = double(value);
if isvector(vectors) && numel(vectors) == 3
    vectors = reshape(vectors, 1, 3);
elseif ndims(vectors) ~= 2 || size(vectors, 2) ~= 3 %#ok<ISMAT>
    error("radia:force:VectorShape", "%s must be a three-vector or N-by-3", name);
end
end
