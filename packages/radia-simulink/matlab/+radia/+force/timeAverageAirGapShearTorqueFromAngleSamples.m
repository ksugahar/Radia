function summary = timeAverageAirGapShearTorqueFromAngleSamples( ...
    angles_rad, magneticFluxDensityRadialPhasor_T, ...
    magneticFluxDensityTangentialPhasor_T, radius_m, axialLength_m, ...
    periodic, period_rad, permeability_H_per_m, amplitude)
%TIMEAVERAGEAIRGAPSHEARTORQUEFROMANGLESAMPLES Integrate phasor gap torque.

if nargin < 5 || isempty(axialLength_m), axialLength_m = 1.0; end
if nargin < 6 || isempty(periodic), periodic = true; end
if nargin < 7 || isempty(period_rad), period_rad = 2*pi; end
if nargin < 8 || isempty(permeability_H_per_m), permeability_H_per_m = 4*pi*1e-7; end
if nargin < 9 || isempty(amplitude), amplitude = "peak"; end
angles = localRealTable(angles_rad, "angles_rad");
radial = localComplexTable(magneticFluxDensityRadialPhasor_T, "magneticFluxDensityRadialPhasor_T");
tangential = localComplexTable(magneticFluxDensityTangentialPhasor_T, "magneticFluxDensityTangentialPhasor_T");
if numel(angles) ~= numel(radial) || numel(angles) ~= numel(tangential)
    error("radia:force:SampleCount", "angles and air-gap fields must have the same length");
end
key = lower(strtrim(string(amplitude)));
if key == "peak"
    factor = 0.5;
elseif key == "rms"
    factor = 1.0;
else
    error("radia:force:PhasorAmplitude", "amplitude must be peak or rms");
end
summary = radia.force.airGapShearTorqueFromAngleSamples( ...
    angles, ones(size(radial)), factor * real(radial .* conj(tangential)), ...
    radius_m, axialLength_m, periodic, period_rad, permeability_H_per_m);
summary.phasor_amplitude = key;
end

function values = localRealTable(value, name)
if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || numel(value) < 2 || any(~isfinite(value), "all")
    error("radia:force:Table", "%s must be a finite real vector with at least two samples", name);
end
values = double(value(:));
if name == "angles_rad" && any(diff(values) <= 0)
    error("radia:force:Angles", "angles_rad must be strictly increasing");
end
end

function values = localComplexTable(value, name)
if ~isnumeric(value) || ~isvector(value) || numel(value) < 2 || ...
        any(~isfinite(real(value)), "all") || any(~isfinite(imag(value)), "all")
    error("radia:force:Table", "%s must be a finite phasor vector with at least two samples", name);
end
values = double(value(:));
end
