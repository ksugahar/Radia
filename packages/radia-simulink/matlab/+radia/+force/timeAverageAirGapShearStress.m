function stress_Pa = timeAverageAirGapShearStress( ...
    magneticFluxDensityRadialPhasor_T, magneticFluxDensityTangentialPhasor_T, ...
    permeability_H_per_m, amplitude)
%TIMEAVERAGEAIRGAPSHEARSTRESS Cycle-averaged Re(Br*conj(Bt))/mu.

if nargin < 3 || isempty(permeability_H_per_m)
    permeability_H_per_m = 4*pi*1e-7;
end
if nargin < 4 || isempty(amplitude)
    amplitude = "peak";
end
values = [magneticFluxDensityRadialPhasor_T, magneticFluxDensityTangentialPhasor_T];
if ~isnumeric(values) || numel(values) ~= 2 || any(~isfinite(real(values))) || any(~isfinite(imag(values)))
    error("radia:force:Field", "air-gap flux-density phasors must be finite scalars");
end
if ~isscalar(permeability_H_per_m) || ~isfinite(permeability_H_per_m) || permeability_H_per_m <= 0
    error("radia:force:Permeability", "permeability_H_per_m must be finite and positive");
end
key = lower(strtrim(string(amplitude)));
if key == "peak"
    factor = 0.5;
elseif key == "rms"
    factor = 1.0;
else
    error("radia:force:PhasorAmplitude", "amplitude must be 'peak' or 'rms'");
end
stress_Pa = factor * real(magneticFluxDensityRadialPhasor_T * ...
    conj(magneticFluxDensityTangentialPhasor_T)) / permeability_H_per_m;
end
