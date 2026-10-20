function pressure_Pa = airGapMaxwellPressure( ...
    magneticFluxDensityNormal_T, permeability_H_per_m)
%AIRGAPMAXWELLPRESSURE Return Bn^2/(2*mu) in Pa.

if nargin < 2 || isempty(permeability_H_per_m)
    permeability_H_per_m = 4*pi*1e-7;
end
if ~isscalar(magneticFluxDensityNormal_T) || ~isreal(magneticFluxDensityNormal_T) || ~isfinite(magneticFluxDensityNormal_T)
    error("radia:force:Field", "magneticFluxDensityNormal_T must be one finite real scalar");
end
if ~isscalar(permeability_H_per_m) || ~isfinite(permeability_H_per_m) || permeability_H_per_m <= 0
    error("radia:force:Permeability", "permeability_H_per_m must be finite and positive");
end
pressure_Pa = magneticFluxDensityNormal_T^2 / (2*permeability_H_per_m);
end
