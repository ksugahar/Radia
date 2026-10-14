function stress_Pa = airGapShearStress( ...
    magneticFluxDensityRadial_T, magneticFluxDensityTangential_T, permeability_H_per_m)
%AIRGAPSHEARSTRESS Static cylindrical air-gap shear Br*Bt/mu in Pa.

if nargin < 3 || isempty(permeability_H_per_m)
    permeability_H_per_m = 4*pi*1e-7;
end
if ~isscalar(magneticFluxDensityRadial_T) || ~isreal(magneticFluxDensityRadial_T) || ~isfinite(magneticFluxDensityRadial_T) || ...
        ~isscalar(magneticFluxDensityTangential_T) || ~isreal(magneticFluxDensityTangential_T) || ~isfinite(magneticFluxDensityTangential_T)
    error("radia:force:Field", "air-gap flux-density components must be finite real scalars");
end
if ~isscalar(permeability_H_per_m) || ~isfinite(permeability_H_per_m) || permeability_H_per_m <= 0
    error("radia:force:Permeability", "permeability_H_per_m must be finite and positive");
end
stress_Pa = magneticFluxDensityRadial_T * magneticFluxDensityTangential_T / permeability_H_per_m;
end
