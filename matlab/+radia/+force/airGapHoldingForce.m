function force_N = airGapHoldingForce( ...
    magneticFluxDensityNormal_T, activeArea_m2, faces, permeability_H_per_m)
%AIRGAPHOLDINGFORCE Uniform-gap Maxwell-pressure holding-force estimate.

if nargin < 3 || isempty(faces)
    faces = 1;
end
if nargin < 4
    permeability_H_per_m = [];
end
if ~isscalar(activeArea_m2) || ~isreal(activeArea_m2) || ~isfinite(activeArea_m2) || activeArea_m2 < 0
    error("radia:force:Area", "activeArea_m2 must be finite and nonnegative");
end
if ~isscalar(faces) || ~isreal(faces) || ~isfinite(faces) || faces < 1 || faces ~= fix(faces)
    error("radia:force:Faces", "faces must be a positive integer");
end
force_N = radia.force.airGapMaxwellPressure( ...
    magneticFluxDensityNormal_T, permeability_H_per_m) * activeArea_m2 * faces;
end
