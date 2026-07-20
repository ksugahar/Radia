function [force_N, upwardLift_N, forceSlope_N_per_m] = ...
    evaluateTeam28CLNForce(lut, heightOffset_m, coilCurrent_A)
%EVALUATETEAM28CLNFORCE Evaluate the validated TEAM 28 CLN force LUT.
%   The returned force follows the repository convention: negative is
%   upward.  UPWARDLIFT_N is positive for upward lift.  Current scaling is
%   quadratic relative to lut.reference_coil_current_A; the LUT frequency is
%   fixed by the source artifact and is not extrapolated here.

arguments
    lut (1,1) struct
    heightOffset_m double {mustBeFinite}
    coilCurrent_A double {mustBeFinite}
end

validateLUT(lut);
if ~(isscalar(coilCurrent_A) || isequal(size(coilCurrent_A), size(heightOffset_m)))
    error("radia:simulink:Team28InputSize", ...
        "coilCurrent_A must be scalar or match heightOffset_m.");
end

height = double(heightOffset_m);
grid = lut.height_offset_m(:);
if lut.extrapolation == "clip"
    query = min(max(height, grid(1)), grid(end));
elseif any(height(:) < grid(1) | height(:) > grid(end))
    error("radia:simulink:Team28Extrapolation", ...
        "heightOffset_m is outside the validated TEAM 28 CLN range.");
else
    query = height;
end

forceBase = interp1(grid, lut.force_N(:), query, "pchip");
slopeBase = interp1(grid, lut.force_slope_N_per_m(:), query, "linear");
scale = (double(coilCurrent_A) / lut.reference_coil_current_A).^2;
force_N = forceBase .* scale;
upwardLift_N = -force_N;
forceSlope_N_per_m = slopeBase .* scale;
end

function validateLUT(lut)
if ~isfield(lut, "schema") || lut.schema ~= "radia.team28.cln_lut.v1"
    error("radia:simulink:Team28LUT", "unsupported TEAM 28 CLN LUT schema.");
end
required = ["height_offset_m", "force_N", "force_slope_N_per_m", ...
    "reference_coil_current_A", "extrapolation"];
if ~all(isfield(lut, cellstr(required)))
    error("radia:simulink:Team28LUT", "TEAM 28 CLN LUT fields are incomplete.");
end
end
