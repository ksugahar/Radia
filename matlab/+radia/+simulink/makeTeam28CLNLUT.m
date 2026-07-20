function lut = makeTeam28CLNLUT(source, options)
%MAKETEAM28CLNLUT Load the validated TEAM 28 six-stage CLN force curve.
%   LUT = radia.simulink.makeTeam28CLNLUT() loads the canonical result
%   artifact under docs/maglev/demos/team28.  The stored TEAM 28 integral is
%   twice the physical time-average Lorentz force, so force_N is corrected by
%   PhysicalForceFactor (default 0.5).  HeightOffset_m is dZ relative to the
%   10.8 mm reference disk position.

arguments
    source = ""
    options.PhysicalForceFactor (1,1) double {mustBePositive} = 0.5
    options.ReferenceCoilCurrent_A (1,1) double {mustBePositive} = 20.0
    options.Extrapolation (1,1) string = "clip"
end

if ~ismember(options.Extrapolation, ["clip", "error"])
    error("radia:simulink:Team28Extrapolation", ...
        "Extrapolation must be 'clip' or 'error'.");
end

if isstruct(source)
    data = source;
    sourceFile = "";
else
    source = string(source);
    if strlength(source) == 0
        repoRoot = fileparts(fileparts(fileparts(fileparts(mfilename("fullpath")))));
        source = fullfile(repoRoot, "docs", "maglev", "demos", "team28", ...
            "team28_cln_sweep_results.json");
    end
    if ~isfile(source)
        error("radia:simulink:Team28Source", ...
            "TEAM 28 CLN result file does not exist: %s", source);
    end
    data = jsondecode(fileread(source));
    sourceFile = source;
end

required = ["dZ_mm", "fz_cln_N", "cln_stages", "geometry"];
if ~all(isfield(data, cellstr(required)))
    error("radia:simulink:Team28Source", ...
        "TEAM 28 source must contain dZ_mm, fz_cln_N, cln_stages, and geometry.");
end
if ~isfield(data.geometry, "freq_Hz")
    error("radia:simulink:Team28Source", ...
        "TEAM 28 source geometry must contain freq_Hz.");
end

heightOffset = double(data.dZ_mm(:)) * 1.0e-3;
storedForce = double(data.fz_cln_N(:));
if numel(heightOffset) < 2 || numel(storedForce) ~= numel(heightOffset) || ...
        any(~isfinite(heightOffset)) || any(~isfinite(storedForce))
    error("radia:simulink:Team28Source", ...
        "TEAM 28 height and force vectors must be finite and have matching sizes.");
end
[heightOffset, order] = sort(heightOffset);
storedForce = storedForce(order);
if any(diff(heightOffset) <= 0)
    error("radia:simulink:Team28Source", ...
        "TEAM 28 heights must be strictly increasing.");
end

force = options.PhysicalForceFactor * storedForce;
forceSlope = gradient(force, heightOffset);
lut = struct( ...
    "schema", "radia.team28.cln_lut.v1", ...
    "height_offset_m", heightOffset, ...
    "force_N", force, ...
    "upward_lift_N", -force, ...
    "force_slope_N_per_m", forceSlope, ...
    "frequency_Hz", double(data.geometry.freq_Hz), ...
    "cln_stages", double(data.cln_stages), ...
    "disk_weight_N", double(data.disk_weight_N), ...
    "reference_coil_current_A", options.ReferenceCoilCurrent_A, ...
    "physical_force_factor", options.PhysicalForceFactor, ...
    "extrapolation", options.Extrapolation, ...
    "source_file", string(sourceFile), ...
    "model", "axisymmetric TEAM 28 6-stage CLN force curve");
end
