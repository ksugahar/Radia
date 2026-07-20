function model = interpolateHCurlEddyCLNFamily(family, height_m, options)
%INTERPOLATEHCURLEDDYCLNFAMILY Interpolate one common-basis HCurl snapshot.
%   The default linear interpolation is deliberately explicit.  It blends
%   the state-space matrices directly for fast Simulink use and also blends
%   R/L/P for harmonic diagnostics.  No extrapolation is permitted by
%   default; clamp is an opt-in control-oriented policy.

arguments
    family (1,1) struct
    height_m (1,1) double {mustBeFinite}
    options.Interpolation (1,1) string = ""
    options.Extrapolation (1,1) string = ""
    options.BuildStateSpace (1,1) logical = false
end

if ~isfield(family, "schema") || family.schema ~= ...
        "radia.hcurl.eddy_cln.family.v1"
    error("radia:simulink:HCurlCLNFamily", "unsupported family struct.");
end
if ~isfield(family, "shared_state_basis") || ~family.shared_state_basis
    error("radia:simulink:HCurlCLNFamily", ...
        "family must declare a shared reduced state basis.");
end
method = options.Interpolation;
if strlength(method) == 0
    method = family.interpolation;
end
extrapolation = options.Extrapolation;
if strlength(extrapolation) == 0
    extrapolation = family.extrapolation;
end
if ~ismember(method, ["nearest", "linear", "pchip"]) || ...
        ~ismember(extrapolation, ["error", "clamp"])
    error("radia:simulink:HCurlCLNFamily", "unsupported interpolation policy.");
end

positions = double(family.positions_m(:));
count = numel(positions);
if count < 1 || any(~isfinite(positions)) || any(diff(positions) <= 0)
    error("radia:simulink:HCurlCLNFamily", ...
        "family positions must be finite and strictly increasing.");
end
query = height_m;
if query < positions(1) || query > positions(end)
    if extrapolation == "error"
        error("radia:simulink:HCurlCLNExtrapolation", ...
            "height_m is outside the HCurl family range.");
    end
    query = min(max(query, positions(1)), positions(end));
end

if method == "nearest" || count == 1
    [~, index] = min(abs(positions - query));
    model = family.models{index};
    model.height_m = query;
    model.family_source_height_m = positions(index);
    model.interpolation = "nearest";
    return
end

if method == "linear"
    upper = find(positions >= query, 1, "first");
    if upper == 1
        lower = 1;
        upper = 1;
        alpha = 0;
    elseif upper > count
        lower = count;
        upper = count;
        alpha = 0;
    else
        lower = upper - 1;
        alpha = (query - positions(lower)) / (positions(upper) - positions(lower));
    end
    blend = @(name) localLinearField(family.models{lower}, family.models{upper}, name, alpha);
else
    blend = @(name) localPchipField(family.models, positions, query, name);
    lower = NaN;
    upper = NaN;
    alpha = NaN;
end

base = family.models{1};
model = base;
model.height_m = query;
if method == "pchip"
    model.family_source_height_m = [positions(1), positions(end)];
else
    model.family_source_height_m = [positions(lower), positions(upper)];
end
model.interpolation = method;
model.resistance = blend("resistance");
model.inductance = blend("inductance");
model.port_rhs = blend("port_rhs");
model.surface_mass = zeros(size(model.resistance));
model.A = blend("A");
model.B = blend("B");
model.C = blend("C");
model.D = blend("D");
model.Ad = blend("Ad");
model.Bd = blend("Bd");
model.Cd = blend("Cd");
model.Dd = blend("Dd");
model.force_operator = blend("force_operator");
if options.BuildStateSpace
    rebuilt = radia.simulink.makeHCurlEddyCLNModel( ...
        model.resistance, model.inductance, model.port_rhs, ...
        SampleTime_s=family.sample_time_s);
    rebuilt.height_m = query;
    rebuilt.family_source_height_m = model.family_source_height_m;
    rebuilt.interpolation = method;
    rebuilt.force_operator = model.force_operator;
    model = rebuilt;
end
end

function value = localLinearField(left, right, fieldName, alpha)
value = (1 - alpha) * double(left.(fieldName)) + alpha * double(right.(fieldName));
end

function value = localPchipField(models, positions, query, fieldName)
first = double(models{1}.(fieldName));
shape = size(first);
table = zeros([numel(positions), numel(first)]);
for index = 1:numel(models)
    table(index, :) = double(models{index}.(fieldName)(:)).';
end
value = reshape(interp1(positions, table, query, "pchip"), shape);
end
