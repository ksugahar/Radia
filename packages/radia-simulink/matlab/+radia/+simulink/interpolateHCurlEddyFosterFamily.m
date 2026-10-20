function model = interpolateHCurlEddyFosterFamily(family, height_m, options)
%INTERPOLATEHCURLEDDYFOSTERFAMILY Foster model of a family at one height.
%   The decay rates and discrete state matrix are shared by every height;
%   only the modal input matrix and the modal force operator are
%   interpolated.  No extrapolation is permitted by default; clamp is an
%   opt-in control-oriented policy.

arguments
    family (1,1) struct
    height_m (1,1) double {mustBeFinite}
    options.Interpolation (1,1) string = ""
    options.Extrapolation (1,1) string = ""
end

if ~isfield(family, "schema") || family.schema ~= "radia.hcurl.eddy_foster.family.v1" || ...
        ~isfield(family, "shared_modes") || ~family.shared_modes
    error("radia:simulink:HCurlFosterFamily", ...
        "family must come from loadHCurlEddyFosterFamily or makeHCurlEddyFosterFamily.");
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
    error("radia:simulink:HCurlFosterFamily", "unsupported interpolation policy.");
end

positions = family.positions_m;
query = height_m;
if query < positions(1) || query > positions(end)
    if extrapolation == "error"
        error("radia:simulink:HCurlFosterExtrapolation", ...
            "height_m is outside the Foster family range.");
    end
    query = min(max(query, positions(1)), positions(end));
end

[Bm, K, source] = localOperators(family, positions, query, method);
C = Bm.';
model = struct( ...
    "schema", "radia.hcurl.eddy_foster.state_space.v1", ...
    "decay_rates", family.decay_rates, "modal_port_rhs", Bm, ...
    "A", -diag(family.decay_rates), "B", Bm, "C", C, ...
    "D", zeros(family.port_count), ...
    "Ad", family.Ad, "Bd", family.bd_gain .* Bm, "Cd", C, ...
    "Dd", zeros(family.port_count), ...
    "force_operator", K, ...
    "x0", family.x0, "sample_time_s", family.sample_time_s, ...
    "state_order", family.state_order, "port_count", family.port_count, ...
    "passive", true, ...
    "input_convention", "u=-d(coil_current)/dt", ...
    "state_convention", "z_dot=-diag(decay_rates)*z+modal_port_rhs*u, y=modal_port_rhs.'*z", ...
    "has_sibc_termination", false, ...
    "height_m", query, ...
    "family_source_height_m", source, ...
    "interpolation", method);
end

function [Bm, K, source] = localOperators(family, positions, query, method)
count = numel(positions);
if method == "nearest" || count == 1
    [~, index] = min(abs(positions - query));
    Bm = family.modal_port_rhs(:, :, index);
    K = family.force_operator(:, :, :, index);
    source = positions(index);
    return
end
if method == "linear"
    upper = find(positions >= query, 1, "first");
    if upper == 1
        lower = 1;
        alpha = 0.0;
    else
        lower = upper - 1;
        alpha = (query - positions(lower)) / (positions(upper) - positions(lower));
    end
    Bm = (1 - alpha) * family.modal_port_rhs(:, :, lower) + ...
        alpha * family.modal_port_rhs(:, :, upper);
    K = (1 - alpha) * family.force_operator(:, :, :, lower) + ...
        alpha * family.force_operator(:, :, :, upper);
    source = [positions(lower), positions(upper)];
    return
end
nState = family.state_order;
nPort = family.port_count;
portTable = reshape(family.modal_port_rhs, nState * nPort, count).';
forceTable = reshape(family.force_operator, 3 * nState * nPort, count).';
Bm = reshape(interp1(positions, portTable, query, "pchip"), nState, nPort);
K = reshape(interp1(positions, forceTable, query, "pchip"), 3, nState, nPort);
source = [positions(1), positions(end)];
end
