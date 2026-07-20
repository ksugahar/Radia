function family = loadHCurlEddyCLNFamily(fileName, options)
%LOADHCURLEDDYCLNFAMILY Load a height-indexed HCurl-VIM exchange family.
%   Each snapshot must use the same reduced state coordinates.  The family
%   can then be evaluated at a mechanical height without an implicit state
%   transfer or a hidden change of basis.

arguments
    fileName (1,1) string
    options.SampleTime_s double = []
    options.Interpolation (1,1) string {mustBeMember(options.Interpolation, ...
        ["nearest", "linear", "pchip"])} = "linear"
    options.Extrapolation (1,1) string {mustBeMember(options.Extrapolation, ...
        ["error", "clamp"])} = "error"
end

if ~isfile(fileName)
    error("radia:simulink:HCurlCLNFamily", ...
        "HCurl CLN family file does not exist: %s", fileName);
end
payload = jsondecode(fileread(fileName));
if ~isfield(payload, "schema") || string(payload.schema) ~= ...
        "radia.hcurl.eddy_cln.family.v1"
    error("radia:simulink:HCurlCLNFamily", ...
        "unsupported HCurl CLN family schema.");
end
if ~isfield(payload, "shared_state_basis") || ~payload.shared_state_basis
    error("radia:simulink:HCurlCLNFamily", ...
        "a family must declare a shared reduced state basis.");
end

snapshots = payload.snapshots;
count = numel(snapshots);
if count < 1
    error("radia:simulink:HCurlCLNFamily", ...
        "the family must contain at least one snapshot.");
end
nState = double(payload.state_order);
nPort = double(payload.port_count);
positions = zeros(count, 1);
models = cell(1, count);
if isempty(options.SampleTime_s)
    sampleTime = double(payload.sample_time_s);
else
    if ~isscalar(options.SampleTime_s) || ...
            ~isfinite(options.SampleTime_s) || options.SampleTime_s <= 0
        error("radia:simulink:HCurlCLNFamily", ...
            "SampleTime_s must be a positive finite scalar.");
    end
    sampleTime = options.SampleTime_s;
end

for index = 1:count
    snapshot = snapshots(index);
    if ~isfield(snapshot, "height_m") || ~isfield(snapshot, "arrays")
        error("radia:simulink:HCurlCLNFamily", ...
            "every snapshot needs height_m and arrays.");
    end
    positions(index) = double(snapshot.height_m);
    R = localDecodeMatrix(snapshot.arrays.resistance, [nState, nState]);
    L = localDecodeMatrix(snapshot.arrays.inductance, [nState, nState]);
    P = localDecodeMatrix(snapshot.arrays.port_rhs, [nState, nPort]);
    model = radia.simulink.makeHCurlEddyCLNModel( ...
        R, L, P, SampleTime_s=sampleTime);
    model.height_m = positions(index);
    model.exchange_schema = "radia.hcurl.eddy_cln.exchange.v1";
    if isfield(snapshot.arrays, "force_operator")
        model.force_operator = localDecodeTensor( ...
            snapshot.arrays.force_operator, [3, nState, nPort]);
    else
        model.force_operator = zeros(3, nState, nPort);
    end
    if isfield(snapshot, "metadata")
        model.metadata = snapshot.metadata;
    else
        model.metadata = struct();
    end
    models{index} = model;
end

if any(~isfinite(positions)) || any(diff(positions) <= 0)
    error("radia:simulink:HCurlCLNFamily", ...
        "snapshot heights must be finite and strictly increasing.");
end
family = struct( ...
    "schema", "radia.hcurl.eddy_cln.family.v1", ...
    "source_file", fileName, ...
    "shared_state_basis", true, ...
    "positions_m", positions, ...
    "models", {models}, ...
    "snapshot_count", count, ...
    "state_order", nState, ...
    "port_count", nPort, ...
    "sample_time_s", sampleTime, ...
    "interpolation", options.Interpolation, ...
    "extrapolation", options.Extrapolation, ...
    "metadata", localFieldOrDefault(payload, "metadata", struct()));
end

function value = localFieldOrDefault(data, fieldName, defaultValue)
if isfield(data, fieldName)
    value = data.(fieldName);
else
    value = defaultValue;
end
end

function matrix = localDecodeMatrix(encoded, expectedShape)
values = double(encoded.values(:));
if ~isfield(encoded, "shape") || ~isequal(double(encoded.shape(:)).', expectedShape)
    error("radia:simulink:HCurlCLNFamily", ...
        "exchange matrix shape does not match the family dimensions.");
end
matrix = reshape(values, [expectedShape(2), expectedShape(1)]).';
end

function tensor = localDecodeTensor(encoded, expectedShape)
values = double(encoded.values(:));
if ~isfield(encoded, "shape") || ~isequal(double(encoded.shape(:)).', expectedShape)
    error("radia:simulink:HCurlCLNFamily", ...
        "exchange tensor shape does not match the family dimensions.");
end
raw = reshape(values, [expectedShape(3), expectedShape(2), expectedShape(1)]);
tensor = permute(raw, [3, 2, 1]);
end
