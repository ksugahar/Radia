function family = loadHCurlEddyFosterFamily(fileName, options)
%LOADHCURLEDDYFOSTERFAMILY Load a height-indexed Foster HCurl family.
%   The file is produced by radia.vim.ExportHCurlEddyFosterFamilyJSON.  All
%   snapshots share one reduced R/L pair (a conductor meshed once and moved
%   relative to its source), so one set of decay rates and one discrete
%   state matrix serve every height; only the modal input matrix and the
%   modal force operator vary with height.

arguments
    fileName (1,1) string
    options.SampleTime_s double = []
    options.Interpolation (1,1) string {mustBeMember(options.Interpolation, ...
        ["nearest", "linear", "pchip"])} = "linear"
    options.Extrapolation (1,1) string {mustBeMember(options.Extrapolation, ...
        ["error", "clamp"])} = "error"
end

if ~isfile(fileName)
    error("radia:simulink:HCurlFosterFamily", ...
        "Foster family file does not exist: %s", fileName);
end
payload = jsondecode(fileread(fileName));
if ~isfield(payload, "schema") || string(payload.schema) ~= ...
        "radia.hcurl.eddy_foster.family.v1"
    error("radia:simulink:HCurlFosterFamily", ...
        "unsupported Foster family schema.");
end
if ~isfield(payload, "shared_modes") || ~payload.shared_modes
    error("radia:simulink:HCurlFosterFamily", ...
        "a Foster family must declare one shared mode set.");
end
nState = double(payload.state_order);
nPort = double(payload.port_count);
lambda = radia.internal.decodeHCurlEddyFosterArray(payload.decay_rates, nState);
snapshots = payload.snapshots;
count = numel(snapshots);
if count < 1
    error("radia:simulink:HCurlFosterFamily", ...
        "the family must contain at least one snapshot.");
end
positions = zeros(count, 1);
modalPort = zeros(nState, nPort, count);
forceOperator = zeros(3, nState, nPort, count);
for index = 1:count
    snapshot = snapshots(index);
    if ~isfield(snapshot, "height_m") || ~isfield(snapshot, "arrays")
        error("radia:simulink:HCurlFosterFamily", ...
            "every snapshot needs height_m and arrays.");
    end
    positions(index) = double(snapshot.height_m);
    modalPort(:, :, index) = radia.internal.decodeHCurlEddyFosterArray( ...
        snapshot.arrays.modal_port_rhs, [nState, nPort]);
    if isfield(snapshot.arrays, "modal_force_operator")
        forceOperator(:, :, :, index) = radia.internal.decodeHCurlEddyFosterArray( ...
            snapshot.arrays.modal_force_operator, [3, nState, nPort]);
    end
end
if isempty(options.SampleTime_s)
    sampleTime = double(payload.sample_time_s);
else
    sampleTime = double(options.SampleTime_s);
end
family = radia.simulink.makeHCurlEddyFosterFamily(positions, lambda, ...
    modalPort, forceOperator, SampleTime_s=sampleTime, ...
    Interpolation=options.Interpolation, Extrapolation=options.Extrapolation);
family.source_file = fileName;
if isfield(payload, "metadata")
    family.metadata = payload.metadata;
end
end
