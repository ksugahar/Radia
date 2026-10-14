function model = loadHCurlEddyFosterModel(fileName, options)
%LOADHCURLEDDYFOSTERMODEL Load a Foster HCurl exchange written by Python.
%   The file is produced by radia.vim.ExportHCurlEddyFosterJSON and holds
%   the decay rates, the modal input matrix and, optionally, the modal
%   Lorentz force operator.

arguments
    fileName (1,1) string
    options.SampleTime_s double = []
end

if ~isfile(fileName)
    error("radia:simulink:HCurlFosterExchange", ...
        "Foster exchange file does not exist: %s", fileName);
end
payload = jsondecode(fileread(fileName));
if ~isfield(payload, "schema") || string(payload.schema) ~= ...
        "radia.hcurl.eddy_foster.exchange.v1"
    error("radia:simulink:HCurlFosterExchange", ...
        "unsupported Foster exchange schema.");
end
if isfield(payload, "has_sibc_termination") && payload.has_sibc_termination
    error("radia:simulink:HCurlFosterSIBC", ...
        "SIBC termination must be rationalized before MATLAB state-space export.");
end
nState = double(payload.state_order);
nPort = double(payload.port_count);
lambda = radia.internal.decodeHCurlEddyFosterArray(payload.decay_rates, nState);
Bm = radia.internal.decodeHCurlEddyFosterArray( ...
    payload.arrays.modal_port_rhs, [nState, nPort]);
K = [];
if isfield(payload.arrays, "modal_force_operator")
    K = radia.internal.decodeHCurlEddyFosterArray( ...
        payload.arrays.modal_force_operator, [3, nState, nPort]);
end
model = radia.simulink.makeHCurlEddyFosterModel(lambda, Bm, ...
    SampleTime_s=localSampleTime(payload, options.SampleTime_s), ...
    ForceOperator=K);
model.exchange_schema = string(payload.schema);
model.source_file = fileName;
model.metadata = localField(payload, "metadata", struct());
if isfield(payload, "basis_names")
    model.basis_names = string(payload.basis_names(:));
end
if isfield(payload, "blocks")
    model.blocks = payload.blocks;
end
end

function sampleTime = localSampleTime(payload, requested)
if isempty(requested)
    sampleTime = double(payload.sample_time_s);
elseif isscalar(requested) && isfinite(requested) && requested > 0
    sampleTime = double(requested);
else
    error("radia:simulink:HCurlFosterExchange", ...
        "SampleTime_s must be a positive finite scalar.");
end
end

function value = localField(data, name, default)
if isfield(data, name)
    value = data.(name);
else
    value = default;
end
end
