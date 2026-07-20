function model = loadHCurlEddyCLNModel(fileName, options)
%LOADHCURLEDDYCLNMODEL Load an NGSolve HCurl-VIM exchange JSON file.
%   The exporter stores row-major flattened arrays with explicit shapes so
%   MATLAB does not depend on nested JSON array orientation.

arguments
    fileName (1,1) string
    options.SampleTime_s double = []
end

if ~isfile(fileName)
    error("radia:simulink:HCurlCLNExchange", ...
        "HCurl CLN exchange file does not exist: %s", fileName);
end
payload = jsondecode(fileread(fileName));
if ~isfield(payload, "schema") || string(payload.schema) ~= ...
        "radia.hcurl.eddy_cln.exchange.v1"
    error("radia:simulink:HCurlCLNExchange", ...
        "unsupported HCurl CLN exchange schema.");
end
if isfield(payload, "has_sibc_termination") && payload.has_sibc_termination
    error("radia:simulink:HCurlCLNSIBC", ...
        "SIBC termination must be rationalized before MATLAB state-space export.");
end

nState = double(payload.state_order);
nPort = double(payload.port_count);
R = localDecodeMatrix(payload.arrays.resistance, [nState, nState]);
L = localDecodeMatrix(payload.arrays.inductance, [nState, nState]);
P = localDecodeMatrix(payload.arrays.port_rhs, [nState, nPort]);
if isempty(options.SampleTime_s)
    sampleTime = double(payload.sample_time_s);
else
    if ~isscalar(options.SampleTime_s) || ...
            ~isfinite(options.SampleTime_s) || options.SampleTime_s <= 0
        error("radia:simulink:HCurlCLNExchange", ...
            "SampleTime_s must be a positive finite scalar.");
    end
    sampleTime = options.SampleTime_s;
end
model = radia.simulink.makeHCurlEddyCLNModel( ...
    R, L, P, SampleTime_s=sampleTime);
model.exchange_schema = string(payload.schema);
model.source_file = fileName;
if isfield(payload, "metadata")
    model.metadata = payload.metadata;
else
    model.metadata = struct();
end
if isfield(payload.arrays, "force_operator")
    model.force_operator = localDecodeTensor( ...
        payload.arrays.force_operator, [3, nState, nPort]);
else
    model.force_operator = zeros(3, nState, nPort);
end
if isfield(payload, "basis_names")
    model.basis_names = string(payload.basis_names(:));
end
if isfield(payload, "blocks")
    model.blocks = payload.blocks;
end
end

function matrix = localDecodeMatrix(encoded, expectedShape)
values = double(encoded.values(:));
if ~isfield(encoded, "shape") || ~isequal(double(encoded.shape(:)).', expectedShape)
    error("radia:simulink:HCurlCLNExchange", ...
        "exchange matrix shape does not match the state-space dimensions.");
end
% Python's C-order [row][column] sequence is MATLAB's column-major reshape
% after swapping the two dimensions, followed by a transpose.
matrix = reshape(values, [expectedShape(2), expectedShape(1)]).';
end

function tensor = localDecodeTensor(encoded, expectedShape)
values = double(encoded.values(:));
if ~isfield(encoded, "shape") || ~isequal(double(encoded.shape(:)).', expectedShape)
    error("radia:simulink:HCurlCLNExchange", ...
        "exchange tensor shape does not match the state-space dimensions.");
end
% Python C-order [force][state][port] is reconstructed as MATLAB
% [port][state][force], then permuted to the public [force][state][port].
raw = reshape(values, [expectedShape(3), expectedShape(2), expectedShape(1)]);
tensor = permute(raw, [3, 2, 1]);
end
