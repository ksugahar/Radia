function model = loadPrimaPortModel(fileName)
%LOADPRIMAPORTMODEL Load a PRIMA port model written by Python as a dss object.
%   model = radia.simulink.loadPrimaPortModel(fileName) reads the exchange
%   written by radia.prima_export.export_prima_lti_json and returns a struct
%   whose sys field is the descriptor model
%
%       E x' = A x + B u,   y = C x + D u
%
%   as a Control System Toolbox dss object.  u is the port current [A]
%   (impedance orientation) or the port voltage [V] (admittance orientation);
%   y is the other quantity.  The file's check response is reproduced with
%   freqresp before the model is returned; a mismatch is an error.
%   Use radia.simulink.addPrimaLTIBlock to place it in a Simulink model.

arguments
    fileName (1,1) string
end

if ~isfile(fileName)
    error("radia:simulink:PrimaExchange", ...
        "PRIMA exchange file does not exist: %s", fileName);
end
payload = jsondecode(fileread(fileName));
if ~isfield(payload, "schema") || string(payload.schema) ~= "radia.prima.port_model.v1"
    error("radia:simulink:PrimaExchange", "unsupported PRIMA exchange schema.");
end
nState = double(payload.state_count);
nPort = double(payload.port_count);
names = string(payload.port_names(:)).';
if numel(names) ~= nPort
    error("radia:simulink:PrimaExchange", "port_names does not match port_count.");
end
arrays = payload.arrays;
E = localDecode(arrays.E, [nState, nState]);
A = localDecode(arrays.A, [nState, nState]);
B = localDecode(arrays.B, [nState, nPort]);
C = localDecode(arrays.C, [nPort, nState]);
D = localDecode(arrays.D, [nPort, nPort]);

sys = dss(A, B, C, D, E);
sys.InputName = names + "." + string(payload.input_quantity);
sys.OutputName = names + "." + string(payload.output_quantity);
sys.InputUnit = repmat(string(payload.input_unit), 1, nPort);
sys.OutputUnit = repmat(string(payload.output_unit), 1, nPort);

% The file may tighten the check limit, never loosen it past the exporter's.
reconstructionLimit = 1e-6;
check = payload.check;
frequency = double(check.frequency_hz(:));
if isempty(frequency) || any(~isfinite(frequency)) || any(frequency <= 0)
    error("radia:simulink:PrimaExchange", ...
        "check frequencies must be a nonempty list of finite positive values [Hz].");
end
limit = double(check.relative_limit);
if ~isscalar(limit) || ~isfinite(limit) || limit <= 0 || limit > reconstructionLimit
    error("radia:simulink:PrimaExchange", ...
        "check relative_limit must be a finite positive scalar <= %g.", reconstructionLimit);
end
nFrequency = numel(frequency);
expected = localDecode(check.response_real, [nFrequency, nPort, nPort]) + ...
    1i * localDecode(check.response_imag, [nFrequency, nPort, nPort]);
actual = freqresp(sys, 2 * pi * frequency);   % [nPort, nPort, nFrequency]
if any(~isfinite(actual(:)))
    error("radia:simulink:PrimaExchange", ...
        "dss model response is non-finite at a check frequency.");
end
worst = 0;
for k = 1:nFrequency
    reference = squeeze(expected(k, :, :));
    if nPort == 1
        reference = reshape(reference, 1, 1);
    end
    worst = max(worst, norm(actual(:, :, k) - reference) / ...
        max(norm(reference), realmin));
end
if ~isfinite(worst) || worst > limit
    error("radia:simulink:PrimaExchange", ...
        "dss model does not reproduce the exported response (relative error %.3g).", worst);
end

model = struct( ...
    "schema", string(payload.schema), ...
    "orientation", string(payload.orientation), ...
    "input_quantity", string(payload.input_quantity), ...
    "output_quantity", string(payload.output_quantity), ...
    "port_names", names, ...
    "E", E, "A", A, "B", B, "C", C, "D", D, ...
    "sys", sys, ...
    "check_relative_error", worst, ...
    "metadata", payload.metadata, ...
    "source_file", fileName);
end

function value = localDecode(encoded, expectedShape)
if ~isfield(encoded, "shape") || ~isfield(encoded, "values") || ...
        ~isequal(double(encoded.shape(:)).', expectedShape)
    error("radia:simulink:PrimaExchange", ...
        "exchange array shape does not match the declared dimensions.");
end
values = double(encoded.values(:));
if numel(values) ~= prod(expectedShape) || any(~isfinite(values))
    error("radia:simulink:PrimaExchange", ...
        "exchange array values are incomplete or non-finite.");
end
% Row-major (C order) on the Python side.
value = permute(reshape(values, fliplr(expectedShape)), numel(expectedShape):-1:1);
end
