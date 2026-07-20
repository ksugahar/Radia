function coefficients = solveHCurlEddyCLNHarmonic(model, frequency_Hz, coilCurrent)
%SOLVEHCURLEDDYCLNHARMONIC Solve a reduced HCurl-VIM model through MEX.
%   This uses the production hybrid_vim.solve kernel for
%       (R + j*omega*L)c = -j*omega*P*i.

arguments
    model (1,1) struct
    frequency_Hz (1,1) double {mustBePositive, mustBeFinite}
    coilCurrent double {mustBeFinite}
end

if ~isfield(model, "schema") || model.schema ~= "radia.hcurl.eddy_cln.state_space.v1"
    error("radia:simulink:HCurlCLNModel", "unsupported HCurl CLN model schema.");
end
if ~isscalar(coilCurrent) && ~isequal(size(coilCurrent), [model.port_count, 1])
    error("radia:simulink:HCurlCLNDrive", ...
        "coilCurrent must be scalar or a port_count-by-1 vector.");
end
if isscalar(coilCurrent)
    drive = repmat(double(coilCurrent), model.port_count, 1);
else
    drive = double(coilCurrent(:));
end
s = 1i * 2.0 * pi * frequency_Hz;
rhs = -s * (model.port_rhs * drive);
coefficients = radia.internal.callMex( ...
    'hybrid_vim.solve', model.resistance + s * model.inductance, rhs);
coefficients = coefficients(:);
end
