function modalCoefficients = solveHCurlEddyFosterHarmonic(model, frequency_Hz, coilCurrent)
%SOLVEHCURLEDDYFOSTERHARMONIC Harmonic response of a Foster HCurl model.
%   With s = j*omega and the coil-current phasor i, the reduced system
%   (R + s*L)c = -s*P*i is diagonal in the Foster modes:
%       z_j = -s*(Bm*i)_j / (s + lambda_j),   c = V*z.
%   The returned modal coefficients z feed evaluateHCurlEddyFosterForce.

arguments
    model (1,1) struct
    frequency_Hz (1,1) double {mustBePositive, mustBeFinite}
    coilCurrent double {mustBeFinite}
end

radia.internal.validateHCurlEddyFosterModel(model);
drive = localDrive(coilCurrent, model.port_count);
s = 1i * 2.0 * pi * frequency_Hz;
modalCoefficients = -s * (model.modal_port_rhs * drive) ./ ...
    (s + model.decay_rates(:));
end

function drive = localDrive(coilCurrent, nPort)
if isscalar(coilCurrent)
    drive = repmat(double(coilCurrent), nPort, 1);
elseif numel(coilCurrent) == nPort
    drive = double(coilCurrent(:));
else
    error("radia:simulink:HCurlFosterDrive", ...
        "coilCurrent must be scalar or contain port_count values.");
end
end
