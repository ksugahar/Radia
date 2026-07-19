function result = simulateIHDrive(plant, time_s, drive_A, ambient_temperature_K, lut, options)
%SIMULATEIHDRIVE Couple motion, drive, EM loss LUT, and thermal state.
%   The LUT inputs may be any subset of:
%   position_rad, speed_rad_s, drive_A, temperature_K, and
%   ambient_temperature_K. The current thermal state is fed back to the LUT
%   before each step, which makes a temperature-dependent ESIM table usable
%   in a control-oriented simulation.

arguments
    plant (1,1) struct
    time_s (:,1) double {mustBeFinite}
    drive_A (:,1) double {mustBeFinite}
    ambient_temperature_K (:,1) double {mustBeFinite}
    lut (1,1) struct
    options.Position_rad (:,1) double {mustBeFinite} = zeros(size(time_s))
    options.Speed_rad_s (:,1) double {mustBeFinite} = zeros(size(time_s))
end

n = numel(time_s);
if n < 1 || numel(drive_A) ~= n || numel(ambient_temperature_K) ~= n || ...
        numel(options.Position_rad) ~= n || numel(options.Speed_rad_s) ~= n
    error("radia:simulink:DriveInputSize", ...
        "drive, ambient, position, and speed vectors must match time_s.");
end
if abs(time_s(1)) > 100 * eps(max(1, abs(time_s(1))))
    error("radia:simulink:TimeOrigin", "time_s must start at zero.");
end
if n > 1
    dt = diff(time_s);
    if any(dt <= 0) || any(abs(dt - plant.sample_time_s) > ...
            1e-8 * max(1, plant.sample_time_s))
        error("radia:simulink:TimeGrid", ...
            "time_s must be strictly increasing with plant.sample_time_s spacing.");
    end
end

state = plant.x0(:);
y = zeros(n, size(plant.C, 1));
xHistory = zeros(n, numel(state));
power_W = zeros(n, 1);
for k = 1:n
    lutInput = zeros(1, numel(lut.input_names));
    for j = 1:numel(lut.input_names)
        switch lut.input_names(j)
            case "position_rad"
                lutInput(j) = options.Position_rad(k);
            case "speed_rad_s"
                lutInput(j) = options.Speed_rad_s(k);
            case "drive_A"
                lutInput(j) = drive_A(k);
            case "temperature_K"
                lutInput(j) = state(1);
            case "ambient_temperature_K"
                lutInput(j) = ambient_temperature_K(k);
            otherwise
                error("radia:simulink:UnknownLUTInput", ...
                    "unsupported motion/LUT input name '%s'.", lut.input_names(j));
        end
    end
    power_W(k) = radia.simulink.evaluateIHPowerLUT(lut, lutInput);
    xHistory(k, :) = state.';
    u = [power_W(k); ambient_temperature_K(k)];
    y(k, :) = (plant.C * state + plant.D * u).';
    state = plant.A * state + plant.B * u;
end

result = struct( ...
    "schema", "radia.ih.simulink.drive_waveform.v1", ...
    "time_s", time_s, "drive_A", drive_A, "power_in_W", power_W, ...
    "ambient_temperature_K", ambient_temperature_K, ...
    "temperature_K", y(:, 1), "heat_loss_W", y(:, 2), ...
    "energy_input_J", y(:, 3), ...
    "temperature_rate_K_per_s", y(:, 4), ...
    "position_rad", options.Position_rad, ...
    "speed_rad_s", options.Speed_rad_s, ...
    "state_history", xHistory, "final_state", state);
end
