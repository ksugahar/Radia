function result = simulateIHWaveform(plant, time_s, power_W, ambient_temperature_K, options)
%SIMULATEIHWAVEFORM Evaluate the Simulink-compatible IH plant offline.
%   RESULT = radia.simulink.simulateIHWaveform(PLANT,T,P,TA) uses the same
%   state update as the generated Simulink Discrete State-Space block. The
%   returned samples are the block outputs before each sample update; the
%   final state after the last input is available as RESULT.final_state.

arguments
    plant (1,1) struct
    time_s (:,1) double {mustBeFinite}
    power_W (:,1) double {mustBeFinite}
    ambient_temperature_K (:,1) double {mustBeFinite}
    options.Position_rad (:,1) double {mustBeFinite} = zeros(size(time_s))
    options.Speed_rad_s (:,1) double {mustBeFinite} = zeros(size(time_s))
end

required = ["A", "B", "C", "D", "x0", "sample_time_s"];
if ~all(isfield(plant, cellstr(required)))
    error("radia:simulink:InvalidPlant", ...
        "plant is missing one or more state-space fields.");
end
if numel(time_s) < 1 || numel(power_W) ~= numel(time_s) || ...
        numel(ambient_temperature_K) ~= numel(time_s)
    error("radia:simulink:InputSize", ...
        "time, power, and ambient-temperature vectors must have equal length.");
end
if numel(options.Position_rad) ~= numel(time_s) || ...
        numel(options.Speed_rad_s) ~= numel(time_s)
    error("radia:simulink:MotionSize", ...
        "position and speed vectors must match time_s.");
end
if abs(time_s(1)) > 100 * eps(max(1, abs(time_s(1))))
    error("radia:simulink:TimeOrigin", "time_s must start at zero.");
end
if numel(time_s) > 1
    dt = diff(time_s);
    if any(dt <= 0) || any(abs(dt - plant.sample_time_s) > ...
            1e-8 * max(1, plant.sample_time_s))
        error("radia:simulink:TimeGrid", ...
            "time_s must be strictly increasing with plant.sample_time_s spacing.");
    end
end

n = numel(time_s);
state = plant.x0(:);
y = zeros(n, size(plant.C, 1));
xHistory = zeros(n, numel(state));
u = [power_W, ambient_temperature_K];
for k = 1:n
    xHistory(k, :) = state.';
    y(k, :) = (plant.C * state + plant.D * u(k, :).').';
    state = plant.A * state + plant.B * u(k, :).';
end

result = struct( ...
    "schema", "radia.ih.simulink.waveform.v1", ...
    "time_s", time_s, ...
    "power_in_W", power_W, ...
    "ambient_temperature_K", ambient_temperature_K, ...
    "temperature_K", y(:, 1), ...
    "heat_loss_W", y(:, 2), ...
    "energy_input_J", y(:, 3), ...
    "temperature_rate_K_per_s", y(:, 4), ...
    "position_rad", options.Position_rad, ...
    "speed_rad_s", options.Speed_rad_s, ...
    "state_history", xHistory, ...
    "final_state", state);
end
