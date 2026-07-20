function environment = makeIHEnvironment(plant, options)
%MAKEIHENVIRONMENT Create a fast discrete induction-heating RL environment.
%   The plant matrices come from radia.simulink.makeIHPlant.  One RL action is
%   the supplied power in watts; one step is exactly one plant sample period.

arguments
    plant (1,1) struct
    options.TargetTemperature_K (1,1) double = NaN
    options.AmbientTemperature_K (1,1) double = NaN
    options.PowerLowerBound_W (1,1) double = 0
    options.PowerUpperBound_W (1,1) double = 1e6
    options.TemperaturePenalty (1,1) double = 1
    options.PowerPenalty (1,1) double = 1e-8
    options.MaxSteps (1,1) double = 1000
    options.InitialState double = []
    options.RewardFcn = []
end

if ~isfield(plant, "schema") || plant.schema ~= "radia.ih.simulink.plant.v1"
    error("radia:rl:InvalidPlant", "plant must come from radia.simulink.makeIHPlant.");
end
if ~(isfinite(options.PowerLowerBound_W) && isfinite(options.PowerUpperBound_W) && ...
        0 <= options.PowerLowerBound_W && options.PowerLowerBound_W < options.PowerUpperBound_W)
    error("radia:rl:PowerBounds", "Power bounds must be finite with lower < upper.");
end
if options.TemperaturePenalty < 0 || options.PowerPenalty < 0
    error("radia:rl:Penalty", "Reward penalties must be nonnegative.");
end
if ~isempty(options.RewardFcn) && ~isa(options.RewardFcn, "function_handle")
    error("radia:rl:RewardFcn", "RewardFcn must be a function handle or empty.");
end

if isnan(options.AmbientTemperature_K)
    ambient = plant.ambient_temperature_K;
else
    ambient = options.AmbientTemperature_K;
end
if isnan(options.TargetTemperature_K)
    target = ambient + 10;
else
    target = options.TargetTemperature_K;
end
if isempty(options.InitialState)
    state = plant.x0;
else
    state = options.InitialState(:);
end
if ~isequal(size(state), size(plant.x0)) || any(~isfinite(state))
    error("radia:rl:InitialState", "InitialState must match plant.x0 and be finite.");
end

initialState = state;
stepIndex = 0;
environment = radia.rl.Environment( ...
    ResetFcn=@resetEpisode, StepFcn=@stepEpisode, MaxSteps=options.MaxSteps);

    function [observation, info] = resetEpisode()
        state = initialState;
        stepIndex = 0;
        [observation, info] = observe(state, 0);
    end

    function [observation, reward, isDone, info] = stepEpisode(action)
        if ~isnumeric(action) || ~isscalar(action) || ~isfinite(action)
            error("radia:rl:Action", "IH action must be one finite numeric power value.");
        end
        power = min(max(double(action), options.PowerLowerBound_W), ...
            options.PowerUpperBound_W);
        input = [power; ambient];
        state = plant.A * state + plant.B * input;
        stepIndex = stepIndex + 1;
        [observation, info] = observe(state, power);
        if isempty(options.RewardFcn)
            reward = -options.TemperaturePenalty * ...
                (observation(1) - target)^2 - options.PowerPenalty * power^2;
        else
            reward = options.RewardFcn(state, observation, power, stepIndex);
        end
        isDone = false;
    end

    function [observation, info] = observe(currentState, power)
        input = [power; ambient];
        observation = plant.C * currentState + plant.D * input;
        info = struct( ...
            "state", currentState, ...
            "power_W", power, ...
            "ambient_temperature_K", ambient, ...
            "target_temperature_K", target, ...
            "temperature_K", observation(1), ...
            "heat_loss_W", observation(2), ...
            "energy_input_J", observation(3), ...
            "temperature_rate_K_per_s", observation(4), ...
            "step", stepIndex);
    end
end
