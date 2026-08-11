function trajectory = track(species,initialState,field,startValue,stopValue,maximumStep,options)
%TRACK Build a transparent fixed-step native beam trajectory.
arguments
    species (1,1) struct
    initialState (1,1) struct
    field (1,1) struct
    startValue (1,1) double {mustBeReal,mustBeFinite}
    stopValue (1,1) double {mustBeReal,mustBeFinite}
    maximumStep (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
    options.Independent (1,1) string = "time"
    options.Stepper (1,1) string {mustBeMember(options.Stepper, ...
        ["boris2","classical-rk4"])} = "boris2"
    options.MaximumSteps (1,1) double {mustBeInteger,mustBePositive} = 1e6
end
config = trackingConfig(species,initialState,field,options.Independent);
config.start = double(startValue);
config.stop = double(stopValue);
config.maximum_step = double(maximumStep);
config.maximum_steps = double(options.MaximumSteps);
config.stepper = char(options.Stepper);
trajectory = radia.internal.callMex('beam.track',config);
end
