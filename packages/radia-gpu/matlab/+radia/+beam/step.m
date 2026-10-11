function result = step(species,state,field,stepSize,options)
%STEP Apply one named native beam integrator step.
arguments
    species (1,1) struct
    state (1,1) struct
    field (1,1) struct
    stepSize (1,1) double {mustBeReal,mustBeFinite}
    options.Independent (1,1) string = "time"
    options.IndependentValue (1,1) double {mustBeReal,mustBeFinite} = 0
    options.Stepper (1,1) string {mustBeMember(options.Stepper, ...
        ["boris2","classical-rk4"])} = "boris2"
end
if stepSize == 0
    error("radia:beam:InvalidStep","Step size must be nonzero.");
end
config = trackingConfig(species,state,field,options.Independent);
config.independent_value = double(options.IndependentValue);
config.step = double(stepSize);
config.stepper = char(options.Stepper);
result = radia.internal.callMex('beam.step',config);
end
