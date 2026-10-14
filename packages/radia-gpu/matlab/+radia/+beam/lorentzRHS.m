function rhs = lorentzRHS(species,state,field,options)
%LORENTZRHS Evaluate the native relativistic Lorentz equation once.
arguments
    species (1,1) struct
    state (1,1) struct
    field (1,1) struct
    options.Independent (1,1) string = "time"
    options.IndependentValue (1,1) double {mustBeReal,mustBeFinite} = 0
end
config = trackingConfig(species,state,field,options.Independent);
config.independent_value = double(options.IndependentValue);
rhs = radia.internal.callMex('beam.equation.rhs',config);
end
