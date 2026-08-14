function trajectory = trackGridFunction( ...
    species, initialState, field, startValue, stopValue, maximumStep, options)
%TRACKGRIDFUNCTION Track directly through a live NGSolve GridFunction.
%   TRAJECTORY = radia.beam.trackGridFunction(SPECIES, STATE, FIELD, START,
%   STOP, MAXSTEP) keeps FIELD inside NGSolve and evaluates it at every
%   native integration stage. No regular-grid field map or MATLAB callback is
%   introduced. FIELD must be a real three-component magnetic GridFunction.
%
%   This direct Lorentz path is the independent check for the local
%   multipole-map approximation returned by propagateGridFunctionMultipoleMap.

arguments
    species (1,1) struct
    initialState (1,1) struct
    field (1,1) radia.ngsolve.GridFunction
    startValue (1,1) double {mustBeReal,mustBeFinite}
    stopValue (1,1) double {mustBeReal,mustBeFinite}
    maximumStep (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
    options.Independent (1,1) string {mustBeMember(options.Independent, ...
        ["time","path_length","azimuth"])} = "path_length"
    options.Stepper (1,1) string {mustBeMember(options.Stepper, ...
        ["boris2","classical-rk4"])} = "classical-rk4"
    options.MaximumSteps (1,1) double {mustBeInteger,mustBePositive} = 1e6
end

% Reuse the public beam species/state validation, then replace its ordinary
% field descriptor with the checked native GridFunction handle at the MEX
% boundary.
config = trackingConfig( ...
    species,initialState,struct(),options.Independent);
config = rmfield(config,'field');
config.start = double(startValue);
config.stop = double(stopValue);
config.maximum_step = double(maximumStep);
config.maximum_steps = double(options.MaximumSteps);
config.stepper = char(options.Stepper);
trajectory = radia.internal.callMex( ...
    'beam.track.grid_function',field.nativeHandle(),config);
end
