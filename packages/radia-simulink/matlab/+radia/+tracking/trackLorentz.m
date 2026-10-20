function result = trackLorentz(chargeC,massKg,initialPositionM, ...
    initialVelocityMPerS,timesS,options)
%TRACKLORENTZ Integrate a charged particle in E and B fields in SI units.
%   Momentum, rather than velocity, is the integrated mechanical state.
%   Magnetic-only tracking therefore conserves gamma and kinetic energy,
%   while the same implementation remains valid when an electric field
%   changes the energy.
arguments
    chargeC (1,1) double {mustBeFinite,mustBeNonzero}
    massKg (1,1) double {mustBeFinite,mustBePositive}
    initialPositionM (3,1) double {mustBeFinite}
    initialVelocityMPerS (3,1) double {mustBeFinite}
    timesS (:,1) double {mustBeFinite}
    options.ElectricField (1,1) function_handle = @localZeroField
    options.MagneticField (1,1) function_handle = @localZeroField
    options.Relativistic (1,1) logical = true
    options.ExitPlane = []
    options.StopBox = []
    options.RelativeTolerance (1,1) double {mustBeFinite,mustBePositive} = 1e-10
    options.AbsoluteTolerance (1,1) double {mustBeFinite,mustBePositive} = 1e-13
end

if numel(timesS)<2 || any(diff(timesS)<=0)
    error("radia:tracking:Times", ...
        "timesS must contain at least two strictly increasing values");
end
speedOfLight = 299792458.0;
speed0 = norm(initialVelocityMPerS);
if options.Relativistic
    if speed0>=speedOfLight
        error("radia:tracking:Speed", ...
            "relativistic initial speed must be below the speed of light");
    end
    gamma0 = 1.0/sqrt(1.0-(speed0/speedOfLight)^2);
else
    gamma0 = 1.0;
end
normalizedMomentum0 = gamma0*initialVelocityMPerS/speedOfLight;
state0 = [initialPositionM;normalizedMomentum0];
rhs = @(time,state) radia.tracking.relativisticLorentzRhs( ...
    time,state,chargeC,massKg,options.ElectricField, ...
    options.MagneticField,options.Relativistic);
odeOptions = odeset("RelTol",options.RelativeTolerance, ...
    "AbsTol",options.AbsoluteTolerance);
plane = localValidatePlane(options.ExitPlane);
box = localValidateBox(options.StopBox);
if ~isempty(plane) || ~isempty(box)
    eventNames = localEventNames(box,plane);
    odeOptions = odeset(odeOptions,"Events", ...
        @(time,state) localStopEvents(time,state,box,plane));
    [time,state,eventTime,eventState,eventIndex] = ode113( ...
        rhs,timesS,state0,odeOptions);
else
    [time,state] = ode113(rhs,timesS,state0,odeOptions);
    eventTime = [];
    eventState = [];
    eventIndex = [];
end

normalizedMomentum = state(:,4:6);
if options.Relativistic
    gamma = sqrt(1.0+sum(normalizedMomentum.^2,2));
    velocity = speedOfLight*normalizedMomentum./gamma;
    kineticEnergy = (gamma-1.0)*massKg*speedOfLight^2;
else
    gamma = ones(size(time));
    velocity = speedOfLight*normalizedMomentum;
    kineticEnergy = 0.5*massKg*speedOfLight^2* ...
        sum(normalizedMomentum.^2,2);
end
momentum = massKg*speedOfLight*normalizedMomentum;
energyScale = max(abs(kineticEnergy(1)),realmin("double"));
result = struct( ...
    "schema","radia-time-domain-lorentz-track/v1", ...
    "backend","MATLAB ode113", ...
    "units",struct( ...
        "position","m", ...
        "velocity","m/s", ...
        "time","s", ...
        "electric_field","V/m", ...
        "magnetic_flux_density","T", ...
        "kinetic_energy","J"), ...
    "relativistic",options.Relativistic, ...
    "success",true, ...
    "message","", ...
    "time_s",time, ...
    "position_m",state(:,1:3), ...
    "velocity_m_s",velocity, ...
    "momentum_kg_m_per_s",momentum, ...
    "gamma",gamma, ...
    "kinetic_energy_j",kineticEnergy, ...
    "maximum_relative_kinetic_energy_drift", ...
        max(abs(kineticEnergy-kineticEnergy(1)))/energyScale, ...
    "stop_event",[]);
if ~isempty(eventTime)
    eventNormalizedMomentum = eventState(1,4:6).';
    if options.Relativistic
        eventGamma = sqrt(1.0+dot( ...
            eventNormalizedMomentum,eventNormalizedMomentum));
    else
        eventGamma = 1.0;
    end
    result.stop_event = struct( ...
        "face",eventNames(eventIndex(1)), ...
        "time_s",eventTime(1), ...
        "position_m",eventState(1,1:3).', ...
        "velocity_m_s",speedOfLight*eventNormalizedMomentum/eventGamma);
end
end

function value = localZeroField(~,~,~)
value = zeros(3,1);
end

function plane = localValidatePlane(value)
if isempty(value)
    plane = [];
    return
end
if ~isstruct(value) || ~isfield(value,"point_m") || ~isfield(value,"normal")
    error("radia:tracking:ExitPlane", ...
        "ExitPlane must contain point_m and normal");
end
point = double(value.point_m(:));
normal = double(value.normal(:));
if ~isequal(size(point),[3,1]) || ~isequal(size(normal),[3,1]) || ...
        any(~isfinite([point;normal])) || norm(normal)==0
    error("radia:tracking:ExitPlane", ...
        "ExitPlane point and normal must be finite three-vectors");
end
direction = 0;
if isfield(value,"direction")
    direction = double(value.direction);
end
if ~isscalar(direction) || ~ismember(direction,[-1,0,1])
    error("radia:tracking:ExitPlane", ...
        "ExitPlane direction must be -1, 0, or +1");
end
plane = struct("point_m",point,"normal",normal/norm(normal), ...
    "direction",direction);
end

function box = localValidateBox(value)
if isempty(value)
    box = [];
    return
end
if ~isstruct(value) || ~isfield(value,"minimum_m") || ...
        ~isfield(value,"maximum_m")
    error("radia:tracking:StopBox", ...
        "StopBox must contain minimum_m and maximum_m");
end
minimum = double(value.minimum_m(:));
maximum = double(value.maximum_m(:));
if ~isequal(size(minimum),[3,1]) || ...
        ~isequal(size(maximum),[3,1]) || ...
        any(~isfinite([minimum;maximum])) || any(minimum>=maximum)
    error("radia:tracking:StopBox", ...
        "StopBox bounds must be finite three-vectors with minimum below maximum");
end
box = struct("minimum_m",minimum,"maximum_m",maximum);
end

function names = localEventNames(box,plane)
names = strings(0,1);
if ~isempty(box)
    names = ["x_minimum";"x_maximum"; ...
        "y_minimum";"y_maximum";"z_minimum";"z_maximum"];
end
if ~isempty(plane)
    names(end+1,1) = "plane";
end
end

function [value,isTerminal,direction] = localStopEvents(~,state,box,plane)
value = zeros(0,1);
direction = zeros(0,1);
if ~isempty(box)
    value = [state(1)-box.minimum_m(1);box.maximum_m(1)-state(1); ...
        state(2)-box.minimum_m(2);box.maximum_m(2)-state(2); ...
        state(3)-box.minimum_m(3);box.maximum_m(3)-state(3)];
    direction = -ones(6,1);
end
if ~isempty(plane)
    value(end+1,1) = dot(state(1:3)-plane.point_m,plane.normal);
    direction(end+1,1) = plane.direction;
end
isTerminal = ones(size(value));
end
