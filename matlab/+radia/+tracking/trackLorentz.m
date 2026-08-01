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
if ~isempty(plane)
    odeOptions = odeset(odeOptions,"Events", ...
        @(time,state) localPlaneEvent(time,state,plane));
    [time,state,eventTime,eventState] = ode113( ...
        rhs,timesS,state0,odeOptions);
else
    [time,state] = ode113(rhs,timesS,state0,odeOptions);
    eventTime = [];
    eventState = [];
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
    eventNormalizedMomentum = eventState(end,4:6).';
    if options.Relativistic
        eventGamma = sqrt(1.0+dot( ...
            eventNormalizedMomentum,eventNormalizedMomentum));
    else
        eventGamma = 1.0;
    end
    result.stop_event = struct( ...
        "face","plane", ...
        "time_s",eventTime(end), ...
        "position_m",eventState(end,1:3).', ...
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

function [value,isTerminal,direction] = localPlaneEvent(~,state,plane)
value = dot(state(1:3)-plane.point_m,plane.normal);
isTerminal = 1;
direction = plane.direction;
end
