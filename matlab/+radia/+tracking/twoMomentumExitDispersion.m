function result = twoMomentumExitDispersion(chargeC,massKg, ...
    initialPositionM,nominalVelocityMPerS,timesS,exitPlane,options)
%TWOMOMENTUMEXITDISPERSION Compare p0*(1-delta) and p0*(1+delta).
arguments
    chargeC (1,1) double {mustBeFinite,mustBeNonzero}
    massKg (1,1) double {mustBeFinite,mustBePositive}
    initialPositionM (3,1) double {mustBeFinite}
    nominalVelocityMPerS (3,1) double {mustBeFinite}
    timesS (:,1) double {mustBeFinite}
    exitPlane (1,1) struct
    options.RelativeMomentumOffset (1,1) double {mustBeFinite,mustBePositive}
    options.TransverseDirection (3,1) double {mustBeFinite}
    options.ElectricField (1,1) function_handle = @localZeroField
    options.MagneticField (1,1) function_handle = @localZeroField
    options.Relativistic (1,1) logical = true
    options.RelativeTolerance (1,1) double {mustBeFinite,mustBePositive} = 1e-10
    options.AbsoluteTolerance (1,1) double {mustBeFinite,mustBePositive} = 1e-13
end

delta = options.RelativeMomentumOffset;
if delta>=1.0
    error("radia:tracking:MomentumOffset", ...
        "RelativeMomentumOffset must be less than one");
end
if ~isfield(exitPlane,"normal")
    error("radia:tracking:ExitPlane","exitPlane must contain normal");
end
normal = double(exitPlane.normal(:));
normal = normal/norm(normal);
transverse = options.TransverseDirection- ...
    dot(options.TransverseDirection,normal)*normal;
if norm(transverse)==0
    error("radia:tracking:TransverseDirection", ...
        "TransverseDirection must not be parallel to the plane normal");
end
transverse = transverse/norm(transverse);

minusVelocity = localMomentumScaledVelocity( ...
    nominalVelocityMPerS,1.0-delta,options.Relativistic);
plusVelocity = localMomentumScaledVelocity( ...
    nominalVelocityMPerS,1.0+delta,options.Relativistic);
common = {"ElectricField",options.ElectricField, ...
    "MagneticField",options.MagneticField, ...
    "Relativistic",options.Relativistic, ...
    "ExitPlane",exitPlane, ...
    "RelativeTolerance",options.RelativeTolerance, ...
    "AbsoluteTolerance",options.AbsoluteTolerance};
minusTrack = radia.tracking.trackLorentz(chargeC,massKg, ...
    initialPositionM,minusVelocity,timesS,common{:});
plusTrack = radia.tracking.trackLorentz(chargeC,massKg, ...
    initialPositionM,plusVelocity,timesS,common{:});
if isempty(minusTrack.stop_event) || isempty(plusTrack.stop_event)
    error("radia:tracking:ExitNotReached", ...
        "both momentum trajectories must reach the exit plane");
end
minusPosition = minusTrack.stop_event.position_m;
plusPosition = plusTrack.stop_event.position_m;
difference = plusPosition-minusPosition;
dispersion = difference/(2.0*delta);
result = struct( ...
    "schema","radia-two-momentum-exit-dispersion/v1", ...
    "relative_momentum_offset",delta, ...
    "minus_position_m",minusPosition, ...
    "plus_position_m",plusPosition, ...
    "position_difference_m",difference, ...
    "dispersion_vector_m",dispersion, ...
    "transverse_direction",transverse, ...
    "eta_m",dot(dispersion,transverse), ...
    "coincident_exit_error_m",norm(difference), ...
    "minus_track",minusTrack, ...
    "plus_track",plusTrack);
end

function velocity = localMomentumScaledVelocity(velocity0,scale,relativistic)
if ~relativistic
    velocity = scale*velocity0;
    return
end
speedOfLight = 299792458.0;
speed = norm(velocity0);
if speed>=speedOfLight
    error("radia:tracking:Speed", ...
        "relativistic initial speed must be below the speed of light");
end
betaGamma = velocity0/speedOfLight/sqrt(1.0-(speed/speedOfLight)^2);
scaled = scale*betaGamma;
gamma = sqrt(1.0+dot(scaled,scaled));
velocity = speedOfLight*scaled/gamma;
end

function value = localZeroField(~,~,~)
value = zeros(3,1);
end
