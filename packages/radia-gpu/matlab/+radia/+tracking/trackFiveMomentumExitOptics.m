function result = trackFiveMomentumExitOptics(chargeC,massKg, ...
    initialPositionM,nominalVelocityMPerS,timesS,exitPlane,options)
%TRACKFIVEMOMENTUMEXITOPTICS Track five momenta to one exit plane.
%   Particles share an initial position and direction. Their momenta are
%   scaled by 1+delta, and transverse position and angle are evaluated at
%   the common exit plane before fitting the five-particle optics metrics.
arguments
    chargeC (1,1) double {mustBeFinite,mustBeNonzero}
    massKg (1,1) double {mustBeFinite,mustBePositive}
    initialPositionM (3,1) double {mustBeFinite}
    nominalVelocityMPerS (3,1) double {mustBeFinite}
    timesS (:,1) double {mustBeFinite}
    exitPlane (1,1) struct
    options.ReferenceExitPointM {mustBeNumeric,mustBeReal}
    options.TransverseDirection {mustBeNumeric,mustBeReal}
    options.LongitudinalDirection {mustBeNumeric,mustBeReal} = []
    options.RelativeMomentumOffsets {mustBeNumeric,mustBeReal} = ...
        [-1.0e-3;-5.0e-4;0.0;5.0e-4;1.0e-3]
    options.ElectricField (1,1) function_handle = @localZeroField
    options.MagneticField (1,1) function_handle = @localZeroField
    options.Relativistic (1,1) logical = false
    options.StopBox = []
    options.X0LimitM (1,1) double = 1.0e-3
    options.Psi0LimitRad (1,1) double = 1.0e-3
    options.EtaLimitM (1,1) double = 5.0e-3
    options.EtaPrimeLimitRad (1,1) double = 1.0e-3
    options.RelativeTolerance (1,1) double {mustBeFinite,mustBePositive} = 1e-11
    options.AbsoluteTolerance (1,1) double {mustBeFinite,mustBePositive} = 1e-14
end

offsets=double(options.RelativeMomentumOffsets(:));
if numel(offsets)~=5 || any(~isfinite(offsets))
    error("radia:tracking:MomentumOffsets", ...
        "RelativeMomentumOffsets must contain five finite values");
end
if any(1.0+offsets<=0)
    error("radia:tracking:MomentumScale", ...
        "every relative momentum scale must be positive");
end
if any(diff(offsets)<=0)
    error("radia:tracking:MomentumOffsets", ...
        "relative momentum offsets must be strictly increasing");
end
if nnz(offsets==0.0)~=1
    error("radia:tracking:MomentumOffsets", ...
        "five-particle offsets must contain exactly one zero");
end
limits=[options.X0LimitM;options.Psi0LimitRad; ...
    options.EtaLimitM;options.EtaPrimeLimitRad];
if any(~isfinite(limits)) || any(limits<=0)
    error("radia:tracking:AcceptanceLimits", ...
        "five-particle acceptance limits must be finite and positive");
end

reference=localThreeVector(options.ReferenceExitPointM, ...
    "ReferenceExitPointM");
transverse=localThreeVector(options.TransverseDirection, ...
    "TransverseDirection");
planeNormal=localExitPlaneNormal(exitPlane);
if isempty(options.LongitudinalDirection)
    longitudinal=planeNormal;
else
    longitudinal=localThreeVector(options.LongitudinalDirection, ...
        "LongitudinalDirection");
end
longitudinalNorm=norm(longitudinal);
if longitudinalNorm==0
    error("radia:tracking:LongitudinalDirection", ...
        "LongitudinalDirection must be nonzero");
end
longitudinal=longitudinal/longitudinalNorm;
transverse=transverse-dot(transverse,longitudinal)*longitudinal;
transverseNorm=norm(transverse);
if transverseNorm==0
    error("radia:tracking:TransverseDirection", ...
        "TransverseDirection and LongitudinalDirection must not be parallel");
end
transverse=transverse/transverseNorm;

common={"ElectricField",options.ElectricField, ...
    "MagneticField",options.MagneticField, ...
    "Relativistic",options.Relativistic, ...
    "ExitPlane",exitPlane, ...
    "StopBox",options.StopBox, ...
    "RelativeTolerance",options.RelativeTolerance, ...
    "AbsoluteTolerance",options.AbsoluteTolerance};
tracks=struct([]);
positions=zeros(5,1);
angles=zeros(5,1);
for index=1:5
    velocity=localMomentumScaledVelocity(nominalVelocityMPerS, ...
        1.0+offsets(index),options.Relativistic);
    track=radia.tracking.trackLorentz(chargeC,massKg, ...
        initialPositionM,velocity,timesS,common{:});
    event=track.stop_event;
    if ~track.success || isempty(event) || string(event.face)~="plane"
        error("radia:tracking:ExitNotReached", ...
            "momentum offset %+.6g did not reach the exit plane: %s", ...
            offsets(index),string(track.message));
    end
    eventPosition=double(event.position_m(:));
    eventVelocity=double(event.velocity_m_s(:));
    positions(index)=dot(eventPosition-reference,transverse);
    angles(index)=atan2(dot(eventVelocity,transverse), ...
        dot(eventVelocity,longitudinal));
    if index==1
        tracks=repmat(track,5,1);
    else
        tracks(index,1)=track;
    end
end

result=radia.tracking.fitFiveMomentumExitOptics(offsets,positions,angles, ...
    X0LimitM=limits(1),Psi0LimitRad=limits(2), ...
    EtaLimitM=limits(3),EtaPrimeLimitRad=limits(4));
result.transverse_positions_m=positions;
result.transverse_angles_rad=angles;
result.reference_exit_point_m=reference;
result.transverse_direction=transverse;
result.longitudinal_direction=longitudinal;
result.tracks=tracks;
end

function value = localThreeVector(value,name)
value=double(value(:));
if ~isequal(size(value),[3,1]) || any(~isfinite(value))
    error("radia:tracking:ExitCoordinates", ...
        "%s must be a finite three-vector",name);
end
end

function normal = localExitPlaneNormal(exitPlane)
if ~isfield(exitPlane,"point_m") || ~isfield(exitPlane,"normal")
    error("radia:tracking:ExitPlane", ...
        "exitPlane must contain point_m and normal");
end
point=double(exitPlane.point_m(:));
normal=double(exitPlane.normal(:));
if ~isequal(size(point),[3,1]) || ~isequal(size(normal),[3,1]) || ...
        any(~isfinite([point;normal])) || norm(normal)==0
    error("radia:tracking:ExitPlane", ...
        "exitPlane point and normal must be finite three-vectors with a nonzero normal");
end
normal=normal/norm(normal);
end

function velocity = localMomentumScaledVelocity(velocity0,scale,relativistic)
if ~relativistic
    velocity=scale*velocity0;
    return
end
speedOfLight=299792458.0;
speed=norm(velocity0);
if speed>=speedOfLight
    error("radia:tracking:Speed", ...
        "relativistic initial speed must be below the speed of light");
end
betaGamma=velocity0/speedOfLight/sqrt(1.0-(speed/speedOfLight)^2);
scaled=scale*betaGamma;
gamma=sqrt(1.0+dot(scaled,scaled));
velocity=speedOfLight*scaled/gamma;
end

function value = localZeroField(~,~,~)
value=zeros(3,1);
end
