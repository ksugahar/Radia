function result = trackReferenceOrbitToPlane(fieldEvaluator,rigidity, ...
        entrancePoint,entranceDirection,exitPlaneNormal,exitPlaneOffset,options)
%TRACKREFERENCEORBITTOPLANE Track an HDiv-MMM orbit to an arbitrary plane.
%   The fixed-step 3-D RK4 integration, composite-field evaluation, cubic-
%   Hermite plane crossing, station sampling, and midpoint curvature all use
%   the same rad_orbit_tracker.cpp kernel as Python.  The native field is
%
%     IronScale * fieldEvaluator + optional Radia object + ConstantField.
%
%   RIGIDITY is the positive physical |B rho|.  CurvatureSign=+1 follows the
%   accelerator-topology convention in which positive Bz gives positive
%   signed curvature; the wrapper performs the C++ charge-sign conversion.
%
%   result fields: positions_m/tangents (3 x N), stations_m (N),
%   signed_curvature_per_m (N-1), length_m, out_of_plane_m,
%   out_of_plane_slope, magnetic_rigidity_t_m, curvature_sign.
arguments
    fieldEvaluator (1,1) radia.HDivFieldEvaluator
    rigidity (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
    entrancePoint (1,3) double {mustBeReal,mustBeFinite}
    entranceDirection (1,3) double {mustBeReal,mustBeFinite}
    exitPlaneNormal (1,3) double {mustBeReal,mustBeFinite}
    exitPlaneOffset (1,1) double {mustBeReal,mustBeFinite}
    options.IronScale (1,1) double {mustBeReal,mustBeFinite} = 1.0e-7
    options.IronAlgorithm (1,1) string = "auto"
    options.RadiaObject (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.MirrorZ (1,1) logical = false
    options.ConstantField (1,3) double {mustBeReal,mustBeFinite} = [0,0,0]
    options.CurvatureSign (1,1) double {mustBeReal,mustBeFinite} = 1.0
    options.Step (1,1) double {mustBePositive} = 1.0e-3
    options.MaximumPath (1,1) double {mustBePositive} = 0.14
    options.PlanarityTolerance (1,1) double {mustBePositive} = 1.0e-6
    options.StationCount (1,1) double {mustBeInteger,mustBePositive} = 65
end
if abs(options.CurvatureSign) ~= 1.0
    error("radia:beam:CurvatureSign", ...
        "CurvatureSign must be +1 or -1.");
end
if ~ismember(options.IronAlgorithm,["auto","direct","tree"])
    error("radia:beam:IronAlgorithm", ...
        "IronAlgorithm must be auto, direct, or tree.");
end
config = struct();
config.iron_scale = double(options.IronScale);
config.iron_algorithm = char(options.IronAlgorithm);
if options.RadiaObject == 0
    config.radia_object = -1.0;
else
    config.radia_object = double(options.RadiaObject);
end
config.mirror_z = logical(options.MirrorZ);
config.constant_field_t = double(options.ConstantField);
config.magnetic_rigidity_t_m = -double(rigidity) / options.CurvatureSign;
config.entrance_point_m = double(entrancePoint);
config.entrance_direction = double(entranceDirection);
config.exit_plane_normal = double(exitPlaneNormal);
config.exit_plane_offset_m = double(exitPlaneOffset);
config.step_m = double(options.Step);
config.maximum_path_m = double(options.MaximumPath);
config.planarity_tolerance_m = double(options.PlanarityTolerance);
config.station_count = double(options.StationCount);
result = radia.internal.callMex( ...
    'beam.orbit.track_reference_to_plane',fieldEvaluator.nativeHandle(),config);
result.magnetic_rigidity_t_m = double(rigidity);
result.curvature_sign = double(options.CurvatureSign);
end
