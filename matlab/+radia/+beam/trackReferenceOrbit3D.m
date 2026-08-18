function result = trackReferenceOrbit3D(radiaObject,rigidity,entrancePoint,entranceDirection,exitX,options)
%TRACKREFERENCEORBIT3D Track the full-3D design orbit of a Radia object.
%   Native fixed-step RK4 tracker (rad_orbit_tracker.cpp) on the full 3D
%   Lorentz force with a cubic-Hermite exit-plane crossing, station
%   sampling, one batched midpoint-curvature evaluation, and the MEASURED
%   planarity gate (tracking is 3D -- a field that leaves the bend plane
%   fails loudly instead of being silently flattened).
%
%   The MEX route drives Radia-object sources (coils, magnet blocks,
%   containers); the HDiv iron-evaluator term is a pybind-owned handle
%   and joins once an evaluator handle exists in the MEX registry.
%
%   result fields: positions_m/tangents (3 x N), stations_m (N),
%   signed_curvature_per_m (N-1), length_m, out_of_plane_m,
%   out_of_plane_slope.
arguments
    radiaObject (1,1) double {mustBeInteger,mustBePositive}
    rigidity (1,1) double {mustBeReal,mustBeFinite}
    entrancePoint (1,3) double {mustBeReal,mustBeFinite}
    entranceDirection (1,3) double {mustBeReal,mustBeFinite}
    exitX (1,1) double {mustBeReal,mustBeFinite}
    options.MirrorZ (1,1) logical = false
    options.Step (1,1) double {mustBePositive} = 1.0e-3
    options.MaximumPath (1,1) double {mustBePositive} = 0.14
    options.PlanarityTolerance (1,1) double {mustBePositive} = 1.0e-6
    options.StationCount (1,1) double {mustBeInteger,mustBePositive} = 65
end
config = struct();
config.radia_object = double(radiaObject);
config.mirror_z = logical(options.MirrorZ);
config.magnetic_rigidity_t_m = double(rigidity);
config.entrance_point_m = double(entrancePoint);
config.entrance_direction = double(entranceDirection);
config.exit_x_m = double(exitX);
config.step_m = double(options.Step);
config.maximum_path_m = double(options.MaximumPath);
config.planarity_tolerance_m = double(options.PlanarityTolerance);
config.station_count = double(options.StationCount);
result = radia.internal.callMex('beam.orbit.track_reference_3d',config);
end
