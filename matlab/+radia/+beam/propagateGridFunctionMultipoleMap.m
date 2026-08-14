function result = propagateGridFunctionMultipoleMap( ...
    field, lengthsM, referencePositionsM, referenceTangents, ...
    magneticRigidityTM, options)
%PROPAGATEGRIDFUNCTIONMULTIPOLEMAP Moving-frame nonlinear beam map.
%   RESULT = radia.beam.propagateGridFunctionMultipoleMap(FIELD, LENGTHS,
%   POSITIONS, TANGENTS, BRHO) evaluates the live real three-component
%   NGSolve FIELD at the center and eight angles of a transverse ring at each
%   reference station. The shared C++ kernel fits
%
%       By + i*Bx = sum(Cn*(x + i*y)^n)
%
%   through MultipoleOrder, builds the paraxial A/F2/F3 jet in
%   (x,px/p0,y,py/p0,sigma,delta), and propagates the region-attributed R/T/U
%   map through MaximumMapOrder. Raw samples, transported frames, multipole
%   coefficients, fit residuals, and local jets remain in RESULT.
%
%   This is a local source-free transverse body-field approximation with a
%   piecewise-constant jet per supplied segment. It is not a complete curved-
%   coordinate Hamiltonian and does not add longitudinal fringe/edge
%   derivatives, closed-orbit finding, or a general symplectic Lie map. Use
%   radia.beam.trackGridFunction as the independent point-evaluation check.

arguments
    field (1,1) radia.ngsolve.GridFunction
    lengthsM double {mustBeReal,mustBeFinite,mustBeNonempty}
    referencePositionsM double {mustBeReal,mustBeFinite,mustBeNonempty}
    referenceTangents double {mustBeReal,mustBeFinite,mustBeNonempty}
    magneticRigidityTM (1,1) double {mustBeReal,mustBeFinite}
    options.SampleRadiusM (1,1) double {mustBeFinite,mustBePositive} = 1e-3
    options.InitialHorizontal double {mustBeReal,mustBeFinite} = [1 0 0]
    options.Names = strings(0,1)
    options.CurvatureSign (1,1) double {mustBeFinite} = 1
    options.GradientSign (1,1) double {mustBeFinite} = 1
    options.MultipoleOrder (1,1) double {mustBeInteger, ...
        mustBeMember(options.MultipoleOrder,[1 2 3])} = 3
    options.MaximumMapOrder (1,1) double {mustBeInteger, ...
        mustBeMember(options.MaximumMapOrder,[1 2 3])} = 3
    options.MaximumStepM (1,1) double {mustBeFinite,mustBePositive} = 1e-3
    options.MaximumSteps (1,1) double {mustBeInteger,mustBePositive} = 1e6
end

lengthsM = double(lengthsM(:));
segmentCount = numel(lengthsM);
if any(lengthsM <= 0)
    error("radia:beam:InvalidLength", ...
        "Every segment length must be positive.");
end
if magneticRigidityTM == 0
    error("radia:beam:InvalidRigidity", ...
        "magneticRigidityTM must be nonzero.");
end
if ~isequal(size(referencePositionsM),[segmentCount 3])
    error("radia:beam:InvalidPositionShape", ...
        "referencePositionsM must have shape N-by-3.");
end
if ~isequal(size(referenceTangents),[segmentCount 3])
    error("radia:beam:InvalidTangentShape", ...
        "referenceTangents must have shape N-by-3.");
end
initialHorizontal = double(options.InitialHorizontal(:));
if numel(initialHorizontal) ~= 3 || norm(initialHorizontal) == 0
    error("radia:beam:InvalidHorizontalAxis", ...
        "InitialHorizontal must contain a nonzero three-vector.");
end

config.schema = 'radia.beam.grid-function-multipole-map.v1';
config.lengths_m = lengthsM;
config.reference_positions_m = double(referencePositionsM);
config.reference_tangents = double(referenceTangents);
config.magnetic_rigidity_t_m = magneticRigidityTM;
config.sample_radius_m = options.SampleRadiusM;
config.initial_horizontal = initialHorizontal;
if ~isempty(options.Names)
    names = string(options.Names(:));
    if numel(names) ~= segmentCount || any(ismissing(names))
        error("radia:beam:InvalidNames", ...
            "Names must contain one nonmissing entry per segment.");
    end
    config.names = cellstr(names);
end
config.curvature_sign = options.CurvatureSign;
config.gradient_sign = options.GradientSign;
config.multipole_order = options.MultipoleOrder;
config.maximum_map_order = options.MaximumMapOrder;
config.maximum_step_m = options.MaximumStepM;
config.maximum_steps = options.MaximumSteps;

result = radia.internal.callMex( ...
    'beam.transfer.multipole_from_grid_function', ...
    field.nativeHandle(),config);
end
