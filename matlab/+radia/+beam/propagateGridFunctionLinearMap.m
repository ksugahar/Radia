function result = propagateGridFunctionLinearMap( ...
    field, lengthsM, referencePositionsM, referenceTangents, ...
    magneticRigidityTM, options)
%PROPAGATEGRIDFUNCTIONLINEARMAP Transfer map from an NGSolve GridFunction.
%   RESULT = radia.beam.propagateGridFunctionLinearMap(FIELD, LENGTHS,
%   POSITIONS, TANGENTS, BRHO) evaluates the real three-component FIELD
%   directly through NGSolve at nine transverse points per reference station.
%   POSITIONS and TANGENTS are N-by-3, with one row per positive segment
%   length. No regular-grid field map is created.
%
%   The native C++ result exposes the transported local frame, center field,
%   fitted field gradient, curvature, normal/skew quadrupole strengths,
%   Maxwell-residual diagnostics, local generators, and accumulated 6-by-6 R
%   map. This entry point is first order; T and U are identically zero.

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

config.schema = 'radia.beam.grid-function-linear-map.v1';
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
config.maximum_step_m = options.MaximumStepM;
config.maximum_steps = options.MaximumSteps;

result = radia.internal.callMex( ...
    'beam.transfer.from_grid_function',field.nativeHandle(),config);
end
