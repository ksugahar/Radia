function result = buildFFAGFixedDesignOrbitTargetFamily( ...
    designOrbits, targetTransferMatrices, options)
%BUILDFFAGFIXEDDESIGNORBITTARGETFAMILY Build a one-pass FFAG target.
%   The design orbits are Python PlanarDesignOrbit objects returned by a
%   Radia Python fallback call. This batch-only boundary does not run during
%   a Simulink time step.
arguments
    designOrbits (1,:) cell
    targetTransferMatrices {mustBeNumeric}
    options.TransferMatrixBand = 1e-3
    options.BendFieldBand = 1e-3
    options.ResponseEntries = []
    options.ControlledComponents = []
    options.RequireSymplectic (1,1) logical = true
    options.SymplecticTolerance (1,1) double {mustBeNonnegative} = 1e-9
    options.CurvatureSign (1,1) double = 1
    options.GradientSign (1,1) double = 1
end
keywords = struct( ...
    "transfer_matrix_band", options.TransferMatrixBand, ...
    "bend_field_band", options.BendFieldBand, ...
    "require_symplectic", options.RequireSymplectic, ...
    "symplectic_tolerance", options.SymplecticTolerance, ...
    "curvature_sign", options.CurvatureSign, ...
    "gradient_sign", options.GradientSign);
if ~isempty(options.ResponseEntries) && ~isempty(options.ControlledComponents)
    error("radia:python:MutuallyExclusiveOptions", ...
        "ResponseEntries and ControlledComponents are mutually exclusive.");
end
if ~isempty(options.ResponseEntries)
    keywords.response_entries = options.ResponseEntries;
end
if ~isempty(options.ControlledComponents)
    keywords.controlled_components = options.ControlledComponents;
end
result = radia.python.ffagTopopt( ...
    "build_ffag_fixed_design_orbit_target_family", ...
    {designOrbits, targetTransferMatrices}, Keywords=keywords);
end
