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
    options.CurvatureSign (1,1) double = 1
    options.GradientSign (1,1) double = 1
end
keywords = struct( ...
    "transfer_matrix_band", options.TransferMatrixBand, ...
    "bend_field_band", options.BendFieldBand, ...
    "curvature_sign", options.CurvatureSign, ...
    "gradient_sign", options.GradientSign);
if ~isempty(options.ResponseEntries)
    keywords.response_entries = options.ResponseEntries;
end
result = radia.python.ffagTopopt( ...
    "build_ffag_fixed_design_orbit_target_family", ...
    {designOrbits, targetTransferMatrices}, Keywords=keywords);
end
