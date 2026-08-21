function result = optimizeFFAGHDivMMMFromDesignOrbits( ...
    designOrbits, targetTransferMatrices, options)
%OPTIMIZEFFAGHDIVMMMFROMDESIGNORBITS Run the one-pass HDiv-MMM inverse.
%   Keywords contains the native/Python solver objects and optimization
%   settings accepted by optimize_ffag_hdiv_mmm_from_design_orbits. This is
%   an explicit setup or batch-solve boundary, never a step-time backend.
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
    options.Keywords (1,1) struct = struct()
end
keywords = options.Keywords;
keywords.transfer_matrix_band = options.TransferMatrixBand;
keywords.bend_field_band = options.BendFieldBand;
keywords.require_symplectic = options.RequireSymplectic;
keywords.symplectic_tolerance = options.SymplecticTolerance;
keywords.curvature_sign = options.CurvatureSign;
keywords.gradient_sign = options.GradientSign;
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
result = radia.internal.callPython( ...
    "radia.ffag_topopt", ...
    "optimize_ffag_hdiv_mmm_from_design_orbits", ...
    {designOrbits, targetTransferMatrices}, Keywords=keywords);
end
