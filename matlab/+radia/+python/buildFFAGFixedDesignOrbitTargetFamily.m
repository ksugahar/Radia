function result = buildFFAGFixedDesignOrbitTargetFamily( ...
    designOrbits, targetTransferMatrices, options)
%BUILDFFAGFIXEDDESIGNORBITTARGETFAMILY Build a one-pass FFAG target.
%   The design orbits are Python PlanarDesignOrbit objects returned by a
%   Radia Python fallback call. This batch-only boundary does not run during
%   a Simulink time step. RESULT.VALUE is a MATLAB inspection view; the
%   complete Python object remains available as RESULT.VALUE.PYTHON_OBJECT.
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
result = radia.internal.callPython( ...
    "radia.ffag_topopt", ...
    "build_ffag_fixed_design_orbit_target_family", ...
    {designOrbits, targetTransferMatrices}, Keywords=keywords);
pythonFamily = result.value;
pythonObjective = pythonFamily.objective;
objective = struct();
objective.python_object = pythonObjective;
objective.orbits = cell(pythonObjective.orbits);
objective.target_matrices = double(pythonObjective.target_matrices);
objective.transfer_matrix_band = double( ...
    pythonObjective.transfer_matrix_band);
bendBands = cell(pythonObjective.bend_field_band);
for index = 1:numel(bendBands)
    bendBands{index} = double(bendBands{index});
end
objective.bend_field_band = bendBands;
objective.response_entries = numericPairMatrix( ...
    pythonObjective.response_entries);
objective.curvature_sign = double(pythonObjective.curvature_sign);
objective.gradient_sign = double(pythonObjective.gradient_sign);

components = cell(pythonFamily.controlled_components);
componentNames = strings(size(components));
for index = 1:numel(components)
    componentNames(index) = string(components{index});
end
value = struct();
value.python_object = pythonFamily;
value.objective = objective;
value.controlled_components = componentNames;
value.target_symplectic_residuals = double( ...
    pythonFamily.target_symplectic_residuals);
value.design_orbits = cell(pythonFamily.design_orbits);
value.magnetic_rigidities_tm = double( ...
    pythonFamily.magnetic_rigidities_tm);
result.value = value;
end

function entries = numericPairMatrix(entries)
if isa(entries, "py.list") || isa(entries, "py.tuple")
    entries = cell(entries);
end
if iscell(entries)
    rows = cell(size(entries));
    for index = 1:numel(entries)
        row = entries{index};
        if isa(row, "py.list") || isa(row, "py.tuple")
            row = double(row);
        elseif iscell(row)
            row = cellfun(@double, row);
        else
            row = double(row);
        end
        rows{index} = reshape(row, 1, []);
    end
    entries = vertcat(rows{:});
end
if ~isnumeric(entries) || size(entries, 2) ~= 2
    error("radia:python:FFAGResponseEntries", ...
        "Python response_entries must be a numeric N-by-2 array.");
end
entries = double(entries);
end
