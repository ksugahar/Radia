function config = validateIHNativeConfig(config)
%VALIDATEIHNATIVECONFIG Validate and normalize the distributed-field IH contract.

arguments
    config (1,1) struct
end

required = [ ...
    "schema", "n_eddy_unknown", "n_heat", "n_temperature", "bh_mode", ...
    "eddy_matrix_real", "eddy_matrix_imag", "eddy_rhs_real", ...
    "eddy_rhs_imag", "heat_projection", "heat_cell_weights", ...
    "heat_to_temperature_projection", ...
    "mass_row_ptr", "mass_col", "mass_value", "stiffness_row_ptr", ...
    "stiffness_col", "stiffness_value", "temperature_cell_weights", ...
    "initial_temperature_K", "sample_time_s", "eddy_solver", ...
    "thermal_solver", "rotation_mode", "angle_origin_rad"];
for name = required
    if ~isfield(config, name) || isempty(config.(name))
        error("radia:simulink:IHConfigMissing", ...
            "Native IH configuration is missing required field '%s'.", name);
    end
end
schema = string(config.schema);
if ~isscalar(schema) || schema ~= "radia.ih.simulink.native_sfunction.v1"
    error("radia:simulink:IHConfigSchema", ...
        "Native IH configuration schema must be radia.ih.simulink.native_sfunction.v1.");
end
config.schema = char(schema);

nUnknown = positiveInteger(config.n_eddy_unknown, "n_eddy_unknown");
nHeat = positiveInteger(config.n_heat, "n_heat");
nTemperature = positiveInteger(config.n_temperature, "n_temperature");
config.n_eddy_unknown = nUnknown;
config.n_heat = nHeat;
config.n_temperature = nTemperature;

mode = lower(string(config.bh_mode));
if mode ~= "linear"
    error("radia:simulink:IHConfigBHMode", ...
        "The IH native preview supports only bh_mode='linear'.");
end
if isfield(config, "backend") && ...
        string(config.backend) ~= "matlab-level2+radia-mex-handles"
    error("radia:simulink:IHConfigBackend", ...
        "IH configuration backend must be matlab-level2+radia-mex-handles.");
end
if isfield(config, "python_fallback") && logical(config.python_fallback)
    error("radia:simulink:IHConfigBackend", ...
        "IH production execution does not permit Python fallback.");
end
config.bh_mode = char(mode);
config.eddy_solver = char(checkedChoice(config.eddy_solver, ...
    ["fem", "peec", "bem-a", "bim"], "eddy_solver"));
config.thermal_solver = char(checkedChoice(config.thermal_solver, ...
    "fem", "thermal_solver"));
config.sample_time_s = positiveScalar(config.sample_time_s, "sample_time_s");
config.rotation_mode = char(checkedChoice(config.rotation_mode, ...
    ["none", "periodic-uniform"], "rotation_mode"));
config.angle_origin_rad = finiteScalar(config.angle_origin_rad, ...
    "angle_origin_rad");

config.eddy_matrix_real = rowMajor(config.eddy_matrix_real, ...
    nUnknown, nUnknown, "eddy_matrix_real");
config.eddy_matrix_imag = rowMajor(config.eddy_matrix_imag, ...
    nUnknown, nUnknown, "eddy_matrix_imag");
config.eddy_rhs_real = finiteVector(config.eddy_rhs_real, ...
    nUnknown, "eddy_rhs_real");
config.eddy_rhs_imag = finiteVector(config.eddy_rhs_imag, ...
    nUnknown, "eddy_rhs_imag");
config.heat_projection = rowMajor(config.heat_projection, ...
    nHeat, nUnknown, "heat_projection");
config.heat_cell_weights = positiveVector(config.heat_cell_weights, ...
    nHeat, "heat_cell_weights");
config.heat_to_temperature_projection = rowMajor( ...
    config.heat_to_temperature_projection, nTemperature, nHeat, ...
    "heat_to_temperature_projection");
config.temperature_cell_weights = positiveVector( ...
    config.temperature_cell_weights, nTemperature, ...
    "temperature_cell_weights");
config.initial_temperature_K = positiveVector(config.initial_temperature_K, ...
    nTemperature, "initial_temperature_K");

[config.mass_row_ptr, config.mass_col, config.mass_value] = checkCSR( ...
    config.mass_row_ptr, config.mass_col, config.mass_value, ...
    nTemperature, "mass");
[config.stiffness_row_ptr, config.stiffness_col, config.stiffness_value] = ...
    checkCSR(config.stiffness_row_ptr, config.stiffness_col, ...
    config.stiffness_value, nTemperature, "stiffness");
if ~isequal(config.mass_row_ptr, config.stiffness_row_ptr) || ...
        ~isequal(config.mass_col, config.stiffness_col)
    error("radia:simulink:IHConfigThermalSparsity", ...
        "Mass and stiffness matrices must use identical CSR sparsity.");
end

convectionNames = ["convection_row_ptr", "convection_col", "convection_value"];
hasConvection = arrayfun(@(name) isfield(config, name), convectionNames);
if any(hasConvection) && ~all(hasConvection)
    error("radia:simulink:IHConfigConvection", ...
        "Convection CSR requires row_ptr, col, and value together.");
end
if all(hasConvection)
    [config.convection_row_ptr, config.convection_col, ...
        config.convection_value] = checkCSR(config.convection_row_ptr, ...
        config.convection_col, config.convection_value, ...
        nTemperature, "convection");
    if ~isequal(config.mass_row_ptr, config.convection_row_ptr) || ...
            ~isequal(config.mass_col, config.convection_col)
        error("radia:simulink:IHConfigThermalSparsity", ...
            "Convection must use the same CSR sparsity as mass and stiffness.");
    end
end

slopeNames = ["eddy_matrix_temperature_slope_real", ...
    "eddy_matrix_temperature_slope_imag"];
hasSlope = arrayfun(@(name) isfield(config, name), slopeNames);
if any(hasSlope) && ~all(hasSlope)
    error("radia:simulink:IHConfigBHSlope", ...
        "Temperature-dependent Eddy slopes require real and imaginary arrays together.");
end
if all(hasSlope)
    slopeSize = nTemperature * nUnknown * nUnknown;
    config.eddy_matrix_temperature_slope_real = finiteVector( ...
        config.eddy_matrix_temperature_slope_real, slopeSize, ...
        "eddy_matrix_temperature_slope_real");
    config.eddy_matrix_temperature_slope_imag = finiteVector( ...
        config.eddy_matrix_temperature_slope_imag, slopeSize, ...
        "eddy_matrix_temperature_slope_imag");
    if ~isfield(config, "bh_reference_temperature_K")
        error("radia:simulink:IHConfigBHReference", ...
            "Temperature-dependent Eddy slopes require bh_reference_temperature_K.");
    end
    config.bh_reference_temperature_K = finiteScalar( ...
        config.bh_reference_temperature_K, "bh_reference_temperature_K");
end

config = optionalPositive(config, "thermal_tolerance", 1e-10, false);
config = optionalPositive(config, "thermal_max_iterations", 500, true);
config = optionalNonnegative(config, "convection_W_per_m2K", 0);
config.configured = true;
config.backend = 'matlab-level2+radia-mex-handles';
config.python_fallback = false;
config.distributed_field = true;
config.surrogate = false;
config.current_change_recomputes_eddy = false;
config.temperature_change_recomputes_eddy = all(hasSlope);
config.temperature_coordinate_system = "workpiece";
if string(config.rotation_mode) == "periodic-uniform"
    config.rotation_transport = "conservative-periodic";
else
    config.rotation_transport = "none";
end
config.dt_order = "eddy;transport(theta_prev,theta_now);thermal";
end

function value = positiveInteger(value, name)
value = finiteScalar(value, name);
if value <= 0 || value ~= fix(value)
    error("radia:simulink:IHConfigInteger", ...
        "%s must be a positive integer.", name);
end
end

function value = positiveScalar(value, name)
value = finiteScalar(value, name);
if value <= 0
    error("radia:simulink:IHConfigPositive", "%s must be positive.", name);
end
end

function value = finiteScalar(value, name)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ~isfinite(value)
    error("radia:simulink:IHConfigScalar", ...
        "%s must be one finite real scalar.", name);
end
value = double(value);
end

function values = finiteVector(values, count, name)
if ~isnumeric(values) || ~isreal(values) || numel(values) ~= count || ...
        any(~isfinite(values), "all")
    error("radia:simulink:IHConfigVector", ...
        "%s must contain %d finite real values.", name, count);
end
values = double(values(:));
end

function values = positiveVector(values, count, name)
values = finiteVector(values, count, name);
if any(values <= 0)
    error("radia:simulink:IHConfigWeights", ...
        "%s must contain only positive values.", name);
end
end

function values = rowMajor(values, rows, columns, name)
if ~isnumeric(values) || ~isreal(values) || numel(values) ~= rows * columns || ...
        any(~isfinite(values), "all")
    error("radia:simulink:IHConfigMatrix", ...
        "%s must be a finite %d-by-%d real matrix.", name, rows, columns);
end
if ~isvector(values) && ~isequal(size(values), [rows, columns])
    error("radia:simulink:IHConfigMatrix", ...
        "%s must be a %d-by-%d matrix or a row-major vector.", ...
        name, rows, columns);
end
if isequal(size(values), [rows, columns]) && rows > 1 && columns > 1
    values = reshape(double(values).', 1, []);
else
    values = reshape(double(values), 1, []);
end
end

function [row, column, value] = checkCSR(row, column, value, n, name)
row = finiteVector(row, n + 1, name + "_row_ptr");
column = finiteVector(column, numel(column), name + "_col");
value = finiteVector(value, numel(value), name + "_value");
if numel(column) ~= numel(value) || any(row ~= fix(row)) || ...
        any(column ~= fix(column)) || row(1) ~= 0 || row(end) ~= numel(column) || ...
        any(diff(row) < 0) || any(column < 0 | column >= n)
    error("radia:simulink:IHConfigCSR", ...
        "%s must be valid zero-based CSR with %d rows.", name, n);
end
end

function value = checkedChoice(value, choices, name)
value = lower(string(value));
if ~isscalar(value) || ~ismember(value, choices)
    error("radia:simulink:IHConfigChoice", ...
        "%s must be one of: %s.", name, strjoin(choices, ", "));
end
end

function config = optionalPositive(config, name, fallback, integerRequired)
if ~isfield(config, name)
    config.(name) = fallback;
end
if integerRequired
    config.(name) = positiveInteger(config.(name), name);
else
    config.(name) = positiveScalar(config.(name), name);
end
end

function config = optionalFinite(config, name, fallback)
if ~isfield(config, name)
    config.(name) = fallback;
end
config.(name) = finiteScalar(config.(name), name);
end

function config = optionalNonnegative(config, name, fallback)
config = optionalFinite(config, name, fallback);
if config.(name) < 0
    error("radia:simulink:IHConfigNonnegative", ...
        "%s must be nonnegative.", name);
end
end
