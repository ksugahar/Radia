function config = makeIHNativeConfig(spec, options)
%MAKEIHNATIVECONFIG Validate the native IH S-Function contract.
%   This function marshals an assembled IHDesignSpec into the checked
%   row-major numeric contract consumed by the native S-Functions.

arguments
    spec (1,1) struct
    options.NHeat (1,1) double {mustBeInteger,mustBePositive} = 1
    options.NTemperature (1,1) double {mustBeInteger,mustBePositive} = 1
    options.CellWeights (:,1) double = []
    options.HeatCellWeights (:,1) double = []
    options.HeatToTemperatureProjection double = []
    options.RotationMode (1,1) string {mustBeMember(options.RotationMode, ...
        ["none","periodic-uniform"])} = "none"
    options.AngleOrigin_rad (1,1) double {mustBeFinite} = 0
    options.SampleTime_s (1,1) double {mustBePositive,mustBeFinite} = 1e-3
    options.WorkpieceVolLabelContract (1,1) string = ""
    options.CoilVolLabelContract (1,1) string = ""
end

required = ["frequency","current","wp_vol"];
missing = required(~isfield(spec, required));
if ~isempty(missing)
    error("radia:simulink:IHConfigMissing", ...
        "IHDesignSpec is missing required fields: %s", strjoin(missing, ", "));
end
for name = required
    if isempty(spec.(name))
        error("radia:simulink:IHConfigMissing", ...
            "IHDesignSpec is missing required field '%s'.", name);
    end
end
hasCoilVol = isfield(spec, "coil_vol") && strlength(string(spec.coil_vol)) > 0;
hasPeecStep = isfield(spec, "peec_step") && strlength(string(spec.peec_step)) > 0;
if ~(hasCoilVol || hasPeecStep)
    error("radia:simulink:IHConfigMissing", ...
        "IHDesignSpec requires coil_vol or peec_step.");
end
if ~isfield(spec, "method") || ~isfield(spec, "solver") || ...
        ~isfield(spec, "thermal_mesh_type")
    error("radia:simulink:IHConfigSolverMissing", ...
        "IHDesignSpec must select method, solver, and thermal_mesh_type.");
end
if ~isfield(spec, "bh_mode") || ~strcmpi(string(spec.bh_mode), "linear")
    error("radia:simulink:IHConfigBHMode", ...
        "The IH native preview supports only bh_mode='linear'.");
end
if isempty(options.CellWeights) || ...
        numel(options.CellWeights) ~= options.NTemperature || ...
        any(~isfinite(options.CellWeights)) || any(options.CellWeights <= 0)
    error("radia:simulink:IHConfigWeights", ...
        "CellWeights must contain one finite positive value per temperature DOF.");
end
if isempty(options.HeatCellWeights) || ...
        numel(options.HeatCellWeights) ~= options.NHeat || ...
        any(~isfinite(options.HeatCellWeights)) || ...
        any(options.HeatCellWeights <= 0)
    error("radia:simulink:IHConfigHeatWeights", ...
        "HeatCellWeights must contain one finite positive value per heat region.");
end

nativeRequired = ["n_eddy_unknown","eddy_matrix_real", ...
    "eddy_matrix_imag","eddy_rhs_real","eddy_rhs_imag", ...
    "heat_projection","mass_row_ptr","mass_col","mass_value", ...
    "stiffness_row_ptr","stiffness_col","stiffness_value", ...
    "initial_temperature_K"];
missingNative = nativeRequired(~isfield(spec, nativeRequired));
if ~isempty(missingNative)
    error("radia:simulink:IHConfigNativeOperator", ...
        "IHDesignSpec is missing native operator data: %s", ...
        strjoin(missingNative, ", "));
end
nUnknown = double(spec.n_eddy_unknown);
if ~(isscalar(nUnknown) && isfinite(nUnknown) && nUnknown > 0 && ...
        fix(nUnknown) == nUnknown)
    error("radia:simulink:IHConfigNativeOperator", ...
        "n_eddy_unknown must be a positive integer.");
end

config = spec;
config.schema = "radia.ih.simulink.native_sfunction.v1";
config.n_heat = options.NHeat;
config.n_temperature = options.NTemperature;
config.temperature_cell_weights = options.CellWeights;
config.heat_cell_weights = options.HeatCellWeights;
config.eddy_matrix_real = rowMajor( ...
    spec.eddy_matrix_real, nUnknown, nUnknown, "eddy_matrix_real");
config.eddy_matrix_imag = rowMajor( ...
    spec.eddy_matrix_imag, nUnknown, nUnknown, "eddy_matrix_imag");
config.eddy_rhs_real = finiteVector( ...
    spec.eddy_rhs_real, nUnknown, "eddy_rhs_real");
config.eddy_rhs_imag = finiteVector( ...
    spec.eddy_rhs_imag, nUnknown, "eddy_rhs_imag");
config.heat_projection = rowMajor( ...
    spec.heat_projection, options.NHeat, nUnknown, "heat_projection");

projection = options.HeatToTemperatureProjection;
if isempty(projection)
    if options.NHeat ~= options.NTemperature
        error("radia:simulink:IHConfigHeatProjection", ...
            "HeatToTemperatureProjection is required when NHeat differs from NTemperature.");
    end
    projection = eye(options.NTemperature);
end
if ~isequal(size(projection), [options.NTemperature, options.NHeat]) || ...
        any(~isfinite(projection), "all")
    error("radia:simulink:IHConfigHeatProjection", ...
        "HeatToTemperatureProjection must be an NTemperature-by-NHeat finite matrix.");
end
config.heat_to_temperature_projection = reshape(projection.', [], 1);
config.sample_time_s = options.SampleTime_s;
config.rotation_mode = char(options.RotationMode);
config.angle_origin_rad = options.AngleOrigin_rad;
config.backend = "native-mex-sfunction";
config.python_fallback = false;
config.release_channel = "preview";
config.operator_assembly = "preassembled";
config.eddy_solver = physicalEddySolver(spec.method);
config.eddy_method = string(spec.method);
config.linear_solver = string(spec.solver);
config.thermal_solver = "fem";
config.thermal_mesh_type = string(spec.thermal_mesh_type);
config.bh_mode = char(lower(string(spec.bh_mode)));
config.current_change_recomputes_eddy = false;
config.temperature_coordinate_system = "workpiece";
config.rotation_transport = options.RotationMode;
config.dt_order = "eddy;transport(theta_prev,theta_now);thermal";

workpieceContract = options.WorkpieceVolLabelContract;
if strlength(workpieceContract) == 0 && isfield(spec, "wp_vol_label_contract")
    workpieceContract = string(spec.wp_vol_label_contract);
end
if strlength(workpieceContract) == 0 || ~isfile(workpieceContract)
    error("radia:simulink:IHConfigVolContract", ...
        "A versioned workpiece .vol label contract is required.");
end
reports = radia.simulink.validateVolFiles(string(spec.wp_vol), ...
    Contract=workpieceContract);
config.workpiece_vol_label_contract = workpieceContract;
if hasCoilVol
    coilContract = options.CoilVolLabelContract;
    if strlength(coilContract) == 0 && isfield(spec, "coil_vol_label_contract")
        coilContract = string(spec.coil_vol_label_contract);
    end
    if strlength(coilContract) == 0 || ~isfile(coilContract)
        error("radia:simulink:IHConfigVolContract", ...
            "A versioned coil .vol label contract is required for coil_vol.");
    end
    reports(end + 1, 1) = radia.simulink.validateVolFiles( ...
        string(spec.coil_vol), Contract=coilContract);
    config.coil_vol_label_contract = coilContract;
end
config.vol_check_reports = reports;
config.vol_check_required = true;

if isfield(spec, "bh_reference_temperature_K")
    config.bh_reference_temperature_K = spec.bh_reference_temperature_K;
end
if xor(isfield(spec, "eddy_matrix_temperature_slope_real"), ...
       isfield(spec, "eddy_matrix_temperature_slope_imag"))
    error("radia:simulink:IHConfigBHSlope", ...
        "Temperature-dependent Eddy matrix slopes require real and imaginary parts together.");
end
if isfield(spec, "eddy_matrix_temperature_slope_real")
    config.eddy_matrix_temperature_slope_real = matrixStackRowMajor( ...
        spec.eddy_matrix_temperature_slope_real, nUnknown, ...
        options.NTemperature, "eddy_matrix_temperature_slope_real");
    config.eddy_matrix_temperature_slope_imag = matrixStackRowMajor( ...
        spec.eddy_matrix_temperature_slope_imag, nUnknown, ...
        options.NTemperature, "eddy_matrix_temperature_slope_imag");
end
config = radia.simulink.validateIHNativeConfig(config);
end

function solver = physicalEddySolver(method)
method = lower(string(method));
if contains(method, "bem-a")
    solver = "bem-a";
elseif contains(method, "bim")
    solver = "bim";
elseif contains(method, "peec")
    solver = "peec";
elseif contains(method, "fem") || contains(method, "full simulation") || ...
        startsWith(method, "thermal:")
    solver = "fem";
else
    error("radia:simulink:IHConfigMethod", ...
        "IHDesignSpec.method cannot be mapped to FEM, PEEC, BEM-A, or BIM: %s", ...
        method);
end
end

function values = finiteVector(value, count, name)
if ~isreal(value) || numel(value) ~= count || any(~isfinite(value), "all")
    error("radia:simulink:IHConfigNativeOperator", ...
        "%s must contain %d finite real values.", name, count);
end
values = double(value(:));
end

function values = rowMajor(value, rows, columns, name)
if ~isreal(value) || numel(value) ~= rows * columns || ...
        any(~isfinite(value), "all")
    error("radia:simulink:IHConfigNativeOperator", ...
        "%s must be a finite %d-by-%d real matrix.", name, rows, columns);
end
if rows > 1 && columns > 1 && ~isequal(size(value), [rows, columns])
    error("radia:simulink:IHConfigNativeOperator", ...
        "%s must have size %d-by-%d before row-major conversion.", ...
        name, rows, columns);
end
values = reshape(double(value).', [], 1);
end

function values = matrixStackRowMajor(value, order, count, name)
if ~isreal(value) || numel(value) ~= count * order * order || ...
        any(~isfinite(value), "all")
    error("radia:simulink:IHConfigBHSlopeShape", ...
        "%s has the wrong size or contains non-finite values.", name);
end
if isvector(value)
    values = double(value(:));
    return
end
if size(value,1) ~= order || size(value,2) ~= order || ...
        size(value,3) ~= count || ndims(value) > 3
    error("radia:simulink:IHConfigBHSlopeShape", ...
        "%s must be order-by-order-by-n_temperature.", name);
end
values = zeros(numel(value), 1);
block = order * order;
for index = 1:count
    offset = (index - 1) * block;
    values(offset + (1:block)) = reshape(double(value(:,:,index)).', [], 1);
end
end
