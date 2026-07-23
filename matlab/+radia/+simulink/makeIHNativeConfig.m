function config = makeIHNativeConfig(spec, options)
%MAKEIHNATIVECONFIG Validate the native IH S-Function contract.
%   This function only marshals an IHDesignSpec into a checked MATLAB struct;
%   it does not solve the electromagnetic or thermal problem.

arguments
    spec (1,1) struct
    options.NHeat (1,1) double {mustBeInteger,mustBePositive} = 1
    options.NTemperature (1,1) double {mustBeInteger,mustBePositive} = 1
    options.CellWeights (:,1) double = []
    options.SampleTime_s (1,1) double {mustBePositive,mustBeFinite} = 1e-3
end

required = ["frequency", "current", "coil_vol", "wp_vol"];
for name = required
    if ~isfield(spec, name) || isempty(spec.(name))
        error("radia:simulink:IHConfigMissing", ...
            "IHDesignSpec is missing required field '%s'.", name);
    end
end
if ~isfield(spec, "method") || ~isfield(spec, "solver") || ...
        ~isfield(spec, "thermal_mesh_type")
    error("radia:simulink:IHConfigSolverMissing", ...
        "IHDesignSpec must select method, solver, and thermal_mesh_type.");
end
if ~isfield(spec, "bh_mode") || ~(strcmpi(string(spec.bh_mode), "linear") || ...
        strcmpi(string(spec.bh_mode), "nonlinear"))
    error("radia:simulink:IHConfigBHMode", ...
        "IHDesignSpec.bh_mode must be 'linear' or 'nonlinear'.");
end
if isempty(options.CellWeights) || numel(options.CellWeights) ~= options.NTemperature || ...
        any(~isfinite(options.CellWeights)) || any(options.CellWeights <= 0)
    error("radia:simulink:IHConfigWeights", ...
        "CellWeights must contain one finite positive value per temperature DOF.");
end

config = spec;
config.schema = "radia.ih.simulink.native_sfunction.v1";
config.n_heat = options.NHeat;
config.n_temperature = options.NTemperature;
config.temperature_cell_weights = options.CellWeights;
config.sample_time_s = options.SampleTime_s;
config.backend = "native-mex-sfunction";
config.python_fallback = false;
config.eddy_solver = string(spec.solver);
config.thermal_solver = string(spec.method);
config.bh_mode = char(lower(string(spec.bh_mode)));
config.current_change_recomputes_eddy = true;
config.temperature_coordinate_system = "workpiece";
config.rotation_transport = "conservative-periodic";
config.dt_order = "eddy;transport(theta_prev,theta_now);thermal";
vol_files = string({spec.coil_vol,spec.wp_vol});
if any(isfile(vol_files))
    existing = vol_files(isfile(vol_files));
    config.vol_check_reports = radia.simulink.validateVolFiles(existing);
    config.vol_check_required = true;
else
    config.vol_check_required = true;
end
if isfield(spec, "bh_reference_temperature_K")
    config.bh_reference_temperature_K = spec.bh_reference_temperature_K;
end
if xor(isfield(spec, "eddy_matrix_temperature_slope_real"), ...
       isfield(spec, "eddy_matrix_temperature_slope_imag"))
    error("radia:simulink:IHConfigBHSlope", ...
        "Temperature-dependent Eddy matrix slopes require real and imaginary parts together.");
end
if isfield(spec, "eddy_matrix_temperature_slope_real")
    expected = options.NTemperature * spec.n_eddy_unknown^2;
    if numel(spec.eddy_matrix_temperature_slope_real) ~= expected || ...
            numel(spec.eddy_matrix_temperature_slope_imag) ~= expected
        error("radia:simulink:IHConfigBHSlopeShape", ...
            "Temperature-dependent Eddy matrix slopes must have n_temperature*n_eddy_unknown^2 entries.");
    end
    config.eddy_matrix_temperature_slope_real = spec.eddy_matrix_temperature_slope_real;
    config.eddy_matrix_temperature_slope_imag = spec.eddy_matrix_temperature_slope_imag;
end
end
