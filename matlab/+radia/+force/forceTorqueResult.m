function result = forceTorqueResult(force_N, torque_Nm, options)
%FORCETORQUERESULT Build the shared radia.force-result/v1 struct.

arguments
    force_N = []
    torque_Nm = []
    options.Method (1,1) string
    options.Frame (1,1) string = "global_cartesian"
    options.Pivot_m = []
    options.FieldConvention (1,1) string = "static"
    options.Amplitude = []
    options.Dimensionality (1,1) string = "3d"
    options.PerUnitDepth (1,1) logical = false
end

force = localOptionalVector(force_N, "force_N");
torque = localOptionalVector(torque_Nm, "torque_Nm");
if isempty(force) && isempty(torque)
    error("radia:force:Result", "at least one of force_N or torque_Nm must be provided");
end
pivot = localOptionalVector(options.Pivot_m, "Pivot_m");
convention = lower(strtrim(options.FieldConvention));
if ~any(convention == ["static", "time_average_phasor"])
    error("radia:force:FieldConvention", "FieldConvention must be static or time_average_phasor");
end
amplitude = options.Amplitude;
if convention == "time_average_phasor"
    amplitude = lower(strtrim(string(amplitude)));
    if ~any(amplitude == ["peak", "rms"])
        error("radia:force:PhasorAmplitude", "Amplitude must be peak or rms for phasor results");
    end
elseif ~isempty(amplitude)
    error("radia:force:PhasorAmplitude", "Amplitude applies only to time_average_phasor results");
else
    amplitude = [];
end
dimension = lower(strtrim(options.Dimensionality));
if ~any(dimension == ["3d", "2d_planar", "axisymmetric"])
    error("radia:force:Dimensionality", "Dimensionality must be 3d, 2d_planar, or axisymmetric");
end
if options.PerUnitDepth
    forceUnit = "N/m";
    torqueUnit = "N";
else
    forceUnit = "N";
    torqueUnit = "N m";
end
result = struct( ...
    "schema", "radia.force-result/v1", ...
    "method", options.Method, ...
    "frame", options.Frame, ...
    "dimensionality", dimension, ...
    "per_unit_depth", options.PerUnitDepth, ...
    "field_convention", convention, ...
    "phasor_amplitude", amplitude, ...
    "pivot_m", pivot, ...
    "force_N", force, ...
    "torque_Nm", torque, ...
    "units", struct("force", forceUnit, "torque", torqueUnit, "pivot", "m"));
end

function vector = localOptionalVector(value, name)
if isempty(value)
    vector = [];
    return
end
if ~isnumeric(value) || ~isreal(value) || numel(value) ~= 3 || any(~isfinite(value), "all")
    error("radia:force:Vector", "%s must be one finite real three-vector", name);
end
vector = reshape(double(value), 1, 3);
end
