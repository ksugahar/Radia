function state = makeFieldState(temperature, magneticField, options)
%MAKEFIELDSTATE Create a workpiece-coordinate field-state contract.
arguments
    temperature (:,1) double
    magneticField (:,1) double
    options.MeshId (1,1) string = ""
    options.CoordinateSystem (1,1) string = "workpiece"
    options.Time_s (1,1) double = 0
end
if numel(temperature) ~= numel(magneticField), error("radia:material:FieldStateSize", "T(x) and H(x) must have equal lengths."); end
if options.CoordinateSystem ~= "workpiece", error("radia:material:FieldStateCoordinates", "FieldState must use workpiece coordinates."); end
if any(~isfinite(temperature)) || any(~isfinite(magneticField)) || ~isfinite(options.Time_s), error("radia:material:FieldStateFinite", "FieldState values must be finite."); end
state = struct("schema","radia.material.field_state.v1", "temperature_K",temperature, ...
    "H_A_per_m",magneticField,"mesh_id",options.MeshId, ...
    "coordinate_system",options.CoordinateSystem,"time_s",options.Time_s);
end
