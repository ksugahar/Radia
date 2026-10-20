function database = validateDatabase(database)
%VALIDATEDATABASE Validate the immutable spatial material-law contract.
arguments
    database (1,1) struct
end
required = ["schema","material_id","conductivity_S_per_m", ...
    "bh_temperature_K","bh_H_A_per_m","bh_B_T","bh_dBdH", ...
    "coordinate_system"];
for name = required
    if ~isfield(database,name), error("radia:material:MissingField", ...
            "MaterialDatabase is missing '%s'.",name); end
end
if string(database.schema) ~= "radia.material.database.v1"
    error("radia:material:Schema", "Unsupported MaterialDatabase schema.");
end
if string(database.coordinate_system) ~= "workpiece"
    error("radia:material:CoordinateSystem", "MaterialDatabase must use workpiece coordinates.");
end
T = database.bh_temperature_K(:); H = database.bh_H_A_per_m(:);
B = database.bh_B_T; dBdH = database.bh_dBdH;
if any(~isfinite(T)) || any(diff(T) <= 0) || any(~isfinite(H)) || any(diff(H) <= 0)
    error("radia:material:Axes", "BH axes must be finite and increasing.");
end
expected = [numel(T),numel(H)];
if ~isequal(size(B),expected) || ~isequal(size(dBdH),expected) || ...
        any(~isfinite(B(:))) || any(~isfinite(dBdH(:))) || any(dBdH(:) <= 0)
    error("radia:material:BH", "BH tables have invalid dimensions or values.");
end
if ~(isscalar(database.conductivity_S_per_m) || isequal(size(database.conductivity_S_per_m),expected))
    error("radia:material:Conductivity", "Conductivity must be scalar or temperature-by-H sized.");
end
end
