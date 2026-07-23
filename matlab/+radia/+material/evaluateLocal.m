function local = evaluateLocal(database, fieldState)
%EVALUATELOCAL Evaluate local material laws at T(x), H(x).
database = radia.material.validateDatabase(database);
if ~isstruct(fieldState) || string(fieldState.schema) ~= "radia.material.field_state.v1"
    error("radia:material:FieldStateSchema", "Invalid FieldState schema.");
end
T = fieldState.temperature_K(:); H = fieldState.H_A_per_m(:);
Ta = database.bh_temperature_K(:); Ha = database.bh_H_A_per_m(:);
n = numel(T); B = zeros(n,1); dBdH = zeros(n,1); sigma = zeros(n,1);
for k = 1:n
    B(k) = interp2(Ha,Ta,database.bh_B_T,H(k),T(k),'linear');
    dBdH(k) = interp2(Ha,Ta,database.bh_dBdH,H(k),T(k),'linear');
    if isscalar(database.conductivity_S_per_m), sigma(k) = database.conductivity_S_per_m;
    else, sigma(k) = interp2(Ha,Ta,database.conductivity_S_per_m,H(k),T(k),'linear'); end
end
if any(~isfinite([B;dBdH;sigma])) || any(dBdH <= 0) || any(sigma <= 0)
    error("radia:material:LocalCoefficients", "Local material evaluation is invalid.");
end
local = struct("schema","radia.material.local_coefficients.v1", ...
    "B_T",B,"dBdH_T_per_Apm",dBdH,"conductivity_S_per_m",sigma, ...
    "coordinate_system","workpiece","mesh_id",fieldState.mesh_id);
end
