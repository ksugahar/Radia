function payload = writeIHPanelImpedance(templateFile, outputFile, Zs)
%WRITEIHPANELIMPEDANCE Set one fixed complex impedance [ohm] per BND triangle.
% Generate templateFile with python -m radia.surface_impedance first. Select
% panels using payload.centroids_m; preserve its mesh identity and frequency.
% Set Geometry Update's panel_zs_file to outputFile before rebuilding.
arguments
    templateFile (1,1) string
    outputFile (1,1) string
    Zs (:,1) double
end
payload = jsondecode(fileread(templateFile));
if ~isfield(payload, "schema") || string(payload.schema) ~= "radia.panel_surface_impedance.v1" || ...
        ~isfield(payload, "layout") || string(payload.layout) ~= "BND-element-order" || ...
        ~isfield(payload, "unit") || string(payload.unit) ~= "ohm" || ...
        ~isfield(payload, "surface_sha256") || ~isfield(payload, "centroids_m") || ...
        ~isfield(payload, "frequency_hz") || ~isscalar(payload.frequency_hz) || ...
        ~isfinite(payload.frequency_hz) || payload.frequency_hz <= 0
    error("radia:simulink:IHPanelZsTemplate", "Invalid mesh-bound panel Zs template.");
end
if size(payload.centroids_m,2) ~= 3 || any(~isfinite(payload.centroids_m), "all") || ...
        numel(Zs) ~= size(payload.centroids_m,1) || isempty(Zs) || ...
        any(~isfinite(Zs)) || any(real(Zs) < 0)
    error("radia:simulink:IHPanelZsValues", "Provide one finite passive Zs per template triangle.");
end
payload.real_ohm = real(Zs);
payload.imag_ohm = imag(Zs);
file = fopen(outputFile, "w", "n", "UTF-8");
if file < 0
    error("radia:simulink:IHPanelZsOutput", "Cannot write %s", outputFile);
end
cleanup = onCleanup(@() fclose(file));
fwrite(file, jsonencode(payload), "char");
end
