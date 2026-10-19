function extensions = ihGeometryExtensions()
%IHGEOMETRYEXTENSIONS Single source of the IH geometry file extensions.
%   extensions.vol / .step / .sol are the accepted (case-insensitive)
%   file suffixes for mesh, CAD, and field-solution inputs.  Every
%   consumer (normalizeIHGeometryRoles, updateIHGeometry's pair
%   classifier, browseIHGeometryFile) reads this table so the accepted
%   formats cannot drift apart between the mask dialogs, the browse
%   filters, and the role-repair logic.

extensions = struct( ...
    "vol", [".vol", ".vol.gz"], ...
    "step", [".step", ".stp"], ...
    "sol", ".sol");
end
