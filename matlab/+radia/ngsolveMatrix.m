function [A, info] = ngsolveMatrix(volPath, space, order, form, options)
%NGSOLVEMATRIX Assemble a native NGSolve volume matrix as MATLAB sparse data.
%
%   A = radia.ngsolveMatrix(volPath, "hcurl", 3, "mass")
%   [A, info] = radia.ngsolveMatrix(volPath, "hcurl", 3, "stiffness")
%
% The finite-element space, element mappings, quadrature, and global DoF
% ordering are owned by NGSolve.  The returned matrix uses that exact DoF
% ordering.  INFO contains the mesh/form metadata and the structural and
% numeric nonzero counts returned by the native gateway.

arguments
    volPath (1,1) string
    space (1,1) string
    order (1,1) double {mustBeInteger, mustBePositive}
    form (1,1) string
    options.NoGrads (1,1) logical = true
end

space = lower(space);
form = lower(form);
if ~ismember(space, ["h1", "hcurl", "hdiv"])
    error("radia:ngsolveMatrix:Space", ...
        "space must be ""h1"", ""hcurl"", or ""hdiv"".");
end
if ~ismember(form, ["mass", "stiffness", "curlcurl", "curl_curl", ...
        "divdiv", "div_div"])
    error("radia:ngsolveMatrix:Form", "Unsupported NGSolve form ""%s"".", form);
end

raw = radia.internal.callMex( ...
    'ngsolve.matrix_dump', char(volPath), char(space), order, char(form), ...
    options.NoGrads);
A = sparse(raw.row, raw.col, raw.values, raw.shape(1), raw.shape(2));
if nargout > 1
    info = rmfield(raw, {'row', 'col', 'values'});
end
end
