function out = hcurl_eddy_cln_native_basis(vol_path, order, ports, steps, options)
%HCURL_EDDY_CLN_NATIVE_BASIS Build a native C++ HCurl response basis.
%   OUT = radia.ngsolve.hcurl_eddy_cln_native_basis(VOL_PATH, ORDER, PORTS,
%   STEPS) calls the radia_mex gateway directly. NGSolve owns mesh loading,
%   HCurl orientation, sparse assembly, free-DoF handling, and factorization;
%   the MEX gateway performs the response compression without Python. The
%   result also contains `mass_gram = V'*M*V`, `curlcurl_gram = V'*K*V`, and
%   `port_rhs = V'*ports`.

arguments
    vol_path (1,1) string
    order (1,1) double {mustBeInteger, mustBePositive}
    ports double {mustBeReal, mustBeFinite, mustBeNonempty}
    steps (1,1) double {mustBeInteger, mustBePositive}
    options.no_grads (1,1) logical = true
    options.rtol (1,1) double {mustBePositive, mustBeFinite} = 1.0e-12
end

if ndims(ports) ~= 2
    error("radia:ngsolve:HCurlPorts", ...
        "ports must be a two-dimensional ndof-by-nports real matrix.");
end

out = radia.internal.callMex( ...
    'hcurl.eddy_cln.native_basis', char(vol_path), order, ports, steps, ...
    options.no_grads, options.rtol);
end
