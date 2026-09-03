function delta = bdfDelta(zeta, method)
%BDFDELTA BDF generating function for convolution quadrature.
arguments
    zeta double {mustBeFinite}
    method (1,1) string {mustBeMember(method,["BDF1","BDF2"])} = "BDF2"
end
if method == "BDF1"
    delta = 1 - zeta;
else
    delta = 1.5 - 2 * zeta + 0.5 * zeta.^2;
end
end
