function delta = bdfDelta(zeta, method)
%BDFDELTA Native BDF generating function for convolution quadrature.
arguments
    zeta double {mustBeFinite}
    method (1,1) string {mustBeMember(method,["BDF1","BDF2"])} = "BDF2"
end
delta = radia.internal.callMex("acoustic.bdf_delta", zeta, char(method));
end
