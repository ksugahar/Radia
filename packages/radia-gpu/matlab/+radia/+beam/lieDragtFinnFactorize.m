function result = lieDragtFinnFactorize(R,T,U,V)
%LIEDRAGTFINNFACTORIZE Factor a factorial R/T/U/V map into Lie generators.
%   Native value path of the fourth-order Dragt-Finn factorization
%   (rad_lie_map_kernel.cpp): f3/f4/f5 generator extraction with
%   symmetrization, symplectic reconstruction of T/U/V, and the formal
%   symplectic-residual audit.
%
%   result fields: f3/f4/f5 generators, reconstructed T/U/V, the
%   generator symmetry defects, relative_reconstruction_error, and the
%   raw/reconstructed residual coefficient rows
%   [constant linear quadratic cubic].
arguments
    R (6,6) double {mustBeReal,mustBeFinite}
    T (6,6,6) double {mustBeReal,mustBeFinite}
    U (6,6,6,6) double {mustBeReal,mustBeFinite}
    V (6,6,6,6,6) double {mustBeReal,mustBeFinite}
end
config = struct();
config.R = double(R);
config.T = double(T);
config.U = double(U);
config.V = double(V);
result = radia.internal.callMex('beam.lie.dragt_finn_factorize',config);
end
