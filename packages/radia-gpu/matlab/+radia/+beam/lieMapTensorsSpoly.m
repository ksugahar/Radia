function result = lieMapTensorsSpoly(Ay,As,lengths,curvature,rigidity,options)
%LIEMAPTENSORSSPOLY Integrate the fourth-order factorial map tensors.
%   Native value path of the CanonicalHCurl Lie-map construction
%   (rad_lie_map_kernel.cpp): per-segment degree-5 Hamiltonian jets from
%   the (Ay, As) zeta-polynomial coefficient arrays, the nonautonomous
%   stage RK4 flow, and the sequential segment composition.
%
%   Ay, As    (d+1, d+1, k_s+1, n) coefficient arrays; entry
%             (x+1, y+1, k+1, s) multiplies zeta^k x^x y^y on segment s.
%   lengths   (n) segment lengths in metres.
%   curvature (n) constant curvature per segment, or (k_c+1, n)
%             zeta-polynomial columns.
%   rigidity  magnetic rigidity B*rho in T*m.
%
%   result fields: R (6x6), T/U/V (rank 3..5 tensors), hamiltonian_linear
%   (6 x n worst-stage H1 per segment), worst_hamiltonian_linear.
arguments
    Ay double {mustBeReal,mustBeFinite}
    As double {mustBeReal,mustBeFinite}
    lengths (1,:) double {mustBeReal,mustBeFinite,mustBePositive}
    curvature double {mustBeReal,mustBeFinite}
    rigidity (1,1) double {mustBeReal,mustBeFinite}
    options.CurvatureSign (1,1) double = 1.0
    options.ReferenceBeta (1,1) double = 1.0
    options.LongitudinalCovariant (1,1) logical = true
    options.MaximumStep (1,1) double {mustBePositive} = 1.0e-3
    options.MaximumSteps (1,1) double {mustBeInteger,mustBePositive} = 1e6
    options.ReferenceOrbitTolerance (1,1) double {mustBePositive} = 1.0e-8
end
config = struct();
config.schema = 'radia.beam.lie-map-spoly.v1';
config.Ay_t_m = double(Ay);
config.As_t_m = double(As);
config.lengths_m = double(lengths);
config.curvature_per_m = double(curvature);
config.magnetic_rigidity_t_m = double(rigidity);
config.curvature_sign = double(options.CurvatureSign);
config.reference_beta = double(options.ReferenceBeta);
config.longitudinal_covariant = logical(options.LongitudinalCovariant);
config.maximum_step_m = double(options.MaximumStep);
config.maximum_steps = double(options.MaximumSteps);
config.reference_orbit_tolerance = double(options.ReferenceOrbitTolerance);
result = radia.internal.callMex('beam.lie.map_tensors_spoly',config);
end
