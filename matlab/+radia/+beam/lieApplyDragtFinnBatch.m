function states = lieApplyDragtFinnBatch(factorization,states,options)
%LIEAPPLYDRAGTFINNBATCH Advance a (6,N) ensemble through Dragt-Finn factors.
%   Native tracking hot path (rad_lie_map_batch.cpp): the ensemble flows
%   through exp(:f5:) -> exp(:f4:) -> exp(:f3:) -> R with particle-parallel
%   implicit-midpoint generator flows that mirror the single-state
%   reference to roundoff (~0.7 us/particle measured).
%
%   factorization is the struct returned by lieDragtFinnFactorize (fields
%   R, f3, f4, and optionally f5); states is (6, n_states).
arguments
    factorization (1,1) struct
    states (6,:) double {mustBeReal,mustBeFinite}
    options.GeneratorSubsteps (1,1) double ...
        {mustBeInteger,mustBePositive} = 1
    options.NewtonTolerance (1,1) double {mustBePositive} = 1.0e-13
    options.MaximumNewtonIterations (1,1) double ...
        {mustBeInteger,mustBePositive} = 20
end
config = struct();
config.R = double(factorization.R);
config.f3 = double(factorization.f3);
config.f4 = double(factorization.f4);
if isfield(factorization,'f5')
    config.f5 = double(factorization.f5);
end
config.states = double(states);
config.generator_substeps = double(options.GeneratorSubsteps);
config.newton_tolerance = double(options.NewtonTolerance);
config.maximum_newton_iterations = ...
    double(options.MaximumNewtonIterations);
states = radia.internal.callMex('beam.lie.apply_dragt_finn_batch',config);
end
