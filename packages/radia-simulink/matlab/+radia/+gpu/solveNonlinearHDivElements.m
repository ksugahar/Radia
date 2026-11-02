function result=solveNonlinearHDivElements(setup,appliedField,bhTable,options)
% One-shot wrapper includes preparation, all transfers, solve and destruction.
arguments
    setup (1,1) struct
    appliedField (1,3) double {mustBeReal,mustBeFinite}
    bhTable (:,2) double {mustBeReal,mustBeFinite}
    options.NonlinearTolerance (1,1) double {mustBePositive}=5e-14
    options.MaxIterations (1,1) double {mustBeInteger,mustBePositive}=50
    options.LinearTolerance (1,1) double {mustBePositive}=5e-14
    options.LinearMaxIterations (1,1) double {mustBeInteger,mustBePositive}=10000
    options.Preconditioner (1,:) char {mustBeMember(options.Preconditioner,{'reference-diagonal','reference-chebyshev','reference-chebyshev-graph'})}='reference-diagonal'
    options.ChebyshevDegree (1,1) double {mustBeInteger,mustBePositive}=32
    options.ChebyshevLowerFraction (1,1) double {mustBePositive,mustBeFinite}= .01
end
whole=tic;values=struct2cell(options);names=fieldnames(options);pairs=reshape([names';values'],1,[]);
obj=radia.gpu.PreparedHDivNonlinear(setup,appliedField,bhTable,pairs{:});
owner=onCleanup(@()delete(obj));info=obj.stats();result=obj.solve();
clear owner obj
wait(gpuDevice);result.solve_seconds=result.whole_seconds;result.prepare_seconds=info.prepare_seconds;result.prepare_phase_seconds=info.prepare_phase_seconds;result.whole_seconds=toc(whole);result.live_count=radia_hdiv_device('live_count');result.prepared_reuse=false;
result.scope='One-shot device nonlinear API including preparation and destruction';
end
