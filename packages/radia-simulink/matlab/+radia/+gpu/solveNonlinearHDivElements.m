function result=solveNonlinearHDivElements(meshPath,order,appliedField,bhTable,options)
% One-shot wrapper includes preparation, all transfers, solve and destruction.
arguments
    meshPath (1,:) char
    order (1,1) double {mustBeMember(order,[1 2])}
    appliedField (1,3) double {mustBeReal,mustBeFinite}
    bhTable (:,2) double {mustBeReal,mustBeFinite}
    options.NonlinearTolerance (1,1) double {mustBePositive}=1e-6
    options.MaxIterations (1,1) double {mustBeInteger,mustBePositive}=50
    options.LinearTolerance (1,1) double {mustBePositive}=1e-10
    options.LinearMaxIterations (1,1) double {mustBeInteger,mustBePositive}=1000
    options.MaterialBudgetBytes (1,1) double {mustBePositive,mustBeInteger}=512*2^20
    options.CacheBytes (1,1) double {mustBePositive,mustBeInteger}=32*2^20
    options.MaterialChunkElements (1,1) double {mustBeInteger,mustBeNonnegative}=0
    options.InexactNewton (1,1) logical=false
    options.FinalResidualOnly (1,1) logical=false
    options.Preconditioner (1,:) char {mustBeMember(options.Preconditioner,{'mass-diagonal','mass-cholesky','mass-polynomial','mass-polynomial-graph','mass-chebyshev','mass-chebyshev-graph'})}='mass-cholesky'
    options.ChebyshevDegree (1,1) double {mustBeInteger,mustBePositive}=8
    options.ChebyshevLowerFraction (1,1) double {mustBePositive,mustBeFinite}=.1
end
whole=tic;values=struct2cell(options);names=fieldnames(options);pairs=reshape([names';values'],1,[]);
obj=radia.gpu.PreparedHDivNonlinear(meshPath,order,appliedField,bhTable,pairs{:});
owner=onCleanup(@()delete(obj));info=obj.stats();result=obj.solve();
clear owner obj
wait(gpuDevice);result.solve_seconds=result.whole_seconds;result.prepare_seconds=info.prepare_seconds;result.prepare_phase_seconds=info.prepare_phase_seconds;result.whole_seconds=toc(whole);result.live_count=radia_cuda_hmatrix('live_count');result.prepared_reuse=false;
result.scope='One-shot bounded nonlinear API including preparation and destruction';
end
