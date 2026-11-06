function result=solveHDivElements(setup,appliedField,invChi,options)
% Device HDiv linear isotropic material solve over a radia.hdiv-device-setup.v1
% (affine TET RT0, BDM1 or BDM2; the setup selects the order and Gram compression).
% Solves (B'*G*B + invChi*M) x = f(H0) with the device-built Gram, charge map
% and reference mass; only vectors cross. Includes creation and destruction.
arguments
    setup (1,1) struct
    appliedField (1,3) double {mustBeReal,mustBeFinite}
    invChi (1,1) double {mustBePositive,mustBeFinite}
    options.Tolerance (1,1) double {mustBePositive,mustBeFinite}=1e-10
    options.MaxIterations (1,1) double {mustBeInteger,mustBePositive}=1000
    options.Preconditioner (1,:) char {mustBeMember(options.Preconditioner,{'mass-diagonal','mass-cholesky','mass-polynomial','mass-polynomial-graph','mass-chebyshev','mass-chebyshev-graph'})}='mass-diagonal'
    options.ChebyshevDegree double {mustBeScalarOrEmpty,mustBeInteger,mustBePositive}=[]
    options.ChebyshevLowerFraction double {mustBeScalarOrEmpty,mustBePositive,mustBeFinite}=[]
end
assert(options.Tolerance<=1e-6,'radia:gpu:HDivLinear','True residual tolerance must be at most 1e-6');
setup=radia.gpu.readHDivDeviceSetup(setup);
timer=tic;device=gpuDevice;
stage=tic;mass=radia.gpu.HDivDeviceMass(setup);massOwner=onCleanup(@()delete(mass));
op=radia.gpu.HDivDeviceOperator(setup);opOwner=onCleanup(@()delete(op));wait(device);createSeconds=toc(stage);
stage=tic;solution=op.solveLinear(mass,setup.rhs_operator*appliedField(:),invChi,options.Tolerance, ...
    options.MaxIterations,options.Preconditioner,options.ChebyshevDegree,options.ChebyshevLowerFraction);
solveSeconds=toc(stage);assembly=op.assembly();stats=op.stats();
result=struct('coefficients',solution.x,'mean_magnetization',setup.mean_operator*solution.x, ...
    'true_relative_residual',solution.relative_residual,'iterations',solution.iterations, ...
    'preconditioner',solution.preconditioner,'preconditioner_graph_launches',solution.preconditioner_graph_launches, ...
    'operator_owned_device_bytes',stats.owned_device_bytes,'solver_owned_device_bytes',solution.solver_owned_device_bytes, ...
    'setup_peak_device_bytes',assembly.setup_peak_device_bytes,'assembly',assembly, ...
    'order',setup.order,'dofs',setup.dofs,'mesh_sha256',setup.mesh_sha256,'backend','gpu', ...
    'gram_compression',assembly.compression,'family',setup.family,'scope','Device HDiv orders 0/1/2, linear isotropic material');
stage=tic;clear opOwner massOwner
result.destroy_seconds=toc(stage);result.create_seconds=createSeconds;result.solve_transfer_seconds=solveSeconds;
result.whole_seconds=toc(timer);
result.timing_limitations='Stage wall times include native transfers; transfer-only and kernel-only times are not separated.';
end
