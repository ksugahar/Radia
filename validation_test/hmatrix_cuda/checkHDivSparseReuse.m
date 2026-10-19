function result=checkHDivSparseReuse(meshPath,referencePreconditioner,cacheBytes)
if nargin<2,referencePreconditioner='mass-cholesky';end
if nargin<3,cacheBytes=32*2^20;end
assert(ismember(referencePreconditioner,{'mass-diagonal','mass-cholesky'}));
% Bounded resident reuse/SPD action gate; independent CPU acceptance separate.
assert(radia_cuda_hmatrix('live_count')==0);
setup=radia_cuda_projection('assemble_hdiv_elements_sparse',meshPath,2,[1 .2 .1]);
if setup.curved_cells>0,cfg=radia.gpu.curvedHDivQuadrature(2);cfg.geometry_mode='curved';
else,cfg=radia.gpu.affineHDivQuadrature(2);cfg.geometry_mode='affine';end
cfg.schema_version=1;
for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo','charge_map'},cfg.(name{1})=setup.(name{1});end
config=struct('assembly','hdiv_hacapk_device_dense','operator',cfg,'epsilon',1e-10,'leaf_size',64,'eta',2.,'max_rank',64,'cache_bytes',cacheBytes);
t=tic;op=radia.gpu.HMatrixOperator(config);creation=toc(t);owner=onCleanup(@()delete(op));
baseline=op.solveMaterial(setup.mass,setup.rhs,2.5,1e-11,1000,referencePreconditioner);
scales=[1 .5 -1 0];errors=zeros(1,4);residuals=errors;times=errors;iterations=errors;bytes=errors;
for k=1:4
    rhs=setup.rhs*scales(k);t=tic;s=op.solveMaterial(setup.mass,rhs,2.5,1e-11,1000,'mass-polynomial');times(k)=toc(t);
    errors(k)=norm(s.x-baseline.x*scales(k))/max(norm(baseline.x*scales(k)),realmin);
    op.upload(s.x);op.apply();ax=op.download()+2.5*setup.mass*s.x;
    residuals(k)=norm(rhs-ax)/max(norm(rhs),realmin);iterations(k)=s.iterations;bytes(k)=s.solver_owned_device_bytes;
end
assert(all(errors<1e-8) && all(residuals<=1e-11));
assert(all(bytes<=baseline.solver_owned_device_bytes));
if strcmp(referencePreconditioner,'mass-cholesky'),assert(all(bytes<baseline.solver_owned_device_bytes));end
rejected=false;try,op.solveMaterial(-setup.mass,setup.rhs,2.5,1e-11,1000,'mass-polynomial');catch,rejected=true;end
assert(rejected && radia_cuda_hmatrix('live_count')==1);
s=op.solveMaterial(setup.mass,setup.rhs,2.5,1e-11,1000,'mass-polynomial');
assert(norm(s.x-baseline.x)/max(norm(baseline.x),realmin)<1e-8);
clear owner
assert(radia_cuda_hmatrix('live_count')==0);
result=struct('accepted',true,'dofs',numel(setup.rhs),'relative_errors',errors,'true_residuals',residuals,'iterations',iterations, ...
    'solver_owned_bytes',bytes,'reference_solver_owned_bytes',baseline.solver_owned_device_bytes,'reference_preconditioner',referencePreconditioner, ...
    'creation_seconds',creation,'reused_solve_seconds',times,'invalid_mass_rejected',rejected,'live_count',0, ...
    'scope','One frozen geometry/operator, four explicit RHS scalings; no solver-workspace reuse or independent CPU certification; solve timings exclude FE/setup/destruction and are diagnostic only');
end
