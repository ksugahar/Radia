function result=checkFinalResidual(mesh,cpu,large)
% Explicit fixed Chebyshev inverse action; independent CPU whole-solve gate.
[reference,~]=radia.gpu.readHDivH5(cpu);field=reference.applied_field(:)';
if large,chunk=128;cache=ceil(1.6*2^30);else,chunk=7;cache=32*2^20;end
budgetRejected=false;
try
    failed=radia.gpu.PreparedHDivNonlinear(mesh,2,field,reference.table,'MaterialBudgetBytes',1,'MaterialChunkElements',chunk,'CacheBytes',cache,'Preconditioner','mass-chebyshev');
    delete(failed);
catch ex
    assert(strcmp(ex.identifier,'radia:gpu:MaterialChunkBudget'),'Unexpected budget rejection');budgetRejected=true;
end
assert(budgetRejected && radia_cuda_hmatrix('live_count')==0);
opts={'NonlinearTolerance',1e-12,'LinearTolerance',1e-12,'LinearMaxIterations',6000,'MaxIterations',100,'MaterialBudgetBytes',2^30,'CacheBytes',cache,'MaterialChunkElements',chunk,'InexactNewton',true};
baseline=radia.gpu.solveNonlinearHDivElements(mesh,2,field,reference.table,opts{:},'Preconditioner','mass-chebyshev');
candidate=radia.gpu.solveNonlinearHDivElements(mesh,2,field,reference.table,opts{:},'Preconditioner','mass-chebyshev','FinalResidualOnly',true);
cpuError=norm(candidate.coefficients-reference.coefficients(:))/norm(reference.coefficients(:));
assert(cpuError<=1e-8 && candidate.nonlinear_relative_residual<=1e-12);
assert(candidate.residual_only_outer_evaluations==1 && candidate.outer_tangent_evaluations==candidate.nonlinear_iterations);
assert(baseline.outer_tangent_evaluations==baseline.nonlinear_iterations+1 && baseline.residual_only_outer_evaluations==0);
obj=radia.gpu.PreparedHDivNonlinear(mesh,2,field,reference.table,opts{:},'Preconditioner','mass-chebyshev','FinalResidualOnly',true);owner=onCleanup(@()delete(obj));
repeated=obj.solve();assert(norm(repeated.coefficients-candidate.coefficients)/norm(candidate.coefficients)<=1e-11);
rejected=false;try,obj.solve(1,ones(numel(reference.coefficients)+1,1));catch,rejected=true;end;assert(rejected);
recovered=obj.solve(1,repeated.coefficients);assert(recovered.nonlinear_iterations==0);
clear owner obj;wait(gpuDevice);assert(radia_cuda_hmatrix('live_count')==0);
result=struct('accepted',true,'baseline',baseline,'candidate',candidate,'repeated',repeated,'cpu_coefficient_relative_error',cpuError,'invalid_input_rejected',rejected,'early_material_budget_rejected',budgetRejected,'live_count',0,'scope','Full final tangent versus residual-only final recertification; same Chebyshev preconditioner; CPU agreement and lifecycle, single-call timings diagnostic only');
end
