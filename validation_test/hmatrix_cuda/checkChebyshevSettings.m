function result=checkChebyshevSettings(fixture,mesh,cpu)
inputs=load(fixture);settings=[8 .1;16 .1;16 .01;32 .01;32 .01];modes={'mass-chebyshev','mass-chebyshev','mass-chebyshev','mass-chebyshev','mass-chebyshev-graph'};linear={};nonlinear={};polynomialGraph={};
for k=1:numel(inputs.solver_cases)
    t=inputs.solver_cases(k);command='create';if isfield(t.config,'assembly'),command='create_analytic_constant_block';end
    h=radia_cuda_hmatrix(command,t.config);owner=onCleanup(@()radia_cuda_hmatrix('destroy',h));
    radia_cuda_hmatrix('upload',h,t.rhs);before=radia_cuda_hmatrix('download',h);
    for j=1:size(settings,1)
        degree=settings(j,1);lower=settings(j,2);
        s=radia_cuda_hmatrix('solve_linear_material',h,t.mass,t.rhs,t.inv_chi,1e-11,200,modes{j},degree,lower);
        e=norm(s.x-t.expected)/norm(t.expected);r=norm(t.operator*s.x-t.rhs)/norm(t.rhs);
        assert(e<=1e-10 && r<=1e-10 && s.true_relative_residual<=1e-11);
        assert(s.chebyshev_degree==degree && s.chebyshev_lower_fraction==lower);
        if j==5,assert(s.preconditioner_graph_launches>0);else,assert(s.preconditioner_graph_launches==0);end
        linear{end+1}=struct('case_index',k,'preconditioner',modes{j},'graph_launches',s.preconditioner_graph_launches,'degree',degree,'lower_fraction',lower,'cpu_error',e,'independent_true_residual',r,'native_true_residual',s.true_relative_residual,'iterations',s.iterations);
    end
    for invalid=[0 .1;65 .1;8 0;8 1;8 NaN]'
        rejected=false;try,unused=radia_cuda_hmatrix('solve_linear_material',h,t.mass,t.rhs,t.inv_chi,1e-11,200,'mass-chebyshev',invalid(1),invalid(2));catch,rejected=true;end %#ok<NASGU>
        assert(rejected && isequal(before,radia_cuda_hmatrix('download',h)));
    end
    q=radia_cuda_hmatrix('solve_linear_material',h,t.mass,t.rhs,t.inv_chi,1e-11,200,'mass-polynomial-graph');
    qe=norm(q.x-t.expected)/norm(t.expected);qr=norm(t.operator*q.x-t.rhs)/norm(t.rhs);
    assert(qe<=1e-10 && qr<=1e-10 && q.true_relative_residual<=1e-11 && q.preconditioner_graph_launches>0);
    polynomialGraph{end+1}=struct('cpu_error',qe,'independent_true_residual',qr,'native_true_residual',q.true_relative_residual,'graph_launches',q.preconditioner_graph_launches);
    clear owner
end
[ref,~]=radia.gpu.readHDivH5(cpu);
for j=1:size(settings,1)
    s=radia.gpu.solveNonlinearHDivElements(mesh,2,ref.applied_field(:)',ref.table,'NonlinearTolerance',1e-12,'LinearTolerance',1e-12,'LinearMaxIterations',6000,'MaxIterations',100,'MaterialBudgetBytes',2^30,'CacheBytes',32*2^20,'MaterialChunkElements',7,'InexactNewton',true,'FinalResidualOnly',true,'Preconditioner',modes{j},'ChebyshevDegree',settings(j,1),'ChebyshevLowerFraction',settings(j,2));
    e=norm(s.coefficients-ref.coefficients(:))/norm(ref.coefficients(:));assert(e<=1e-8 && s.nonlinear_relative_residual<=1e-12 && s.live_count==0);
    if j==5,assert(s.preconditioner_graph_launches>0);else,assert(s.preconditioner_graph_launches==0);end
    nonlinear{end+1}=struct('preconditioner',modes{j},'graph_launches',s.preconditioner_graph_launches,'degree',settings(j,1),'lower_fraction',settings(j,2),'cpu_error',e,'native_true_residual',s.nonlinear_relative_residual,'updates',s.nonlinear_iterations,'inner_iterations',sum([s.history.inner_iterations]),'whole_seconds_diagnostic',s.whole_seconds,'combined_array_bound_bytes',s.combined_array_bound_bytes);
end
assert(radia_cuda_hmatrix('live_count')==0 && radia_cuda_projection('live_count')==0);
result=struct('accepted',true,'linear',[linear{:}],'polynomial_graph',[polynomialGraph{:}],'nonlinear',[nonlinear{:}],'live_count',0,'scope','RTX bounded algebra CPU original residual and 288 DoF strong nonlinear CPU coefficients; timings diagnostic');
end
