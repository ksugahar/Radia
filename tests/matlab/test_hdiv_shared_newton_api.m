function report=test_hdiv_shared_newton_api()
% Small shared-kernel API/lifecycle regression, not FE physics certification.
G=[4 .2 .1;.2 3 -.1;.1 -.1 5];
leaf=struct('row',1,'column',1,'rows',3,'columns',3,'rank',0, ...
    'a1',reshape(G.',[],1),'a2',[]);
config=struct('schema_version',1,'dimension',3,'permutation',(1:3)', ...
    'scale',ones(3,1),'leaves',leaf,'operator_kind','demag','charge_map',speye(3));
op=radia.gpu.HMatrixOperator(config);
opCleanup=onCleanup(@()delete(op)); %#ok<NASGU>
law=struct('knots',[0;1;2],'coefficients',[1 0 1 0;1 3 4 2], ...
    'integral_base',[0;.75;6],'cap',1.5,'barrier_slope',7, ...
    'last_field',4.875,'initial_reluctivity',1);
material=radia.gpu.HDivBHMaterial(speye(3),1,law);
materialCleanup=onCleanup(@()delete(material)); %#ok<NASGU>
target=[.2;-.3;.1];rhs=G*target+(1+target'*target)*target;
options=struct('nonlinear_tolerance',1e-11,'linear_tolerance',5e-14, ...
    'maximum_iterations',30,'maximum_linear_iterations',1000, ...
    'preconditioner','reference-diagonal','degree',32,'lower_fraction',.01);
cases=cell(1,2);before=material.stats();
for k=1:2
    names={'reference-diagonal','reference-chebyshev-graph'};
    options.preconditioner=names{k};
    result=op.solveNonlinearMaterial(material,speye(3),rhs,zeros(3,1),options);
    residual=G*result.x+(1+result.x'*result.x)*result.x-rhs;
    coefficientError=norm(result.x-target)/norm(target);
    originalResidual=norm(residual)/norm(rhs);
    assert(coefficientError<=1e-10 && originalResidual<=1e-11);
    cases{k}=struct('preconditioner',names{k},'coefficient_error',coefficientError, ...
        'independent_original_residual',originalResidual);
end
after=material.stats();assert(before.setup_upload_bytes==after.setup_upload_bytes);
options.preconditioner='mass-cholesky';
failed=false;
try,op.solveNonlinearMaterial(material,speye(3),rhs,zeros(3,1),options);catch,failed=true;end
assert(failed,'Unsupported preconditioner must not silently change route');
options.preconditioner='reference-diagonal';
recovered=op.solveNonlinearMaterial(material,speye(3),rhs,zeros(3,1),options);
assert(norm(recovered.x-target)<=1e-10*norm(target));
delete(material);delete(op);
assert(radia_cuda_hmatrix('live_count')==0);
report=struct('accepted',true,'scope','small MATLAB shared Newton API and lifecycle', ...
    'cases',[cases{:}],'unsupported_preconditioner_rejected',true, ...
    'failure_recovery',true,'setup_uploaded_once',true);
end
