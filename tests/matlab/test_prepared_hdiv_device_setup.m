function report=test_prepared_hdiv_device_setup(setup,bhTable,outPath)
% Guards: the MATLAB device workflow (radia_hdiv_device) over a checked
% radia.hdiv-device-setup.v1 solves, reuses prepared state, uploads the
% setup once, rejects malformed setups before allocation and releases every
% handle; the frozen tangent equals the evaluated tangent action and refuses a
% stale epoch. Assumes: independent CPU certification is external (validation_test).
modes={'reference-diagonal','reference-chebyshev-graph'};rows=cell(1,2);solutions=cell(1,2);
for k=1:2
    % A Newton step whose data is near rounding reaches a true residual of about
    % 1e-12 at best on this fixture, so the inner tolerance stays well above that.
    obj=radia.gpu.PreparedHDivNonlinear(setup,[1000 -2000 5000],bhTable,Preconditioner=modes{k}, ...
        NonlinearTolerance=5e-14,LinearTolerance=1e-8);
    cleanup=onCleanup(@()delete(obj));
    before=obj.stats();first=obj.solve();second=obj.solve(1,first.coefficients);after=obj.stats();
    assert(strcmp(first.mesh_sha256,setup.mesh_sha256) && strcmp(first.backend,'gpu'));
    assert(first.nonlinear_relative_residual<=1e-11 && second.nonlinear_relative_residual<=1e-11, ...
        'residuals %.3e %.3e',first.nonlinear_relative_residual,second.nonlinear_relative_residual);
    assert(norm(first.coefficients-second.coefficients)<=1e-10*norm(first.coefficients));
    assert(before.material.setup_upload_bytes==after.material.setup_upload_bytes);
    rows{k}=struct('preconditioner',modes{k},'reported_residual',first.nonlinear_relative_residual, ...
        'warm_reuse_residual',second.nonlinear_relative_residual,'setup_uploaded_once',true, ...
        'mean_magnetization',first.mean_magnetization);
    solutions{k}=first.coefficients;
    clear cleanup obj
    assert(radia_hdiv_device('live_count')==0);
end
linear=radia.gpu.solveHDivElements(setup,[0 0 1e5],.1,Tolerance=1e-11,MaxIterations=2000,Preconditioner='mass-cholesky');
assert(linear.true_relative_residual<=1e-11 && radia_hdiv_device('live_count')==0);
law=radia.gpu.makeHDivBHMaterial(bhTable);material=radia.gpu.HDivBHMaterial(setup,law);
state=1e3*sin((1:setup.dofs)');direction=cos((1:setup.dofs)');
evaluated=material.evaluate(state,direction);frozen=material.linearize(state);
action=material.applyTangent(direction,frozen.linearization_epoch);
assert(norm(action-evaluated.tangent_action)<=1e-11*norm(evaluated.tangent_action), ...
    'tangent %.3e',norm(action-evaluated.tangent_action)/norm(evaluated.tangent_action));
material.evaluate(state+direction);failed=false;
try
    material.applyTangent(direction,frozen.linearization_epoch);
catch
    failed=true;
end
assert(failed,'A stale linearization epoch must be refused');
samples=radia.gpu.HDivBHMaterial.inverseSamples(law,[1;1e3;1e5]);
assert(numel(samples.secant)==3 && all(samples.secant>0) && all(samples.energy>=0));
clear material
assert(radia_hdiv_device('live_count')==0);
bad=setup;bad.order=3;reject(bad,bhTable);
bad=setup;bad.evaluation(1).element_dofs=bad.evaluation(1).element_dofs(1:end-1);reject(bad,bhTable);
bad=setup;bad.charges.charge_host=double(bad.charges.charge_host)+.5;reject(bad,bhTable);
report=struct('accepted',true,'scope','MATLAB device prepared HDiv workflow; independent certification external', ...
    'cases',[rows{:}],'linear_relative_residual',linear.true_relative_residual, ...
    'final_live_count',radia_hdiv_device('live_count'));
radia.gpu.writeHDivH5(string(outPath),struct('coefficients',[solutions{:}],'linear_coefficients',linear.coefficients), ...
    struct('units','SI','labels',{{'body','exterior'}},'report',report));
end
function reject(setup,bh)
failed=false;
try
    obj=radia.gpu.PreparedHDivNonlinear(setup,[1 0 0],bh);delete(obj);
catch
    failed=true;
end
assert(failed && radia_hdiv_device('live_count')==0,'Malformed device setup must fail before allocation');
end
