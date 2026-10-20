function report=test_prepared_hdiv_host_setup(hostSetup,bhTable,outPath)
% Numeric host setup only; independent CPU certification belongs to validation.
solutions=cell(1,2);rows=cell(1,2);modes={'reference-diagonal','reference-chebyshev-graph'};
for k=1:2
    obj=radia.gpu.PreparedHDivNonlinear('',hostSetup.order,[1000 -2000 5000],bhTable, ...
        Assembly='ngsolve-host-normalized',HostSetup=hostSetup,Preconditioner=modes{k}, ...
        NonlinearTolerance=5e-14,LinearTolerance=5e-14);
    cleanup=onCleanup(@()delete(obj));before=obj.stats();first=obj.solve();
    second=obj.solve(1,first.coefficients);after=obj.stats();
    assert(strcmp(first.assembly_route,'ngsolve-host-normalized') && strcmp(first.mesh_sha256,hostSetup.mesh_sha256));
    assert(first.nonlinear_relative_residual<=1e-11 && second.nonlinear_relative_residual<=1e-11);
    assert(norm(first.coefficients-second.coefficients)<=1e-10*norm(first.coefficients));
    assert(before.material.setup_upload_bytes==after.material.setup_upload_bytes);
    solutions{k}=first.coefficients;
    rows{k}=struct('preconditioner',modes{k},'reported_residual',first.nonlinear_relative_residual, ...
        'setup_uploaded_once',true,'mean_magnetization',first.mean_magnetization);
    clear cleanup obj
    assert(radia_cuda_hmatrix('live_count')==0);
end
bad=hostSetup;bad.order=3;reject(bad,bhTable);
bad=hostSetup;bad.evaluation=bad.evaluation(:,1:end-1);reject(bad,bhTable);
report=struct('accepted',true,'scope','Shared MATLAB normalized host FE setup; CPU certificate external', ...
    'cases',[rows{:}],'final_live_count',radia_cuda_hmatrix('live_count'));
radia.gpu.writeHDivH5(string(outPath),struct('coefficients',[solutions{:}]), ...
    struct('units','SI','labels',{{'body','exterior'}},'report',report));
end
function reject(h,bh)
failed=false;try
    obj=radia.gpu.PreparedHDivNonlinear('',h.order,[1 0 0],bh,Assembly='ngsolve-host-normalized',HostSetup=h);delete(obj);
catch,failed=true;end
assert(failed && radia_cuda_hmatrix('live_count')==0);
end
