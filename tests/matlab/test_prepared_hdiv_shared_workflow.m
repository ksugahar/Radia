function report=test_prepared_hdiv_shared_workflow(meshPath,order,bhTable,outPath)
% Small checked-FE workflow test; independent CPU certification is external.
rows=cell(1,2);solutions=cell(1,2);
for k=1:2
    modes={'reference-diagonal','reference-chebyshev-graph'};
    obj=radia.gpu.PreparedHDivNonlinear(meshPath,order,[1000 -2000 5000],bhTable, ...
        Preconditioner=modes{k},NonlinearTolerance=5e-14,LinearTolerance=5e-14);
    cleanup=onCleanup(@()delete(obj));
    before=obj.stats();first=obj.solve();second=obj.solve(1,first.coefficients);after=obj.stats();
    assert(strcmp(first.implementation,'shared-native') && first.nonlinear_relative_residual<=1e-11);
    assert(norm(first.coefficients-second.coefficients)<=1e-10*norm(first.coefficients));
    assert(before.material.setup_upload_bytes==after.material.setup_upload_bytes);
    assert(before.material.epoch==after.material.epoch);
    rows{k}=struct('preconditioner',modes{k},'reported_residual',first.nonlinear_relative_residual, ...
        'warm_reuse_residual',second.nonlinear_relative_residual,'setup_uploaded_once',true);
    solutions{k}=first.coefficients;
    clear cleanup obj
    assert(radia_cuda_hmatrix('live_count')==0);
end
failed=false;
try
    bad=radia.gpu.PreparedHDivNonlinear(meshPath,order,[1000 -2000 5000],bhTable,MaterialChunkElements=1);
    delete(bad);
catch
    failed=true;
end
assert(failed && radia_cuda_hmatrix('live_count')==0);
report=struct('accepted',true,'scope','MATLAB shared native prepared FE workflow; independent certification external', ...
    'cases',[rows{:}],'unsupported_chunk_rejected',true,'final_live_count',0);
radia.gpu.writeHDivH5(string(outPath),struct('coefficients',[solutions{:}]), ...
    struct('units','SI','labels',{{'body','exterior'}},'report',report));
end
