function result=checkChebyshevLinear(meshPath,cacheBytes)
% Same frozen G, identical mass/RHS: compare explicit fixed inverse actions.
setup=radia_cuda_projection('assemble_hdiv_elements_sparse',meshPath,2,[1 .2 .1]);
if setup.curved_cells>0,cfg=radia.gpu.curvedHDivQuadrature(2);cfg.geometry_mode='curved';else,cfg=radia.gpu.affineHDivQuadrature(2);cfg.geometry_mode='affine';end
cfg.schema_version=1;
for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo','charge_map'},cfg.(name{1})=setup.(name{1});end
config=struct('assembly','hdiv_hacapk_device_dense','operator',cfg,'epsilon',1e-10,'leaf_size',64,'eta',2.,'max_rank',64,'cache_bytes',cacheBytes);
op=radia.gpu.HMatrixOperator(config);owner=onCleanup(@()delete(op));rows={};
for scale=[1 .5 -1 0]
    rhs=scale*setup.rhs;
    a=op.solveMaterial(setup.mass,rhs,2.5,1e-11,4000,'mass-polynomial');
    b=op.solveMaterial(setup.mass,rhs,2.5,1e-11,4000,'mass-chebyshev');
    err=norm(a.x-b.x)/max(norm(a.x),realmin);assert(err<=1e-9);
    op.upload(b.x);op.apply();r=norm(op.download()+2.5*setup.mass*b.x-rhs)/max(norm(rhs),realmin);
    assert(r<=1e-11 && a.preconditioner_graph_launches==0);
    if scale~=0,assert(b.preconditioner_graph_launches==0);end
    rows{end+1}=struct('scale',scale,'relative_error',err,'true_residual',r,'original_iterations',a.iterations,'chebyshev_iterations',b.iterations,'unused_graph_launches',b.preconditioner_graph_launches);
end
rejected=false;try,op.solveMaterial(-setup.mass,setup.rhs,2.5,1e-11,4000,'mass-chebyshev');catch,rejected=true;end;assert(rejected);
recovered=op.solveMaterial(setup.mass,setup.rhs,2.5,1e-11,4000,'mass-chebyshev');assert(recovered.true_relative_residual<=1e-11);
clear owner op;assert(radia_cuda_hmatrix('live_count')==0);
result=struct('accepted',true,'rows',[rows{:}],'invalid_mass_rejected',rejected,'live_count',0,'scope','Chebyshev versus Jacobi polynomial on identical G/mass/RHS; independent CPU acceptance is a separate whole-solve gate');
end
