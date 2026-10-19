function result=checkNonlinearChunks(meshPath,cpuPath,chunkElements)
% Native point identity and GPU chunk accumulation versus independent CPU material.
validateattributes(chunkElements,{'double'},{'scalar','integer','positive','<=',128});
[cpu,~]=radia.gpu.readHDivH5(cpuPath);setup=radia_cuda_projection('assemble_hdiv_elements_sparse',meshPath,2,[0 0 100]);ne=size(setup.cell_nodes,1);
chunks={};seen=false(numel(cpu.weights),1);errors=struct('weights',0.,'points',0.,'values',0.);entries=0;tangentBound=0;
for first=0:chunkElements:ne-1
    q=radia_cuda_projection('hdiv_material_quadrature_chunk',meshPath,2,first,min(chunkElements,ne-first));
    ids=q.global_point_indices(:)+1;nq=numel(ids);assert(q.total_points==numel(seen) && q.first_element==first && all(ids>=1 & ids<=numel(seen)) && numel(unique(ids))==nq && ~any(seen(ids)));seen(ids)=true;
    assert(nq<=65536 && q.potential_shape_entries<=4000000);
    values=reshape(q.evaluation*cpu.m(:),nq,3);
    errors.weights=max(errors.weights,norm(q.weights-cpu.weights(ids))/norm(cpu.weights(ids)));
    errors.points=max(errors.points,norm(q.points-cpu.points(ids,:),'fro')/max(norm(cpu.points(ids,:),'fro'),realmin));
    errors.values=max(errors.values,norm(values-cpu.values(ids,:),'fro')/max(norm(cpu.values(ids,:),'fro'),realmin));
    item=struct('E',{{}},'weights',gpuArray(q.weights));for c=1:3,item.E{c}=gpuArray(q.evaluation((c-1)*nq+(1:nq),:));end
    chunks{end+1}=item;entries=entries+q.potential_shape_entries;tangentBound=tangentBound+q.tangent_nonzero_bound;
end
assert(all(seen) && max([errors.weights errors.points errors.values])<=1e-11);
law=radia.gpu.makeHDivBHMaterial(cpu.table);law.device=struct('knots',gpuArray(law.knots),'coefficients',gpuArray(law.coefficients),'integral_base',gpuArray(law.integral_base));
out=radia.gpu.evaluateHDivBHChunks(chunks,law,cpu.m(:),true);
errors.energy=abs(out.energy-cpu.energy)/max(abs(cpu.energy),realmin);errors.gradient=norm(out.gradient-cpu.gradient(:))/norm(cpu.gradient(:));errors.tangent=norm(out.tangent-cpu.tangent,'fro')/norm(cpu.tangent,'fro');assert(max([errors.energy errors.gradient errors.tangent])<=1e-8);
invalid=0;for spec={[0 0],[ne 1],[-1 1],[0 ne+1],[0 .5]}
    s=spec{1};try,ignored=radia_cuda_projection('hdiv_material_quadrature_chunk',meshPath,2,s(1),s(2));clear ignored;catch,invalid=invalid+1;end
end
assert(invalid==5);result=struct('accepted',true,'chunks',numel(chunks),'chunk_elements',chunkElements,'dofs',numel(cpu.m),'points',numel(seen),'shape_entry_bound',entries,'tangent_nonzero_bound',tangentBound,'errors',errors,'invalid_ranges_rejected',invalid,'scope','Material chunks only, no complete nonlinear solve or performance acceptance');
end
