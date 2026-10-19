function result=checkProjectionOwnerPhases(mesh)
% Requires the isolated RADIA_PROJECTION_LIFECYCLE_TEST build.
% Visual phases are negative controls: they deliberately register mesh owners.
phases={'mesh-only','targets-cf','fespace-construct','fespace-update', ...
    'gridfunction-create','gridfunction-update','interpolate', ...
    'gridfunction-novisual-create','gridfunction-novisual-update','interpolate-novisual'};
expectedMesh=[0 0 0 0 1 1 1 0 0 0];
counts=cell(size(phases));
for p=1:numel(phases)
    rows=zeros(10,8);
    for k=1:10
        rows(k,:)=radia_cuda_projection('diagnostic_owner_phase',mesh,phases{p});
        assert(rows(k,5)==expectedMesh(p) && all(rows(k,6:8)==0), ...
            'Projection ownership phase changed: %s',phases{p});
        assert(radia_cuda_projection('live_count')==0);
    end
    counts{p}=rows;
end
result=struct('accepted',true,'phases',{phases},'counts',{counts}, ...
    'expected_released_mesh_count',expectedMesh, ...
    'scope','Construction ownership regression with deliberate visualization negative controls');
end
