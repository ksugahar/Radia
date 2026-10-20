function result=checkHDivFusedMoments(meshPath,outputPath,profileAssembly)
% Kernel operator parity only; this does not replace an independent CPU solve.
arguments
    meshPath (1,:) char
    outputPath (1,:) char=''
    profileAssembly (1,1) logical=false
end
setup=radia_cuda_projection('assemble_hdiv_elements_sparse',meshPath,2,[1 .2 .1]);
config=radia.gpu.curvedHDivQuadrature(2);config.geometry_mode='curved';config.schema_version=1;config.host_vector_cache=true;config.profile_assembly=profileAssembly;
for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo'}
    config.(name{1})=setup.(name{1});
end
n=numel(config.host);config.charge_map=speye(n);
a=radia_cuda_hmatrix('create_hdiv_elements',config);ca=onCleanup(@()radia_cuda_hmatrix('destroy',a));
config.host_moment_cache=true;
b=radia_cuda_hmatrix('create_hdiv_elements',config);cb=onCleanup(@()radia_cuda_hmatrix('destroy',b));
A=zeros(n);B=zeros(n);
for column=1:n
    vector=zeros(n,1);vector(column)=1;
    radia_cuda_hmatrix('upload',a,vector);radia_cuda_hmatrix('apply',a,1);A(:,column)=radia_cuda_hmatrix('download',a);
    radia_cuda_hmatrix('upload',b,vector);radia_cuda_hmatrix('apply',b,1);B(:,column)=radia_cuda_hmatrix('download',b);
end
scalarProfile=radia_cuda_hmatrix('stats',a);fusedProfile=radia_cuda_hmatrix('stats',b);
if profileAssembly
    assert(isfield(scalarProfile,'assembly_profile_enabled') && scalarProfile.assembly_profile_enabled);
    assert(isfield(fusedProfile,'assembly_profile_enabled') && fusedProfile.assembly_profile_enabled);
    names=fieldnames(fusedProfile);duration=0;
    for i=1:numel(names)
        if startsWith(names{i},'profile_'),value=fusedProfile.(names{i});assert(isfinite(value) && value>=0);duration=duration+value;end
    end
    assert(duration>0,'No assembly event durations were recorded');
end
error=norm(A-B,'fro')/norm(A,'fro');assert(error<=1e-12,'Fused source moment operator mismatch');
bad=config;bad.host=[config.host;config.host(1:3)];bad.kind=[config.kind;config.kind(1:3)];bad.expo=[config.expo;config.expo(1:3,:)];bad.charge_map=speye(n+3);
rejected=false;
try
    unexpected=radia_cuda_hmatrix('create_hdiv_elements',bad);
catch exception
    rejected=contains(exception.message,'at most six');
end
assert(rejected,'Oversized source group must be rejected');assert(radia_cuda_hmatrix('live_count')==2);
if ~isempty(outputPath)
    radia.gpu.writeHDivH5(outputPath,struct('scalar_operator',A,'fused_operator',B),struct('units','SI','labels',{{'magnetic_body','exterior'}}));
end
result=struct('accepted',true,'scope','GPU operator parity only; independent CPU solver gate remains required', ...
    'charges',n,'relative_frobenius_error',error,'oversized_group_rejected',rejected,'operation_kind','operator','residual_applicability','not_applicable_operator');
clear ca cb
assert(radia_cuda_hmatrix('live_count')==0);result.live_count=0;result.scalar_profile=scalarProfile;result.fused_profile=fusedProfile;
end
