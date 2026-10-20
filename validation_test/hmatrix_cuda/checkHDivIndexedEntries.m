function result=checkHDivIndexedEntries(meshPath,outputPath)
% Small full-operator comparison before connecting the oracle to compression.
arguments
    meshPath (1,:) char
    outputPath (1,:) char=''
end
setup=radia_cuda_projection('assemble_hdiv_elements_sparse',meshPath,2,[1 .2 .1]);
if setup.curved_cells>0
    config=radia.gpu.curvedHDivQuadrature(2);config.geometry_mode='curved';
else
    config=radia.gpu.affineHDivQuadrature(2);config.geometry_mode='affine';
end
config.schema_version=1;
for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo'}
    config.(name{1})=setup.(name{1});
end
n=numel(config.host);assert(n<=256,'Diagnostic full-operator budget exceeded');
config.charge_map=speye(n);
h=radia_cuda_hmatrix('create_hdiv_elements',config);closer=onCleanup(@()radia_cuda_hmatrix('destroy',h));
A=zeros(n);
for column=1:n
    vector=zeros(n,1);vector(column)=1;
    radia_cuda_hmatrix('upload',h,vector);radia_cuda_hmatrix('apply',h,1);A(:,column)=radia_cuda_hmatrix('download',h);
end
[rows,columns]=ndgrid(1:n,1:n);
values=radia_cuda_hmatrix('curved_entry_batch',config,rows(:),columns(:));
B=reshape(values,n,n);error=norm(A-B,'fro')/norm(A,'fro');assert(error<=1e-12,'Indexed curved entry mismatch');
assert(isempty(radia_cuda_hmatrix('curved_entry_batch',config,zeros(0,1),zeros(0,1))));
for invalid={0,n+1,1.5,NaN}
    rejected=false;
    try
        unexpected=radia_cuda_hmatrix('curved_entry_batch',config,invalid{1},1); %#ok<NASGU>
    catch
        rejected=true;
    end
    assert(rejected,'Invalid entry index was accepted');
end
% Repeat a diagonal request across the 256-entry workspace boundary.
repeated=radia_cuda_hmatrix('curved_entry_batch',config,ones(257,1),ones(257,1));
assert(all(repeated==B(1,1)),'Entry batch boundary/repeat mismatch');
clear closer
assert(radia_cuda_hmatrix('live_count')==0);
result=struct('accepted',true,'charges',n,'relative_frobenius_error',error,'live_count',0, ...
    'scope','Indexed GPU entries versus dense GPU operator only; no CPU physics or HACApK compression acceptance', ...
    'operation_kind','operator','true_relative_residual',[],'residual_applicability','not_applicable_operator');
if ~isempty(outputPath)
    radia.gpu.writeHDivH5(outputPath,struct('dense_operator',A,'indexed_operator',B),struct('units','SI','labels',{{'magnetic_body','exterior'}}));
end
end
