function [action,description]=makeGatewayMemoryWorkload(kind,inputsDir,outputDir)
% Caller retains fixtures; each step releases outputs and native handles.
% The common driver initializes gpuDevice and waits BEFORE calling this.
description=struct('kind',kind,'persistent_owner','Caller-owned frozen fixture inputs', ...
    'temporary_owner','Step-local MATLAB outputs and native contexts');
switch kind
    case 'csr-pcg'
        data=load(fullfile(inputsDir,'ic-reference.mat'));f=struct('op',csr(data.sparse_A), ...
            'rhs',gpuArray(data.sparse_b),'inverse',gpuArray(1./diag(data.sparse_A)), ...
            'matrix',data.sparse_A,'expected',data.sparse_expected);
        action=@()csrStep(f);
    case 'ic'
        data=load(fullfile(inputsDir,'ic-reference.mat'));f=struct('op',csr(data.sparse_A), ...
            'rhs',gpuArray(data.sparse_b),'matrix',data.sparse_A,'expected',data.sparse_expected);
        action=@()icStep(f);
    case 'ic-harness'
        % Exact original validation path, including file IO and rejection tests.
        action=@()validate_ic(fullfile(inputsDir,'ic-reference.mat'), ...
            fileparts(which('radia_cuda_ic')),fileparts(fileparts(fileparts(which('radia.gpu.ICContext')))), ...
            fullfile(outputDir,'ic-harness-result.mat'));
        description.temporary_owner='Original validation harness, exceptions and MAT-file IO; not a Simulink step';
    case 'field'
        f=load(fullfile(inputsDir,'field-reference.mat'));action=@()fieldStep(f);
    case 'fmm'
        f=load(fullfile(inputsDir,'fmm-reference.mat'));f.charge=gpuArray(f.charge);f.dipole=gpuArray(f.dipole);
        action=@()fmmStep(f);
    case 'projection'
        f=jsondecode(fileread(fullfile(inputsDir,'projection-reference.json')));
        f.mesh=fullfile(inputsDir,'curved-sphere.vol');action=@()projectionStep(f);
    case 'ams'
        f=jsondecode(fileread(fullfile(inputsDir,'ams-independent-cpu.json')));m=f.model;
        m.A.rowptr=int32(m.A.rowptr);m.A.columns=int32(m.A.columns);
        m.G.rowptr=int32(m.G.rowptr);m.G.columns=int32(m.G.columns);m.free=int32(m.free);
        f.hierarchy=radia.gpu.prepareAMS(m,SetupThreads=1);f.gpu_rhs=gpuArray(f.rhs);
        row=repelem((1:m.A.rows)',diff(double(m.A.rowptr)));
        f.matrix=sparse(row,double(m.A.columns)+1,m.A.values,m.A.rows,m.A.cols);
        action=@()amsStep(f);
    case 'ih'
        f=struct('config',makeResidentAlgebraConfig(), ...
            'reference',jsondecode(fileread(fullfile(inputsDir,'real-thermal-reference.json'))));
        action=@()ihStep(f);
    case {'hmatrix','hmatrix-profile'}
        f=struct('setup',radia_cuda_projection('assemble_hdiv_elements_sparse', ...
            fullfile(inputsDir,'curved-sphere.vol'),1,[1 .2 .1]), ...
            'reference',load(fullfile(inputsDir,'curved-cpu.mat')));
        f.config=radia.gpu.curvedHDivQuadrature(1);f.config.geometry_mode='curved';f.config.schema_version=1;
        for name={'cell_nodes','face_nodes','cell_vertices','face_vertices','host','kind','expo','charge_map'}
            f.config.(name{1})=f.setup.(name{1});
        end
        f.profile=strcmp(kind,'hmatrix-profile');
        if f.profile
            f.config.host_vector_cache=true;f.config.host_moment_cache=true;
            f.config.profile_assembly=true;
        end
        action=@()hmatrixStep(f);
    otherwise,error('radia:test:memoryKind','Unsupported memory workload');
end
end

function op=csr(a)
n=size(a,1);[columns,rows,values]=find(sparse(a.'));counts=accumarray(rows,1,[n 1]);
op=struct('rows',double(n),'cols',double(n),'rowptr',gpuArray(int32([0;cumsum(counts)])), ...
    'columns',gpuArray(int32(columns-1)),'values',gpuArray(values));
end
function csrStep(f)
[x,info]=radia_cuda_csr_pcg(f.op.rowptr,f.op.columns,f.op.values,f.rhs,f.inverse,1e-12,1000);
x=gather(x);assert(norm(x-f.expected)/norm(f.expected)<1e-11);
assert(norm(f.matrix*x-gather(f.rhs))/norm(gather(f.rhs))<1e-11 && info.true_relative_residual<=1e-11);
end
function icStep(f)
[h,~]=radia_cuda_ic('create',f.op,1.05);owner=onCleanup(@()radia_cuda_ic('release',h));
[x,info]=radia_cuda_ic('solve',h,f.rhs,1e-12,1000);x=gather(x);
assert(norm(x-f.expected)/norm(f.expected)<1e-10);
assert(norm(f.matrix*x-gather(f.rhs))/norm(gather(f.rhs))<1e-10 && info.true_relative_residual<=1e-10);
clear owner;assert(radia_cuda_ic('live_count')==0);
end
function fieldStep(f)
empty=zeros(0,1);[h,~]=radia_cuda_field('create_cloud',f.sources,f.charge,empty,empty,empty,1,0.05);
owner=onCleanup(@()radia_cuda_field('destroy',h));[value,~]=radia_cuda_field('field',h,f.targets,'direct');
assert(norm(value-f.expected,'fro')/norm(f.expected,'fro')<1e-12);
clear owner;assert(radia_cuda_field('live_count')==0);
end
function fmmStep(f)
h=radia_cuda_fmm('create',f.sources,f.targets,f.order,f.leaf,f.separation,false);
owner=onCleanup(@()radia_cuda_fmm('destroy',h));[potential,value,info]=radia_cuda_fmm('apply',h,f.charge,f.dipole);
assert(norm(gather(potential)-f.expectedP,'fro')/norm(f.expectedP,'fro')<=1e-5);
assert(norm(gather(value)-f.expectedE,'fro')/norm(f.expectedE,'fro')<=1e-4 && info.far_pairs>0);
clear owner;assert(radia_cuda_fmm('live_count')==0);
end
function projectionStep(f)
h=radia_cuda_projection('create',f.mesh,'hdiv',1);owner=onCleanup(@()radia_cuda_projection('destroy',h));
points=radia_cuda_projection('targets',h);
[value,~]=radia_cuda_projection('project',h,repmat([1 .2 .1],size(points,1),1));
assert(norm(value(:)-f.coefficients(:))/norm(f.coefficients(:))<=1e-11);
clear owner;assert(radia_cuda_projection('live_count')==0);
end
function amsStep(f)
h=radia_cuda_ams('create',f.hierarchy);owner=onCleanup(@()radia_cuda_ams('release',h));
[x,info]=radia_cuda_ams('solve',h,f.hierarchy.A,f.gpu_rhs,1e-10,4000,true);x=gather(x);free=logical(f.model.free);
assert(norm(x(free)-f.reference(free))/norm(f.reference(free))<=1e-8);
res=f.matrix*x-f.rhs;assert(norm(res(free))/norm(f.rhs(free))<=1e-10 && info.true_relative_residual<=1e-10);
clear owner;assert(all(structfun(@(v)v==0,radia_cuda_ams('resident_status'))));
end
function ihStep(f)
h=radia_cuda_ih('create',f.config);owner=onCleanup(@()radia_cuda_ih('destroy',h));
for step=1:5,radia_cuda_ih('update',h,1,300);end
[temperature,heat,info]=radia_cuda_ih('outputs',h);
assert(norm(temperature-f.reference.temperature(end,:)')/norm(f.reference.temperature(end,:))<=1e-8);
assert(norm(heat-f.reference.heat)/norm(f.reference.heat)<=1e-9);
assert(info.eddy_true_residual<=1e-10 && info.thermal_true_residual<=1e-10);
clear owner;assert(radia_cuda_ih('live_count')==0);
end
function hmatrixStep(f)
if f.profile
    h=radia_cuda_hmatrix('create_curved_hdiv',f.config);
else
    h=radia_cuda_hmatrix('create_hdiv_hacapk_device_dense',f.config,1e-12,64,2.,64,32*2^20);
end
owner=onCleanup(@()radia_cuda_hmatrix('destroy',h));
value=radia_cuda_hmatrix('solve_linear_material',h,f.setup.mass,f.setup.rhs,2.5,1e-11,1000);
assert(norm(value.x-f.reference.expected)/norm(f.reference.expected)<1e-8);
assert(norm(f.reference.rhs-f.reference.system*value.x)/norm(f.reference.rhs)<1e-10);
assert(value.true_relative_residual<=1e-11);
clear owner;assert(radia_cuda_hmatrix('live_count')==0);
end
