function result=runBoundaryLawGPU(caseFile,artifactDirectory,outputDirectory,refinementLevel)
%RUNBOUNDARYLAWGPU Prescribed patch-Zs solve; not nonlinear BH acceptance.
% Execute only after candidate/source and compute admission are agreed.
arguments
    caseFile (1,:) char
    artifactDirectory (1,:) char
    outputDirectory (1,:) char
    refinementLevel (1,1) double {mustBeMember(refinementLevel,[1,2])}
end
assert(~isfolder(outputDirectory),'Refusing to overwrite previous evidence');
mkdir(outputDirectory);
caseData=jsondecode(fileread(caseFile));
assert(strcmp(caseData.schema,'radia.ih-element-impedance-case.v1'));
expected=fullfile(artifactDirectory,'radia_cuda_ih.mexw64');
assert(isfile(expected),'Exact accepted candidate artifact is required');
[~,loaded]=inmem('-completenames');
for k=1:numel(loaded)
    if endsWith(loaded{k},'radia_cuda_ih.mexw64','IgnoreCase',true)
        assert(strcmpi(canonical(loaded{k}),canonical(expected)),...
            'A different IH artifact is loaded; preserve its owner session');
    end
end
oldPath=path;pathCleanup=onCleanup(@()path(oldPath));
addpath(artifactDirectory);
assert(strcmpi(canonical(which('radia_cuda_ih')),canonical(expected)),...
    'Loaded IH path differs from the selected artifact');
vertices=caseData.vertices_m;triangles=caseData.triangles_one_based;
parent=(1:size(triangles,1))';
for level=1:refinementLevel
    [vertices,triangles,parent]=refine(vertices,triangles,parent);
end
assert(4*size(vertices,1)+5<=512,'Complete block exceeds the declared bound');
coil=caseData.coil_vertices_m;
geometry=struct('vertices_m',vertices,'triangles',triangles,...
    'segments_m',[coil,coil([2:end,1],:)],...
    'segment_coil',ones(size(coil,1),1),...
    'regular_quad_degree',9,'singular_order',6);
nv=size(vertices,1);nf=size(triangles,1);
thermal=struct('schema_version',1,'sample_time_s',1,...
    'initial_temperature_K',300*ones(nv,1),'constant_coefficients',ones(nv,1),...
    'conductivity_scale',1,'convection_scale',0,...
    'stiffness',sparse(nv,nv),'convection',sparse(nv,nv),...
    'thermal_system',speye(nv),'heat_projection',sparse(nv,nf));
% This is an EM-only prescribed-boundary experiment. The thermal identity
% preserves state; it does not represent physical electromagnetic heat transfer.
material=struct('esim_cell_parameters',repmat(...
    [32,100,0,.001,2*pi*caseData.frequency_Hz,1e-8,.5],nf,1));
coilZ=.1+1i*.01;
[cfg,assembly]=radia.simulink.makeCompletePEECBEMConfig(geometry,coilZ,material,thermal);
adapter=which('radia.simulink.makeCompletePEECBEMConfig');
patchZ=complex(caseData.patch_zs_ohm_real,caseData.patch_zs_ohm_imag);
zs=patchZ(parent);assert(all(real(zs)>0)&&all(isfinite(zs)));
n=size(cfg.eddy,1);
cfg.eddy=cfg.eddy+reshape(cfg.operator_impedance_map*zs,n,n).';
cfg.loss_weights=cfg.loss_weights.*repelem(real(zs),3);
cfg.material_mode='fixed';cfg.solver_relative_tolerance=1e-13;
h=radia_cuda_ih('create',cfg);
handleCleanup=onCleanup(@()radia_cuda_ih('destroy',h));
current=caseData.current_peak_A;
radia_cuda_ih('update',h,current,300);
[~,facePower,info]=radia_cuda_ih('outputs',h);
em=radia_cuda_ih('electromagnetic_state',h);
ht=reshape(cfg.loss_projection*em,3,nf).';
reaction=em(assembly.unknown_indices.port)/current-coilZ;
points=zeros(3*nf,3);weights=zeros(3*nf,1);patches=zeros(3*nf,1);
sampleH=complex(zeros(3*nf,3));sampleQ=zeros(3*nf,1);
bary=[2/3,1/6,1/6;1/6,2/3,1/6;1/6,1/6,2/3];
for face=1:nf
    rows=(face-1)*3+(1:3);
    points(rows,:)=bary*vertices(triangles(face,:),:);
    weights(rows)=assembly.face_area_m2(face)/3;
    patches(rows)=parent(face);
    sampleH(rows,:)=repmat(ht(face,:),3,1);
    sampleQ(rows)=facePower(face)/assembly.face_area_m2(face);
end
result=struct('schema','radia.ih-element-impedance-bem.v1',...
    'status','computed-not-yet-crosschecked',...
    'scope','Prescribed heterogeneous boundary law, independent FEM pending',...
    'case_file_sha256',fileHash(caseFile),...
    'mex_sha256',fileHash(expected),'actual_mex_file',canonical(expected),...
    'assembly_adapter_sha256',fileHash(adapter),'actual_assembly_adapter_file',canonical(adapter),...
    'refinement_level',refinementLevel,'vertices_m',vertices,...
    'triangles_one_based',triangles,'parent_patch_one_based',parent,...
    'ndof',n,'frequency_Hz',caseData.frequency_Hz,'current_peak_A',current,...
    'total_power_W',sum(facePower),...
    'reaction_impedance_ohm_real',real(reaction),...
    'reaction_impedance_ohm_imag',imag(reaction),...
    'power_identity_relative',abs(sum(facePower)-.5*real(reaction)*current^2)/max(sum(facePower),1e-18),...
    'linear_true_relative_residual',info.eddy_true_residual,...
    'points_m',points,'patch_one_based',patches,'area_weights_m2',weights,...
    'H_t_real',real(sampleH),'H_t_imag',imag(sampleH),'q_surf_W_per_m2',sampleQ,...
    'nonlinear_volumetric_validation','not-performed',...
    'experimental_validation','not-performed');
save(fullfile(outputDirectory,'input-and-state.mat'),'cfg','geometry','assembly','em','zs');
f=fopen(fullfile(outputDirectory,'result.json'),'w');assert(f>=0);
writer=onCleanup(@()fclose(f));fprintf(f,'%s',jsonencode(result));clear writer handleCleanup pathCleanup;
end

function [vertices,triangles,parent]=refine(vertices,triangles,parent)
edges=containers.Map('KeyType','char','ValueType','double');
nextVertex=size(vertices,1);vertices=[vertices;zeros(3*size(triangles,1),3)];
children=zeros(4*size(triangles,1),3);childParent=zeros(size(children,1),1);
for k=1:size(triangles,1)
    tri=triangles(k,:);mid=zeros(1,3);
    for j=1:3
        edge=sort([tri(j),tri(mod(j,3)+1)]);key=sprintf('%d:%d',edge);
        if ~isKey(edges,key)
            nextVertex=nextVertex+1;
            vertices(nextVertex,:)=(vertices(edge(1),:)+vertices(edge(2),:))/2;
            edges(key)=nextVertex;
        end
        mid(j)=edges(key);
    end
    rows=(k-1)*4+(1:4);
    children(rows,:)=[tri(1),mid(1),mid(3);mid(1),tri(2),mid(2);...
        mid(3),mid(2),tri(3);mid(1),mid(2),mid(3)];
    childParent(rows)=parent(k);
end
vertices=vertices(1:nextVertex,:);triangles=children;parent=childParent;
end

function value=fileHash(file)
f=fopen(file,'rb');assert(f>=0);closer=onCleanup(@()fclose(f));
md=java.security.MessageDigest.getInstance('SHA-256');md.update(fread(f,Inf,'*uint8'));
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end

function value=canonical(file)
value=char(java.io.File(file).getCanonicalPath());
end
