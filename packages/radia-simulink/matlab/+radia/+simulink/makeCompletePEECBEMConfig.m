function [cfg,assembly]=makeCompletePEECBEMConfig(geometry,coilImpedance,material,thermal)
%MAKECOMPLETEPEECBEMCONFIG Bounded flat-P1 complete reaction block.
% Preserve the body BIE, mass action, surface-Poisson projection and its
% bilinear transpose without computing a CPU inverse or a frozen body field.
assembly=radia_cuda_ih('assemble_peec_bem',geometry);
m=assembly.mass;p=assembly.surface_stiffness;s=assembly.single_layer;d=assembly.double_layer;
f=assembly.field_gradient;q=assembly.incident_weak_map;area=assembly.face_area_m2;
v=size(m,1);coils=size(q,2);faces=numel(area);
assert(isequal(size(coilImpedance),[coils,coils]) && all(isfinite(coilImpedance),'all'),'Invalid coil impedance');
assert(size(material.esim_cell_parameters,1)==faces && size(material.esim_cell_parameters,2)==7,'One ESIM cell per face required');
omega=material.esim_cell_parameters(1,5);
assert(all(material.esim_cell_parameters(:,5)==omega) && isfinite(omega) && omega>0,'One positive angular frequency required');
coil=1:coils;port=coils+1;psi=port+(1:v);psiGauge=psi(end)+1;
phi=psiGauge+(1:v);bodyGauge=phi(end)+1;flux=bodyGauge+(1:v);
reaction=flux(end)+(1:v);reactionGauge=reaction(end)+1;n=reactionGauge;
a=complex(zeros(n));a(coil,coil)=coilImpedance;a(coil,port)=-1;a(port,coil)=1;
a(psi,psi)=p;a(psi,coil)=-q;a(psi,psiGauge)=1;a(psiGauge,psi)=1;
a(phi,phi)=.5*m-d;a(phi,flux)=s;a(phi,psi)=-m;
gauge=m*ones(v,1);a(phi,bodyGauge)=gauge;a(bodyGauge,phi)=gauge.';
a(flux,flux)=m;a(reaction,reaction)=p;a(reaction,reactionGauge)=1;a(reactionGauge,reaction)=1;
magnetic=zeros(coils,v);
for face=1:faces
 rows=(face-1)*3+(1:3);normal=assembly.face_normals(face,:);
 crossNormal=[0,-normal(3),normal(2);normal(3),0,-normal(1);-normal(2),normal(1),0];
 magnetic=magnetic+assembly.unit_a(rows,:).'*crossNormal*f(rows,:)*area(face);
end
a(coil,phi)=1i*omega*magnetic;a(coil,reaction)=q.';
entryMap=sparse(n*n,faces);mu0=4*pi*1e-7;
for face=1:faces
 rows=(face-1)*3+(1:3);kz=area(face)*(f(rows,:).'*f(rows,:));
 derivative=complex(zeros(n));derivative(flux,phi)=-kz/(1i*omega*mu0);derivative(reaction,phi)=-kz;
 entryMap(:,face)=reshape(derivative.',[],1); % native dense CSR is row-major
end
cfg=thermal;cfg.schema_version=1;cfg.material_mode='element-esim-operators';
cfg.eddy=a;cfg.operator_impedance_map=entryMap;
cfg.preconditioner_precision='fp64';cfg.preconditioner_block_ends=[port;n];
cfg.loss_projection=sparse(size(f,1),n);cfg.loss_projection(:,phi)=f;
cfg.loss_components=3;cfg.loss_weights=repelem(.5*area,3);
cfg.rhs_per_amp=zeros(n,1);cfg.rhs_per_amp(port)=1;
fields=fieldnames(material);for k=1:numel(fields),cfg.(fields{k})=material.(fields{k});end
assembly.unknown_indices=struct('coil',coil,'port',port,'incident_potential',psi,'body_potential',phi,'mass_flux',flux,'reaction_potential',reaction);
assembly.complete_reaction='magnetic plus surface-Poisson transpose work, bilinear phasors';
end
