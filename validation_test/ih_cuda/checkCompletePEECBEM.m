function result=checkCompletePEECBEM(artifact,out)
if ~isfolder(out),mkdir(out);end
result=struct('accepted',false,'research_only',true,'release_accepted',false);
try
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(artifact);addpath(fullfile(root,'packages','radia-gpu','matlab'));addpath(fileparts(mfilename('fullpath')));
g=struct('vertices_m',.01*[0,0,0;1,0,0;0,1,0;0,0,1],'triangles',[1,3,2;1,2,4;1,4,3;2,3,4],'regular_quad_degree',9,'singular_order',6);
p=[-.02,-.02,.03;.03,-.02,.03;.03,.03,.03;-.02,.03,.03];g.segments_m=[p,p([2,3,4,1],:)];g.segment_coil=ones(4,1);
t=makeResidentAlgebraConfig();t.thermal_system=speye(4);t.stiffness=sparse(4,4);t.convection=sparse(4,4);t.heat_projection=speye(4);t.initial_temperature_K=[300;310;320;330];t.sample_time_s=.001;t.constant_coefficients=ones(4,1);
m=struct('temperature_evaluation',speye(4),'temperature_axis_K',[290;600],'resistivity_ohm_m',repmat([1e-7,4e-7],4,1),'initial_permeability_H_per_m',repmat([4*pi*1e-7,4*pi*1e-7],4,1),'saturation_magnetization_A_per_m',zeros(4,2),'esim_cell_parameters',repmat([32,100,0,.0001,2*pi*1000,1e-8,.5],4,1),'initial_field_A_per_m',ones(4,1),'max_constitutive_iterations',100,'constitutive_tolerance',1e-8,'constitutive_relaxation',.8,'field_floor_A_per_m',1e-6);
[cfg,assembly]=radia.simulink.makeCompletePEECBEMConfig(g,.1+1i*.01,m,t);
save(fullfile(out,'input.mat'),'g','cfg','assembly');
h=radia_cuda_ih('create',cfg);c=onCleanup(@()radia_cuda_ih('destroy',h));
for step=1:3
radia_cuda_ih('update',h,1,300);[temperature,heat,info]=radia_cuda_ih('outputs',h);
end
result.accepted=true;result.temperature=temperature;result.heat=heat;result.info=info;result.em=radia_cuda_ih('electromagnetic_state',h);clear c;
catch e
result.error=getReport(e,'extended','hyperlinks','off');
end
result.matlab=version;f=fopen(fullfile(out,'result.json'),'w');fprintf(f,'%s',jsonencode(result));fclose(f);disp(result);
end
