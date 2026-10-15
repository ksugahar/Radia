function [cfg,assembly]=makeTemperatureBHConfig()
% Original public tetrahedron and external square filament, SI units.
% Synthetic anhysteretic B-H tables; never experimental material data.
g=struct('vertices_m',.01*[0,0,0;1,0,0;0,1,0;0,0,1],...
 'triangles',[1,3,2;1,2,4;1,4,3;2,3,4],...
 'regular_quad_degree',9,'singular_order',6);
p=[-.02,-.02,.03;.03,-.02,.03;.03,.03,.03;-.02,.03,.03];
g.segments_m=[p,p([2,3,4,1],:)];g.segment_coil=ones(4,1);
t=struct('schema_version',1,'sample_time_s',.01);
t.initial_temperature_K=[300;310;320;330];
t.constant_coefficients=ones(4,1);t.conductivity_scale=1;t.convection_scale=1;
m=struct('temperature_axis_K',[290;600],...
 'resistivity_ohm_m',repmat([1e-7,4e-7],4,1),...
 'esim_cell_parameters',repmat([64,400,0,.001,2*pi*1000,1e-10,.3],4,1),...
 'initial_field_A_per_m',ones(4,1),'max_constitutive_iterations',300,...
 'constitutive_tolerance',1e-11,'constitutive_relaxation',.5,'field_floor_A_per_m',1e-6);
m.bh_knots=[0;100;1000;10000];m.bh_values_T=zeros(4,2,4);
mu0=4*pi*1e-7;
for e=1:4
 for j=1:2
  susceptibility=(1+.1*e)*(3-j);
  b=mu0*(m.bh_knots+susceptibility*100*tanh(m.bh_knots/100));
  m.bh_values_T(e,j,:)=b;
 end
end
m.initial_permeability_H_per_m=(m.bh_values_T(:,:,2)-m.bh_values_T(:,:,1))/100;
m.saturation_magnetization_A_per_m=m.bh_values_T(:,:,end)/mu0-m.bh_knots(end);
% First assemble geometry to obtain consistent surface thermal operators.
a=radia_cuda_ih('assemble_peec_bem',g);
% Thin surface layer: independent synthetic rho*c=3.8e6 J/(m3 K),
% thickness 2 mm, conductivity 20 W/(m K), convection 10 W/(m2 K).
capacity=3.8e6*.002*a.mass;
t.stiffness=sparse(20*.002*a.surface_stiffness);
t.convection=sparse(10*a.mass);
t.thermal_system=sparse(capacity+t.sample_time_s*(t.stiffness+t.convection));
t.heat_projection=sparse(4,4);m.temperature_evaluation=sparse(4,4);
for e=1:4
 t.heat_projection(g.triangles(e,:),e)=1/3;
 m.temperature_evaluation(e,g.triangles(e,:))=1/3;
end
[cfg,assembly]=radia.simulink.makeCompletePEECBEMConfig(g,.1+1i*.01,m,t);
cfg.observe_element_impedance=true;cfg.solver_relative_tolerance=1e-13;
assembly.thermal_capacity=capacity;
assembly.public_case='Original tetrahedron, external square coil, synthetic amplitude B-H';
end
