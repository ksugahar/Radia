function result=checkPersistentElementESIM(artifact,outdir)
%CHECKPERSISTENTELEMENTESIM Explicit generic operator ownership diagnostic.
% This does not certify complete coil/body reaction geometry assembly.
arguments
    artifact (1,:) char
    outdir (1,:) char
end
if ~isfolder(outdir), mkdir(outdir); end
addpath(artifact);
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'packages','radia-gpu','matlab'));
addpath(fileparts(mfilename('fullpath')));
cfg=makeResidentAlgebraConfig();
cfg.material_mode='element-esim-operators';
cfg.eddy=complex(diag([2,2,2]),.1*eye(3));
cfg.loss_projection=speye(3); cfg.loss_components=1;
cfg.thermal_system=speye(3);cfg.stiffness=sparse(3,3);cfg.convection=sparse(3,3);
cfg.heat_projection=.001*speye(3);cfg.rhs_per_amp=100*ones(3,1);
cfg.initial_temperature_K=[300;400;500];cfg.loss_weights=.5*ones(3,1);cfg.sample_time_s=.001;
cfg.temperature_evaluation=sparse([.5,.5,0;0,1,0;0,0,1]);
cfg.operator_impedance_map=sparse([1,5,9],[1,2,3],10,9,3);
cfg.temperature_axis_K=[300;600];cfg.resistivity_ohm_m=[1e-7,4e-7;2e-7,8e-7;3e-7,12e-7];
cfg.initial_permeability_H_per_m=repmat([1.2566370614359173e-4,6.283185307179586e-5],3,1);
cfg.saturation_magnetization_A_per_m=zeros(3,2);
cfg.esim_cell_parameters=[32,50,0,.001,6283.185307179586,1e-6,.5;32,50,1,.001,6283.185307179586,1e-6,.5;32,50,0,.001,6283.185307179586,1e-6,.5];
cfg.initial_field_A_per_m=[100;200;300];cfg.max_constitutive_iterations=100;
cfg.constitutive_tolerance=1e-8;cfg.constitutive_relaxation=.8;cfg.field_floor_A_per_m=1e-6;
assert(radia_cuda_ih('live_count')==0);
h=radia_cuda_ih('create',cfg);cleanup=onCleanup(@()destroyIfLive(h));
t=zeros(6,3);q=t;
for step=0:5
 [t(step+1,:),q(step+1,:),info]=radia_cuda_ih('outputs',h);
 if step<5,radia_cuda_ih('update',h,1,300);end
end
assert(info.steps==5 && info.constitutive_true_residual<=1e-8 && info.eddy_true_residual<=1e-10 && info.thermal_true_residual<=1e-10);
em=radia_cuda_ih('electromagnetic_state',h);rejected=false;
try,radia_cuda_ih('update',h,realmax,300);catch,rejected=true;end
[afterT,afterQ,afterInfo]=radia_cuda_ih('outputs',h);
assert(rejected && isequal(afterT,t(end,:)') && isequal(afterQ,q(end,:)') && afterInfo.steps==5 && isequal(em,radia_cuda_ih('electromagnetic_state',h)));
radia_cuda_ih('reset',h);assert(all(radia_cuda_ih('electromagnetic_state',h)==0));clear cleanup;assert(radia_cuda_ih('live_count')==0);
model='ih_element_esim_diagnostic';path=fullfile(outdir,[model '.slx']);
if bdIsLoaded(model),error('Owned diagnostic model already loaded');end
new_system(model);closer=onCleanup(@()close_system(model,0));
set_param(model,'Solver','FixedStepDiscrete','FixedStep','.001','StopTime','.005');
workspace=get_param(model,'ModelWorkspace');assignin(workspace,'ihConfig',cfg);
add_block('simulink/Sources/Constant',[model '/Current'],'Value','1');
add_block('simulink/Sources/Constant',[model '/Ambient'],'Value','300');
add_block('simulink/User-Defined Functions/Level-2 MATLAB S-Function',[model '/IH'],'FunctionName','radia_ih_cuda_sfun','Parameters','ihConfig');
add_block('simulink/Sinks/To Workspace',[model '/Temperature'],'VariableName','temperature','SaveFormat','Timeseries');
add_block('simulink/Sinks/To Workspace',[model '/Heat'],'VariableName','heat','SaveFormat','Timeseries');
set_param(model,'SimulationCommand','update');
add_line(model,'Current/1','IH/1');add_line(model,'Ambient/1','IH/2');add_line(model,'IH/1','Temperature/1');add_line(model,'IH/2','Heat/1');
set_param(model,'SimulationCommand','update');save_system(model,path);close_system(model,0);open_system(path,'loadonly');
one=sim(model,'ReturnWorkspaceOutputs','on');two=sim(model,'ReturnWorkspaceOutputs','on');
assert(isequal(one.temperature.Data,two.temperature.Data)&&isequal(one.heat.Data,two.heat.Data));
assert(norm(one.temperature.Data-t,'fro')<1e-10 && norm(one.heat.Data-q,'fro')<1e-10);assert(radia_cuda_ih('live_count')==0);
result=struct('accepted',true,'scope','Generic temperature-dependent operator MEX/Simulink diagnostic','full_reaction_assembly',false,'release_accepted',false,'temperature',t,'heat',q,'info',info,'repeat_identical',true,'save_close_reopen',true,'contexts_released',true,'failed_candidate_preserved',true,'matlab',version,'toolkit_accepted',false);
f=fopen(fullfile(outdir,'result.json'),'w');assert(f>=0);writer=onCleanup(@()fclose(f));fprintf(f,'%s',jsonencode(result));clear writer closer;
end
function destroyIfLive(h)
radia_cuda_ih('destroy',h);
end

