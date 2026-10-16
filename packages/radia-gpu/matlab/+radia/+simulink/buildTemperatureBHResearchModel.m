function path=buildTemperatureBHResearchModel(config,path,options)
%BUILDTEMPERATUREBHRESEARCHMODEL Native element-Zs research demonstration.
% A single-valued amplitude B(H,T), not hysteresis history or bulk thermal FE.
arguments
    config (1,1) struct
    path (1,:) char
    options.CurrentA (1,1) double {mustBeFinite}=1
    options.AmbientK (1,1) double {mustBeFinite,mustBePositive}=300
    options.Samples (1,1) double {mustBeInteger,mustBePositive}=5
    options.Overwrite (1,1) logical=false
end
assert(options.Samples<=1e6,'Sample budget exceeded');
assert(isfield(config,'bh_knots')&&isfield(config,'bh_values_T'),...
    'Explicit temperature-dependent amplitude B-H material required');
assert(strcmp(config.material_mode,'element-esim-operators'),'Unsupported material mode');
dt=config.sample_time_s;assert(isscalar(dt)&&isfinite(dt)&&dt>0,'Invalid sample time');
[folder,model,ext]=fileparts(path);
assert(strcmpi(ext,'.slx')&&isvarname(model),'Expected a named SLX file');
assert(~bdIsLoaded(model),'Close the existing model explicitly before creating another');
assert(options.Overwrite||~isfile(path),'Refusing to overwrite an existing model');
if ~isempty(folder)&&~isfolder(folder),mkdir(folder);end
config.observe_element_impedance=true;
new_system(model);closer=onCleanup(@()close_system(model,0));
set_param(model,'Solver','FixedStepDiscrete','FixedStep',num2str(dt,17),...
    'StopTime',num2str(options.Samples*dt,17));
impedanceNote=['Initial and reset Zs=0 is a nonphysical placeholder. ' ...
    'Omit sample zero when plotting physical impedance. Later Zs and heat ' ...
    'belong to the last accepted EM update at its starting temperature.'];
set_param(model,'Description',impedanceNote);
workspace=get_param(model,'ModelWorkspace');assignin(workspace,'ihConfig',config);
add_block('simulink/Sources/Constant',[model '/Current'],...
    'Value',num2str(options.CurrentA,17));
add_block('simulink/Sources/Constant',[model '/Ambient'],...
    'Value',num2str(options.AmbientK,17));
add_block('simulink/User-Defined Functions/Level-2 MATLAB S-Function',...
    [model '/IH'],'FunctionName','radia_ih_cuda_sfun','Parameters','ihConfig');
set_param([model '/IH'],'Description',impedanceNote);
names={'Temperature','Heat','ZsReal','ZsImag'};
variables={'temperature','heat','zsReal','zsImag'};
for k=1:4
    add_block('simulink/Sinks/To Workspace',[model '/' names{k}],...
        'VariableName',variables{k},'SaveFormat','Timeseries',...
        'SampleTime',num2str(dt,17));
end
set_param(model,'SimulationCommand','update');
add_line(model,'Current/1','IH/1');add_line(model,'Ambient/1','IH/2');
for k=1:4,add_line(model,sprintf('IH/%d',k),[names{k} '/1']);end
set_param(model,'SimulationCommand','update');
Simulink.BlockDiagram.arrangeSystem(model);
save_system(model,path);clear closer;
end
