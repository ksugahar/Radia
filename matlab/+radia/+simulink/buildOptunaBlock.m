function blockPath=buildOptunaBlock(modelName,options)
%BUILDOPTUNABLOCK Add a triggerable MATLAB Optuna block to a Simulink model.
arguments
 modelName (1,1) string; options.ObjectiveFcn (1,1) string; options.NumTrials (1,1) double {mustBeInteger,mustBePositive}=20
 options.Direction (1,1) string {mustBeMember(options.Direction,["minimize","maximize"])}="minimize"
 options.StoragePath (1,1) string=""; options.SampleTime_s (1,1) double {mustBePositive}=1
 options.LiveVisualization (1,1) logical=true; options.Save (1,1) logical=true
end
if ~bdIsLoaded(modelName), new_system(modelName); end
blockPath=modelName+"/Optuna Optimization";
parameters=sprintf('''%s'',%d,''%s'',''%s'',%.17g,%d',options.ObjectiveFcn,options.NumTrials,options.Direction,options.StoragePath,options.SampleTime_s,options.LiveVisualization);
add_block('simulink/User-Defined Functions/Level-2 MATLAB S-Function',char(blockPath),'FunctionName','radia_optuna_sfun','Parameters',parameters);
set_param(char(blockPath),'Position',[180 90 360 170]);
if options.Save, save_system(modelName); end
end
