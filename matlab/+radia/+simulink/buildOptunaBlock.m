function blockPath=buildOptunaBlock(modelName,options)
%BUILDOPTUNABLOCK Add a triggerable MATLAB Optuna block to a Simulink model.
arguments
 modelName (1,1) string; options.ObjectiveFcn (1,1) string; options.NumTrials (1,1) double {mustBeInteger,mustBePositive}=20
 options.Direction (1,1) string {mustBeMember(options.Direction,["minimize","maximize"])}="minimize"
 options.Directions string=string.empty
 options.StoragePath (1,1) string=""; options.SampleTime_s (1,1) double {mustBePositive}=1
 options.Sampler (1,1) string {mustBeMember(options.Sampler, ...
  ["auto","random","tpe","cmaes","motpe","nsgaii"])}="auto"
 options.LiveVisualization (1,1) logical=false; options.Save (1,1) logical=true
end
if ~bdIsLoaded(modelName), new_system(modelName); end
blockPath=modelName+"/Optuna Optimization";
directions=options.Directions;
if isempty(directions), directions=options.Direction; end
if any(~ismember(directions,["minimize","maximize"]))
 error("radia:simulink:OptunaDirection","Directions must contain only minimize or maximize.");
end
directionExpression=formatDirections(directions);
objectiveExpression=formatObjective(options.ObjectiveFcn);
initialParameters=sprintf('%s,%d,%s,''%s'',%.17g,%d,''%s''', ...
 objectiveExpression,options.NumTrials,directionExpression, ...
 options.StoragePath,options.SampleTime_s,options.LiveVisualization, ...
 options.Sampler);
add_block('simulink/User-Defined Functions/Level-2 MATLAB S-Function', ...
 char(blockPath),'FunctionName','radia_optuna_sfun', ...
 'Parameters',initialParameters);
set_param(char(blockPath),'Position',[180 90 360 170]);
mask=Simulink.Mask.create(blockPath);
mask.Description="Incremental MATLAB Optuna study. One trial executes per sample; connect the live outputs to an Optuna Monitor.";
mask.addParameter(Type="edit",Name="objective_fcn",Prompt="Objective function", ...
 Value=objectiveExpression,Evaluate="on");
mask.addParameter(Type="edit",Name="num_trials",Prompt="Number of trials", ...
 Value=string(options.NumTrials),Evaluate="on");
mask.addParameter(Type="edit",Name="directions",Prompt="Objective directions", ...
 Value=directionExpression,Evaluate="on");
mask.addParameter(Type="edit",Name="storage_path",Prompt="Study MAT file", ...
 Value=quoteString(options.StoragePath),Evaluate="on");
mask.addParameter(Type="edit",Name="sample_time_s",Prompt="Trial sample time (s)", ...
 Value=compose("%.17g",options.SampleTime_s),Evaluate="on");
liveValue="off"; if options.LiveVisualization,liveValue="on";end
mask.addParameter(Type="checkbox",Name="live_visualization", ...
 Prompt="External MATLAB monitor",Value=liveValue);
mask.addParameter(Type="popup",Name="sampler_name",Prompt="Sampler", ...
 TypeOptions={"auto","random","tpe","cmaes","motpe","nsgaii"}, ...
 Value=options.Sampler);
set_param(blockPath,"Parameters", ...
 "objective_fcn,num_trials,directions,storage_path,sample_time_s," + ...
 "strcmp(live_visualization,'on'),sampler_name");
mask.Display="disp('Optuna Optimization');" + ...
 "port_label('input',1,'start');" + ...
 "port_label('output',1,'best');port_label('output',2,'best trial');" + ...
 "port_label('output',3,'status');port_label('output',4,'trials');" + ...
 "port_label('output',5,'last');port_label('output',6,'elapsed');" + ...
 "port_label('output',7,'best update');port_label('output',8,'pareto N');" + ...
 "port_label('output',9,'pareto X');port_label('output',10,'pareto Y');" + ...
 "port_label('output',11,'pareto rev');";
set_param(blockPath,"UserData",struct("execution","one-trial-per-sample", ...
 "visualization","simulink-monitor","browser_required",false),"UserDataPersistent","on");
if options.Save, save_system(modelName); end
end

function expression=formatDirections(directions)
quoted='"'+reshape(string(directions),1,[])+'"';
if isscalar(quoted)
 expression=char(quoted);
else
 expression=char("["+strjoin(quoted,",")+"]");
end
end

function expression=formatObjective(value)
value=strtrim(string(value));
if startsWith(value,"@")
 expression=char(value);
else
 expression=char(quoteString(value));
end
end

function expression=quoteString(value)
expression="'"+replace(string(value),"'","''")+"'";
end
