function blockPath = buildLTspiceBlock(modelName, options)
%BUILDLTSPICEBLOCK Create a readable Level-2 MATLAB LTspice coupling block.
arguments
 modelName (1,1) string; options.Netlist (1,1) string; options.InputNames (:,1) string=strings(0,1); options.OutputTraces (:,1) string=strings(0,1)
 options.InputName (1,1) string="control"; options.OutputTrace (1,1) string=""; options.SampleTime_s (1,1) double=1e-3
 options.MaxStep_s (1,1) double {mustBePositive}=inf; options.Timeout_s (1,1) double {mustBePositive}=300; options.Executable (1,1) string=""; options.Save (1,1) logical=true
end
if ~bdIsLoaded(modelName), new_system(modelName); end
blockPath=modelName+"/LTspice Circuit";
if isempty(options.InputNames),options.InputNames=options.InputName;end
if isempty(options.OutputTraces),options.OutputTraces=options.OutputTrace;end
if any(strlength(options.OutputTraces)==0),error("radia:simulink:LTspiceOutputTrace","OutputTraces is required.");end
parameters=sprintf('''%s'',%s,%s,%.17g,''C:\\temp\\radia_ltspice_block'',%.17g,%.17g,''%s''',options.Netlist,cellLiteral(options.InputNames),cellLiteral(options.OutputTraces),options.SampleTime_s,options.MaxStep_s,options.Timeout_s,options.Executable);
add_block('simulink/User-Defined Functions/Level-2 MATLAB S-Function',char(blockPath),'FunctionName','radia_ltspice_sfun','Parameters',parameters);
if options.Save, save_system(modelName); end
end
function text=cellLiteral(values),text="{"+join("'"+replace(values(:)',"'","''")+"'",",")+"}";end
