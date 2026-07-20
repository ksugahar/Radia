function blockPath=buildHystereticLTspiceBlock(modelName,options)
%BUILDHYSTERETICLTSPICEBLOCK Add the full LTspice/hysteresis coupling block.
arguments
 modelName (1,1) string;options.Netlist (1,1) string {mustBeFile};options.Tables (1,1) cell
 options.K (1,1) double {mustBeInteger,mustBePositive}=1;options.EtaOrChi (1,1) double {mustBePositive}=0.1
 options.HysteresisKind (1,1) string {mustBeMember(options.HysteresisKind,["play","energy"])}="play";options.Epsilon (1,1) double {mustBePositive}=1e-8
 options.CommandName (1,1) string="command";options.BackEmfName (1,1) string="back_emf";options.CurrentTrace (1,1) string
 options.Turns (1,1) double {mustBePositive};options.CoreArea_m2 (1,1) double {mustBePositive};options.MagneticPath_m (1,1) double {mustBePositive};options.CoreVolume_m3 (1,1) double {mustBePositive}
 options.GapPathFactor (1,1) double {mustBeNonnegative}=2;options.SampleTime_s (1,1) double {mustBePositive}
 options.MaxIterations (1,1) double {mustBeInteger,mustBePositive}=12;options.RelativeTolerance (1,1) double {mustBePositive}=1e-3;options.Relaxation (1,1) double {mustBePositive}=0.5
 options.MaxStep_s (1,1) double {mustBePositive}=inf;options.Timeout_s (1,1) double {mustBePositive}=300;options.CouplingSamples (1,1) double {mustBeInteger,mustBePositive}=101
 options.ConfigFile (1,1) string="";options.Save (1,1) logical=true
end
config=options;config=rmfield(config,["ConfigFile","Save"]);
if strlength(options.ConfigFile)==0,folder="C:\temp\radia_hysteretic_ltspice_configs";if ~isfolder(folder),mkdir(folder);end;configFile=fullfile(folder,modelName+"_config.mat");else,configFile=options.ConfigFile;end
save(configFile,"config");if ~bdIsLoaded(modelName),new_system(modelName);end
blockPath=modelName+"/Hysteretic LTspice Plant";params=sprintf('''%s'',%.17g,''C:\\temp\\radia_hysteretic_ltspice_block''',configFile,options.SampleTime_s);
add_block('simulink/User-Defined Functions/Level-2 MATLAB S-Function',char(blockPath),'FunctionName','radia_hysteretic_ltspice_sfun','Parameters',params,'Position',[150 70 350 140]);
if options.Save,save_system(modelName);end
end
