function result=runLTspiceInterval(netlistFile,inputNames,inputValues,state,options)
%RUNLTSPICEINTERVAL Advance an LTspice circuit by one reset-time interval.
arguments
 netlistFile (1,1) string {mustBeFile}; inputNames (:,1) string; inputValues (:,1) double
 state (1,1) struct
 options.Duration_s (1,1) double {mustBePositive}; options.OutputDirectory (1,1) string
 options.MaxStep_s (1,1) double {mustBePositive}=inf; options.Timeout_s (1,1) double {mustBePositive}=300
 options.Executable (1,1) string=""
end
if numel(inputNames)~=numel(inputValues),error("radia:simulink:LTspiceInputCount","Input name/value counts differ.");end
if ~isfolder(options.OutputDirectory),mkdir(options.OutputDirectory);end
stateNetlist=fullfile(options.OutputDirectory,"interval_state.cir");
radia.ltspice.applyTransientState(netlistFile,state,stateNetlist,Duration_s=options.Duration_s,MaxStep_s=options.MaxStep_s);
signals=struct();
for k=1:numel(inputNames)
 field=matlab.lang.makeValidName(inputNames(k));
 if field~=inputNames(k),error("radia:simulink:LTspiceInputName","Input names must be valid MATLAB field names: %s",inputNames(k));end
 signals.(field)=[0,inputValues(k);options.Duration_s,inputValues(k)];
end
simulation=radia.simulink.runLTspice(stateNetlist,InputSignals=signals,Executable=options.Executable, ...
 OutputDirectory=options.OutputDirectory,Timeout_s=options.Timeout_s);
nextState=radia.ltspice.extractTransientState(simulation);
nextState.time_s=state.time_s+options.Duration_s;
result=struct("schema","radia.simulink.ltspice.interval.v1","simulation",simulation, ...
 "state",nextState,"duration_s",options.Duration_s,"input_names",inputNames,"input_values",inputValues);
end
