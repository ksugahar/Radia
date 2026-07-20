function result = runLTspiceClosedLoop(netlistFile, controller, options)
%RUNLTSPICECLOSEDLOOP Causal sampled-data control using cumulative LTspice runs.
arguments
 netlistFile (1,1) string {mustBeFile}; controller (1,1) function_handle
 options.SampleTime_s (1,1) double {mustBePositive}; options.Steps (1,1) double {mustBeInteger,mustBePositive}
 options.InitialInput (1,1) double = 0; options.InputName (1,1) string = "control"
 options.OutputTrace (1,1) string; options.Parameters (1,1) struct = struct()
 options.OutputDirectory (1,1) string = "C:\temp\radia_ltspice_closed_loop"
end
t=(0:options.Steps-1)'*options.SampleTime_s; u=zeros(options.Steps,1); y=zeros(options.Steps,1); u(1)=options.InitialInput;
for k=1:options.Steps
 p=options.Parameters; p.Tstop=max(t(k),eps);
 sig=struct(); sig.(options.InputName)=[t(1:k),u(1:k)];
 r=radia.simulink.runLTspice(netlistFile,InputSignals=sig,Parameters=p, ...
   OutputDirectory=fullfile(options.OutputDirectory,sprintf("step_%06d",k)));
 j=find(r.waveform.names==options.OutputTrace,1);
 if isempty(j), error("radia:simulink:LTspiceTrace","Trace not found: %s",options.OutputTrace); end
 y(k)=r.waveform.values(end,j);
 if k<options.Steps, u(k+1)=controller(t(k),y(k),u(k)); end
end
result=struct("schema","radia.simulink.ltspice.closed_loop.v1","time_s",t,"input",u,"output",y);
end
