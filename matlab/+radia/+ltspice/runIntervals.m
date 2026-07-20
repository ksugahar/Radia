function result=runIntervals(netlistFile,intervals,options)
%RUNINTERVALS Run long transients in reset-time windows with explicit state handoff.
arguments
 netlistFile (1,1) string {mustBeFile}; intervals (:,1) double {mustBePositive}
 options.Parameters (1,1) struct=struct(); options.Executable (1,1) string=""
 options.OutputDirectory (1,1) string=""; options.MaxStep_s (1,1) double {mustBePositive}=inf
end
if strlength(options.OutputDirectory)==0,folder=string(tempname("C:\temp"));else,folder=options.OutputDirectory;end
if ~isfolder(folder),mkdir(folder);end
runs=cell(numel(intervals),1); states=cell(numel(intervals),1); offset=0;
for k=1:numel(intervals)
 runFolder=fullfile(folder,sprintf("interval_%04d",k));
 source=fullfile(folder,sprintf("state_%04d.cir",k));
 if k==1
   initial=struct("schema","radia.ltspice.transient_state.v1","time_s",0, ...
     "node_names",strings(0,1),"node_voltages_V",zeros(0,1), ...
     "inductor_names",strings(0,1),"inductor_currents_A",zeros(0,1));
 else, initial=states{k-1}; end
 radia.ltspice.applyTransientState(netlistFile,initial,source,Duration_s=intervals(k),MaxStep_s=options.MaxStep_s);
 runs{k}=radia.ltspice.run(source,Parameters=options.Parameters,Executable=options.Executable,OutputDirectory=runFolder);
 runs{k}.waveform.values(:,1)=runs{k}.waveform.values(:,1)+offset;
 states{k}=radia.ltspice.extractTransientState(runs{k}); offset=offset+intervals(k);
end
result=struct("schema","radia.ltspice.interval_run.v1","runs",{runs},"states",{states}, ...
 "total_duration_s",sum(intervals),"output_directory",folder);
end
