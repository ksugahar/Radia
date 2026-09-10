function ltspiceSFunction(block)
%LTSPICESFUNCTION Level-2 block for sampled LTspice closed-loop coupling.
setup(block);
end
function setup(b)
b.NumDialogPrms=8; b.NumInputPorts=1; b.NumOutputPorts=1;
b.SetPreCompInpPortInfoToDynamic; b.SetPreCompOutPortInfoToDynamic;
b.InputPort(1).Dimensions=numel(names(b.DialogPrm(2).Data)); b.InputPort(1).DirectFeedthrough=true;
b.OutputPort(1).Dimensions=numel(names(b.DialogPrm(3).Data)); b.SampleTimes=[b.DialogPrm(4).Data 0];
b.SimStateCompliance='CustomSimState';
b.RegBlockMethod('PostPropagationSetup',@postSetup); b.RegBlockMethod('Start',@start); b.RegBlockMethod('Outputs',@outputs); b.RegBlockMethod('Update',@update); b.RegBlockMethod('GetSimState',@getSimState);b.RegBlockMethod('SetSimState',@setSimState);b.RegBlockMethod('Terminate',@terminate);
end
function postSetup(b)
b.NumDworks=1; b.Dwork(1).Name='step'; b.Dwork(1).Dimensions=1; b.Dwork(1).DatatypeID=0; b.Dwork(1).Complexity='Real'; b.Dwork(1).UsedAsDiscState=true;
end
function start(b)
b.Dwork(1).Data=0; root=string(b.DialogPrm(5).Data);if ~isfolder(root),mkdir(root);end
logPath=radia.simulink.internal.runArtifacts("logpath",root,b.BlockHandle);
entry=struct("state",emptyState(),"pending_state",[],"pending_output",[], ...
 "has_pending",false,"folder",string(tempname(root)), ...
 "log_path",logPath,"kept",strings(0,1),"keep_count",defaultKeptStepFolders(),"failed",false);
radia.simulink.internal.runArtifacts("begin",logPath,sprintf( ...
 'schema=radia.simulink.ltspice.run.v1 block=%s netlist=%s sample_time_s=%.17g started=%s artifacts=%s kept_steps=%d', ...
 getfullname(b.BlockHandle),string(b.DialogPrm(1).Data),b.DialogPrm(4).Data, ...
 string(datetime("now","TimeZone","local"),"yyyy-MM-dd'T'HH:mm:ssXXX"),entry.folder,entry.keep_count));
mkdir(entry.folder); storage("set",key(b),entry);
end
function outputs(b)
k=b.Dwork(1).Data+1; entry=storage("get",key(b));
if entry.has_pending,b.OutputPort(1).Data=entry.pending_output;return;end
inames=names(b.DialogPrm(2).Data); traces=names(b.DialogPrm(3).Data);
stepFolder=fullfile(entry.folder,sprintf("step_%06d",k));
try
 r=radia.simulink.runLTspiceInterval(string(b.DialogPrm(1).Data),inames,b.InputPort(1).Data(:),entry.state, ...
   Duration_s=b.DialogPrm(4).Data,OutputDirectory=stepFolder,MaxStep_s=b.DialogPrm(6).Data, ...
   Timeout_s=b.DialogPrm(7).Data,Executable=string(b.DialogPrm(8).Data));
 y=zeros(numel(traces),1); for n=1:numel(traces),j=find(r.simulation.waveform.names==traces(n),1);if isempty(j),error("radia:simulink:LTspiceTrace","Trace not found: %s",traces(n));end;y(n)=real(r.simulation.waveform.values(end,j));end
 b.OutputPort(1).Data=y;entry.pending_state=r.state;entry.pending_output=y;
 entry.has_pending=true;entry.kept(end+1,1)=string(stepFolder);storage("set",key(b),entry);
 radia.simulink.internal.runArtifacts("append",entry.log_path,sprintf( ...
  'step=%d t=%.17g inputs=[%s] outputs=[%s] artifacts=%s',k,b.CurrentTime, ...
  join(compose("%.17g",b.InputPort(1).Data(:).'),","),join(compose("%.17g",y.'),","),stepFolder));
catch cause
 entry.failed=true;storage("set",key(b),entry);
 radia.simulink.internal.runArtifacts("append",entry.log_path,sprintf( ...
  'step=%d t=%.17g FAILED identifier=%s message=%s',k,b.CurrentTime,cause.identifier,cause.message));
 error("radia:simulink:LTspiceStepFailed","LTspice step %d at Simulink time %.17g failed: %s",k,b.CurrentTime,cause.message);
end
end
function update(b)
entry=storage("get",key(b));if ~entry.has_pending,return;end
entry.state=entry.pending_state;entry.pending_state=[];entry.pending_output=[];
entry.has_pending=false;entry.kept=radia.simulink.internal.runArtifacts("prune",entry.kept,entry.keep_count);
storage("set",key(b),entry);b.Dwork(1).Data=b.Dwork(1).Data+1;
end
function state=getSimState(b),state=struct("step",b.Dwork(1).Data,"entry",storage("get",key(b)));end
function setSimState(b,state)
if ~isstruct(state)||~all(isfield(state,["step","entry"])),error("radia:simulink:LTspiceSimState","Invalid LTspice block SimState.");end
b.Dwork(1).Data=state.step;entry=state.entry;if ~isfolder(entry.folder),mkdir(entry.folder);end;storage("set",key(b),entry)
end
function terminate(b)
k=key(b);folder="";failed=false;logPath="";
try entry=storage("get",k);folder=entry.folder;failed=entry.failed;logPath=entry.log_path;catch,end
storage("remove",k);
if strlength(logPath)>0,radia.simulink.internal.runArtifacts("append",logPath,sprintf('finished failed=%d artifacts_retained=%d',failed,failed));end
radia.simulink.internal.runArtifacts("finish",folder,failed);
end
function answer=names(value)
answer=string(value(:));
answer=answer(strlength(answer)>0);if isempty(answer),error("radia:simulink:LTspiceNames","At least one name is required.");end
end
function answer=key(b),answer=sprintf('%.0f',b.BlockHandle);end
function state=emptyState(),state=struct("schema","radia.ltspice.transient_state.v1","time_s",0,"node_names",strings(0,1),"node_voltages_V",zeros(0,1),"inductor_names",strings(0,1),"inductor_currents_A",zeros(0,1));end
function value=storage(action,k,value)
persistent entries;if isempty(entries),entries=containers.Map('KeyType','char','ValueType','any');end
switch action
 case "set",entries(k)=value;
 case "get",if ~isKey(entries,k),error("radia:simulink:LTspiceState","Block state was not initialized.");end;value=entries(k);
 case "remove",if isKey(entries,k),remove(entries,k);end;value=[];
end
end
function count=defaultKeptStepFolders(),count=3;end
