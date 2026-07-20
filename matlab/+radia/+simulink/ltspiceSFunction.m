function ltspiceSFunction(block)
%LTSPICESFUNCTION Level-2 block for sampled LTspice closed-loop coupling.
setup(block);
end
function setup(b)
b.NumDialogPrms=8; b.NumInputPorts=1; b.NumOutputPorts=1;
b.SetPreCompInpPortInfoToDynamic; b.SetPreCompOutPortInfoToDynamic;
b.InputPort(1).Dimensions=numel(names(b.DialogPrm(2).Data)); b.InputPort(1).DirectFeedthrough=true;
b.OutputPort(1).Dimensions=numel(names(b.DialogPrm(3).Data)); b.SampleTimes=[b.DialogPrm(4).Data 0];
b.RegBlockMethod('PostPropagationSetup',@postSetup); b.RegBlockMethod('Start',@start); b.RegBlockMethod('Outputs',@outputs); b.RegBlockMethod('Terminate',@terminate);
end
function postSetup(b)
b.NumDworks=1; b.Dwork(1).Name='step'; b.Dwork(1).Dimensions=1; b.Dwork(1).DatatypeID=0; b.Dwork(1).Complexity='Real'; b.Dwork(1).UsedAsDiscState=true;
end
function start(b)
b.Dwork(1).Data=0; root=string(b.DialogPrm(5).Data);if ~isfolder(root),mkdir(root);end
entry=struct("state",emptyState(),"folder",string(tempname(root)));
mkdir(entry.folder); storage("set",key(b),entry);
end
function outputs(b)
k=b.Dwork(1).Data+1; entry=storage("get",key(b)); inames=names(b.DialogPrm(2).Data); traces=names(b.DialogPrm(3).Data);
stepFolder=fullfile(entry.folder,sprintf("step_%06d",k));
try
 r=radia.simulink.runLTspiceInterval(string(b.DialogPrm(1).Data),inames,b.InputPort(1).Data(:),entry.state, ...
   Duration_s=b.DialogPrm(4).Data,OutputDirectory=stepFolder,MaxStep_s=b.DialogPrm(6).Data, ...
   Timeout_s=b.DialogPrm(7).Data,Executable=string(b.DialogPrm(8).Data));
 y=zeros(numel(traces),1); for n=1:numel(traces),j=find(r.simulation.waveform.names==traces(n),1);if isempty(j),error("radia:simulink:LTspiceTrace","Trace not found: %s",traces(n));end;y(n)=real(r.simulation.waveform.values(end,j));end
 b.OutputPort(1).Data=y; entry.state=r.state;storage("set",key(b),entry);b.Dwork(1).Data=k;
catch cause
 error("radia:simulink:LTspiceStepFailed","LTspice step %d at Simulink time %.17g failed: %s",k,b.CurrentTime,cause.message);
end
end
function terminate(b),storage("remove",key(b));end
function answer=names(value)
if iscell(value),answer=string(value(:));else,answer=string(value(:));end
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
