function hystereticLTspiceSFunction(block)
%HYSTERETICLTSPICESFUNCTION Full circuit/hysteresis interval coupling block.
% Input 1 is [command; mechanical_gap_displacement_m]. GapPathFactor maps
% that displacement to total series air-gap length; it is not extra iron.
setup(block);
end
function setup(b)
b.NumDialogPrms=3;b.NumInputPorts=1;b.NumOutputPorts=1;b.SetPreCompInpPortInfoToDynamic;b.SetPreCompOutPortInfoToDynamic;
b.InputPort(1).Dimensions=2;b.InputPort(1).DirectFeedthrough=true;b.OutputPort(1).Dimensions=6;b.SampleTimes=[b.DialogPrm(2).Data 0];
b.SimStateCompliance='CustomSimState';
b.RegBlockMethod('PostPropagationSetup',@postSetup);b.RegBlockMethod('Start',@start);b.RegBlockMethod('Outputs',@outputs);b.RegBlockMethod('Update',@update);b.RegBlockMethod('GetSimState',@getSimState);b.RegBlockMethod('SetSimState',@setSimState);b.RegBlockMethod('Terminate',@terminate);
end
function postSetup(b),b.NumDworks=1;b.Dwork(1).Name='step';b.Dwork(1).Dimensions=1;b.Dwork(1).DatatypeID=0;b.Dwork(1).Complexity='Real';b.Dwork(1).UsedAsDiscState=true;end
function start(b)
loaded=load(string(b.DialogPrm(1).Data),"config");c=loaded.config;
if c.HysteresisKind~="play",error("radia:simulink:UnsupportedHysteresisKind","The LTspice coupling block currently supports only Play hysteresis materials.");end
material=[];
try
 material=radia.MatPlayHysteresis(c.K,c.EtaOrChi,c.Tables);root=string(b.DialogPrm(3).Data);if ~isfolder(root),mkdir(root);end
 keepCount=defaultKeptStepFolders();if isfield(c,"KeepStepArtifacts"),keepCount=c.KeepStepArtifacts;end
 logPath=radia.simulink.internal.runArtifacts("logpath",root,b.BlockHandle);
 entry=struct("config",c,"material",material,"hysteresis_state",radia.MatHysSaveState(material), ...
  "circuit_state",emptyState(),"previous_flux_Wb",0,"pending_hysteresis_state",[],"pending_circuit_state",[], ...
  "pending_flux_Wb",[],"pending_output",[],"has_pending",false,"folder",string(tempname(root)), ...
  "log_path",logPath,"kept",strings(0,1),"keep_count",keepCount,"failed",false);
 radia.simulink.internal.runArtifacts("begin",logPath,sprintf( ...
  'schema=radia.simulink.ltspice_hysteresis.run.v1 block=%s netlist=%s sample_time_s=%.17g started=%s artifacts=%s kept_steps=%d', ...
  getfullname(b.BlockHandle),c.Netlist,b.DialogPrm(2).Data,string(datetime("now","TimeZone","local"),"yyyy-MM-dd'T'HH:mm:ssXXX"),entry.folder,keepCount));
 mkdir(entry.folder);store("set",key(b),entry);b.Dwork(1).Data=0;
catch cause
 if ~isempty(material)
  try
   radia.UtiDel(material);
  catch
  end
 end
 rethrow(cause)
end
end
function outputs(b)
entry=store("get",key(b));c=entry.config;step=b.Dwork(1).Data+1;u=b.InputPort(1).Data(:);airGap=c.GapPathFactor*max(u(2),0);
if entry.has_pending,b.OutputPort(1).Data=entry.pending_output;return,end
try
 r=radia.simulink.runHystereticLTspiceInterval(c.Netlist,entry.material,entry.hysteresis_state,entry.circuit_state, ...
   CommandName=c.CommandName,CommandValue=u(1),BackEmfName=c.BackEmfName,CurrentTrace=c.CurrentTrace, ...
   Duration_s=b.DialogPrm(2).Data,Turns=c.Turns,CoreArea_m2=c.CoreArea_m2,MagneticPath_m=c.MagneticPath_m,AirGap_m=airGap,CoreVolume_m3=c.CoreVolume_m3, ...
   PreviousFlux_Wb=entry.previous_flux_Wb,OutputDirectory=fullfile(entry.folder,sprintf("step_%06d",step)), ...
   MaxIterations=c.MaxIterations,RelativeTolerance=c.RelativeTolerance,Relaxation=c.Relaxation,MaxStep_s=c.MaxStep_s,Timeout_s=c.Timeout_s,CouplingSamples=c.CouplingSamples);
 mu0=4*pi*1e-7;force=r.B_T(end)^2*c.CoreArea_m2/(2*mu0);
 y=[r.current_A(end);r.B_T(end);r.flux_Wb(end);r.back_emf_V(end);force;r.hysteresis_energy_J];b.OutputPort(1).Data=y;
 entry.pending_hysteresis_state=r.hysteresis_state;entry.pending_circuit_state=r.circuit_state;entry.pending_flux_Wb=r.flux_Wb(end);entry.pending_output=y;entry.has_pending=true;
 entry.kept(end+1,1)=string(r.output_directory);store("set",key(b),entry);
 radia.simulink.internal.runArtifacts("append",entry.log_path,sprintf( ...
  'step=%d t=%.17g command=%.17g position_m=%.17g iron_path_m=%.17g air_gap_m=%.17g iterations=%d residual=%.6g converged=%d current_A=%.17g B_T=%.17g flux_Wb=%.17g back_emf_V=%.17g force_N=%.17g energy_J=%.17g energy_balance_residual_J=%.17g artifacts=%s', ...
  step,b.CurrentTime,u(1),u(2),c.MagneticPath_m,airGap,r.iterations,r.relative_residual,r.converged,y(1),y(2),y(3),y(4),y(5),y(6),r.energy_balance_residual_J,r.output_directory));
catch cause
 entry.failed=true;store("set",key(b),entry);
 radia.simulink.internal.runArtifacts("append",entry.log_path,sprintf( ...
  'step=%d t=%.17g command=%.17g position=%.17g FAILED identifier=%s message=%s',step,b.CurrentTime,u(1),u(2),cause.identifier,cause.message));
 error("radia:simulink:HystereticLTspiceStepFailed","Circuit/hysteresis step %d failed: %s",step,cause.message);
end
end
function update(b)
entry=store("get",key(b));if ~entry.has_pending,return,end
entry.hysteresis_state=entry.pending_hysteresis_state;entry.circuit_state=entry.pending_circuit_state;entry.previous_flux_Wb=entry.pending_flux_Wb;
entry.pending_hysteresis_state=[];entry.pending_circuit_state=[];entry.pending_flux_Wb=[];entry.pending_output=[];entry.has_pending=false;
entry.kept=radia.simulink.internal.runArtifacts("prune",entry.kept,entry.keep_count);
store("set",key(b),entry);b.Dwork(1).Data=b.Dwork(1).Data+1;
end
function state=getSimState(b)
entry=store("get",key(b));state=struct("step",b.Dwork(1).Data,"hysteresis_state",entry.hysteresis_state,"circuit_state",entry.circuit_state,"previous_flux_Wb",entry.previous_flux_Wb, ...
 "pending_hysteresis_state",entry.pending_hysteresis_state,"pending_circuit_state",entry.pending_circuit_state,"pending_flux_Wb",entry.pending_flux_Wb,"pending_output",entry.pending_output,"has_pending",entry.has_pending);
end
function setSimState(b,state)
required=["step","hysteresis_state","circuit_state","previous_flux_Wb","pending_hysteresis_state","pending_circuit_state","pending_flux_Wb","pending_output","has_pending"];
if ~isstruct(state)||~all(isfield(state,required)),error("radia:simulink:HystereticLTspiceSimState","Invalid hysteretic LTspice block SimState.");end
entry=store("get",key(b));for name=required(2:end),entry.(name)=state.(name);end;radia.MatHysRestoreState(entry.material,entry.hysteresis_state);store("set",key(b),entry);b.Dwork(1).Data=state.step;
end
function terminate(b)
k=key(b);folder="";failed=false;logPath="";
try
 entry=store("get",k);folder=entry.folder;failed=entry.failed;logPath=entry.log_path;radia.UtiDel(entry.material);
catch
end
store("remove",k);
if strlength(logPath)>0
 radia.simulink.internal.runArtifacts("append",logPath,sprintf('finished failed=%d artifacts_retained=%d',failed,failed));
end
radia.simulink.internal.runArtifacts("finish",folder,failed);
end
function count=defaultKeptStepFolders(),count=3;end
function state=emptyState(),state=struct("schema","radia.ltspice.transient_state.v1","time_s",0,"node_names",strings(0,1),"node_voltages_V",zeros(0,1),"inductor_names",strings(0,1),"inductor_currents_A",zeros(0,1));end
function answer=key(b),answer=sprintf('%.0f',b.BlockHandle);end
function value=store(action,k,value)
persistent entries;if isempty(entries),entries=containers.Map('KeyType','char','ValueType','any');end
switch action,case "set",entries(k)=value;case "get",if ~isKey(entries,k),error("radia:simulink:HystereticState","Block state is unavailable.");end;value=entries(k);case "remove",if isKey(entries,k),remove(entries,k);end;value=[];end
end
