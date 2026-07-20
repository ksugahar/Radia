function result=runHystereticLTspiceInterval(netlistFile,material,hysteresisState,circuitState,options)
%RUNHYSTERETICLTSPICEINTERVAL Waveform-relaxed LTspice/hysteresis coupling.
arguments
 netlistFile (1,1) string {mustBeFile}; material (1,1) double {mustBePositive}
 hysteresisState (:,1) double; circuitState (1,1) struct
 options.CommandName (1,1) string="command"; options.CommandValue (1,1) double
 options.BackEmfName (1,1) string="back_emf"; options.CurrentTrace (1,1) string
 options.Duration_s (1,1) double {mustBePositive}; options.Turns (1,1) double {mustBePositive}
 options.CoreArea_m2 (1,1) double {mustBePositive}; options.MagneticPath_m (1,1) double {mustBePositive}
 options.CoreVolume_m3 (1,1) double {mustBePositive}; options.PreviousFlux_Wb (1,1) double=0
 options.OutputDirectory (1,1) string=""; options.MaxIterations (1,1) double {mustBeInteger,mustBePositive}=12
 options.RelativeTolerance (1,1) double {mustBePositive}=1e-4; options.Relaxation (1,1) double {mustBeGreaterThan(options.Relaxation,0),mustBeLessThanOrEqual(options.Relaxation,1)}=0.5
 options.MaxStep_s (1,1) double {mustBePositive}=inf; options.Timeout_s (1,1) double {mustBePositive}=300
 options.Executable (1,1) string=""; options.MaxAbsB_T (1,1) double {mustBePositive}=5
 options.CouplingSamples (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(options.CouplingSamples,3)}=101
end
root=options.OutputDirectory;if strlength(root)==0,root=string(tempname("C:\temp"));end;if ~isfolder(root),mkdir(root);end
saved=radia.MatHysSaveState(material);cleanup=onCleanup(@()radia.MatHysRestoreState(material,saved));
tOld=[0;options.Duration_s];eOld=zeros(2,1);converged=false;history=zeros(options.MaxIterations,1);
for iteration=1:options.MaxIterations
 folder=fullfile(root,sprintf("iteration_%03d",iteration));if ~isfolder(folder),mkdir(folder);end
 stateNetlist=fullfile(folder,"coupled_state.cir");radia.ltspice.applyTransientState(netlistFile,circuitState,stateNetlist,Duration_s=options.Duration_s,MaxStep_s=options.MaxStep_s);
 signals=struct();signals.(char(options.CommandName))=[0,options.CommandValue;options.Duration_s,options.CommandValue];signals.(char(options.BackEmfName))=[tOld,eOld];
 simulation=radia.simulink.runLTspice(stateNetlist,InputSignals=signals,Executable=options.Executable,OutputDirectory=folder,Timeout_s=options.Timeout_s);
 names=simulation.waveform.names;j=find(names==options.CurrentTrace,1);if isempty(j),error("radia:simulink:HystereticCurrentTrace","Current trace not found: %s",options.CurrentTrace);end
 rawTime=real(simulation.waveform.values(:,1));rawCurrent=real(simulation.waveform.values(:,j));
 t=linspace(0,options.Duration_s,options.CouplingSamples).';current=interp1(rawTime,rawCurrent,t,"linear","extrap");
 [B,H,states,flux,emf,energy]=hysteresisWaveform(material,hysteresisState,current,t,options);
 previous=interp1(tOld,eOld,t,"linear","extrap");scale=max([max(abs(emf)),max(abs(previous)),1e-12]);history(iteration)=max(abs(emf-previous))/scale;
 if history(iteration)<=options.RelativeTolerance,converged=true;break,end
 eOld=(1-options.Relaxation)*previous+options.Relaxation*emf;tOld=t;
end
if ~converged,error("radia:simulink:HystereticCouplingNotConverged","LTspice/hysteresis waveform iteration did not converge in %d iterations; residual=%g.",options.MaxIterations,history(options.MaxIterations));end
nextCircuit=radia.ltspice.extractTransientState(simulation);nextCircuit.time_s=circuitState.time_s+options.Duration_s;
result=struct("schema","radia.simulink.ltspice_hysteresis.interval.v1","simulation",simulation, ...
 "circuit_state",nextCircuit,"hysteresis_state",states(end,:).',"time_s",t,"current_A",current, ...
 "B_T",B,"H_A_per_m",H,"flux_Wb",flux,"back_emf_V",emf,"hysteresis_energy_J",energy, ...
 "iterations",iteration,"relative_residual",history(iteration),"converged",converged,"output_directory",root);
clear cleanup;radia.MatHysRestoreState(material,saved);
end

function [B,H,states,flux,emf,energy]=hysteresisWaveform(material,state0,current,t,o)
n=numel(t);B=zeros(n,1);H=zeros(n,1);states=zeros(n,numel(state0));state=state0(:).';guess=0;energy=0;
for k=1:n
 target=o.Turns*current(k)/o.MagneticPath_m;
 fun=@(b)localH(material,b,state)-target;lo=-o.MaxAbsB_T;hi=o.MaxAbsB_T;
 if fun(lo)>0||fun(hi)<0,error("radia:simulink:HystereticFieldRange","Required field %.6g A/m exceeds the B search range.",target);end
 if abs(fun(guess))<1e-12,b=guess;else,b=fzero(fun,[lo hi]);end
 h=localH(material,b,state);next=radia.MatHysCommitBatch(material,[b,0,0],state);
 B(k)=b;H(k)=h;states(k,:)=next; if k>1,energy=energy+0.5*(H(k-1)+H(k))*(B(k)-B(k-1))*o.CoreVolume_m3;end
 state=next;guess=b;
end
flux=o.Turns*o.CoreArea_m2*B;emf=zeros(n,1);
if n>1,dt=max(t(2)-t(1),eps);else,dt=o.Duration_s;end
emf(1)=(flux(1)-o.PreviousFlux_Wb)/dt;if n>1,emf(2:end)=diff(flux)./max(diff(t),eps);end
end
function h=localH(material,b,state)
value=radia.MatHysForwardBatch(material,[b,0,0],state);h=value(1);
end
