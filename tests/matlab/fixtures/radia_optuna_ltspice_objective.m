function value=radia_optuna_ltspice_objective(trial)
%RADIA_OPTUNA_LTSPICE_OBJECTIVE Exercise the complete MATLAB/LTspice path.
folder=fileparts(mfilename("fullpath")); netlist=fullfile(folder,"ltspice_rc.cir");
r=trial.suggestFloat("Rval",500,2500); result=radia.ltspice.run(netlist,Parameters=struct("Rval",r));
t=real(result.waveform.values(:,1)); j=find(result.waveform.names=="V(out)",1); [~,k]=min(abs(t-1e-3));
value=abs(result.waveform.values(k,j)-0.5);
end
