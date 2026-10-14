function result=solveLPUpdate(density,objectiveGradient,cellVolumes,volumeMax,options)
%SOLVELPUPDATE Solve a volume-constrained material update using LINPROG.
arguments
 density (:,1) double {mustBeBetween(density,0,1)}
 objectiveGradient (:,1) double {mustBeFinite}
 cellVolumes (:,1) double {mustBePositive}
 volumeMax (1,1) double {mustBePositive}
 options.MoveLimit (1,1) double {mustBePositive}=0.2
 options.Aineq double=double.empty; options.bineq double=double.empty
end
n=numel(density); if numel(objectiveGradient)~=n||numel(cellVolumes)~=n, error("radia:topopt:Shape","Cell vectors must have equal length."); end
if options.MoveLimit>1, error("radia:topopt:MoveLimit","MoveLimit must not exceed one."); end
Aineq=cellVolumes'; bineq=volumeMax;
if ~isempty(options.Aineq), if size(options.Aineq,2)~=n||size(options.Aineq,1)~=numel(options.bineq), error("radia:topopt:Shape","Aineq/bineq mismatch."); end, Aineq=[Aineq;options.Aineq]; bineq=[bineq;options.bineq(:)]; end
lower=max(0,density-options.MoveLimit); upper=min(1,density+options.MoveLimit);
settings=optimoptions("linprog","Display","none","Algorithm","dual-simplex-highs");
[next,value,exitflag,output]=linprog(objectiveGradient,Aineq,bineq,[],[],lower,upper,settings);
if exitflag<=0, error("radia:topopt:LPFailed","Topology LP failed: %s",output.message); end
result=struct("schema","radia.topopt.lp-update/v1","density",next,"delta",next-density, ...
 "predicted_objective",value,"exitflag",exitflag,"iterations",output.iterations);
end
