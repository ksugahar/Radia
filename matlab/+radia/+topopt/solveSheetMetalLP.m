function result=solveSheetMetalLP(normalDisplacement,thickness,activation,objectiveGradient,cellAreas,options)
%SOLVESHEETMETALLP Manufacturing-aware LP update for a sheet-metal design.
arguments
 normalDisplacement (:,1) double {mustBeFinite}
 thickness (:,1) double {mustBePositive}
 activation (:,1) double {mustBeBetween(activation,0,1)}
 objectiveGradient (:,1) double {mustBeFinite}
 cellAreas (:,1) double {mustBePositive}
 options.VolumeMax (1,1) double {mustBePositive}
 options.DisplacementMove double {mustBePositive}
 options.ThicknessMove (1,1) double {mustBePositive}
 options.ActivationMove (1,1) double {mustBePositive}=0.2
 options.ThicknessBounds (1,2) double {mustBePositive}
 options.Laplacian double=double.empty
 options.CurvatureLimit double=double.empty
 options.Aineq double=double.empty
 options.bineq double=double.empty
end
n=numel(normalDisplacement);
if n==0||numel(thickness)~=n||numel(activation)~=n||numel(cellAreas)~=n||numel(objectiveGradient)~=3*n
 error("radia:topopt:Shape","Sheet vectors must match and objectiveGradient must contain 3*n values.");
end
tmin=options.ThicknessBounds(1); tmax=options.ThicknessBounds(2);
if tmin>tmax, error("radia:topopt:ThicknessBounds","Thickness bounds are reversed."); end
move=options.DisplacementMove;
if isscalar(move), move=repmat(move,n,1); else, move=move(:); end
if numel(move)~=n, error("radia:topopt:Shape","DisplacementMove must be scalar or n-by-1."); end
current=[normalDisplacement;thickness;activation];
lower=[normalDisplacement-move;max(tmin,thickness-options.ThicknessMove);max(0,activation-options.ActivationMove)];
upper=[normalDisplacement+move;min(tmax,thickness+options.ThicknessMove);min(1,activation+options.ActivationMove)];
volume0=sum(cellAreas.*thickness.*activation);
Aineq=[zeros(1,n),(cellAreas.*activation)',(cellAreas.*thickness)']; bineq=options.VolumeMax+volume0;
if ~isempty(options.Laplacian)
 L=options.Laplacian;
 if size(L,2)~=n||isempty(options.CurvatureLimit), error("radia:topopt:Curvature","Laplacian requires n columns and CurvatureLimit."); end
 limit=options.CurvatureLimit;
 if isscalar(limit), limit=repmat(limit,size(L,1),1); else, limit=limit(:); end
 if numel(limit)~=size(L,1), error("radia:topopt:Curvature","CurvatureLimit size mismatch."); end
 pad=zeros(size(L,1),2*n); Aineq=[Aineq;L,pad;-L,pad]; bineq=[bineq;limit;limit];
end
if ~isempty(options.Aineq)
 if size(options.Aineq,2)~=3*n||size(options.Aineq,1)~=numel(options.bineq), error("radia:topopt:Shape","Aineq/bineq mismatch."); end
 Aineq=[Aineq;options.Aineq]; bineq=[bineq;options.bineq(:)];
end
settings=optimoptions("linprog","Display","none","Algorithm","dual-simplex-highs");
[next,value,exitflag,output]=linprog(objectiveGradient,Aineq,bineq,[],[],lower,upper,settings);
if exitflag<=0, error("radia:topopt:LPFailed","Sheet-metal LP failed: %s",output.message); end
result=struct("schema","radia.topopt.sheet-metal-lp/v1","normal_displacement",next(1:n), ...
 "thickness",next(n+1:2*n),"activation",next(2*n+1:end),"delta",next-current, ...
 "predicted_objective",value,"exitflag",exitflag,"iterations",output.iterations);
end
