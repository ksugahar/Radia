function result=acceptTrafoStep(qualityFcn,relativeDisplacements,options)
%ACCEPTTRAFOSTEP Backtrack a candidate until its NGSolve Trafo quality is safe.
arguments
 qualityFcn (1,1) function_handle
 relativeDisplacements (:,1) double {mustBeNonnegative}
 options.MinimumScale (1,1) double {mustBePositive}=1/64
 options.Contraction (1,1) double {mustBePositive}=0.5
 options.MinimumJacobianRatio (1,1) double=0.2
 options.MaximumCondition (1,1) double {mustBePositive}=20
 options.RefineThreshold (1,1) double {mustBePositive}=0.25
 options.RebuildThreshold (1,1) double {mustBePositive}=0.5
end
if options.MinimumScale>1||options.Contraction>=1, error("radia:topopt:Backtracking","Invalid backtracking parameters."); end
scale=1; trials=0; accepted=false; route=struct.empty;
while scale>=options.MinimumScale
 trials=trials+1; quality=qualityFcn(scale);
 route=radia.topopt.routeMeshUpdate(quality.jacobian_determinants,quality.jacobian_conditions, ...
  scale*relativeDisplacements,MinimumJacobian=options.MinimumJacobianRatio,MaximumCondition=options.MaximumCondition, ...
  RefineThreshold=options.RefineThreshold,RebuildThreshold=options.RebuildThreshold);
 unsafe=any(quality.jacobian_determinants<=options.MinimumJacobianRatio)|| ...
  any(quality.jacobian_conditions>=2*options.MaximumCondition);
 if ~unsafe, accepted=true; break, end
 scale=scale*options.Contraction;
end
if ~accepted, scale=0; end
result=struct("schema","radia.topopt.trafo-acceptance/v1","accepted",accepted,"scale",scale,"route",route,"trials",trials);
end
