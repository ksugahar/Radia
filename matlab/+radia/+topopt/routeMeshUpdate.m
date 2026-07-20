function result=routeMeshUpdate(jacobianDeterminants,jacobianConditions,relativeDisplacements,options)
%ROUTEMESHUPDATE Route to NGSolve deformation/refinement or Cubit rebuild.
arguments
 jacobianDeterminants (:,1) double {mustBeFinite}
 jacobianConditions (:,1) double {mustBeFinite}
 relativeDisplacements (:,1) double {mustBeFinite}
 options.RefineThreshold (1,1) double {mustBePositive}=0.25
 options.RebuildThreshold (1,1) double {mustBePositive}=0.5
 options.MinimumJacobian (1,1) double=0.2
 options.MaximumCondition (1,1) double {mustBePositive}=20
 options.TopologyChanged (1,1) logical=false
end
n=numel(jacobianDeterminants);
if n==0||numel(jacobianConditions)~=n||numel(relativeDisplacements)~=n, error("radia:topopt:Shape","Quality vectors must be non-empty and equal length."); end
bad=jacobianDeterminants<=options.MinimumJacobian|jacobianConditions>=options.MaximumCondition|relativeDisplacements>=options.RefineThreshold;
reasons=strings(0,1);
if options.TopologyChanged, reasons(end+1)="material topology changed"; end
if any(jacobianDeterminants<=0), reasons(end+1)="inverted element"; end
if any(relativeDisplacements>=options.RebuildThreshold), reasons(end+1)="deformation exceeded rebuild threshold"; end
if any(jacobianConditions>=2*options.MaximumCondition), reasons(end+1)="severe Jacobian distortion"; end
if ~isempty(reasons), route="cubit_rebuild";
elseif any(bad), route="ngsolve_refine"; reasons(end+1)="local deformation quality threshold exceeded";
else, route="ngsolve_deform"; reasons(end+1)="deformation remains inside quality limits";
end
result=struct("schema","radia.topopt.mesh-update-route/v1","route",route,"refine_elements",find(bad)-1, ...
 "reasons",reasons,"minimum_jacobian",min(jacobianDeterminants),"maximum_condition",max(jacobianConditions), ...
 "maximum_relative_displacement",max(relativeDisplacements));
end
