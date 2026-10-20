function result=backtrackTrafoDeformation(mesh,deformationFactory,currentNormal,targetNormal,relativeDisplacements,options)
%BACKTRACKTRAFODEFORMATION Backtrack an absolute VectorH1 mesh deformation.
%   NGSolve owns the active deformation and GetTrafo quality evaluation.
arguments
    mesh (1,1) radia.ngsolve.Mesh
    deformationFactory (1,1) function_handle
    currentNormal (:,1) double {mustBeFinite}
    targetNormal (:,1) double {mustBeFinite}
    relativeDisplacements (:,1) double {mustBeNonnegative,mustBeFinite}
    options.MinimumScale (1,1) double {mustBePositive}=1/64
    options.Contraction (1,1) double {mustBePositive}=0.5
    options.MinimumJacobianRatio (1,1) double {mustBePositive}=0.2
    options.MaximumCondition (1,1) double {mustBePositive}=20
    options.RefineThreshold (1,1) double {mustBePositive}=0.25
    options.RebuildThreshold (1,1) double {mustBePositive}=0.5
    options.IntegrationOrder (1,1) double {mustBeInteger,mustBePositive}=2
    options.TopologyChanged (1,1) logical=false
end
if options.MinimumScale>1||options.Contraction>=1
    error("radia:topopt:Backtracking","Invalid backtracking parameters.");
end
if ~isequal(size(currentNormal),size(targetNormal))
    error("radia:topopt:Shape","Current and target displacement shapes must match.");
end

mesh.unsetDeformation();
reference=mesh.trafoQuality(IntegrationOrder=options.IntegrationOrder);
scale=1; trials=0; accepted=false; decision=struct.empty;
quality=struct.empty; deformation=[]; candidate=currentNormal;
while scale>=options.MinimumScale
    trials=trials+1;
    candidate=currentNormal+scale*(targetNormal-currentNormal);
    deformation=deformationFactory(mesh,candidate);
    if ~isa(deformation,"radia.ngsolve.GridFunction")||deformation.Space~="vectorh1"
        error("radia:topopt:Deformation", ...
            "DeformationFactory must return a VectorH1 radia.ngsolve.GridFunction.");
    end
    mesh.setDeformation(deformation);
    quality=mesh.trafoQuality(IntegrationOrder=options.IntegrationOrder, ...
        ReferenceDeterminants=reference.raw_jacobian_determinants);
    decision=radia.topopt.routeMeshUpdate(quality.jacobian_determinants, ...
        quality.jacobian_conditions,scale*relativeDisplacements, ...
        MinimumJacobian=options.MinimumJacobianRatio, ...
        MaximumCondition=options.MaximumCondition, ...
        RefineThreshold=options.RefineThreshold, ...
        RebuildThreshold=options.RebuildThreshold, ...
        TopologyChanged=options.TopologyChanged);
    unsafe=any(quality.jacobian_determinants<=options.MinimumJacobianRatio)|| ...
        any(quality.jacobian_conditions>=2*options.MaximumCondition);
    if ~unsafe
        accepted=true;
        if decision.route~="ngsolve_deform"
            mesh.unsetDeformation();
        end
        break
    end
    mesh.unsetDeformation();
    deformation=[];
    scale=scale*options.Contraction;
end
if ~accepted
    scale=0;
    mesh.unsetDeformation();
end
result=struct("schema","radia.topopt.trafo-deformation/v1", ...
    "accepted",accepted,"scale",scale,"decision",decision, ...
    "quality",quality,"trials",trials,"candidate",candidate, ...
    "deformation",deformation);
end
