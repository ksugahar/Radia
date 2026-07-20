function result=linearizeAndSolveHCurlActivationSheetJouleLP(state,operator,cellCurlGrams,frequenciesHz,rhs,conductivity,area,options)
%LINEARIZEANDSOLVEHCURLACTIVATIONSHEETJOULELP Activation-only HCurl LP step.
arguments
    state (1,1) struct
    operator (1,1) radia.topopt.HCurlTopologyOperator
    cellCurlGrams double {mustBeFinite}
    frequenciesHz (:,1) double {mustBePositive,mustBeFinite}
    rhs double
    conductivity (1,1) struct
    area (:,1) double {mustBePositive,mustBeFinite}
    options.VolumeMax (1,1) double {mustBePositive}
    options.ActivationMove (1,1) double {mustBePositive}=0.2
    options.Weights (:,1) double {mustBeNonnegative,mustBeFinite}=ones(numel(frequenciesHz),1)
    options.RHSJacobian double=double.empty
    options.InductancePower (1,1) double {mustBeGreaterThanOrEqual(options.InductancePower,1)}=1
end
linearization=radia.topopt.linearizeHCurlActivationMultifrequencyJoule( ...
    operator,cellCurlGrams,state.activation(:),frequenciesHz,rhs,conductivity, ...
    Weights=options.Weights,RHSJacobian=options.RHSJacobian, ...
    InductancePower=options.InductancePower);
n=numel(state.activation);
if numel(linearization.gradient)~=n
    error("radia:topopt:Shape", ...
        "Activation gradient and sheet-cell count differ.");
end
gradient=[zeros(2*n,1);linearization.gradient(:)];
frozenMove=eps;
update=radia.topopt.solveSheetMetalLP( ...
    state.normal_displacement(:),state.thickness(:),state.activation(:), ...
    gradient,area,VolumeMax=options.VolumeMax, ...
    DisplacementMove=frozenMove,ThicknessMove=frozenMove, ...
    ActivationMove=options.ActivationMove, ...
    ThicknessBounds=[min(state.thickness),max(state.thickness)]);
result=struct("schema","radia.hcurl.topopt.activation-sheet-joule-lp/v1", ...
    "linearization",linearization,"update",update);
end
