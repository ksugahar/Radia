function result=linearizeAndSolveHCurlSheetJouleLP(state,operator,resistance,cellVertexVelocities,frequenciesHz,rhs,designModeJacobian,area,options)
%LINEARIZEANDSOLVEHCURLSHEETJOULELP Close HCurl adjoints into the sheet LP.
arguments
    state (1,1) struct
    operator (1,1) radia.topopt.HCurlTopologyOperator
    resistance (1,1) struct
    cellVertexVelocities double {mustBeFinite}
    frequenciesHz (:,1) double {mustBePositive,mustBeFinite}
    rhs double
    designModeJacobian double {mustBeFinite}
    area (:,1) double {mustBePositive,mustBeFinite}
    options.VolumeMax (1,1) double {mustBePositive}
    options.DisplacementMove double {mustBePositive}
    options.ThicknessMove (1,1) double {mustBePositive}
    options.ActivationMove (1,1) double {mustBePositive}=0.2
    options.ThicknessBounds (1,2) double {mustBePositive}
    options.Weights (:,1) double {mustBeNonnegative,mustBeFinite}=ones(numel(frequenciesHz),1)
    options.RHSJacobian double=double.empty
    options.AdditionalSheetGradient double=double.empty
    options.Laplacian double=double.empty
    options.CurvatureLimit double=double.empty
end
linearization=radia.topopt.linearizeHCurlMultifrequencyJoule( ...
    operator,resistance,cellVertexVelocities,frequenciesHz,rhs, ...
    Weights=options.Weights,RHSJacobian=options.RHSJacobian);
n=numel(state.normal_displacement);
if ~isequal(size(designModeJacobian),[numel(linearization.gradient),3*n])
    error("radia:topopt:Shape", ...
        "designModeJacobian must have q rows and 3*n columns.");
end
gradient=designModeJacobian'*linearization.gradient(:);
if ~isempty(options.AdditionalSheetGradient)
    extra=options.AdditionalSheetGradient(:);
    if numel(extra)~=3*n
        error("radia:topopt:Shape", ...
            "AdditionalSheetGradient must have 3*n entries.");
    end
    gradient=gradient+extra;
end
update=radia.topopt.solveSheetMetalLP( ...
    state.normal_displacement(:),state.thickness(:),state.activation(:), ...
    gradient,area,VolumeMax=options.VolumeMax, ...
    DisplacementMove=options.DisplacementMove, ...
    ThicknessMove=options.ThicknessMove,ActivationMove=options.ActivationMove, ...
    ThicknessBounds=options.ThicknessBounds,Laplacian=options.Laplacian, ...
    CurvatureLimit=options.CurvatureLimit);
result=struct("schema","radia.hcurl.topopt.sheet-joule-lp/v1", ...
    "linearization",linearization,"update",update);
end
