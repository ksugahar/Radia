function result=magneticShieldLPStep(normalDisplacement,thickness,activation,cellAreas,parallelField,analyticResponseJacobian,options)
%MAGNETICSHIELDLPSTEP One finite-difference-free thin-shield LP update.
arguments
 normalDisplacement (:,1) double
 thickness (:,1) double {mustBePositive}
 activation (:,1) double {mustBeBetween(activation,0,1)}
 cellAreas (:,1) double {mustBePositive}
 parallelField (:,1) double
 analyticResponseJacobian (:,:) double
 options.VolumeMax (1,1) double {mustBePositive}
 options.DisplacementMove double {mustBePositive}
 options.ThicknessMove (1,1) double {mustBePositive}
 options.ThicknessBounds (1,2) double {mustBePositive}
 options.Laplacian double=double.empty
 options.CurvatureLimit double=double.empty
end
gradient=radia.topopt.magneticShieldRMSGradient(parallelField,analyticResponseJacobian);
result=radia.topopt.solveSheetMetalLP(normalDisplacement,thickness,activation,gradient,cellAreas, ...
 VolumeMax=options.VolumeMax,DisplacementMove=options.DisplacementMove,ThicknessMove=options.ThicknessMove, ...
 ThicknessBounds=options.ThicknessBounds,Laplacian=options.Laplacian,CurvatureLimit=options.CurvatureLimit);
end
