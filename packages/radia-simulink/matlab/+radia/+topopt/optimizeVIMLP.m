function result=optimizeVIMLP(initialDensity,cellVolumes,volumeFraction,linearizeFcn,options)
%OPTIMIZEVIMLP Sequential Radia-VIM linearization and LP material updates.
arguments
 initialDensity (:,1) double {mustBeBetween(initialDensity,0,1)}
 cellVolumes (:,1) double {mustBePositive}
 volumeFraction (1,1) double {mustBePositive}
 linearizeFcn (1,1) function_handle
 options.ObjectiveWeights (:,1) double
 options.MoveLimit (1,1) double {mustBePositive}=0.2
 options.MaxIterations (1,1) double {mustBeInteger,mustBePositive}=30
 options.DensityTolerance (1,1) double {mustBePositive}=1e-3
 options.ProgressFcn=[]
end
if volumeFraction>1, error("radia:topopt:VolumeFraction","VolumeFraction must not exceed one."); end
if numel(initialDensity)~=numel(cellVolumes), error("radia:topopt:Shape","Density and cell volumes must match."); end
density=initialDensity; history=table('Size',[0,4],'VariableTypes',{'double','double','double','double'}, ...
 'VariableNames',{'Iteration','Objective','Volume','MaxDensityChange'}); converged=false;
for iteration=0:options.MaxIterations-1
 model=linearizeFcn(density); jacobian=model.response_jacobian; response=model.response(:);
 if ~isequal(size(jacobian),[numel(options.ObjectiveWeights),numel(density)]), error("radia:topopt:Shape","Linearized response shape mismatch."); end
 gradient=real(options.ObjectiveWeights(:)'*jacobian)';
 update=radia.topopt.solveLPUpdate(density,gradient,cellVolumes,volumeFraction*sum(cellVolumes),MoveLimit=options.MoveLimit);
 change=max(abs(update.delta)); objective=real(options.ObjectiveWeights(:)'*response);
 history(end+1,:)={iteration,objective,cellVolumes'*update.density,change}; %#ok<AGROW>
 density=update.density;
 if ~isempty(options.ProgressFcn), options.ProgressFcn(struct("iteration",iteration,"density",density,"objective",objective,"change",change)); end
 if change<=options.DensityTolerance, converged=true; break, end
end
result=struct("schema","radia.topopt.vim-lp/v1","density",density,"history",history,"converged",converged);
end
