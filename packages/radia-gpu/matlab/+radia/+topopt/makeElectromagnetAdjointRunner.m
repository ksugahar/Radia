function runner = makeElectromagnetAdjointRunner( ...
        baseMatrix, sourceVector, responseMatrix, cellMatrix, target, options)
%MAKEELECTROMAGNETADJOINTRUNNER Build a VIM density-topology optimizer.
%   The state equation is A(rho)*m=b(rho). One state solve and one adjoint
%   solve provide the complete objective gradient for MMA or SQP.
arguments
    baseMatrix (:,:) double {mustBeReal,mustBeFinite}
    sourceVector (:,1) double {mustBeReal,mustBeFinite}
    responseMatrix (:,:) double {mustBeReal,mustBeFinite}
    cellMatrix double {mustBeReal,mustBeFinite}
    target (:,1) double {mustBeReal,mustBeFinite}
    options.CellSource double {mustBeReal,mustBeFinite} = double.empty
    options.CellResponse double {mustBeReal,mustBeFinite} = double.empty
    options.ResponseWeights double {mustBeReal,mustBeFinite} = ...
        ones(size(responseMatrix,1),1)
    options.ElementVolumes double {mustBeReal,mustBeFinite} = ...
        ones(size(cellMatrix,1),1)
    options.VolumeFraction (1,1) double {mustBePositive,mustBeFinite} = 0.5
    options.PenalizationPower (1,1) double ...
        {mustBeGreaterThanOrEqual(options.PenalizationPower,1),mustBeFinite} = 1
    options.DensityRegularization (1,1) double ...
        {mustBeNonnegative,mustBeFinite} = 0
    options.InitialDesign double {mustBeReal,mustBeFinite} = double.empty
    options.LowerBounds double {mustBeReal,mustBeFinite} = double.empty
    options.UpperBounds double {mustBeReal,mustBeFinite} = double.empty
    options.MoveLimit double {mustBeReal,mustBeFinite} = 0.15
    options.MaxIterations (1,1) double {mustBeInteger,mustBePositive} = 60
    options.Solver (1,1) string ...
        {mustBeMember(options.Solver,["mma","sqp"])} = "mma"
    options.OptimizerOptions (1,1) struct = struct
end

stateCount = size(baseMatrix,1);
if size(baseMatrix,2) ~= stateCount || numel(sourceVector) ~= stateCount
    error("radia:topopt:ElectromagnetShape", ...
        "baseMatrix must be square and sourceVector must match it.");
end
responseCount = size(responseMatrix,1);
if size(responseMatrix,2) ~= stateCount || numel(target) ~= responseCount
    error("radia:topopt:ElectromagnetShape", ...
        "responseMatrix and target must match the state equation.");
end
if ndims(cellMatrix) ~= 3 || size(cellMatrix,2) ~= stateCount || ...
        size(cellMatrix,3) ~= stateCount || size(cellMatrix,1) < 1
    error("radia:topopt:ElectromagnetShape", ...
        "cellMatrix must be design_count-by-state_count-by-state_count.");
end
designCount = size(cellMatrix,1);
cellSource = double(options.CellSource);
if isempty(cellSource)
    cellSource = zeros(designCount,stateCount);
end
cellResponse = double(options.CellResponse);
if isempty(cellResponse)
    cellResponse = zeros(designCount,responseCount,stateCount);
end
if ~isequal(size(cellSource),[designCount,stateCount]) || ...
        ~isequal(size(cellResponse),[designCount,responseCount,stateCount])
    error("radia:topopt:ElectromagnetShape", ...
        "CellSource or CellResponse has an incompatible shape.");
end
weights = localExpand(options.ResponseWeights,responseCount, ...
    "ResponseWeights");
volumes = localExpand(options.ElementVolumes,designCount,"ElementVolumes");
if any(weights <= 0) || any(volumes <= 0)
    error("radia:topopt:ElectromagnetWeights", ...
        "ResponseWeights and ElementVolumes must be positive.");
end
if options.VolumeFraction > 1
    error("radia:topopt:ElectromagnetVolumeFraction", ...
        "VolumeFraction must not exceed one.");
end

lower = localDefaultExpand(options.LowerBounds,zeros(designCount,1), ...
    designCount,"LowerBounds");
upper = localDefaultExpand(options.UpperBounds,ones(designCount,1), ...
    designCount,"UpperBounds");
initial = localDefaultExpand(options.InitialDesign, ...
    options.VolumeFraction*ones(designCount,1),designCount,"InitialDesign");
if any(lower < 0) || any(upper > 1) || any(lower >= upper) || ...
        any(initial < lower) || any(initial > upper)
    error("radia:topopt:ElectromagnetBounds", ...
        "Density bounds must lie in [0,1] and bracket InitialDesign.");
end
moveLimit = localExpand(options.MoveLimit,designCount,"MoveLimit");
if any(moveLimit <= 0)
    error("radia:topopt:ElectromagnetMoveLimit", ...
        "MoveLimit must be positive.");
end

optimizerOptions = options.OptimizerOptions;
optimizerOptions.LowerBounds = lower;
optimizerOptions.UpperBounds = upper;
optimizerOptions.MoveLimit = moveLimit;
optimizerOptions.MaxIterations = options.MaxIterations;
metadata = struct( ...
    "domain","electromagnet-topology", ...
    "formulation","vim-density-adjoint", ...
    "state_equation","A(rho)*m=b(rho)", ...
    "gradient","one-state-one-adjoint", ...
    "design_count",designCount, ...
    "response_count",responseCount, ...
    "volume_fraction_limit",options.VolumeFraction, ...
    "element_volumes",volumes, ...
    "penalization_power",options.PenalizationPower, ...
    "python_per_step",false);
runner = radia.topopt.AdjointRunner(initial,@evaluate, ...
    Solver=options.Solver,OptimizerOptions=optimizerOptions,Metadata=metadata);

    function evaluation = evaluate(density)
        interpolation = density.^options.PenalizationPower;
        interpolationDerivative = options.PenalizationPower * ...
            density.^(options.PenalizationPower-1);
        matrix = baseMatrix;
        rhs = sourceVector;
        observation = responseMatrix;
        for cellIndex = 1:designCount
            matrix = matrix + interpolation(cellIndex) * ...
                reshape(cellMatrix(cellIndex,:,:),stateCount,stateCount);
            rhs = rhs + interpolation(cellIndex) * cellSource(cellIndex,:).';
            observation = observation + interpolation(cellIndex) * ...
                reshape(cellResponse(cellIndex,:,:),responseCount,stateCount);
        end
        state = matrix\rhs;
        if any(~isfinite(state))
            error("radia:topopt:ElectromagnetState", ...
                "The VIM state solve returned a non-finite value.");
        end
        response = observation*state;
        residual = response-target;
        weightedResidual = weights.*residual;
        adjoint = matrix.'\(observation.'*weightedResidual);
        if any(~isfinite(adjoint))
            error("radia:topopt:ElectromagnetAdjoint", ...
                "The VIM adjoint solve returned a non-finite value.");
        end
        gradient = zeros(designCount,1);
        for cellIndex = 1:designCount
            scale = interpolationDerivative(cellIndex);
            dMatrix = scale * ...
                reshape(cellMatrix(cellIndex,:,:),stateCount,stateCount);
            dSource = scale * cellSource(cellIndex,:).';
            dObservation = scale * ...
                reshape(cellResponse(cellIndex,:,:),responseCount,stateCount);
            gradient(cellIndex) = weightedResidual.' * ...
                (dObservation*state) + adjoint.'*(dSource-dMatrix*state);
        end
        regularization = 0.5*options.DensityRegularization * ...
            sum(volumes.*density.^2);
        gradient = gradient + options.DensityRegularization*volumes.*density;
        volumeFraction = (volumes.'*density)/sum(volumes);
        evaluation = struct( ...
            "objective",0.5*(residual.'*weightedResidual)+regularization, ...
            "gradient",gradient, ...
            "constraints",volumeFraction-options.VolumeFraction, ...
            "constraint_jacobian",volumes/sum(volumes), ...
            "state",state,"adjoint",adjoint,"response",response, ...
            "target",target,"volume_fraction",volumeFraction, ...
            "field_residual_norm",norm(residual));
    end
end

function value = localDefaultExpand(value,defaultValue,count,name)
if isempty(value)
    value = defaultValue;
else
    value = localExpand(value,count,name);
end
end

function value = localExpand(value,count,name)
value = double(value(:));
if isscalar(value)
    value = repmat(value,count,1);
end
if numel(value) ~= count
    error("radia:topopt:ElectromagnetShape", ...
        "%s must be scalar or contain %d values.",name,count);
end
end
