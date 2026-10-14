function runner = makeStreamFunctionAdjointRunner( ...
        responseMatrix, target, options)
%MAKESTREAMFUNCTIONADJOINTRUNNER Build an analytic stream-function optimizer.
%   The design is the vector of stream-function coefficients. The objective
%   combines weighted field residual and quadratic regularization. Optional
%   current samples impose two-sided linear inequalities.
arguments
    responseMatrix (:,:) double {mustBeReal,mustBeFinite}
    target (:,1) double {mustBeReal,mustBeFinite}
    options.RegularizationMatrix double {mustBeReal,mustBeFinite} = ...
        zeros(0,size(responseMatrix,2))
    options.Alpha (1,1) double {mustBeNonnegative,mustBeFinite} = 0
    options.ResponseWeights double {mustBeReal,mustBeFinite} = ...
        ones(size(responseMatrix,1),1)
    options.CurrentOperator double {mustBeReal,mustBeFinite} = ...
        zeros(0,size(responseMatrix,2))
    options.CurrentLimit double {mustBeReal,mustBeFinite} = double.empty
    options.InitialDesign double {mustBeReal,mustBeFinite} = ...
        zeros(size(responseMatrix,2),1)
    options.LowerBounds double {mustBeReal,mustBeFinite} = ...
        -ones(size(responseMatrix,2),1)
    options.UpperBounds double {mustBeReal,mustBeFinite} = ...
        ones(size(responseMatrix,2),1)
    options.MoveLimit double {mustBeReal,mustBeFinite} = double.empty
    options.MaxIterations (1,1) double {mustBeInteger,mustBePositive} = 60
    options.Solver (1,1) string ...
        {mustBeMember(options.Solver,["mma","sqp"])} = "mma"
    options.OptimizerOptions (1,1) struct = struct
end

[measurementCount,designCount] = size(responseMatrix);
if designCount == 0 || measurementCount == 0
    error("radia:topopt:StreamFunctionShape", ...
        "responseMatrix must have at least one row and one column.");
end
if numel(target) ~= measurementCount
    error("radia:topopt:StreamFunctionShape", ...
        "target must contain one value per response row.");
end
regularization = double(options.RegularizationMatrix);
if size(regularization,2) ~= designCount
    error("radia:topopt:StreamFunctionShape", ...
        "RegularizationMatrix must have one column per design coefficient.");
end
weights = localExpand(options.ResponseWeights,measurementCount, ...
    "ResponseWeights");
if any(weights <= 0)
    error("radia:topopt:StreamFunctionWeights", ...
        "ResponseWeights must be positive.");
end
currentOperator = double(options.CurrentOperator);
if size(currentOperator,2) ~= designCount
    error("radia:topopt:StreamFunctionShape", ...
        "CurrentOperator must have one column per design coefficient.");
end
if isempty(currentOperator)
    if ~isempty(options.CurrentLimit)
        error("radia:topopt:StreamFunctionCurrentLimit", ...
            "CurrentLimit requires a nonempty CurrentOperator.");
    end
    currentLimit = zeros(0,1);
else
    currentLimit = localExpand(options.CurrentLimit, ...
        size(currentOperator,1),"CurrentLimit");
    if any(currentLimit <= 0)
        error("radia:topopt:StreamFunctionCurrentLimit", ...
            "CurrentLimit must be positive.");
    end
end
initialDesign = localExpand(options.InitialDesign,designCount,"InitialDesign");
lower = localExpand(options.LowerBounds,designCount,"LowerBounds");
upper = localExpand(options.UpperBounds,designCount,"UpperBounds");
if any(lower >= upper) || any(initialDesign < lower) || ...
        any(initialDesign > upper)
    error("radia:topopt:StreamFunctionBounds", ...
        "Bounds must strictly bracket the initial stream-function design.");
end

optimizerOptions = options.OptimizerOptions;
optimizerOptions.LowerBounds = lower;
optimizerOptions.UpperBounds = upper;
optimizerOptions.MaxIterations = options.MaxIterations;
if ~isempty(options.MoveLimit)
    optimizerOptions.MoveLimit = localExpand( ...
        options.MoveLimit,designCount,"MoveLimit");
end
metadata = struct( ...
    "domain","stream-function", ...
    "parameterization","stream-function-coefficients", ...
    "gradient","analytic-adjoint", ...
    "response_count",measurementCount, ...
    "design_count",designCount, ...
    "current_constraint_count",size(currentOperator,1));
runner = radia.topopt.AdjointRunner(initialDesign,@evaluate, ...
    Solver=options.Solver,OptimizerOptions=optimizerOptions, ...
    Metadata=metadata);

    function evaluation = evaluate(design)
        residual = responseMatrix*design-target;
        regularized = regularization*design;
        evaluation = struct( ...
            "objective",0.5*(residual'*(weights.*residual)) + ...
                0.5*options.Alpha*(regularized'*regularized), ...
            "gradient",responseMatrix'*(weights.*residual) + ...
                options.Alpha*(regularization'*regularized));
        if ~isempty(currentOperator)
            current = currentOperator*design;
            evaluation.constraints = [ ...
                current-currentLimit; -current-currentLimit];
            evaluation.constraint_jacobian = [ ...
                currentOperator',-currentOperator'];
        end
    end
end

function value = localExpand(value,count,name)
value = double(value(:));
if isscalar(value)
    value = repmat(value,count,1);
end
if numel(value) ~= count
    error("radia:topopt:StreamFunctionShape", ...
        "%s must be scalar or contain %d values.",name,count);
end
end
