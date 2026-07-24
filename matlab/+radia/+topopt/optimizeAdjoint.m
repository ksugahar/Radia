function result = optimizeAdjoint(initialDesign, evaluateFcn, options)
%OPTIMIZEADJOINT Minimize a checked adjoint problem with MMA or gradient SQP.
%   EVALUATEFCN(X) returns objective, gradient, and optional constraints <= 0
%   with a design-by-constraint constraint_jacobian. Optional equalities and
%   equality_jacobian are accepted by SQP. MMA requires finite design bounds
%   and inequality constraints only.
arguments
    initialDesign (:,1) double {mustBeReal,mustBeFinite}
    evaluateFcn (1,1) function_handle
    options.Solver (1,1) string ...
        {mustBeMember(options.Solver,["mma","sqp"])}="mma"
    options.LowerBounds double=-inf(size(initialDesign))
    options.UpperBounds double=inf(size(initialDesign))
    options.MaxIterations (1,1) double {mustBeInteger,mustBePositive}=50
    options.MaxFunctionEvaluations (1,1) double ...
        {mustBeInteger,mustBePositive}=2000
    options.OptimalityTolerance (1,1) double {mustBePositive,mustBeFinite}=1e-6
    options.ConstraintTolerance (1,1) double {mustBePositive,mustBeFinite}=1e-6
    options.StepTolerance (1,1) double {mustBePositive,mustBeFinite}=1e-7
    options.MoveLimit double=double.empty
    options.InitialAsymptote (1,1) double {mustBePositive,mustBeFinite}=0.5
    options.AsymptoteIncrease (1,1) double ...
        {mustBeGreaterThan(options.AsymptoteIncrease,1),mustBeFinite}=1.2
    options.AsymptoteDecrease (1,1) double ...
        {mustBePositive,mustBeLessThan(options.AsymptoteDecrease,1)}=0.7
    options.MinimumAsymptote (1,1) double {mustBePositive,mustBeFinite}=0.01
    options.MaximumAsymptote (1,1) double {mustBePositive,mustBeFinite}=10
    options.ConstraintPenalty double {mustBePositive,mustBeFinite}=1000
    options.SlackQuadratic double {mustBePositive,mustBeFinite}=1
    options.CurvatureRegularization (1,1) double {mustBePositive,mustBeFinite}=1e-5
    options.SubproblemMaxIterations (1,1) double ...
        {mustBeInteger,mustBePositive}=200
    options.Display (1,1) string ...
        {mustBeMember(options.Display,["none","iter","final"])}="none"
    options.ProgressFcn=[]
end

n = numel(initialDesign);
lower = localExpand(options.LowerBounds, n, "LowerBounds");
upper = localExpand(options.UpperBounds, n, "UpperBounds");
if any(~isreal(lower)) || any(~isreal(upper)) || any(lower >= upper)
    error("radia:topopt:AdjointBounds", ...
        "LowerBounds and UpperBounds must strictly bracket every variable.");
end
if any(initialDesign < lower) || any(initialDesign > upper)
    error("radia:topopt:AdjointBounds", ...
        "The initial design must lie inside the supplied bounds.");
end
if ~isempty(options.ProgressFcn) && ~isa(options.ProgressFcn,"function_handle")
    error("radia:topopt:AdjointProgress", ...
        "ProgressFcn must be empty or a function handle.");
end

initialEvaluation = radia.topopt.internal.evaluateAdjoint( ...
    evaluateFcn, initialDesign);
switch options.Solver
    case "sqp"
        result = localSQP(initialDesign, evaluateFcn, initialEvaluation, ...
            lower, upper, options);
    case "mma"
        if any(~isfinite(lower)) || any(~isfinite(upper))
            error("radia:topopt:MMABounds", ...
                "MMA requires finite lower and upper bounds.");
        end
        if ~isempty(initialEvaluation.equalities)
            error("radia:topopt:MMAEqualities", ...
                "MMA accepts inequalities only; use SQP for equality constraints.");
        end
        result = localMMA(initialDesign, evaluateFcn, initialEvaluation, ...
            lower, upper, options);
end
end

function result = localSQP(initialDesign, evaluateFcn, initialEvaluation, ...
        lower, upper, options)
constraintCount = numel(initialEvaluation.constraints);
equalityCount = numel(initialEvaluation.equalities);
cacheDesign = initialDesign;
cacheEvaluation = initialEvaluation;
evaluationCount = 1;
history = localHistoryTable();
previousDesign = initialDesign;

settings = optimoptions("fmincon", ...
    "Algorithm", "sqp", ...
    "SpecifyObjectiveGradient", true, ...
    "SpecifyConstraintGradient", true, ...
    "Display", char(options.Display), ...
    "MaxIterations", options.MaxIterations, ...
    "MaxFunctionEvaluations", options.MaxFunctionEvaluations, ...
    "OptimalityTolerance", options.OptimalityTolerance, ...
    "ConstraintTolerance", options.ConstraintTolerance, ...
    "StepTolerance", options.StepTolerance, ...
    "OutputFcn", @recordIteration);
if constraintCount == 0 && equalityCount == 0
    nonlinearConstraints = [];
else
    nonlinearConstraints = @constraintFunction;
end
[design, ~, exitflag, output, multipliers] = fmincon( ...
    @objectiveFunction, initialDesign, [], [], [], [], lower, upper, ...
    nonlinearConstraints, settings);
evaluation = evaluateCached(design);
converged = exitflag > 0 && ...
    localMaxViolation(evaluation.constraints) <= options.ConstraintTolerance && ...
    localMaxEquality(evaluation.equalities) <= options.ConstraintTolerance;
result = localResult("sqp", design, evaluation, history, converged, ...
    exitflag, output, multipliers, evaluationCount);

    function evaluation = evaluateCached(designValue)
        designValue = designValue(:);
        if ~isequaln(designValue, cacheDesign)
            cacheEvaluation = radia.topopt.internal.evaluateAdjoint( ...
                evaluateFcn, designValue, ...
                ConstraintCount=constraintCount, EqualityCount=equalityCount);
            cacheDesign = designValue;
            evaluationCount = evaluationCount + 1;
        end
        evaluation = cacheEvaluation;
    end

    function [value, gradient] = objectiveFunction(designValue)
        evaluation = evaluateCached(designValue);
        value = evaluation.objective;
        gradient = evaluation.gradient;
    end

    function [constraints, equalities, jacobian, equalityJacobian] = ...
            constraintFunction(designValue)
        evaluation = evaluateCached(designValue);
        constraints = evaluation.constraints;
        equalities = evaluation.equalities;
        jacobian = evaluation.constraint_jacobian;
        equalityJacobian = evaluation.equality_jacobian;
    end

    function stop = recordIteration(designValue, optimValues, state)
        stop = false;
        if state ~= "iter" && state ~= "init"
            return
        end
        evaluation = evaluateCached(designValue);
        if isempty(history)
            step = 0;
        else
            step = norm(designValue(:) - previousDesign, inf);
        end
        previousDesign = designValue(:);
        history(end+1,:) = {double(optimValues.iteration), ... %#ok<AGROW>
            evaluation.objective, localMaxViolation(evaluation.constraints), ...
            localMaxEquality(evaluation.equalities), step, NaN, ...
            double(optimValues.firstorderopt), evaluationCount, 0, 0};
        localNotify(options.ProgressFcn, "sqp", designValue(:), ...
            evaluation, history(end,:));
    end
end

function result = localMMA(initialDesign, evaluateFcn, initialEvaluation, ...
        lower, upper, options)
n = numel(initialDesign);
constraintCount = numel(initialEvaluation.constraints);
span = upper - lower;
if isempty(options.MoveLimit)
    move = 0.2 * span;
else
    move = localExpand(options.MoveLimit, n, "MoveLimit");
end
if any(~isfinite(move)) || any(move <= 0)
    error("radia:topopt:MMAMoveLimit", ...
        "MoveLimit must contain finite positive design-unit values.");
end

design = initialDesign;
oldDesign1 = design;
oldDesign2 = design;
low = design - options.InitialAsymptote * span;
upp = design + options.InitialAsymptote * span;
evaluation = initialEvaluation;
evaluationCount = 1;
history = localHistoryTable();
history(end+1,:) = {0,evaluation.objective, ...
    localMaxViolation(evaluation.constraints),0,0,NaN, ...
    norm(evaluation.gradient,inf),evaluationCount,0,0};
localNotify(options.ProgressFcn, "mma", design, evaluation, history(end,:));
converged = false;
exitflag = 0;
message = "Maximum MMA iterations reached.";
multipliers = struct("ineqlin",zeros(0,1),"eqlin",zeros(0,1), ...
    "ineqnonlin",zeros(constraintCount,1),"eqnonlin",zeros(0,1), ...
    "lower",zeros(n,1),"upper",zeros(n,1));

for iteration = 1:options.MaxIterations
    [low, upp] = localUpdateAsymptotes(design, oldDesign1, oldDesign2, ...
        low, upp, lower, upper, iteration, options);
    alpha = max([lower, design - move, ...
        low + 0.1 * (design - low)], [], 2);
    beta = min([upper, design + move, ...
        upp - 0.1 * (upp - design)], [], 2);
    subproblem = radia.topopt.internal.solveMMASubproblem( ...
        design, evaluation, low, upp, alpha, beta, ...
        ConstraintPenalty=options.ConstraintPenalty, ...
        SlackQuadratic=options.SlackQuadratic, ...
        CurvatureRegularization=options.CurvatureRegularization, ...
        MaxIterations=options.SubproblemMaxIterations, ...
        ConstraintTolerance=options.ConstraintTolerance * 0.1);
    nextDesign = subproblem.design;
    nextEvaluation = radia.topopt.internal.evaluateAdjoint( ...
        evaluateFcn, nextDesign, ConstraintCount=constraintCount, EqualityCount=0);
    evaluationCount = evaluationCount + 1;
    step = norm(nextDesign - design, inf);
    objectiveChange = abs(nextEvaluation.objective - evaluation.objective);
    multipliers.ineqnonlin = subproblem.multipliers;
    firstOrder = localKKTResidual(nextDesign, nextEvaluation, ...
        lower, upper, subproblem.multipliers);
    history(end+1,:) = {iteration,nextEvaluation.objective, ...
        localMaxViolation(nextEvaluation.constraints),0,step,objectiveChange, ...
        firstOrder,evaluationCount,subproblem.iterations, ...
        localMaxViolation(subproblem.slack)}; %#ok<AGROW>
    localNotify(options.ProgressFcn, "mma", nextDesign, ...
        nextEvaluation, history(end,:));

    oldDesign2 = oldDesign1;
    oldDesign1 = design;
    design = nextDesign;
    evaluation = nextEvaluation;
    feasible = localMaxViolation(evaluation.constraints) <= ...
        options.ConstraintTolerance;
    objectiveScale = max(1, abs(evaluation.objective));
    if feasible && step <= options.StepTolerance && ...
            (firstOrder <= options.OptimalityTolerance || ...
            objectiveChange <= options.OptimalityTolerance * objectiveScale)
        converged = true;
        exitflag = 1;
        message = "MMA step, feasibility, and first-order tolerances satisfied.";
        break
    end
end
output = struct("iterations",height(history)-1, ...
    "message",message,"algorithm","moving-asymptotes/interior-point-subproblem");
result = localResult("mma", design, evaluation, history, converged, ...
    exitflag, output, multipliers, evaluationCount);
end

function [low, upp] = localUpdateAsymptotes(design, old1, old2, ...
        low, upp, lower, upper, iteration, options)
span = upper - lower;
if iteration <= 2
    low = design - options.InitialAsymptote * span;
    upp = design + options.InitialAsymptote * span;
else
    trend = (design - old1) .* (old1 - old2);
    factor = ones(size(design));
    factor(trend > 0) = options.AsymptoteIncrease;
    factor(trend < 0) = options.AsymptoteDecrease;
    low = design - factor .* (old1 - low);
    upp = design + factor .* (upp - old1);
    low = max(design - options.MaximumAsymptote * span, ...
        min(design - options.MinimumAsymptote * span, low));
    upp = max(design + options.MinimumAsymptote * span, ...
        min(design + options.MaximumAsymptote * span, upp));
end
end

function value = localKKTResidual(design, evaluation, lower, upper, lambda)
lagrangianGradient = evaluation.gradient + ...
    evaluation.constraint_jacobian * lambda(:);
tolerance = sqrt(eps) * max(1, max(abs([lower;upper])));
projected = lagrangianGradient;
atLower = design <= lower + tolerance;
atUpper = design >= upper - tolerance;
projected(atLower & lagrangianGradient > 0) = 0;
projected(atUpper & lagrangianGradient < 0) = 0;
value = norm(projected, inf);
end

function result = localResult(solver, design, evaluation, history, ...
        converged, exitflag, output, multipliers, evaluationCount)
result = struct( ...
    "schema", "radia.topopt.adjoint-optimization/v1", ...
    "solver", solver, ...
    "design", design, ...
    "objective", evaluation.objective, ...
    "gradient", evaluation.gradient, ...
    "constraints", evaluation.constraints, ...
    "constraint_jacobian", evaluation.constraint_jacobian, ...
    "equalities", evaluation.equalities, ...
    "equality_jacobian", evaluation.equality_jacobian, ...
    "evaluation", evaluation, ...
    "history", history, ...
    "converged", logical(converged), ...
    "exitflag", double(exitflag), ...
    "output", output, ...
    "multipliers", multipliers, ...
    "evaluation_count", double(evaluationCount));
end

function history = localHistoryTable()
history = table('Size',[0,10], ...
    'VariableTypes',repmat({'double'},1,10), ...
    'VariableNames',{'Iteration','Objective','MaxConstraint', ...
    'MaxEquality','StepInfinity','ObjectiveChange','FirstOrderOptimality', ...
    'EvaluationCount','SubproblemIterations','MaxSlack'});
end

function value = localMaxViolation(constraints)
if isempty(constraints)
    value = 0;
else
    value = max(0, max(constraints));
end
end

function value = localMaxEquality(equalities)
if isempty(equalities)
    value = 0;
else
    value = max(abs(equalities));
end
end

function localNotify(progressFcn, solver, design, evaluation, historyRow)
if isempty(progressFcn)
    return
end
progressFcn(struct( ...
    "schema", "radia.topopt.adjoint-progress/v1", ...
    "solver", solver, "design", design, ...
    "evaluation", evaluation, "history", historyRow));
end

function value = localExpand(value, count, name)
value = double(value(:));
if isscalar(value)
    value = repmat(value, count, 1);
end
if numel(value) ~= count
    error("radia:topopt:AdjointShape", ...
        "%s must be scalar or match the design size.", name);
end
end
