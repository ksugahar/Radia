function result = checkAdjointGradient(design, evaluateFcn, options)
%CHECKADJOINTGRADIENT Directional finite-difference QA for adjoint derivatives.
%   This diagnostic never participates in production optimization. It checks
%   the complete objective/constraint derivative returned by EVALUATEFCN.
arguments
    design (:,1) double {mustBeReal,mustBeFinite}
    evaluateFcn (1,1) function_handle
    options.Directions double=double.empty
    options.RelativeStep (1,1) double {mustBePositive,mustBeFinite}=1e-6
    options.RelativeTolerance (1,1) double {mustBePositive,mustBeFinite}=2e-4
    options.AbsoluteTolerance (1,1) double {mustBePositive,mustBeFinite}=1e-7
    options.ErrorOnFailure (1,1) logical=true
end

base = radia.topopt.internal.evaluateAdjoint(evaluateFcn, design);
directions = options.Directions;
if isempty(directions)
    directions = eye(numel(design));
end
if size(directions,1) ~= numel(design) || ...
        ~isreal(directions) || any(~isfinite(directions), "all")
    error("radia:topopt:AdjointDirections", ...
        "Directions must be design_count-by-direction_count and finite.");
end
directionCount = size(directions,2);
history = table('Size',[directionCount,7], ...
    'VariableTypes',{'double','double','double','double','double','double','logical'}, ...
    'VariableNames',{'Direction','ObjectiveAnalytic','ObjectiveNumeric', ...
    'ObjectiveError','ConstraintError','Step','Passed'});

scale = max(1, norm(design, inf));
for index = 1:directionCount
    direction = directions(:,index);
    directionNorm = norm(direction);
    if directionNorm == 0
        error("radia:topopt:AdjointDirections", ...
            "Every QA direction must be nonzero.");
    end
    direction = direction / directionNorm;
    step = options.RelativeStep * scale;
    plus = radia.topopt.internal.evaluateAdjoint( ...
        evaluateFcn, design + step * direction, ...
        ConstraintCount=numel(base.constraints), ...
        EqualityCount=numel(base.equalities));
    minus = radia.topopt.internal.evaluateAdjoint( ...
        evaluateFcn, design - step * direction, ...
        ConstraintCount=numel(base.constraints), ...
        EqualityCount=numel(base.equalities));

    objectiveAnalytic = base.gradient.' * direction;
    objectiveNumeric = (plus.objective - minus.objective) / (2 * step);
    objectiveError = abs(objectiveAnalytic - objectiveNumeric);
    objectiveLimit = options.AbsoluteTolerance + options.RelativeTolerance * ...
        max(abs(objectiveAnalytic), abs(objectiveNumeric));

    analyticConstraints = [base.constraint_jacobian, ...
        base.equality_jacobian].' * direction;
    numericConstraints = ([plus.constraints; plus.equalities] - ...
        [minus.constraints; minus.equalities]) / (2 * step);
    if isempty(analyticConstraints)
        constraintError = 0;
        constraintPassed = true;
    else
        errors = abs(analyticConstraints - numericConstraints);
        limits = options.AbsoluteTolerance + options.RelativeTolerance * ...
            max(abs(analyticConstraints), abs(numericConstraints));
        constraintError = max(errors);
        constraintPassed = all(errors <= limits);
    end
    passed = objectiveError <= objectiveLimit && constraintPassed;
    history(index,:) = {index, objectiveAnalytic, objectiveNumeric, ...
        objectiveError, constraintError, step, passed};
end

result = struct( ...
    "schema", "radia.topopt.adjoint-gradient-check/v1", ...
    "passed", all(history.Passed), ...
    "evaluation", base, ...
    "history", history, ...
    "evaluation_count", 1 + 2 * directionCount);
if options.ErrorOnFailure && ~result.passed
    failed = history(~history.Passed,:);
    error("radia:topopt:AdjointGradientCheck", ...
        "Adjoint directional derivative check failed in direction(s): %s", ...
        strjoin(string(failed.Direction), ", "));
end
end
