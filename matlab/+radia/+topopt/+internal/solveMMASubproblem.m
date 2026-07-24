function result = solveMMASubproblem(design, evaluation, low, upp, ...
        alpha, beta, options)
%SOLVEMMASUBPROBLEM Solve one separable convex MMA approximation.
arguments
    design (:,1) double {mustBeReal,mustBeFinite}
    evaluation (1,1) struct
    low (:,1) double {mustBeReal,mustBeFinite}
    upp (:,1) double {mustBeReal,mustBeFinite}
    alpha (:,1) double {mustBeReal,mustBeFinite}
    beta (:,1) double {mustBeReal,mustBeFinite}
    options.ConstraintPenalty double {mustBePositive,mustBeFinite}=1000
    options.SlackQuadratic double {mustBePositive,mustBeFinite}=1
    options.CurvatureRegularization (1,1) double {mustBePositive,mustBeFinite}=1e-5
    options.MaxIterations (1,1) double {mustBeInteger,mustBePositive}=200
    options.ConstraintTolerance (1,1) double {mustBePositive,mustBeFinite}=1e-9
end
n = numel(design);
m = numel(evaluation.constraints);
if any(upp <= design) || any(design <= low) || any(beta <= alpha)
    error("radia:topopt:MMAAsymptote", ...
        "MMA asymptotes and subproblem bounds must strictly bracket the design.");
end

ux = upp - design;
xl = design - low;
span = max(upp - low, sqrt(eps));
regularization = options.CurvatureRegularization ./ span;

positiveObjective = max(evaluation.gradient, 0);
negativeObjective = max(-evaluation.gradient, 0);
objectiveCorrection = 1e-3 * ...
    (positiveObjective + negativeObjective) + regularization;
p0 = (positiveObjective + objectiveCorrection) .* ux.^2;
q0 = (negativeObjective + objectiveCorrection) .* xl.^2;

if m == 0
    next = localUnconstrainedMinimizer(p0, q0, low, upp, alpha, beta);
    result = struct( ...
        "schema", "radia.topopt.mma-subproblem/v1", ...
        "design", next, "slack", zeros(0,1), ...
        "multipliers", zeros(0,1), "exitflag", 1, ...
        "iterations", 0, "message", "closed-form unconstrained MMA step");
    return
end

jacobian = evaluation.constraint_jacobian;
positiveConstraints = max(jacobian.', 0);
negativeConstraints = max(-jacobian.', 0);
constraintCorrection = 1e-3 * ...
    (positiveConstraints + negativeConstraints) + ...
    options.CurvatureRegularization * ones(m,1) * (1 ./ span).';
P = (positiveConstraints + constraintCorrection) .* (ux.^2).';
Q = (negativeConstraints + constraintCorrection) .* (xl.^2).';
b = P * (1 ./ ux) + Q * (1 ./ xl) - evaluation.constraints;

penalty = localExpand(options.ConstraintPenalty, m, "ConstraintPenalty");
quadratic = localExpand(options.SlackQuadratic, m, "SlackQuadratic");
x0 = min(max(design, alpha), beta);
approximateConstraints = localApproximateConstraints(x0, low, upp, P, Q, b);
y0 = max(0, approximateConstraints + 10 * eps);
z0 = [x0; y0];
lowerBounds = [alpha; zeros(m,1)];
upperBounds = [beta; inf(m,1)];

settings = optimoptions("fmincon", ...
    "Algorithm", "interior-point", ...
    "SpecifyObjectiveGradient", true, ...
    "SpecifyConstraintGradient", true, ...
    "Display", "none", ...
    "MaxIterations", options.MaxIterations, ...
    "MaxFunctionEvaluations", max(1000, 20 * (n + m)), ...
    "OptimalityTolerance", 1e-10, ...
    "ConstraintTolerance", options.ConstraintTolerance, ...
    "StepTolerance", 1e-12);
[solution, ~, exitflag, output, lambda] = fmincon( ...
    @objective, z0, [], [], [], [], lowerBounds, upperBounds, ...
    @constraints, settings);
if exitflag <= 0
    error("radia:topopt:MMASubproblemFailed", ...
        "MMA convex subproblem failed: %s", output.message);
end

result = struct( ...
    "schema", "radia.topopt.mma-subproblem/v1", ...
    "design", solution(1:n), ...
    "slack", solution(n+1:end), ...
    "multipliers", lambda.ineqnonlin(:), ...
    "exitflag", exitflag, ...
    "iterations", output.iterations, ...
    "message", string(output.message));

    function [value, gradient] = objective(z)
        x = z(1:n);
        y = z(n+1:end);
        value = sum(p0 ./ (upp - x) + q0 ./ (x - low)) + ...
            penalty.' * y + 0.5 * quadratic.' * (y.^2);
        gradient = [p0 ./ (upp - x).^2 - q0 ./ (x - low).^2; ...
            penalty + quadratic .* y];
    end

    function [c, ceq, gradient, equalityGradient] = constraints(z)
        x = z(1:n);
        y = z(n+1:end);
        c = localApproximateConstraints(x, low, upp, P, Q, b) - y;
        derivative = P ./ ((upp - x).^2).' - ...
            Q ./ ((x - low).^2).';
        gradient = [derivative.'; -eye(m)];
        ceq = zeros(0,1);
        equalityGradient = zeros(n+m,0);
    end
end

function values = localApproximateConstraints(x, low, upp, P, Q, b)
values = P * (1 ./ (upp - x)) + Q * (1 ./ (x - low)) - b;
end

function design = localUnconstrainedMinimizer(p, q, low, upp, alpha, beta)
sqrtP = sqrt(p);
sqrtQ = sqrt(q);
design = (sqrtQ .* upp + sqrtP .* low) ./ (sqrtP + sqrtQ);
design = min(max(design, alpha), beta);
end

function value = localExpand(value, count, name)
value = double(value(:));
if isscalar(value)
    value = repmat(value, count, 1);
end
if numel(value) ~= count
    error("radia:topopt:MMASubproblemOption", ...
        "%s must be scalar or match the inequality-constraint count.", name);
end
end
