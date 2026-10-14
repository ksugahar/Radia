function result = optimizeIH(objective, x0, lower_bound, upper_bound, options)
%OPTIMIZEIH Optimize an IH drive/controller parameter vector.
%   RESULT = radia.simulink.optimizeIH(OBJECTIVE,X0,LB,UB) provides a
%   toolbox-aware outer loop. OBJECTIVE must return one finite scalar and
%   may call sim(Simulink.SimulationInput(...)), a fast waveform model, or a
%   Radia-generated LUT. fmincon is used when available; otherwise a bounded
%   fminsearch fallback keeps the workflow usable without Optimization
%   Toolbox.

arguments
    objective (1,1) function_handle
    x0 double {mustBeVector, mustBeFinite}
    lower_bound double {mustBeVector, mustBeFinite}
    upper_bound double {mustBeVector, mustBeFinite}
    options.UseFmincon (1,1) logical = true
    options.MaxIterations (1,1) double {mustBeInteger, mustBePositive} = 100
    options.FunctionTolerance (1,1) double {mustBePositive} = 1.0e-6
    options.Display (1,1) string = "final"
    options.UsesSimulink (1,1) logical = false
end

x0 = x0(:).';
lower_bound = lower_bound(:).';
upper_bound = upper_bound(:).';
if ~isequal(size(x0), size(lower_bound), size(upper_bound)) || ...
        any(lower_bound >= upper_bound) || any(x0 < lower_bound | x0 > upper_bound)
    error("radia:simulink:OptimizationBounds", ...
        "x0, lower_bound, and upper_bound must have equal sizes and strict bounds.");
end

evaluate = @(x) localObjective(objective, min(max(x(:).', lower_bound), upper_bound));
useFmincon = options.UseFmincon && exist("fmincon", "file") == 2;
if useFmincon
    solverOptions = optimoptions("fmincon", ...
        "Display", options.Display, ...
        "MaxIterations", options.MaxIterations, ...
        "FunctionTolerance", options.FunctionTolerance, ...
        "SpecifyObjectiveGradient", false);
    [x, fval, exitflag, output] = fmincon(evaluate, x0, [], [], [], [], ...
        lower_bound, upper_bound, [], solverOptions);
    solver = "fmincon";
else
    solverOptions = optimset("Display", char(options.Display), ...
        "MaxIter", options.MaxIterations, ...
        "TolFun", options.FunctionTolerance);
    [z, ~, exitflag, output] = fminsearch(evaluate, x0, solverOptions);
    x = min(max(z(:).', lower_bound), upper_bound);
    fval = evaluate(x);
    solver = "bounded-fminsearch";
end

result = struct( ...
    "schema", "radia.ih.simulink.optimization.v1", ...
    "x", x(:).', ...
    "objective_value", fval, ...
    "exitflag", exitflag, ...
    "output", output, ...
    "solver", solver, ...
    "used_simulink", options.UsesSimulink, ...
    "note", "The objective owns Simulink.SimulationInput and parsim policy.");
end

function value = localObjective(objective, x)
value = objective(x);
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
    error("radia:simulink:Objective", ...
        "objective must return one finite numeric scalar.");
end
value = double(value);
end
