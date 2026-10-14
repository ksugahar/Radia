function options = optimoptions(varargin)
%OPTIMOPTIONS Create options for radia.optuna.optimize.
%   OPTIONS = RADIA.OPTUNA.OPTIMOPTIONS(Name=Value, ...) returns an options
%   object, the way Global Optimization Toolbox solvers are configured.
%
%   OPTIONS = RADIA.OPTUNA.OPTIMOPTIONS("optuna", Name=Value, ...) accepts a
%   leading solver name, matching optimoptions(@ga, ...) usage.  "optuna" is
%   the only solver this package provides; any other name is an error rather
%   than a silently ignored argument.
%
%   OPTIONS = RADIA.OPTUNA.OPTIMOPTIONS(OLD, Name=Value, ...) copies OLD and
%   overrides the named options, matching optimoptions(oldopts, ...).
%
%   Example:
%       opts = radia.optuna.optimoptions(MaxTrials=200, Sampler="tpe", ...
%           Display="iter", OutputFcn=@myOutputFcn, UseParallel=true);
%
%   See also radia.optuna.optimize, radia.optuna.OptimizeOptions.

arguments (Repeating)
    varargin
end
base = [];
first = 1;
if ~isempty(varargin)
    leading = varargin{1};
    if isa(leading, "radia.optuna.OptimizeOptions")
        base = leading;
        first = 2;
    elseif isa(leading, "function_handle") || ...
            (isstring(leading) && isscalar(leading)) || ...
            (ischar(leading) && isrow(leading))
        name = string(leading);
        if isa(leading, "function_handle")
            name = string(func2str(leading));
        end
        if mod(numel(varargin) - 1, 2) == 0 && ...
                ismember(lower(name), ["optuna", "@optuna", "radia.optuna.optimize"])
            first = 2;
        elseif mod(numel(varargin), 2) ~= 0
            error("radia:optuna:OptimOptions", ...
                "'%s' is not a solver this package provides. " + ...
                "Use ""optuna"", or pass only Name=Value options.", name);
        end
    end
end

pairs = varargin(first:end);
if mod(numel(pairs), 2) ~= 0
    error("radia:optuna:OptimOptions", ...
        "Options must be supplied as Name=Value pairs.");
end

if isempty(base)
    options = radia.optuna.OptimizeOptions();
else
    options = base;
end

known = string(properties("radia.optuna.OptimizeOptions"));
for index = 1:2:numel(pairs)
    name = string(pairs{index});
    match = known(strcmpi(known, name));
    if isempty(match)
        error("radia:optuna:OptimOptions", ...
            "Unknown option '%s'. Available options: %s.", ...
            name, strjoin(sort(known), ", "));
    end
    options.(match(1)) = pairs{index + 1};
end

% Re-run the constructor so callback normalization and cross-field checks
% apply to values assigned here, not only to those passed positionally.
settings = cell(1, 2 * numel(known));
for index = 1:numel(known)
    settings{2 * index - 1} = known(index);
    settings{2 * index} = options.(known(index));
end
options = radia.optuna.OptimizeOptions(settings{:});
end
