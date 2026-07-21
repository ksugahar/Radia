function result = callPython(moduleName, functionName, positional, options)
%CALLPYTHON Invoke one explicit Radia Python fallback in-process.
% This boundary is for initialization, update, artifact generation, and batch
% solves. It is intentionally unsuitable for per-step Simulink execution.

arguments
    moduleName (1,1) string
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end

radia.internal.pythonProcessPath("capture");
runtime = requireInProcessPython();
module = py.importlib.import_module(char(moduleName));
callable = py.getattr(module, char(functionName));
pythonArguments = cellfun(@toPython, positional, UniformOutput=false);

keywordNames = fieldnames(options.Keywords);
if isempty(keywordNames)
    raw = callable(pythonArguments{:});
else
    keywordArguments = cell(1, 2 * numel(keywordNames));
    for index = 1:numel(keywordNames)
        keywordArguments{2 * index - 1} = keywordNames{index};
        keywordArguments{2 * index} = toPython(options.Keywords.(keywordNames{index}));
    end
    raw = callable(pythonArguments{:}, pyargs(keywordArguments{:}));
end

result = struct( ...
    "backend", "python-fallback", ...
    "module", moduleName, ...
    "function", functionName, ...
    "value", toMatlab(raw), ...
    "python", runtime);
end

function runtime = requireInProcessPython()
environment = pyenv;
if environment.Status == "NotLoaded"
    try
        environment = pyenv(ExecutionMode="InProcess");
    catch exception
        error("radia:python:Unavailable", ...
            ["Radia's MATLAB fallback requires an in-process CPython runtime. " ...
             "Configure pyenv with the supported Python 3.12 executable.\n%s"], ...
            exception.message);
    end
end
if environment.ExecutionMode ~= "InProcess"
    error("radia:python:ExecutionMode", ...
        ["Radia's MATLAB fallback requires pyenv ExecutionMode='InProcess' " ...
         "so Python objects and arrays have a deterministic lifetime. " ...
         "Restart MATLAB and configure pyenv before loading Python."]);
end
if ~startsWith(string(environment.Version), "3.12")
    error("radia:python:Version", ...
        "Radia currently supports Python 3.12 for MATLAB fallback; pyenv reports %s.", ...
        string(environment.Version));
end
runtime = struct( ...
    "status", string(environment.Status), ...
    "version", string(environment.Version), ...
    "executable", string(environment.Executable), ...
    "library", string(environment.Library), ...
    "execution_mode", string(environment.ExecutionMode));
end

function value = toPython(value)
if startsWith(string(class(value)), "py.")
    return
end
if isstring(value)
    if isscalar(value)
        value = py.str(char(value));
    else
        value = py.list(cellstr(value));
    end
    return
end
if ischar(value)
    value = py.str(value);
    return
end
if isnumeric(value) || islogical(value)
    if isscalar(value)
        return
    end
    value = py.numpy.asarray(value);
    return
end
if iscell(value)
    items = cellfun(@toPython, value, UniformOutput=false);
    value = py.list(items);
    return
end
if isstruct(value) && isscalar(value)
    names = fieldnames(value);
    keywordArguments = cell(1, 2 * numel(names));
    for index = 1:numel(names)
        keywordArguments{2 * index - 1} = names{index};
        keywordArguments{2 * index} = toPython(value.(names{index}));
    end
    value = py.dict(pyargs(keywordArguments{:}));
    return
end
error("radia:python:InputType", ...
    "Cannot convert MATLAB class %s to the Radia Python fallback contract.", ...
    class(value));
end

function value = toMatlab(value)
if isa(value, "py.NoneType")
    value = [];
elseif isa(value, "py.numpy.ndarray")
    value = double(value);
elseif isa(value, "py.numpy.generic")
    value = toMatlab(value.item());
elseif isa(value, "py.bool")
    value = logical(value);
elseif isa(value, "py.int") || isa(value, "py.float") || isa(value, "py.complex")
    value = double(value);
elseif isa(value, "py.str")
    value = string(value);
elseif isa(value, "py.dict")
    keys = cell(py.list(value.keys()));
    output = struct();
    for index = 1:numel(keys)
        sourceName = string(keys{index});
        fieldName = matlab.lang.makeValidName(sourceName);
        output.(fieldName) = toMatlab(value{keys{index}});
    end
    value = output;
elseif isa(value, "py.list") || isa(value, "py.tuple")
    items = cell(value);
    value = cellfun(@toMatlab, items, UniformOutput=false);
end
end
