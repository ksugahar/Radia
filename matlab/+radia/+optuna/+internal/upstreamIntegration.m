function value=upstreamIntegration(symbolName,varargin)
%UPSTREAMINTEGRATION Resolve or invoke an Optuna 4.9 integration export.
arguments
    symbolName (1,1) string
end
arguments (Repeating)
    varargin
end
environment=pyenv;
if environment.Status=="NotLoaded"
    environment=pyenv(ExecutionMode="InProcess");
end
if environment.ExecutionMode~="InProcess"
    error("radia:optuna:IntegrationPython", ...
        "Optuna integrations require in-process Python.");
end
try
    optunaModule=py.importlib.import_module("optuna");
    version=string(py.builtins.getattr(optunaModule,"__version__"));
    if version~="4.9.0"
        error("radia:optuna:IntegrationVersion", ...
            "Optuna integrations require optuna==4.9.0, found %s.",version);
    end
    module=py.importlib.import_module("optuna.integration");
    symbol=py.builtins.getattr(module,char(symbolName));
catch cause
    error("radia:optuna:IntegrationUnavailable", ...
        "Optuna integration '%s' is unavailable: %s", ...
        symbolName,cause.message);
end
if isempty(varargin)
    value=symbol;
    return
end
converted=cellfun(@convertValue,varargin,"UniformOutput",false);
if mod(numel(converted),2)==0 && ...
        all(cellfun(@(item)ischar(item)||isstring(item),converted(1:2:end)))
    names=cellfun(@char,converted(1:2:end),"UniformOutput",false);
    keyword=cell(1,numel(converted));
    keyword(1:2:end)=names;
    keyword(2:2:end)=converted(2:2:end);
    value=symbol(pyargs(keyword{:}));
else
    value=symbol(converted{:});
end
end

function value=convertValue(source)
if startsWith(class(source),"py.")
    value=source;
elseif isstring(source) && isscalar(source)
    value=char(source);
elseif iscell(source)
    value=py.list(cellfun(@convertValue,source,"UniformOutput",false));
elseif isstruct(source)
    jsonModule=py.importlib.import_module("json");
    loads=py.builtins.getattr(jsonModule,"loads");
    value=loads(char(jsonencode(source)));
else
    value=source;
end
end
