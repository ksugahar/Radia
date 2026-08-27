function value=upstreamVisualization(functionName,study,varargin)
%UPSTREAMVISUALIZATION Invoke pinned Optuna Plotly or matplotlib output.
arguments
    functionName (1,1) string
    study = []
end
arguments (Repeating)
    varargin
end
backend="plotly";
keep=true(1,numel(varargin));
for index=1:2:numel(varargin)-1
    if isstring(varargin{index}) || ischar(varargin{index})
        if strcmpi(string(varargin{index}),"Backend")
            backend=lower(string(varargin{index+1}));
            keep(index:index+1)=false;
        end
    end
end
varargin=varargin(keep);
if ~ismember(backend,["plotly","matplotlib"])
    error("radia:optuna:VisualizationBackend", ...
        "Backend must be 'plotly' or 'matplotlib'.");
end
moduleName="optuna.visualization";
if backend=="matplotlib"
    matplotlibModule=py.importlib.import_module("matplotlib");
    useBackend=py.builtins.getattr(matplotlibModule,"use");
    useBackend("Agg",pyargs("force",true));
    moduleName=moduleName+".matplotlib";
end
try
    module=py.importlib.import_module(char(moduleName));
    operation=py.builtins.getattr(module,char(functionName));
catch cause
    error("radia:optuna:VisualizationUnavailable", ...
        "Optuna visualization backend '%s' is unavailable: %s", ...
        backend,cause.message);
end
if functionName=="is_available"
    value=logical(operation());
    return
end
[pythonStudy]=radia.optuna.internal.toUpstreamStudy(study);
converted=cellfun(@convertValue,varargin,"UniformOutput",false);
if isempty(converted)
    value=operation(pythonStudy);
elseif mod(numel(converted),2)==0 && ...
        all(cellfun(@(item)ischar(item)||isstring(item),converted(1:2:end)))
    keyword=cell(1,numel(converted));
    keyword(1:2:end)=cellfun(@char,converted(1:2:end), ...
        "UniformOutput",false);
    keyword(2:2:end)=converted(2:2:end);
    value=operation(pythonStudy,pyargs(keyword{:}));
else
    value=operation(pythonStudy,converted{:});
end
end

function value=convertValue(source)
if startsWith(class(source),"py.")
    value=source;
elseif isstring(source)
    if isscalar(source)
        value=char(source);
    else
        value=py.list(cellstr(reshape(source,1,[])));
    end
elseif isnumeric(source) && ~isscalar(source)
    value=py.list(num2cell(reshape(source,1,[])));
elseif iscell(source)
    value=py.list(cellfun(@convertValue,source,"UniformOutput",false));
else
    value=source;
end
end
