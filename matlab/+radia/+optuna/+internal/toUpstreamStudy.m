function [pythonStudy,pythonTrials,optunaModule]=toUpstreamStudy( ...
        trials,studyDirection)
%TOUPSTREAMSTUDY Convert MATLAB FrozenTrial snapshots to Optuna 4.9.
arguments
    trials
    studyDirection = []
end
environment=pyenv;
if environment.Status=="NotLoaded"
    environment=pyenv(ExecutionMode="InProcess");
end
if environment.ExecutionMode~="InProcess"
    error("radia:optuna:UpstreamPython", ...
        "Upstream adapters require in-process Python.");
end
optunaModule=py.importlib.import_module("optuna");
version=string(py.builtins.getattr(optunaModule,"__version__"));
if version~="4.9.0"
    error("radia:optuna:UpstreamVersion", ...
        "Upstream adapters require optuna==4.9.0, found %s.",version);
end
if isa(trials,"radia.optuna.Study")
    studyDirection=trials.Directions;
    trials=trials.get_trials();
end
directions=lower(reshape(string(studyDirection),1,[]));
if isempty(directions), directions="minimize"; end
createStudy=py.builtins.getattr(optunaModule,"create_study");
if numel(directions)==1
    pythonStudy=createStudy(pyargs("direction",char(directions)));
else
    pythonStudy=createStudy(pyargs("directions", ...
        py.list(cellstr(directions))));
end
trialModule=py.importlib.import_module("optuna.trial");
distributionModule=py.importlib.import_module("optuna.distributions");
jsonModule=py.importlib.import_module("json");
loads=py.builtins.getattr(jsonModule,"loads");
createTrial=py.builtins.getattr(trialModule,"create_trial");
stateClass=py.builtins.getattr(trialModule,"TrialState");
toDistribution=py.builtins.getattr( ...
    distributionModule,"json_to_distribution");
pythonTrialCells=cell(1,numel(trials));
for index=1:numel(trials)
    source=trials(index);
    params=loads(char(jsonencode(source.Params)));
    userAttrs=loads(char(jsonencode(source.UserAttrs)));
    systemAttrs=loads(char(jsonencode(source.SystemAttrs)));
    distributions=py.dict;
    names=string(fieldnames(source.Distributions));
    for nameIndex=1:numel(names)
        spec=source.Distributions.(names(nameIndex));
        distributions{char(names(nameIndex))}=toDistribution( ...
            char(radia.optuna.distribution_to_json(spec)));
    end
    intermediate=py.dict;
    for row=1:height(source.IntermediateValues)
        intermediate{int64(source.IntermediateValues.Step(row))}= ...
            double(source.IntermediateValues.Value(row));
    end
    state=py.builtins.getattr(stateClass,char(source.State));
    keyword={"state",state,"params",params, ...
        "distributions",distributions, ...
        "intermediate_values",intermediate, ...
        "user_attrs",userAttrs,"system_attrs",systemAttrs};
    values=reshape(double(source.Values),1,[]);
    if numel(values)==1 && ~isnan(values)
        keyword(end+1:end+2)={"value",values};
    elseif ~all(isnan(values))
        keyword(end+1:end+2)={"values",py.list(num2cell(values))};
    end
    pythonTrialCells{index}=createTrial(pyargs(keyword{:}));
    pythonStudy.add_trial(pythonTrialCells{index});
end
pythonTrials=py.list(pythonTrialCells);
end
