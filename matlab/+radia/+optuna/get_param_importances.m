function result=get_param_importances(study,options)
%GET_PARAM_IMPORTANCES Evaluate parameter importance with Optuna 4.9.0.
arguments
    study (1,1) radia.optuna.Study
    options.evaluator (1,1) string = "fanova"
    options.params string = strings(1,0)
    options.normalize (1,1) logical = true
    options.objective_index (1,1) double {mustBeInteger,mustBePositive} = 1
    options.seed (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.n_trees (1,1) double {mustBeInteger,mustBePositive} = 64
    options.max_depth (1,1) double {mustBeInteger,mustBePositive} = 64
    options.target_quantile (1,1) double = 0.1
    options.region_quantile (1,1) double = 1.0
    options.baseline_quantile (1,1) double = NaN
    options.evaluate_on_local (1,1) logical = true
    options.target = []
end
if options.objective_index>numel(study.Directions)
    error("radia:optuna:ImportanceObjective", ...
        "objective_index exceeds the number of study objectives.");
end
optuna=importPinnedOptuna();
importance=py.importlib.import_module("optuna.importance");
pythonStudy=optuna.create_study(pyargs( ...
    "direction",char(study.Directions(options.objective_index))));
numbers=study.TrialTable.TrialNumber(study.TrialTable.State=="COMPLETE");
if isempty(numbers)
    error("radia:optuna:ImportanceTrials", ...
        "Parameter importance requires at least one COMPLETE trial.");
end
for number=reshape(numbers,1,[])
    parameterRows=find(study.ParamTable.TrialNumber==number)';
    params=py.dict;
    distributions=py.dict;
    for row=parameterRows
        name=char(study.ParamTable.Name(row));
        distribution=radia.optuna.internal.DistributionCodec.decode( ...
            study.ParamTable.Kind(row),study.ParamTable.Distribution(row));
        value=parameterValue(study.ParamTable(row,:));
        params{name}=pythonValue(value);
        distributions{name}=pythonDistribution(optuna,distribution);
    end
    objectiveRows=sortrows(study.ObjectiveTable( ...
        study.ObjectiveTable.TrialNumber==number,:),"ObjectiveIndex");
    if height(objectiveRows)<options.objective_index
        error("radia:optuna:ImportanceObjective", ...
            "Trial %d does not contain objective %d.", ...
            number,options.objective_index);
    end
    objectiveValue=objectiveRows.Value(options.objective_index);
    if ~isempty(options.target)
        trials=study.get_trials();
        numbersFromTrials=reshape([trials.Number],1,[]);
        source=trials(numbersFromTrials==number);
        if ~isscalar(source)
            error("radia:optuna:ImportanceTarget", ...
                "Could not resolve trial %d for the target callback.",number);
        end
        objectiveValue=options.target(source);
        if ~isnumeric(objectiveValue) || ~isscalar(objectiveValue)
            error("radia:optuna:ImportanceTarget", ...
                "The target callback must return a numeric scalar.");
        end
    end
    frozen=optuna.trial.create_trial(pyargs( ...
        "value",double(objectiveValue), ...
        "params",params,"distributions",distributions));
    pythonStudy.add_trial(frozen);
end

switch lower(options.evaluator)
    case "fanova"
        evaluator=importance.FanovaImportanceEvaluator(pyargs( ...
            "n_trees",int64(options.n_trees), ...
            "max_depth",int64(options.max_depth), ...
            "seed",int64(options.seed)));
    case {"mdi","mean_decrease_impurity"}
        evaluator=importance. ...
            MeanDecreaseImpurityImportanceEvaluator(pyargs( ...
            "n_trees",int64(options.n_trees), ...
            "max_depth",int64(options.max_depth), ...
            "seed",int64(options.seed)));
    case {"ped_anova","ped-anova"}
        pairs={"target_quantile",options.target_quantile, ...
            "region_quantile",options.region_quantile, ...
            "evaluate_on_local",options.evaluate_on_local};
        if isfinite(options.baseline_quantile)
            pairs=[pairs,{"baseline_quantile",options.baseline_quantile}];
        end
        evaluator=importance.PedAnovaImportanceEvaluator(pyargs(pairs{:}));
    otherwise
        error("radia:optuna:ImportanceEvaluator", ...
            "evaluator must be fanova, mdi, or ped_anova.");
end
pairs={"evaluator",evaluator,"normalize",options.normalize};
if ~isempty(options.params)
    pairs=[pairs,{"params",py.list(cellstr(options.params))}];
end
raw=importance.get_param_importances( ...
    pythonStudy,pyargs(pairs{:}));
keys=cell(py.list(raw.keys()));
names=strings(numel(keys),1);
values=zeros(numel(keys),1);
for index=1:numel(keys)
    names(index)=string(keys{index});
    values(index)=double(raw{keys{index}});
end
result=table(names,values,'VariableNames',{'Parameter','Importance'});
end

function optuna=importPinnedOptuna()
environment=pyenv;
if environment.Status=="NotLoaded"
    environment=pyenv(ExecutionMode="InProcess");
end
if environment.ExecutionMode~="InProcess" || ...
        ~startsWith(string(environment.Version),"3.12")
    error("radia:optuna:ImportancePython", ...
        "Parameter importance requires in-process Python 3.12.");
end
optuna=py.importlib.import_module("optuna");
if string(py.builtins.getattr(optuna,"__version__"))~="4.9.0"
    error("radia:optuna:ImportancePython", ...
        "Parameter importance requires optuna==4.9.0.");
end
end

function value=parameterValue(row)
if ~isnan(row.ValueNumeric)
    value=row.ValueNumeric;
else
    value=jsondecode(row.ValueText);
end
end

function result=pythonValue(value)
if isstring(value) || ischar(value)
    result=py.builtins.str(char(string(value)));
elseif islogical(value) && isscalar(value)
    result=py.builtins.bool(value);
elseif isnumeric(value) && isscalar(value)
    result=py.builtins.float(double(value));
else
    error("radia:optuna:ImportanceParameter", ...
        "Only scalar string, logical, and numeric parameters are supported.");
end
end

function result=pythonDistribution(optuna,distribution)
switch distribution.kind
    case "float"
        pairs={"log",logical(distribution.log)};
        if isfinite(distribution.step)
            pairs=[pairs,{"step",distribution.step}];
        end
        result=optuna.distributions.FloatDistribution( ...
            distribution.low,distribution.high,pyargs(pairs{:}));
    case "integer"
        result=optuna.distributions.IntDistribution( ...
            int64(distribution.low),int64(distribution.high),pyargs( ...
            "log",logical(distribution.log), ...
            "step",int64(distribution.step)));
    case "categorical"
        choices=cell(1,numel(distribution.choices));
        for index=1:numel(choices)
            choices{index}=pythonValue( ...
                radia.optuna.internal.DistributionCodec.choiceAt( ...
                distribution.choices,index));
        end
        result=optuna.distributions.CategoricalDistribution(py.list(choices));
    otherwise
        error("radia:optuna:ImportanceParameter", ...
            "Unsupported distribution kind '%s'.",distribution.kind);
end
end
