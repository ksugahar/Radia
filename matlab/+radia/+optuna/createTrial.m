function trial=createTrial(options)
%CREATETRIAL Construct a FrozenTrial for Study.addTrial().
arguments
    options.State = "COMPLETE"
    options.Value double = NaN
    options.Values double = NaN
    options.Params (1,1) struct = struct()
    options.Distributions (1,1) struct = struct()
    options.IntermediateValues table = table()
    options.UserAttrs (1,1) struct = struct()
    options.SystemAttrs (1,1) struct = struct()
    options.Constraints = struct()
    options.ConstraintPresent (1,1) logical = false
    options.DatetimeStart datetime = NaT
    options.DatetimeComplete datetime = NaT
    options.ErrorMessage (1,1) string = ""
end
state=radia.optuna.TrialState.toStorage(options.State);
[constraintNames,constraintValues]=normalizeConstraints(options.Constraints);
startTime=options.DatetimeStart;
completeTime=options.DatetimeComplete;
if state=="RUNNING" && isnat(startTime)
    startTime=datetime("now","TimeZone","local");
elseif ismember(state,["COMPLETE","PRUNED","FAIL"])
    if isnat(startTime), startTime=datetime("now","TimeZone","local"); end
    if isnat(completeTime), completeTime=startTime; end
end
trial=radia.optuna.FrozenTrial(State=state,Value=options.Value, ...
    Values=options.Values,Params=options.Params, ...
    Distributions=options.Distributions, ...
    IntermediateValues=options.IntermediateValues, ...
    UserAttrs=options.UserAttrs,SystemAttrs=options.SystemAttrs, ...
    Constraints=constraintValues,ConstraintNames=constraintNames, ...
    ConstraintPresent=options.ConstraintPresent || ~isempty(constraintNames), ...
    DatetimeStart=startTime, ...
    DatetimeComplete=completeTime, ...
    ErrorMessage=options.ErrorMessage);
end

function [names,values]=normalizeConstraints(source)
if isempty(source)
    names=strings(1,0);
    values=zeros(1,0);
elseif isa(source,"dictionary")
    names=reshape(string(keys(source)),1,[]);
    values=zeros(1,numel(names));
    for index=1:numel(names)
        values(index)=double(source(names(index)));
    end
elseif isstruct(source) && isscalar(source)
    names=reshape(string(fieldnames(source)),1,[]);
    values=zeros(1,numel(names));
    for index=1:numel(names)
        values(index)=double(source.(names(index)));
    end
elseif isnumeric(source)
    values=reshape(double(source),1,[]);
    names=string(0:numel(values)-1);
else
    error("radia:optuna:ConstraintType", ...
        "Constraints must be a dictionary, scalar struct, or numeric vector.");
end
if any(isnan(values)) || any(~isreal(values)) || ...
        numel(unique(names))~=numel(names)
    error("radia:optuna:Constraints", ...
        "Constraint names must be unique and values must not contain NaN.");
end
end
