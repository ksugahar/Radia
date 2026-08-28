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
    options.Constraints double = zeros(1,0)
    options.ConstraintPresent (1,1) logical = false
    options.DatetimeStart datetime = NaT
    options.DatetimeComplete datetime = NaT
    options.ErrorMessage (1,1) string = ""
end
state=radia.optuna.TrialState.toStorage(options.State);
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
    Constraints=options.Constraints, ...
    ConstraintPresent=options.ConstraintPresent, ...
    DatetimeStart=startTime, ...
    DatetimeComplete=completeTime, ...
    ErrorMessage=options.ErrorMessage);
end
