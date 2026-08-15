classdef FrozenTrial
    %FROZENTRIAL Immutable snapshot passed to callbacks and addTrial().

    properties (SetAccess=private)
        Number (1,1) double = -1
        State (1,1) string = "COMPLETE"
        Value (1,1) double = NaN
        Values double = NaN
        Params struct = struct()
        Distributions struct = struct()
        IntermediateValues table = table()
        UserAttrs struct = struct()
        SystemAttrs struct = struct()
        Constraints double = zeros(1,0)
        ConstraintPresent (1,1) logical = false
        DatetimeStart datetime = NaT
        DatetimeComplete datetime = NaT
        Duration double = NaN
        ErrorMessage (1,1) string = ""
    end

    methods
        function obj=FrozenTrial(options)
            arguments
                options.Number (1,1) double = -1
                options.State (1,1) string = "COMPLETE"
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
            state=upper(options.State);
            if ~ismember(state,["WAITING","RUNNING","COMPLETE","PRUNED","FAIL"])
                error("radia:optuna:TrialState","Unknown trial state '%s'.",state);
            end
            values=reshape(double(options.Values),1,[]);
            if all(isnan(values)) && ~all(isnan(options.Value))
                values=reshape(double(options.Value),1,[]);
            end
            if isempty(values), values=NaN; end
            obj.Number=options.Number;
            obj.State=state;
            obj.Values=values;
            obj.Value=values(1);
            obj.Params=options.Params;
            obj.Distributions=options.Distributions;
            obj.IntermediateValues=options.IntermediateValues;
            obj.UserAttrs=options.UserAttrs;
            obj.SystemAttrs=options.SystemAttrs;
            obj.Constraints=reshape(double(options.Constraints),1,[]);
            obj.ConstraintPresent=options.ConstraintPresent || ...
                ~isempty(options.Constraints);
            obj.DatetimeStart=options.DatetimeStart;
            obj.DatetimeComplete=options.DatetimeComplete;
            obj.ErrorMessage=options.ErrorMessage;
            if ~isnat(obj.DatetimeStart) && ~isnat(obj.DatetimeComplete)
                obj.Duration=seconds(obj.DatetimeComplete-obj.DatetimeStart);
            end
        end

        function step=last_step(obj)
            if isempty(obj.IntermediateValues)
                step=[];
            else
                step=max(obj.IntermediateValues.Step);
            end
        end
    end
end
