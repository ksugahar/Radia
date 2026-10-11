classdef MaxTrialsCallback
    %MAXTRIALSCALLBACK Stop after a selected number of settled trials.

    properties (SetAccess=private)
        NTrials (1,1) double
        States string
    end

    methods
        function obj=MaxTrialsCallback(nTrials,options)
            arguments
                nTrials (1,1) double {mustBeInteger,mustBeNonnegative}
                options.States string = "COMPLETE"
            end
            states=upper(options.States);
            allowed=["WAITING","RUNNING","COMPLETE","PRUNED","FAIL"];
            if any(~ismember(states,allowed))
                error("radia:optuna:TrialState", ...
                    "States contains an unsupported trial state.");
            end
            obj.NTrials=nTrials;
            obj.States=states;
        end

        function invoke(obj,study,~)
            if isempty(obj.States)
                count=height(study.TrialTable);
            else
                count=sum(ismember(study.TrialTable.State,obj.States));
            end
            if count>=obj.NTrials, study.stop(); end
        end

        function callback=callback(obj)
            callback=@obj.invoke;
        end
    end
end
