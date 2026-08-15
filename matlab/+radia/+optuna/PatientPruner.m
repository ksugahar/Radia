classdef PatientPruner < handle
    %PATIENTPRUNER Prune after a patience window without improvement.

    properties (SetAccess=private)
        WrappedPruner
        Patience (1,1) double
        MinDelta (1,1) double
    end

    methods
        function obj=PatientPruner(wrappedPruner,options)
            arguments
                wrappedPruner = []
                options.Patience (1,1) double = 0
                options.MinDelta (1,1) double = 0
            end
            if options.Patience<0 || options.Patience~=floor(options.Patience)
                error("radia:optuna:PrunerPatience", ...
                    "Patience must be a nonnegative integer.");
            end
            if options.MinDelta<0 || ~isfinite(options.MinDelta)
                error("radia:optuna:PrunerDelta", ...
                    "MinDelta must be finite and nonnegative.");
            end
            if ~isempty(wrappedPruner) && ...
                    ~ismethod(wrappedPruner,"shouldPrune")
                error("radia:optuna:PrunerWrapper", ...
                    "WrappedPruner must implement shouldPrune.");
            end
            obj.WrappedPruner=wrappedPruner;
            obj.Patience=options.Patience;
            obj.MinDelta=options.MinDelta;
        end

        function decision=shouldPrune(obj,study,trial)
            decision=false;
            count=height(trial.IntermediateValues);
            if count<=obj.Patience+1
                return
            end
            ordered=sortrows(trial.IntermediateValues,"Step");
            split=count-obj.Patience-1;
            before=ordered.Value(1:split);
            after=ordered.Value((split+1):end);
            before=before(~isnan(before));
            after=after(~isnan(after));
            if isempty(before) || isempty(after)
                return
            end
            if study.Directions(1)=="minimize"
                stalled=min(before)+obj.MinDelta<min(after);
            else
                stalled=max(before)-obj.MinDelta>max(after);
            end
            if ~stalled
                return
            end
            if isempty(obj.WrappedPruner)
                decision=true;
            else
                decision=obj.WrappedPruner.shouldPrune(study,trial);
            end
        end
    end
end
