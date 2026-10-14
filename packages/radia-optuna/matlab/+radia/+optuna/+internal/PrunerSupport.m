classdef PrunerSupport
    %PRUNERSUPPORT Shared Optuna-compatible intermediate-value helpers.

    methods (Static)
        function result=isFirstInInterval(trial,warmupSteps,intervalSteps)
            if isempty(trial.IntermediateValues)
                result=false;
                return
            end
            step=max(trial.IntermediateValues.Step);
            if step<warmupSteps
                result=false;
                return
            end
            nearest=floor((step-warmupSteps)/intervalSteps)* ...
                intervalSteps+warmupSteps;
            previous=trial.IntermediateValues.Step( ...
                trial.IntermediateValues.Step~=step);
            if isempty(previous)
                secondLast=-1;
            else
                secondLast=max(previous);
            end
            result=secondLast<nearest;
        end

        function value=bestIntermediate(trial,direction)
            values=trial.IntermediateValues.Value;
            values=values(~isnan(values));
            if isempty(values)
                value=NaN;
            elseif direction=="maximize"
                value=max(values);
            else
                value=min(values);
            end
        end

        function value=latestIntermediate(trial)
            if isempty(trial.IntermediateValues)
                value=NaN;
                return
            end
            [~,index]=max(trial.IntermediateValues.Step);
            value=trial.IntermediateValues.Value(index);
        end

        function value=percentile(values,percent)
            % NumPy's default nanpercentile uses linear interpolation.
            values=sort(reshape(double(values(~isnan(values))),[],1));
            if isempty(values)
                value=NaN;
                return
            end
            position=1+(percent/100)*(numel(values)-1);
            lower=floor(position);
            upper=ceil(position);
            if lower==upper
                value=values(lower);
            else
                weight=position-lower;
                value=(1-weight)*values(lower)+weight*values(upper);
            end
        end

        function validateSchedule(warmupSteps,intervalSteps)
            if warmupSteps<0 || warmupSteps~=floor(warmupSteps)
                error("radia:optuna:PrunerWarmup", ...
                    "NWarmupSteps must be a nonnegative integer.");
            end
            if intervalSteps<1 || intervalSteps~=floor(intervalSteps)
                error("radia:optuna:PrunerInterval", ...
                    "IntervalSteps must be a positive integer.");
            end
        end
    end
end
