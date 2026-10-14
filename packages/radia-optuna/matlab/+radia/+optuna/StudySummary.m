classdef StudySummary
    %STUDYSUMMARY Optuna 5.0 study metadata and aggregate results.

    properties (SetAccess=private)
        study_name (1,1) string = ""
        best_trial = []
        user_attrs (1,1) struct = struct()
        n_trials (1,1) double = 0
        datetime_start datetime = NaT
    end

    properties (Dependent, SetAccess=private)
        direction
        directions
    end

    properties (Access=private)
        DirectionValues = []
        StudyId (1,1) double = -1
    end

    methods
        function obj=StudySummary(studyName,direction,bestTrial,userAttrs, ...
                nTrials,datetimeStart,studyId,options)
            arguments
                studyName (1,1) string
                direction = []
                bestTrial = []
                userAttrs (1,1) struct = struct()
                nTrials (1,1) double {mustBeInteger,mustBeNonnegative} = 0
                datetimeStart datetime = NaT
                studyId (1,1) double {mustBeInteger} = -1
                options.directions = []
            end
            suppliedDirection=~isempty(direction);
            suppliedDirections=~isempty(options.directions);
            if suppliedDirection==suppliedDirections
                error("radia:optuna:StudySummaryDirection", ...
                    "Specify exactly one of direction and directions.");
            end
            if suppliedDirections
                values=radia.optuna.StudyDirection.from(options.directions);
            else
                values=radia.optuna.StudyDirection.from(direction);
            end
            obj.study_name=studyName;
            obj.DirectionValues=reshape(values,1,[]);
            obj.best_trial=bestTrial;
            obj.user_attrs=userAttrs;
            obj.n_trials=nTrials;
            obj.datetime_start=datetimeStart;
            obj.StudyId=studyId;
        end

        function value=get.direction(obj)
            if numel(obj.DirectionValues)>1
                error("radia:optuna:MultiObjectiveDirection", ...
                    "This attribute is not available during multi-objective optimization.");
            end
            value=obj.DirectionValues(1);
        end

        function value=get.directions(obj)
            value=obj.DirectionValues;
        end

        function result=eq(left,right)
            if ~isa(right,"radia.optuna.StudySummary")
                result=false;
                return
            end
            result=left.StudyId==right.StudyId && ...
                left.study_name==right.study_name && ...
                isequaln(left.DirectionValues,right.DirectionValues) && ...
                isequaln(left.best_trial,right.best_trial) && ...
                isequaln(left.user_attrs,right.user_attrs) && ...
                left.n_trials==right.n_trials && ...
                isequaln(left.datetime_start,right.datetime_start);
        end

        function result=lt(left,right)
            result=left.StudyId<right.StudyId;
        end

        function result=le(left,right)
            result=left.StudyId<=right.StudyId;
        end
    end
end
