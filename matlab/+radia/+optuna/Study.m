classdef Study < handle
    %STUDY Table-backed implementation of a verified Optuna 4.9 subset.

    properties (SetAccess=private)
        Name (1,1) string
        Directions string
        StoragePath (1,1) string = ""
        Sampler
        Pruner
        AutoSave (1,1) logical = true
        IntermediateTable table
        UserAttrTable table
        SystemAttrTable table
        ConstraintTable table
        ConstraintCountTable table
        SamplerStateTable table
        QueueParamTable table
        ProgressFcn = []
        UserAttrs struct = struct()
        SystemAttrs struct = struct()
        MetricNames string = strings(1,0)
    end

    properties (Dependent, SetAccess=private)
        % Public compatibility views.  The optimizer's hot path uses the
        % column stores below and materializes these tables only on demand.
        TrialTable
        ParamTable
        ObjectiveTable
    end

    properties (Access=private)
        NextTrialNumber (1,1) double = 0
        StopRequested (1,1) logical = false
        InOptimize (1,1) logical = false
        TrialNumberData double = zeros(0,1)
        TrialStateData string = strings(0,1)
        TrialValueData double = zeros(0,1)
        TrialStartTimeData
        TrialEndTimeData
        TrialDurationData double = zeros(0,1)
        TrialParamsData cell = cell(0,1)
        TrialIntermediateData cell = cell(0,1)
        TrialErrorData string = strings(0,1)
        ParamTrialNumberData double = zeros(0,1)
        ParamNameData string = strings(0,1)
        ParamKindData string = strings(0,1)
        ParamValueNumericData double = zeros(0,1)
        ParamValueTextData string = strings(0,1)
        ParamDistributionData string = strings(0,1)
        ObjectiveTrialNumberData double = zeros(0,1)
        ObjectiveIndexData double = zeros(0,1)
        ObjectiveValueData double = zeros(0,1)
        TrialTableCache = table()
        ParamTableCache = table()
        ObjectiveTableCache = table()
        TrialTableDirty (1,1) logical = true
        ParamTableDirty (1,1) logical = true
        ObjectiveTableDirty (1,1) logical = true
        SamplerHasBeforeTrial (1,1) logical = false
        SamplerHasAfterTrial (1,1) logical = false
    end

    properties (Constant, Access=private)
        StorageSchema = "radia.optuna.study"
        StorageVersion = 4
    end

    methods
        function obj = Study(options)
            arguments
                options.Name (1,1) string = ""
                options.StoragePath (1,1) string = ""
                options.Directions = "minimize"
                options.Sampler = []
                options.Pruner = []
                options.AutoSave (1,1) logical = true
                options.ProgressFcn = []
            end
            obj.Name = options.Name;
            if strlength(obj.Name) == 0
                obj.Name = obj.generateStudyName();
            end
            obj.Directions = ...
                radia.optuna.StudyDirection.toStorage(options.Directions);
            if isempty(obj.Directions) || any(~ismember(obj.Directions, ["minimize", "maximize"]))
                error("radia:optuna:Direction", ...
                    "Directions must contain only 'minimize' or 'maximize'.");
            end
            obj.StoragePath = options.StoragePath;
            obj.AutoSave = options.AutoSave;
            if ~isempty(options.ProgressFcn) && ~isa(options.ProgressFcn, "function_handle")
                error("radia:optuna:ProgressFcn", ...
                    "ProgressFcn must be a function handle.");
            end
            obj.ProgressFcn = options.ProgressFcn;
            obj.Sampler = options.Sampler;
            obj.Pruner = options.Pruner;
            if isempty(obj.Sampler)
                if isscalar(obj.Directions)
                    obj.Sampler = radia.optuna.TPESampler( ...
                        NStartupTrials=10);
                else
                    % Optuna 4.9 selects NSGA-II for a multi-objective
                    % study when no sampler is supplied.
                    obj.Sampler = radia.optuna.NSGAIISampler( ...
                        PopulationSize=50);
                end
            end
            if isempty(obj.Pruner)
                obj.Pruner = radia.optuna.MedianPruner();
            end
            obj.SamplerHasBeforeTrial=ismethod(obj.Sampler,"beforeTrial");
            obj.SamplerHasAfterTrial=ismethod(obj.Sampler,"afterTrial");
            obj.initializeTables();
            if strlength(obj.StoragePath) > 0 && ...
                    (isfile(obj.StoragePath) || isfile(obj.backupStoragePath()))
                obj.loadState();
            elseif strlength(obj.StoragePath) > 0
                % Upstream create_study registers even an empty study in
                % its storage immediately.
                obj.persist();
            end
        end

        function trial = ask(obj, fixedDistributions)
            if nargin < 2
                fixedDistributions = struct();
            end
            if ~isstruct(fixedDistributions) || ~isscalar(fixedDistributions)
                error("radia:optuna:FixedDistributions", ...
                    "fixedDistributions must be a scalar struct of distributions.");
            end
            waiting=find(obj.TrialStateData=="WAITING",1);
            if isempty(waiting)
                trial = radia.optuna.Trial(obj, obj.NextTrialNumber);
                obj.NextTrialNumber = obj.NextTrialNumber + 1;
                obj.appendTrialRecord(trial.Number,"RUNNING",NaN, ...
                    trial.startTimeSerial(),NaN,NaN, ...
                    trial.Params,trial.IntermediateValues,"");
            else
                number=obj.TrialNumberData(waiting);
                trial=radia.optuna.Trial(obj,number);
                claimedAt=trial.startTimeSerial();
                trial.restoreSnapshot(obj.freezeTrial(number));
                trial.markRunning(claimedAt);
                rows=obj.QueueParamTable.TrialNumber==number;
                queuedParams=obj.QueueParamTable(rows,:);
                trial.setFixedParameters(queuedParams.Name,queuedParams.Value);
                obj.TrialStateData(waiting)="RUNNING";
                obj.TrialStartTimeData(waiting)=trial.startTimeSerial();
                obj.TrialEndTimeData(waiting)=NaN;
                obj.TrialDurationData(waiting)=NaN;
                obj.TrialErrorData(waiting)="";
                obj.TrialTableDirty=true;
            end
            if obj.SamplerHasBeforeTrial
                obj.Sampler.beforeTrial(obj, trial);
            end
            if ~isempty(fieldnames(fixedDistributions))
                obj.suggestFixedDistributions(trial, fixedDistributions);
            end
            if ~isempty(waiting) || ~isempty(fieldnames(fixedDistributions))
                obj.updateTrialSnapshot(trial);
            end
            if obj.AutoSave && strlength(obj.StoragePath)>0
                obj.save();
            end
            if ~isempty(obj.ProgressFcn)
                obj.reportProgress("TRIAL_STARTED", trial);
            end
        end

        function frozen = tell(obj, trial, value, options)
            arguments
                obj
                trial
                value double = double.empty(1,0)
                options.State = string.empty
                options.SkipIfFinished (1,1) logical = false
            end
            trial=obj.resolveTrial(trial);
            if obj.skipFinishedTrial(trial, options.SkipIfFinished)
                if nargout>0
                    frozen=obj.freezeTrial(trial.Number);
                end
                return
            end
            if isempty(options.State)
                state="";
            else
                state=radia.optuna.TrialState.toStorage(options.State);
            end
            if ~ismember(state, ["", "COMPLETE", "PRUNED", "FAIL"])
                error("radia:optuna:TrialState", ...
                    "tell state must be COMPLETE, PRUNED, FAIL, or empty.");
            end
            value = reshape(double(value), 1, []);
            if state == ""
                if isempty(value) || any(isnan(value))
                    warning("radia:optuna:InvalidObjectiveValue", ...
                        "The objective value is missing or NaN; the trial is marked FAIL.");
                    state = "FAIL";
                else
                    state = "COMPLETE";
                end
            end
            if state == "COMPLETE"
                if numel(value) ~= numel(obj.Directions) || any(isnan(value))
                    error("radia:optuna:Value", ...
                        "A COMPLETE trial requires one non-NaN value per direction.");
                end
            elseif ~isempty(value)
                error("radia:optuna:Value", ...
                    "Values cannot be specified for a PRUNED or FAIL trial.");
            elseif state == "PRUNED" && ~isempty(trial.IntermediateValues)
                ordered = sortrows(trial.IntermediateValues, "Step");
                value = ordered.Value(end);
            else
                value = NaN;
            end
            obj.finishTrial(trial, state, value, "");
            if nargout>0
                frozen=obj.freezeTrial(trial.Number);
            end
        end

        function fail(obj, trial, message, options)
            arguments
                obj
                trial (1,1) radia.optuna.Trial
                message (1,1) string
                options.SkipIfFinished (1,1) logical = false
            end
            if obj.skipFinishedTrial(trial, options.SkipIfFinished)
                return
            end
            obj.finishTrial(trial, "FAIL", NaN, message);
        end

        function enqueueTrial(obj,params,options)
            arguments
                obj
                params (1,1) struct
                options.UserAttrs (1,1) struct = struct()
                options.SkipIfExists (1,1) logical = false
            end
            if options.SkipIfExists && obj.containsParameterSet(params)
                return
            end
            number=obj.NextTrialNumber;
            obj.NextTrialNumber=number+1;
            notATime=NaN;
            obj.appendTrialRecord(number,"WAITING",NaN,notATime, ...
                notATime,NaN,struct(),table(),"");
            names=string(fieldnames(params));
            for index=1:numel(names)
                value=params.(names(index));
                obj.QueueParamTable(end+1,:)={number,names(index),{value}};
            end
            attrNames=string(fieldnames(options.UserAttrs));
            for index=1:numel(attrNames)
                obj.UserAttrTable(end+1,:)={number,attrNames(index), ...
                    string(jsonencode(options.UserAttrs.(attrNames(index))))};
            end
            obj.persist();
        end

        function enqueue_trial(obj,params,options)
            arguments
                obj
                params (1,1) struct
                options.UserAttrs (1,1) struct = struct()
                options.SkipIfExists (1,1) logical = false
            end
            obj.enqueueTrial(params,UserAttrs=options.UserAttrs, ...
                SkipIfExists=options.SkipIfExists);
        end

        function results = optimize(obj, objective, nTrials, options)
            arguments
                obj
                objective (1,1) function_handle
                nTrials (1,1) double {mustBeInteger, mustBeNonnegative}
                options.Timeout (1,1) double = Inf
                options.Callbacks = {}
                options.Catch string = strings(1,0)
            end
            if ~(isinf(options.Timeout) || ...
                    (isfinite(options.Timeout) && options.Timeout>=0))
                error("radia:optuna:Timeout", ...
                    "Timeout must be nonnegative or Inf.");
            end
            if obj.InOptimize
                error("radia:optuna:NestedOptimize", ...
                    "Nested Study.optimize calls are not supported.");
            end
            callbacks=obj.normalizeCallbacks(options.Callbacks);
            obj.InOptimize=true;
            obj.StopRequested=false;
            cleanup=onCleanup(@()obj.finishOptimize());
            started=tic;
            for k = 1:nTrials
                if obj.StopRequested || toc(started)>=options.Timeout
                    break
                end
                trial = obj.ask();
                try
                    value = objective(trial);
                    if trial.State == "RUNNING"
                        obj.tell(trial, value);
                    end
                catch exception
                    if exception.identifier == "radia:optuna:TrialPruned"
                        if trial.State == "RUNNING"
                            obj.tell(trial, State="PRUNED");
                        end
                    else
                        if trial.State == "RUNNING"
                            obj.fail(trial, exception.message);
                        end
                        caught=any(options.Catch==string(exception.identifier)) || ...
                            any(options.Catch=="*");
                        if ~caught
                            rethrow(exception);
                        end
                    end
                end
                if ~isempty(callbacks)
                    frozen=obj.freezeTrial(trial.Number);
                    for callback=callbacks
                        callback{1}(obj,frozen);
                    end
                end
            end
            results = obj.TrialTable;
            clear cleanup
        end

        function stop(obj)
            if ~obj.InOptimize
                error("radia:optuna:StopOutsideOptimize", ...
                    "Study.stop() may only be called from an objective or callback.");
            end
            obj.StopRequested=true;
        end

        function results = optimizeSimulink(obj, runner, nTrials)
            arguments
                obj
                runner (1,1) radia.optuna.SimulinkRunner
                nTrials (1,1) double {mustBeInteger, mustBeNonnegative}
            end
            results = runner.optimize(obj, nTrials);
        end

        function value = sampleFloat(obj, trial, name, low, high, options)
            value = obj.Sampler.sampleFloat(obj, trial, name, low, high, options);
        end

        function value = sampleInteger(obj, trial, name, low, high)
            value = obj.Sampler.sampleInteger(obj, trial, name, low, high);
        end

        function value = sampleCategorical(obj, trial, name, choices)
            value = obj.Sampler.sampleCategorical(obj, trial, name, choices);
        end

        function values = sampleJoint(obj, trial, names, lows, highs, options)
            arguments
                obj
                trial (1,1) radia.optuna.Trial
                names (1,:) string
                lows (1,:) double
                highs (1,:) double
                options.Log (1,:) logical = false(1, numel(names))
            end
            if numel(names) ~= numel(lows) || numel(names) ~= numel(highs)
                error("radia:optuna:JointShape", ...
                    "Joint names and bounds must have the same length.");
            end
            if numel(options.Log) ~= numel(names)
                error("radia:optuna:JointShape", ...
                    "Joint Log flags must have the same length as names.");
            end
            if ismethod(obj.Sampler, "sampleJoint")
                values = obj.Sampler.sampleJoint(obj, trial, names, lows, ...
                    highs, options);
            else
                values = zeros(1, numel(names));
                for index = 1:numel(names)
                    values(index) = obj.sampleFloat(trial, names(index), ...
                        lows(index), highs(index), struct( ...
                        "Log", options.Log(index), "Step", NaN));
                end
            end
            values = reshape(double(values), 1, []);
            if numel(values) ~= numel(names) || any(~isfinite(values))
                error("radia:optuna:JointValue", ...
                    "Joint sampler returned invalid values.");
            end
        end

        function setUserAttr(obj, name, value)
            arguments
                obj
                name (1,1) string
                value
            end
            key = matlab.lang.makeValidName(name);
            obj.UserAttrs.(key) = value;
            obj.persist();
        end

        function set_user_attr(obj, name, value)
            obj.setUserAttr(name, value);
        end

        function value = user_attrs(obj)
            value = obj.UserAttrs;
        end

        function setSystemAttr(obj,name,value)
            arguments
                obj
                name (1,1) string
                value
            end
            obj.SystemAttrs.(matlab.lang.makeValidName(name))=value;
            obj.persist();
        end

        function set_system_attr(obj,name,value)
            obj.setSystemAttr(name,value);
        end

        function value = system_attrs(obj)
            value=obj.SystemAttrs;
        end

        function setMetricNames(obj,names)
            names=reshape(string(names),1,[]);
            if numel(names)~=numel(obj.Directions) || ...
                    any(strlength(names)==0) || numel(unique(names))~=numel(names)
                error("radia:optuna:MetricNames", ...
                    "Metric names must be unique, nonempty, and match Directions.");
            end
            obj.MetricNames=names;
            obj.persist();
        end

        function set_metric_names(obj,names)
            obj.setMetricNames(names);
        end

        function value = metric_names(obj)
            value=obj.MetricNames;
        end

        function result = bestTrial(obj)
            if numel(obj.Directions) ~= 1
                error("radia:optuna:MultiObjectiveBest", ...
                    "Use paretoFront() for a multi-objective study.");
            end
            complete = obj.TrialTable.State == "COMPLETE" & ...
                ~isnan(obj.TrialTable.Value);
            if ~any(complete)
                result = obj.TrialTable([],:);
                return;
            end
            rows = find(complete);
            [feasible, constrained] = obj.feasibleTrials( ...
                obj.TrialTable.TrialNumber(rows));
            rows = rows(feasible);
            if constrained && isempty(rows)
                error("radia:optuna:NoFeasibleTrial", ...
                    "The study has completed trials, but none has a " + ...
                    "complete set of constraints satisfying c <= 0.");
            end
            values = obj.TrialTable.Value(rows);
            if obj.Directions(1) == "minimize"
                [~, index] = min(values);
            else
                [~, index] = max(values);
            end
            result = obj.TrialTable(rows(index),:);
        end

        function value = bestValue(obj)
            %BESTVALUE Return the remembered single-objective best value.
            best = obj.bestTrial();
            if isempty(best)
                value = NaN;
            else
                value = best.Value(1);
            end
        end

        function params = bestParams(obj)
            %BESTPARAMS Return the parameters associated with the best trial.
            best = obj.bestTrial();
            if isempty(best)
                params = struct();
            else
                params = best.Params{1};
                if isempty(params)
                    params = struct();
                end
            end
        end

        function result = bestSolution(obj)
            %BESTSOLUTION Return a reloadable best-solution snapshot.
            best = obj.bestTrial();
            if isempty(best)
                result = struct( ...
                    "available", false, ...
                    "trial_number", NaN, ...
                    "value", NaN, ...
                    "params", struct(), ...
                    "state", "NONE", ...
                    "end_time", NaT);
                return;
            end
            result = struct( ...
                "available", true, ...
                "trial_number", best.TrialNumber(1), ...
                "value", best.Value(1), ...
                "params", best.Params{1}, ...
                "state", best.State(1), ...
                "end_time", best.EndTime(1));
        end

        function result = paretoFront(obj)
            complete = obj.TrialTable.State == "COMPLETE";
            trialNumbers = obj.TrialTable.TrialNumber(complete);
            if isempty(trialNumbers)
                result = obj.emptyParetoFront();
                return;
            end
            feasible = obj.feasibleTrials(trialNumbers);
            trialNumbers = trialNumbers(feasible);
            if isempty(trialNumbers)
                result = obj.emptyParetoFront();
                return;
            end
            values = NaN(numel(trialNumbers), numel(obj.Directions));
            for k = 1:numel(trialNumbers)
                rows = obj.ObjectiveTable.TrialNumber == trialNumbers(k);
                objectiveRows = obj.ObjectiveTable(rows,:);
                for j = 1:height(objectiveRows)
                    values(k, objectiveRows.ObjectiveIndex(j)) = objectiveRows.Value(j);
                end
            end
            valid = all(~isnan(values),2);
            trialNumbers = trialNumbers(valid); values = values(valid,:);
            if isempty(trialNumbers)
                result = obj.emptyParetoFront();
                return;
            end
            dominated = false(size(trialNumbers));
            signs = ones(1,numel(obj.Directions));
            signs(obj.Directions == "maximize") = -1;
            normalized = values .* signs;
            for i = 1:size(normalized,1)
                dominated(i) = any(all(normalized <= normalized(i,:),2) & ...
                    any(normalized < normalized(i,:),2));
            end
            keep = ~dominated;
            params = cell(sum(keep),1);
            keptNumbers = trialNumbers(keep);
            for k = 1:numel(keptNumbers)
                row = obj.TrialTable.TrialNumber == keptNumbers(k);
                params{k} = obj.TrialTable.Params{find(row,1)};
            end
            result = table(keptNumbers, num2cell(values(keep,:),2), params, ...
                'VariableNames',{'TrialNumber','Values','Params'});
        end

        function result = trials(obj, states)
            if nargin < 2 || isempty(states)
                result = obj.TrialTable;
            else
                result = obj.TrialTable(ismember(obj.TrialTable.State, string(states)),:);
            end
        end

        function result = best_trial(obj)
            best = obj.bestTrial();
            if isempty(best)
                error("radia:optuna:NoCompletedTrials", ...
                    "No trials are completed yet.");
            end
            result = obj.freezeTrial(best.TrialNumber(1));
        end

        function value = best_value(obj)
            value = obj.bestValue();
        end

        function params = best_params(obj)
            params = obj.bestParams();
        end

        function result = best_solution(obj)
            result = obj.bestSolution();
        end

        function result = best_trials(obj)
            front=obj.paretoFront();
            result=obj.freezeTrials(front.TrialNumber);
        end

        function result = trials_dataframe(obj,options)
            %TRIALS_DATAFRAME Export upstream-shaped trial columns as a table.
            %   MATLAB tables have no pandas MultiIndex. When multi_index is
            %   true, the exact two-level upstream column labels are retained
            %   in result.Properties.UserData.column_levels while the table
            %   variables keep the flattened Optuna names for direct access.
            arguments
                obj
                options.attrs (1,:) string = ["number","value", ...
                    "datetime_start","datetime_complete","duration", ...
                    "params","user_attrs","system_attrs","state"]
                options.multi_index (1,1) logical = false
            end
            result=obj.buildTrialsDataframe( ...
                options.attrs,options.multi_index);
        end

        function value = direction(obj)
            if numel(obj.Directions)~=1
                error("radia:optuna:MultiObjectiveDirection", ...
                    "Use Directions for a multi-objective study.");
            end
            value=radia.optuna.StudyDirection.from(obj.Directions(1));
        end

        function result = get_trials(obj, states)
            if nargin < 2
                rows = obj.TrialTable;
            else
                rows = obj.trials(states);
            end
            result=obj.freezeTrials(rows.TrialNumber);
        end

        function values = intermediateValuesAtStep(obj, step)
            rows = obj.IntermediateTable.Step == step;
            completed = obj.TrialTable.State == "COMPLETE";
            completedNumbers = obj.TrialTable.TrialNumber(completed);
            completedRows = ismember( ...
                obj.IntermediateTable.TrialNumber, completedNumbers);
            values = obj.IntermediateTable.Value(rows & completedRows);
        end

        function recovered = recoverStaleRunning(obj, timeoutSeconds, options)
            %RECOVERSTALERUNNING Mark abandoned RUNNING trials as FAIL.
            %   This deliberately does not reconstruct a Trial handle or
            %   resubmit parameters. A scheduler may inspect the returned
            %   rows and explicitly ask for replacement work.
            arguments
                obj
                timeoutSeconds (1,1) double {mustBeFinite,mustBeNonnegative}
                options.Message (1,1) string = ""
            end
            nowTime = now; %#ok<TNOW1> hot-path serial timestamp
            ages = (nowTime-obj.TrialStartTimeData)*86400;
            rows = obj.TrialStateData == "RUNNING" & ...
                isfinite(ages) & ages >= timeoutSeconds;
            if ~any(rows)
                recovered = obj.TrialTable([],:);
                return
            end
            message = options.Message;
            if strlength(message) == 0
                message = "Recovered stale RUNNING trial after timeout.";
            end
            trialNumbers = obj.TrialNumberData(rows);
            obj.TrialStateData(rows) = "FAIL";
            obj.TrialValueData(rows) = NaN;
            obj.TrialEndTimeData(rows) = nowTime;
            obj.TrialDurationData(rows) = ages(rows);
            obj.TrialErrorData(rows) = message;
            keep=~ismember(obj.ObjectiveTrialNumberData,trialNumbers);
            obj.ObjectiveTrialNumberData=obj.ObjectiveTrialNumberData(keep);
            obj.ObjectiveIndexData=obj.ObjectiveIndexData(keep);
            obj.ObjectiveValueData=obj.ObjectiveValueData(keep);
            obj.TrialTableDirty=true;
            obj.ObjectiveTableDirty=true;
            obj.persist();
            recovered = obj.TrialTable(rows,:);
        end

        function recovered = recover_stale_running(obj, timeoutSeconds, options)
            arguments
                obj
                timeoutSeconds (1,1) double {mustBeFinite,mustBeNonnegative}
                options.Message (1,1) string = ""
            end
            recovered = obj.recoverStaleRunning( ...
                timeoutSeconds, Message=options.Message);
        end

        function save(obj)
            if strlength(obj.StoragePath) == 0
                error("radia:optuna:Storage", "StoragePath is empty.");
            end
            StudyData = struct( ...
                "Schema", obj.StorageSchema, ...
                "Version", obj.StorageVersion, ...
                "Name", obj.Name, ...
                "Directions", obj.Directions, ...
                "NextTrialNumber", obj.NextTrialNumber, ...
                "TrialTable", obj.TrialTable, ...
                "ParamTable", obj.ParamTable, ...
                "IntermediateTable", obj.IntermediateTable, ...
                "UserAttrTable", obj.UserAttrTable, ...
                "SystemAttrTable", obj.SystemAttrTable, ...
                "ConstraintTable", obj.ConstraintTable, ...
                "ConstraintCountTable", obj.ConstraintCountTable, ...
                "UserAttrs", obj.UserAttrs, ...
                "SystemAttrs", obj.SystemAttrs, ...
                "MetricNames", obj.MetricNames, ...
                "ObjectiveTable", obj.ObjectiveTable, ...
                "SamplerStateTable", obj.SamplerStateTable, ...
                "QueueParamTable",obj.QueueParamTable); %#ok<NASGU>
            folder = fileparts(obj.StoragePath);
            if strlength(folder) > 0 && ~isfolder(folder)
                mkdir(folder);
            end
            [targetPath, folder] = obj.resolvedStoragePath();
            temporaryPath = string(tempname(folder)) + ".mat";
            cleanup = onCleanup(@() obj.removeTemporaryFile(temporaryPath));
            builtin("save", temporaryPath, "StudyData", "-mat");
            obj.readStudyData(temporaryPath);
            obj.replaceStorageFile(temporaryPath, targetPath, ...
                targetPath + ".bak");
            clear cleanup
        end
    end

    methods
        function value=get.TrialTable(obj)
            if obj.TrialTableDirty
                value=table(obj.TrialNumberData,obj.TrialStateData, ...
                    obj.TrialValueData,obj.serialDatetimes(obj.TrialStartTimeData), ...
                    obj.serialDatetimes(obj.TrialEndTimeData),obj.TrialDurationData, ...
                    obj.TrialParamsData,obj.TrialIntermediateData, ...
                    obj.TrialErrorData, ...
                    'VariableNames',{'TrialNumber','State','Value', ...
                    'StartTime','EndTime','Duration_s','Params', ...
                    'IntermediateValues','ErrorMessage'});
                obj.TrialTableCache=value;
                obj.TrialTableDirty=false;
            else
                value=obj.TrialTableCache;
            end
        end

        function set.TrialTable(obj,value)
            obj.TrialNumberData=reshape(double(value.TrialNumber),[],1);
            obj.TrialStateData=reshape(string(value.State),[],1);
            obj.TrialValueData=reshape(double(value.Value),[],1);
            obj.TrialStartTimeData=reshape(datenum(value.StartTime),[],1); %#ok<DATNM>
            obj.TrialEndTimeData=reshape(datenum(value.EndTime),[],1); %#ok<DATNM>
            obj.TrialDurationData=reshape(double(value.Duration_s),[],1);
            obj.TrialParamsData=reshape(value.Params,[],1);
            obj.TrialIntermediateData=reshape(value.IntermediateValues,[],1);
            obj.TrialErrorData=reshape(string(value.ErrorMessage),[],1);
            obj.TrialTableCache=value;
            obj.TrialTableDirty=false;
        end

        function value=get.ParamTable(obj)
            if obj.ParamTableDirty
                value=table(obj.ParamTrialNumberData,obj.ParamNameData, ...
                    obj.ParamKindData,obj.ParamValueNumericData, ...
                    obj.ParamValueTextData,obj.ParamDistributionData, ...
                    'VariableNames',{'TrialNumber','Name','Kind', ...
                    'ValueNumeric','ValueText','Distribution'});
                obj.ParamTableCache=value;
                obj.ParamTableDirty=false;
            else
                value=obj.ParamTableCache;
            end
        end

        function set.ParamTable(obj,value)
            obj.ParamTrialNumberData=reshape(double(value.TrialNumber),[],1);
            obj.ParamNameData=reshape(string(value.Name),[],1);
            obj.ParamKindData=reshape(string(value.Kind),[],1);
            obj.ParamValueNumericData=reshape(double(value.ValueNumeric),[],1);
            obj.ParamValueTextData=reshape(string(value.ValueText),[],1);
            obj.ParamDistributionData=reshape(string(value.Distribution),[],1);
            obj.ParamTableCache=value;
            obj.ParamTableDirty=false;
        end

        function value=get.ObjectiveTable(obj)
            if obj.ObjectiveTableDirty
                value=table(obj.ObjectiveTrialNumberData, ...
                    obj.ObjectiveIndexData,obj.ObjectiveValueData, ...
                    'VariableNames',{'TrialNumber','ObjectiveIndex','Value'});
                obj.ObjectiveTableCache=value;
                obj.ObjectiveTableDirty=false;
            else
                value=obj.ObjectiveTableCache;
            end
        end

        function set.ObjectiveTable(obj,value)
            obj.ObjectiveTrialNumberData=reshape(double(value.TrialNumber),[],1);
            obj.ObjectiveIndexData=reshape(double(value.ObjectiveIndex),[],1);
            obj.ObjectiveValueData=reshape(double(value.Value),[],1);
            obj.ObjectiveTableCache=value;
            obj.ObjectiveTableDirty=false;
        end

    end

    methods (Hidden=true)

        function result=buildTrialsDataframe(obj,attrs,multiIndex)
            attrs=reshape(string(attrs),1,[]);
            if isempty(attrs)
                error("radia:optuna:TrialsDataframeAttrs", ...
                    "attrs must contain at least one FrozenTrial field.");
            end
            allowed=["number","value","values","datetime_start", ...
                "datetime_complete","duration","params","user_attrs", ...
                "system_attrs","state","intermediate_values", ...
                "distributions"];
            unknown=attrs(~ismember(attrs,allowed));
            if ~isempty(unknown)
                error("radia:optuna:TrialsDataframeAttribute", ...
                    "FrozenTrial has no public field '%s'.",unknown(1));
            end

            if isempty(obj.TrialNumberData)
                result=table();
                result.Properties.UserData=obj.dataframeMetadata( ...
                    strings(0,1),strings(0,1),logical(multiIndex));
                return
            end

            if numel(obj.Directions)>1 && any(attrs=="value")
                attrs(attrs=="value")="values";
            end
            columns=struct('Top',{},'Sub',{},'Values',{});
            for attr=attrs
                switch attr
                    case "number"
                        columns=obj.appendDataframeColumn(columns,attr,"", ...
                            obj.TrialNumberData);
                    case "value"
                        if isempty(obj.MetricNames)
                            columns=obj.appendDataframeColumn( ...
                                columns,attr,"",obj.TrialValueData);
                        else
                            columns=obj.dataframeObjectiveColumns( ...
                                columns,"value");
                        end
                    case "values"
                        columns=obj.dataframeObjectiveColumns( ...
                            columns,"values");
                    case "datetime_start"
                        columns=obj.appendDataframeColumn(columns,attr,"", ...
                            obj.serialDatetimes(obj.TrialStartTimeData));
                    case "datetime_complete"
                        columns=obj.appendDataframeColumn(columns,attr,"", ...
                            obj.serialDatetimes(obj.TrialEndTimeData));
                    case "duration"
                        columns=obj.appendDataframeColumn(columns,attr,"", ...
                            seconds(obj.TrialDurationData));
                    case "params"
                        columns=obj.dataframeParameterColumns( ...
                            columns,"params",false);
                    case "distributions"
                        columns=obj.dataframeParameterColumns( ...
                            columns,"distributions",true);
                    case "user_attrs"
                        columns=obj.dataframeAttributeColumns( ...
                            columns,"user_attrs",obj.UserAttrTable);
                    case "system_attrs"
                        columns=obj.dataframeAttributeColumns( ...
                            columns,"system_attrs",obj.SystemAttrTable);
                    case "intermediate_values"
                        columns=obj.dataframeIntermediateColumns(columns);
                    case "state"
                        columns=obj.appendDataframeColumn(columns,attr,"", ...
                            obj.TrialStateData);
                end
            end

            top=reshape(string({columns.Top}),[],1);
            sub=reshape(string({columns.Sub}),[],1);
            metadata=obj.dataframeMetadata(top,sub,logical(multiIndex));
            result=table();
            storageNames=matlab.lang.makeUniqueStrings( ...
                cellstr(metadata.flat_columns));
            for index=1:numel(columns)
                result.(storageNames{index})=columns(index).Values;
            end
            result.Properties.UserData=metadata;
            descriptions=top;
            nested=strlength(sub)>0;
            descriptions(nested)=top(nested)+"."+sub(nested);
            result.Properties.VariableDescriptions=cellstr(descriptions);
        end

        function metadata=dataframeMetadata(~,top,sub,multiIndex)
            flat=top;
            nested=strlength(sub)>0;
            flat(nested)=top(nested)+"_"+sub(nested);
            metadata=struct( ...
                "schema","radia.optuna.trials-dataframe.v1", ...
                "multi_index",logical(multiIndex), ...
                "column_levels",[top,sub], ...
                "flat_columns",flat);
        end

        function columns=appendDataframeColumn(~,columns,top,sub,values)
            entry=struct();
            entry.Top=string(top);
            entry.Sub=string(sub);
            entry.Values=reshape(values,[],1);
            columns(end+1)=entry;
        end

        function columns=dataframeObjectiveColumns(obj,columns,top)
            objectiveCount=numel(obj.Directions);
            if ~isempty(obj.MetricNames)
                labels=reshape(obj.MetricNames,[],1);
                indices=(1:objectiveCount)';
                [labels,order]=sort(labels);
                indices=indices(order);
            else
                labels=string((0:objectiveCount-1)');
                indices=(1:objectiveCount)';
            end
            for item=1:numel(indices)
                values=NaN(numel(obj.TrialNumberData),1);
                sources=find(obj.ObjectiveIndexData==indices(item));
                [present,targets]=obj.dataframeTrialRows( ...
                    obj.ObjectiveTrialNumberData(sources));
                values(targets(present))= ...
                    obj.ObjectiveValueData(sources(present));
                columns=obj.appendDataframeColumn( ...
                    columns,top,labels(item),values);
            end
        end

        function columns=dataframeParameterColumns( ...
                obj,columns,top,decodeDistribution)
            names=sort(unique(obj.ParamNameData));
            trialCount=numel(obj.TrialNumberData);
            for name=reshape(names,1,[])
                raw=cell(trialCount,1);
                present=false(trialCount,1);
                rows=find(obj.ParamNameData==name);
                [known,targets]=obj.dataframeTrialRows( ...
                    obj.ParamTrialNumberData(rows));
                rows=rows(known);
                targets=targets(known);
                if decodeDistribution
                    for item=1:numel(rows)
                        source=rows(item);
                        target=targets(item);
                        raw{target}=radia.optuna.internal.DistributionCodec. ...
                            decode(obj.ParamKindData(source), ...
                            obj.ParamDistributionData(source));
                    end
                else
                    numeric=isfinite(obj.ParamValueNumericData(rows));
                    raw(targets(numeric))=num2cell( ...
                        obj.ParamValueNumericData(rows(numeric)));
                    textRows=~numeric;
                    raw(targets(textRows))=obj.decodeDataframeJsonArray( ...
                        obj.ParamValueTextData(rows(textRows)));
                end
                present(targets)=true;
                values=obj.coerceDataframeColumn(raw,present);
                columns=obj.appendDataframeColumn( ...
                    columns,top,name,values);
            end
        end

        function columns=dataframeAttributeColumns( ...
                obj,columns,top,attributes)
            names=sort(unique(attributes.Name));
            trialCount=numel(obj.TrialNumberData);
            for name=reshape(names,1,[])
                raw=cell(trialCount,1);
                present=false(trialCount,1);
                rows=find(attributes.Name==name);
                [known,targets]=obj.dataframeTrialRows( ...
                    attributes.TrialNumber(rows));
                rows=rows(known);
                targets=targets(known);
                raw(targets)=obj.decodeDataframeJsonArray( ...
                    attributes.ValueJSON(rows));
                present(targets)=true;
                values=obj.coerceDataframeColumn(raw,present);
                columns=obj.appendDataframeColumn( ...
                    columns,top,name,values);
            end
        end

        function columns=dataframeIntermediateColumns(obj,columns)
            steps=sort(unique(obj.IntermediateTable.Step));
            trialCount=numel(obj.TrialNumberData);
            for step=reshape(steps,1,[])
                values=NaN(trialCount,1);
                rows=find(obj.IntermediateTable.Step==step);
                [present,targets]=obj.dataframeTrialRows( ...
                    obj.IntermediateTable.TrialNumber(rows));
                values(targets(present))= ...
                    obj.IntermediateTable.Value(rows(present));
                columns=obj.appendDataframeColumn( ...
                    columns,"intermediate_values",string(step),values);
            end
        end

        function values=coerceDataframeColumn(~,raw,present)
            selected=raw(present);
            numeric=all(cellfun('isclass',selected,'double'));
            text=all(cellfun('isclass',selected,'char')) || ...
                all(cellfun('isclass',selected,'string'));
            if numeric
                values=NaN(numel(raw),1);
                values(present)=reshape(cell2mat(selected),[],1);
            elseif text
                values=strings(numel(raw),1);
                values(:)=missing;
                values(present)=string(selected);
            else
                values=raw;
                values(~present)={[]};
            end
        end

        function [present,rows]=dataframeTrialRows(obj,trialNumbers)
            trialNumbers=reshape(double(trialNumbers),[],1);
            rows=trialNumbers+1;
            present=isfinite(rows) & rows==floor(rows) & ...
                rows>=1 & rows<=numel(obj.TrialNumberData);
            indices=find(present);
            present(indices)=obj.TrialNumberData(rows(indices))== ...
                trialNumbers(indices);
        end

        function values=decodeDataframeJsonArray(~,encoded)
            encoded=reshape(string(encoded),[],1);
            values=cell(numel(encoded),1);
            if isempty(encoded)
                return
            end
            try
                decoded=jsondecode(char("["+strjoin(encoded,",")+"]"));
                if iscell(decoded)
                    values=reshape(decoded,[],1);
                elseif isstruct(decoded)
                    values=reshape(num2cell(decoded),[],1);
                elseif numel(decoded)==numel(encoded)
                    values=num2cell(reshape(decoded,[],1));
                else
                    error("radia:optuna:TrialsDataframeJSON", ...
                        "Decoded JSON array has an unexpected shape.");
                end
                if numel(values)~=numel(encoded)
                    error("radia:optuna:TrialsDataframeJSON", ...
                        "Decoded JSON array has an unexpected length.");
                end
            catch
                for index=1:numel(encoded)
                    values{index}=jsondecode(char(encoded(index)));
                end
            end
        end

        function data=trialData(obj)
            data=struct('TrialNumber',obj.TrialNumberData, ...
                'State',obj.TrialStateData,'Value',obj.TrialValueData, ...
                'StartTime',obj.TrialStartTimeData, ...
                'EndTime',obj.TrialEndTimeData, ...
                'Duration_s',obj.TrialDurationData, ...
                'Params',{obj.TrialParamsData}, ...
                'IntermediateValues',{obj.TrialIntermediateData}, ...
                'ErrorMessage',obj.TrialErrorData);
        end

        function data=parameterData(obj)
            data=struct('TrialNumber',obj.ParamTrialNumberData, ...
                'Name',obj.ParamNameData,'Kind',obj.ParamKindData, ...
                'ValueNumeric',obj.ParamValueNumericData, ...
                'ValueText',obj.ParamValueTextData, ...
                'Distribution',obj.ParamDistributionData);
        end

        function data=objectiveData(obj)
            data=struct('TrialNumber',obj.ObjectiveTrialNumberData, ...
                'ObjectiveIndex',obj.ObjectiveIndexData, ...
                'Value',obj.ObjectiveValueData);
        end

        function appendTrialRecord(obj,number,state,value,startTime, ...
                endTime,duration,params,intermediate,errorMessage)
            row=numel(obj.TrialNumberData)+1;
            obj.TrialNumberData(row,1)=double(number);
            obj.TrialStateData(row,1)=string(state);
            obj.TrialValueData(row,1)=double(value);
            obj.TrialStartTimeData(row,1)=double(startTime);
            obj.TrialEndTimeData(row,1)=double(endTime);
            obj.TrialDurationData(row,1)=double(duration);
            obj.TrialParamsData{row,1}=params;
            obj.TrialIntermediateData{row,1}=intermediate;
            obj.TrialErrorData(row,1)=string(errorMessage);
            obj.TrialTableDirty=true;
        end

        function row=trialRow(obj,number)
            candidate=double(number)+1;
            if candidate>=1 && candidate<=numel(obj.TrialNumberData) && ...
                    obj.TrialNumberData(candidate)==number
                row=candidate;
            else
                row=find(obj.TrialNumberData==number);
            end
        end

        function rows=trialRowsFor(obj,numbers)
            %TRIALROWSFOR Vectorized trialRow; zero marks an unknown trial.
            numbers=reshape(double(numbers),[],1);
            rows=zeros(size(numbers));
            candidate=numbers+1;
            direct=candidate>=1 & candidate<=numel(obj.TrialNumberData);
            direct(direct)= ...
                obj.TrialNumberData(candidate(direct))==numbers(direct);
            rows(direct)=candidate(direct);
            for index=reshape(find(~direct),1,[])
                match=find(obj.TrialNumberData==numbers(index),1);
                if ~isempty(match)
                    rows(index)=match;
                end
            end
        end

        function count=nonRunningTrialCount(obj)
            %NONRUNNINGTRIALCOUNT COMPLETE and PRUNED TPE history trials.
            count=sum(obj.TrialStateData=="COMPLETE" | ...
                obj.TrialStateData=="PRUNED");
        end

        function [steps,values]=lastIntermediateValues(obj,trialNumbers)
            %LASTINTERMEDIATEVALUES Deepest reported step and its value.
            trialNumbers=reshape(double(trialNumbers),[],1);
            steps=NaN(size(trialNumbers));
            values=NaN(size(trialNumbers));
            if height(obj.IntermediateTable)==0
                return
            end
            allNumbers=obj.IntermediateTable.TrialNumber;
            allSteps=obj.IntermediateTable.Step;
            allValues=obj.IntermediateTable.Value;
            for index=1:numel(trialNumbers)
                selected=find(allNumbers==trialNumbers(index));
                if isempty(selected)
                    continue
                end
                [steps(index),position]=max(allSteps(selected));
                values(index)=allValues(selected(position));
            end
        end

        function value=serialDatetimes(~,serial)
            value=datetime(serial,'ConvertFrom','datenum','TimeZone','local');
        end

        function name = generateStudyName(~)
            % Match Optuna's anonymous-study naming shape.
            try
                identifier = string(java.util.UUID.randomUUID());
            catch
                temporary = string(tempname("C:\temp"));
                identifier = extractAfter(temporary, "C:\temp\");
            end
            name = "no-name-" + identifier;
        end

        function suggestFixedDistributions(~, trial, distributions)
            names = string(fieldnames(distributions));
            for index = 1:numel(names)
                name = names(index);
                distribution = distributions.(name);
                if ~radia.optuna.internal.DistributionCodec.isSpec(distribution)
                    error("radia:optuna:FixedDistributions", ...
                        "Distribution '%s' is not a Radia Optuna distribution.", ...
                        name);
                end
                distribution=radia.optuna.internal.DistributionCodec.normalize( ...
                    distribution);
                switch string(distribution.kind)
                    case "float"
                        trial.suggest_float(name, distribution.low, ...
                            distribution.high, Log=distribution.log, ...
                            Step=distribution.step);
                    case "integer"
                        trial.suggest_int(name, distribution.low, ...
                            distribution.high, Step=distribution.step, ...
                            Log=distribution.log);
                    case "categorical"
                        trial.suggest_categorical(name, ...
                            distribution.choices);
                    otherwise
                        error("radia:optuna:FixedDistributions", ...
                            "Unsupported distribution kind '%s'.", ...
                            string(distribution.kind));
                end
            end
        end

        function stopWhenOptimizing(obj)
            if obj.InOptimize
                obj.StopRequested=true;
            end
        end

        function initializeTables(obj)
            persistent templates
            if ~isempty(templates)
                obj.TrialTable=templates.TrialTable;
                obj.ParamTable=templates.ParamTable;
                obj.IntermediateTable=templates.IntermediateTable;
                obj.UserAttrTable=templates.UserAttrTable;
                obj.SystemAttrTable=templates.SystemAttrTable;
                obj.ConstraintTable=templates.ConstraintTable;
                obj.ConstraintCountTable=templates.ConstraintCountTable;
                obj.ObjectiveTable=templates.ObjectiveTable;
                obj.SamplerStateTable=templates.SamplerStateTable;
                obj.QueueParamTable=templates.QueueParamTable;
                return
            end
            trialTable = table('Size', [0, 9], ...
                'VariableTypes', {'double','string','double','datetime','datetime', ...
                'double','cell','cell','string'}, ...
                'VariableNames', {'TrialNumber','State','Value','StartTime','EndTime', ...
                'Duration_s','Params','IntermediateValues','ErrorMessage'});
            trialTable.StartTime.TimeZone = "local";
            trialTable.EndTime.TimeZone = "local";
            obj.TrialTable = trialTable;
            obj.ParamTable = table('Size', [0, 6], ...
                'VariableTypes', {'double','string','string','double','string','string'}, ...
                'VariableNames', {'TrialNumber','Name','Kind','ValueNumeric','ValueText','Distribution'});
            obj.IntermediateTable = table('Size', [0, 4], ...
                'VariableTypes', {'double','double','double','datetime'}, ...
                'VariableNames', {'TrialNumber','Step','Value','Timestamp'});
            obj.UserAttrTable = table('Size', [0, 3], ...
                'VariableTypes', {'double','string','string'}, ...
                'VariableNames', {'TrialNumber','Name','ValueJSON'});
            obj.SystemAttrTable = table('Size', [0, 3], ...
                'VariableTypes', {'double','string','string'}, ...
                'VariableNames', {'TrialNumber','Name','ValueJSON'});
            obj.ConstraintTable = table('Size', [0, 3], ...
                'VariableTypes', {'double','double','double'}, ...
                'VariableNames', {'TrialNumber','ConstraintIndex','Value'});
            obj.ConstraintCountTable = table('Size', [0, 2], ...
                'VariableTypes', {'double','double'}, ...
                'VariableNames', {'TrialNumber','Count'});
            obj.ObjectiveTable = table('Size', [0, 3], ...
                'VariableTypes', {'double','double','double'}, ...
                'VariableNames', {'TrialNumber','ObjectiveIndex','Value'});
            obj.SamplerStateTable = table('Size', [0, 7], ...
                'VariableTypes', ...
                {'string','string','double','double','double','cell','datetime'}, ...
                'VariableNames', {'Sampler','Schema','Revision','TrialNumber', ...
                'Generation','State','Timestamp'});
            obj.QueueParamTable=table('Size',[0,3], ...
                'VariableTypes',{'double','string','cell'}, ...
                'VariableNames',{'TrialNumber','Name','Value'});
            obj.IntermediateTable.Timestamp.TimeZone = "local";
            obj.SamplerStateTable.Timestamp.TimeZone = "local";
            templates=struct('TrialTable',obj.TrialTable, ...
                'ParamTable',obj.ParamTable, ...
                'IntermediateTable',obj.IntermediateTable, ...
                'UserAttrTable',obj.UserAttrTable, ...
                'SystemAttrTable',obj.SystemAttrTable, ...
                'ConstraintTable',obj.ConstraintTable, ...
                'ConstraintCountTable',obj.ConstraintCountTable, ...
                'ObjectiveTable',obj.ObjectiveTable, ...
                'SamplerStateTable',obj.SamplerStateTable, ...
                'QueueParamTable',obj.QueueParamTable);
        end

        function loadState(obj)
            [primaryPath, ~] = obj.resolvedStoragePath();
            backupPath = primaryPath + ".bak";
            primaryFailure = "Primary storage file does not exist.";
            backupFailure = "Backup storage file does not exist.";
            loadedFromBackup = false;
            loadedSuccessfully = false;
            if isfile(primaryPath)
                try
                    data = obj.readStudyData(primaryPath);
                    loadedSuccessfully = true;
                catch exception
                    primaryFailure = string(exception.message);
                end
            end
            if ~loadedSuccessfully && isfile(backupPath)
                try
                    data = obj.readStudyData(backupPath);
                    loadedSuccessfully = true;
                    loadedFromBackup = true;
                catch exception
                    backupFailure = string(exception.message);
                end
            end
            if ~loadedSuccessfully
                error("radia:optuna:StorageRecovery", ...
                    "Unable to load study storage '%s' or its backup. " + ...
                    "Primary: %s Backup: %s", primaryPath, ...
                    primaryFailure, backupFailure);
            end
            if loadedFromBackup
                warning("radia:optuna:RecoveredStorage", ...
                    "Loaded study from backup '%s' after primary load failed: %s", ...
                    backupPath, primaryFailure);
            end
            obj.Name = string(data.Name);
            obj.Directions = string(data.Directions);
            obj.NextTrialNumber = data.NextTrialNumber;
            obj.TrialTable = data.TrialTable;
            obj.ParamTable = data.ParamTable;
            obj.IntermediateTable = data.IntermediateTable;
            obj.UserAttrTable = data.UserAttrTable;
            if isfield(data,"SystemAttrTable")
                obj.SystemAttrTable=data.SystemAttrTable;
            end
            if isfield(data, "ConstraintTable")
                obj.ConstraintTable = data.ConstraintTable;
            end
            if isfield(data, "ConstraintCountTable")
                obj.ConstraintCountTable = data.ConstraintCountTable;
            elseif ~isempty(obj.ConstraintTable)
                numbers = unique(obj.ConstraintTable.TrialNumber,"stable");
                counts = zeros(numel(numbers),1);
                for index = 1:numel(numbers)
                    counts(index) = sum( ...
                        obj.ConstraintTable.TrialNumber == numbers(index));
                end
                obj.ConstraintCountTable = table(numbers,counts, ...
                    'VariableNames',{'TrialNumber','Count'});
            end
            if isfield(data, "UserAttrs")
                obj.UserAttrs = data.UserAttrs;
            end
            if isfield(data,"SystemAttrs")
                obj.SystemAttrs=data.SystemAttrs;
            end
            if isfield(data,"MetricNames")
                obj.MetricNames=string(data.MetricNames);
            end
            if isfield(data, "ObjectiveTable")
                obj.ObjectiveTable = data.ObjectiveTable;
            else
                complete = obj.TrialTable.State == "COMPLETE" & ~isnan(obj.TrialTable.Value);
                obj.ObjectiveTable = table(obj.TrialTable.TrialNumber(complete), ...
                    ones(sum(complete),1), obj.TrialTable.Value(complete), ...
                    'VariableNames',{'TrialNumber','ObjectiveIndex','Value'});
            end
            if isfield(data, "SamplerStateTable")
                obj.SamplerStateTable = data.SamplerStateTable;
            end
            if isfield(data,"QueueParamTable")
                obj.QueueParamTable=data.QueueParamTable;
            end
        end

        function frozen = freezeTrial(obj,trialNumber)
            row=obj.trialRow(trialNumber);
            if numel(row)~=1
                error("radia:optuna:UnknownTrial", ...
                    "Trial %d does not identify exactly one row.",trialNumber);
            end
            values=NaN;
            objectiveRows=obj.ObjectiveTrialNumberData==trialNumber;
            if any(objectiveRows)
                [~,order]=sort(obj.ObjectiveIndexData(objectiveRows));
                objectiveValues=obj.ObjectiveValueData(objectiveRows);
                values=reshape(objectiveValues(order),1,[]);
            elseif ~isnan(obj.TrialValueData(row))
                values=obj.TrialValueData(row);
            end
            distributions=struct();
            parameterRows=find(obj.ParamTrialNumberData==trialNumber)';
            parameterKeys=radia.optuna.Trial.claimKeys( ...
                obj.ParamNameData(parameterRows));
            for index=1:numel(parameterRows)
                parameterRow=parameterRows(index);
                distributions.(char(parameterKeys(index)))= ...
                    radia.optuna.internal.DistributionCodec.decode( ...
                    obj.ParamKindData(parameterRow), ...
                    obj.ParamDistributionData(parameterRow));
            end
            userAttrs=obj.attributesForTrial(obj.UserAttrTable,trialNumber);
            systemAttrs=obj.attributesForTrial(obj.SystemAttrTable,trialNumber);
            intermediateRows=obj.IntermediateTable.TrialNumber==trialNumber;
            intermediate=obj.IntermediateTable(intermediateRows, ...
                ["Step","Value","Timestamp"]);
            [constraintPresent,constraints]=obj.constraintRecord(trialNumber);
            frozen=radia.optuna.FrozenTrial(Number=trialNumber, ...
                State=obj.TrialStateData(row),Values=values, ...
                Params=obj.TrialParamsData{row}, ...
                Distributions=distributions, ...
                IntermediateValues=intermediate,UserAttrs=userAttrs, ...
                SystemAttrs=systemAttrs,Constraints=constraints, ...
                ConstraintPresent=constraintPresent, ...
                DatetimeStart=obj.serialDatetimes(obj.TrialStartTimeData(row)), ...
                DatetimeComplete=obj.serialDatetimes(obj.TrialEndTimeData(row)), ...
                ErrorMessage=obj.TrialErrorData(row));
        end

        function frozen = freezeTrials(obj,trialNumbers)
            trialNumbers=reshape(double(trialNumbers),[],1);
            if isempty(trialNumbers)
                frozen=radia.optuna.FrozenTrial.empty(0,1);
                return
            end
            frozen(numel(trialNumbers),1)=obj.freezeTrial(trialNumbers(end));
            for index=1:numel(trialNumbers)
                frozen(index,1)=obj.freezeTrial(trialNumbers(index));
            end
        end

        function trial=resolveTrial(obj,value)
            if isa(value,"radia.optuna.Trial")
                if ~isscalar(value)
                    error("radia:optuna:Trial", ...
                        "tell requires one Trial or one trial number.");
                end
                trial=value;
                return
            end
            if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
                    value==floor(value) && value>=0)
                error("radia:optuna:Trial", ...
                    "tell requires one Trial or one nonnegative trial number.");
            end
            number=double(value);
            rows=obj.trialRow(number);
            if numel(rows)~=1
                error("radia:optuna:UnknownTrial", ...
                    "Trial %d does not identify exactly one row.",number);
            end
            trial=radia.optuna.Trial(obj,number);
            trial.restoreSnapshot(obj.freezeTrial(number));
        end

        function addTrial(obj,frozen)
            arguments
                obj
                frozen (1,1) radia.optuna.FrozenTrial
            end
            if frozen.State=="COMPLETE" && ...
                    (numel(frozen.Values)~=numel(obj.Directions) || ...
                    any(isnan(frozen.Values)))
                error("radia:optuna:AddTrialValue", ...
                    "A COMPLETE trial must have one non-NaN value per direction.");
            end
            names=string(fieldnames(frozen.Params));
            distributionNames=string(fieldnames(frozen.Distributions));
            if ~isempty(names) && ~all(ismember(names,distributionNames))
                error("radia:optuna:AddTrialDistribution", ...
                    "Every imported parameter requires a distribution.");
            end
            number=obj.NextTrialNumber;
            obj.NextTrialNumber=number+1;
            startTime=frozen.DatetimeStart;
            endTime=frozen.DatetimeComplete;
            if frozen.State=="WAITING"
                startTime=datetime(NaT,"TimeZone","local");
                endTime=datetime(NaT,"TimeZone","local");
                duration=NaN;
            elseif frozen.State=="RUNNING"
                if isnat(startTime)
                    startTime=datetime("now","TimeZone","local");
                elseif strlength(string(startTime.TimeZone))==0
                    startTime.TimeZone="local";
                end
                endTime=datetime(NaT,"TimeZone","local");
                duration=NaN;
            else
                if isnat(startTime)
                    startTime=datetime("now","TimeZone","local");
                elseif strlength(string(startTime.TimeZone))==0
                    startTime.TimeZone="local";
                end
                if isnat(endTime), endTime=startTime; end
                if strlength(string(endTime.TimeZone))==0
                    endTime.TimeZone="local";
                end
                duration=seconds(endTime-startTime);
            end
            scalarValue=NaN;
            if frozen.State=="COMPLETE", scalarValue=frozen.Values(1); end
            startSerial=datenum(startTime); %#ok<DATNM>
            endSerial=datenum(endTime); %#ok<DATNM>
            obj.appendTrialRecord(number,frozen.State,scalarValue, ...
                startSerial,endSerial,duration,frozen.Params, ...
                frozen.IntermediateValues,frozen.ErrorMessage);
            if frozen.State=="COMPLETE"
                count=numel(frozen.Values);
                obj.ObjectiveTrialNumberData(end+(1:count),1)=number;
                obj.ObjectiveIndexData(end+(1:count),1)=(1:count)';
                obj.ObjectiveValueData(end+(1:count),1)= ...
                    reshape(frozen.Values,[],1);
                obj.ObjectiveTableDirty=true;
            end
            for name=reshape(names,1,[])
                distribution=frozen.Distributions.(name);
                obj.appendImportedParameter(number,name,frozen.Params.(name), ...
                    distribution);
            end
            obj.UserAttrTable=obj.appendImportedAttributes( ...
                obj.UserAttrTable,number,frozen.UserAttrs);
            obj.SystemAttrTable=obj.appendImportedAttributes( ...
                obj.SystemAttrTable,number,frozen.SystemAttrs);
            if frozen.ConstraintPresent
                obj.ConstraintCountTable(end+1,:)={number,numel(frozen.Constraints)};
                if ~isempty(frozen.Constraints)
                    obj.ConstraintTable=[obj.ConstraintTable;table( ...
                        repmat(number,numel(frozen.Constraints),1), ...
                        (1:numel(frozen.Constraints))', ...
                        reshape(frozen.Constraints,[],1), ...
                        'VariableNames',obj.ConstraintTable.Properties.VariableNames)];
                end
            end
            if ~isempty(frozen.IntermediateValues)
                count=height(frozen.IntermediateValues);
                timestamps=frozen.IntermediateValues.Timestamp;
                obj.IntermediateTable=[obj.IntermediateTable;table( ...
                    repmat(number,count,1),frozen.IntermediateValues.Step, ...
                    frozen.IntermediateValues.Value,timestamps, ...
                    'VariableNames',obj.IntermediateTable.Properties.VariableNames)];
            end
            obj.persist();
        end

        function add_trial(obj,frozen)
            obj.addTrial(frozen);
        end

        function addTrials(obj,trials)
            for index=1:numel(trials)
                obj.addTrial(trials(index));
            end
        end

        function add_trials(obj,trials)
            obj.addTrials(trials);
        end

        function recordSamplerState(obj, sampler, schema, trialNumber, ...
                generation, state)
            if ~isstruct(state) || ~isscalar(state)
                error("radia:optuna:SamplerState", ...
                    "Sampler state must be one scalar struct.");
            end
            sampler = string(sampler);
            schema = string(schema);
            rows = obj.SamplerStateTable.Sampler == sampler & ...
                obj.SamplerStateTable.Schema == schema;
            revision = 1;
            if any(rows)
                revision = max(obj.SamplerStateTable.Revision(rows)) + 1;
                obj.SamplerStateTable(rows,:) = [];
            end
            obj.SamplerStateTable(end+1,:) = {sampler, schema, revision, ...
                double(trialNumber), double(generation), {state}, ...
                datetime("now", "TimeZone", "local")};
        end

        function state = samplerState(obj, sampler, schema)
            sampler = string(sampler);
            schema = string(schema);
            rows = obj.SamplerStateTable.Sampler == sampler & ...
                obj.SamplerStateTable.Schema == schema;
            if ~any(rows)
                state = [];
                return
            end
            selected = obj.SamplerStateTable(rows,:);
            [~, index] = max(selected.Revision);
            state = selected.State{index};
        end

        function recordParameter(obj, trial, name, kind, value, distribution)
            existing=find(obj.ParamTrialNumberData==trial.Number & ...
                obj.ParamNameData==name);
            numeric = NaN;
            text = "";
            if isnumeric(value) && isscalar(value)
                numeric = double(value);
            else
                text = string(jsonencode(value));
            end
            if isempty(existing)
                row=numel(obj.ParamTrialNumberData)+1;
            else
                row=existing(1);
                if ~isscalar(existing)
                    stale=existing(2:end);
                    obj.ParamTrialNumberData(stale)=[];
                    obj.ParamNameData(stale)=[];
                    obj.ParamKindData(stale)=[];
                    obj.ParamValueNumericData(stale)=[];
                    obj.ParamValueTextData(stale)=[];
                    obj.ParamDistributionData(stale)=[];
                end
            end
            obj.ParamTrialNumberData(row,1)=trial.Number;
            obj.ParamNameData(row,1)=string(name);
            obj.ParamKindData(row,1)=string(kind);
            obj.ParamValueNumericData(row,1)=numeric;
            obj.ParamValueTextData(row,1)=text;
            obj.ParamDistributionData(row,1)=string(distribution);
            obj.ParamTableDirty=true;
            trialRow=obj.trialRow(trial.Number);
            obj.TrialParamsData{trialRow}=trial.Params;
            obj.TrialIntermediateData{trialRow}=trial.IntermediateValues;
            obj.TrialTableDirty=true;
            if obj.AutoSave && strlength(obj.StoragePath)>0
                obj.save();
            end
        end

        function recordIntermediate(obj, trial, value, step)
            rows = obj.IntermediateTable.TrialNumber == trial.Number & ...
                obj.IntermediateTable.Step == step;
            if any(rows)
                obj.IntermediateTable(rows,:) = [];
            end
            obj.IntermediateTable(end+1,:) = {trial.Number, step, value, ...
                datetime("now", "TimeZone", "local")};
            obj.updateTrialSnapshot(trial);
            obj.persist();
        end

        function recordUserAttribute(obj, trial, name, value)
            rows = obj.UserAttrTable.TrialNumber == trial.Number & ...
                obj.UserAttrTable.Name == name;
            if any(rows)
                obj.UserAttrTable(rows,:) = [];
            end
            obj.UserAttrTable(end+1,:) = {trial.Number, name, string(jsonencode(value))};
            obj.persist();
        end

        function recordSystemAttribute(obj,trial,name,value)
            rows=obj.SystemAttrTable.TrialNumber==trial.Number & ...
                obj.SystemAttrTable.Name==name;
            if any(rows)
                obj.SystemAttrTable(rows,:)=[];
            end
            obj.SystemAttrTable(end+1,:)={trial.Number,name, ...
                string(jsonencode(value))};
            obj.persist();
        end

        function recordConstraints(obj, trial, values)
            values = reshape(double(values), [], 1);
            if any(isnan(values))
                error("radia:optuna:Constraints", ...
                    "Constraint values must not contain NaN.");
            end
            if ~isempty(values) && ( ...
                    isa(obj.Sampler, "radia.optuna.CmaEsSampler"))
                error("radia:optuna:SamplerConstraints", ...
                    "The configured sampler does not implement " + ...
                    "constraint-aware ranking. Use TPESampler, " + ...
                    "MOTPESampler, or NSGAIISampler.");
            end
            obj.ConstraintTable( ...
                obj.ConstraintTable.TrialNumber == trial.Number, :) = [];
            obj.ConstraintCountTable( ...
                obj.ConstraintCountTable.TrialNumber == trial.Number, :) = [];
            obj.ConstraintCountTable(end+1,:) = ...
                {trial.Number,double(numel(values))};
            if ~isempty(values)
                obj.ConstraintTable = [obj.ConstraintTable; table( ...
                    repmat(trial.Number, numel(values), 1), ...
                    (1:numel(values))', values, ...
                    'VariableNames', obj.ConstraintTable.Properties.VariableNames)];
            end
            trial.setConstraints(values);
            obj.persist();
        end

        function values = constraintsForTrial(obj, trialNumber)
            [present,values] = obj.constraintRecord(trialNumber);
            if ~present
                values = zeros(1,0);
            end
        end

        function [present, values] = constraintRecord(obj, trialNumber)
            countRows = obj.ConstraintCountTable.TrialNumber == trialNumber;
            if sum(countRows) > 1
                error("radia:optuna:ConstraintShape", ...
                    "Trial %d has duplicate constraint-count records.",trialNumber);
            end
            present = any(countRows);
            if ~present
                values = zeros(1,0);
                return
            end
            count = obj.ConstraintCountTable.Count(countRows);
            if ~isscalar(count) || ~isfinite(count) || count < 0 || ...
                    count ~= floor(count)
                error("radia:optuna:ConstraintShape", ...
                    "Trial %d has an invalid constraint count.",trialNumber);
            end
            rows = obj.ConstraintTable.TrialNumber == trialNumber;
            selected = sortrows(obj.ConstraintTable(rows,:),"ConstraintIndex");
            expectedIndices = (1:count)';
            if height(selected) ~= count || ...
                    ~isequal(selected.ConstraintIndex,expectedIndices)
                error("radia:optuna:ConstraintShape", ...
                    "Trial %d has an incomplete constraint vector.",trialNumber);
            end
            values = reshape(selected.Value,1,[]);
            if any(isnan(values))
                error("radia:optuna:Constraints", ...
                    "Trial %d has NaN constraint values.",trialNumber);
            end
        end

        function result = hasConstraintRecords(obj)
            result = ~isempty(obj.ConstraintCountTable);
        end

        function finishTrial(obj, trial, state, value, message)
            obj.ensureTrialOwnership(trial);
            rows = obj.trialRow(trial.Number);
            if numel(rows) ~= 1
                error("radia:optuna:UnknownTrial", ...
                    "Trial %d does not identify exactly one row in this study.", ...
                    trial.Number);
            end
            storedState = obj.TrialStateData(rows);
            if trial.State ~= "RUNNING" || storedState ~= "RUNNING"
                if trial.State ~= "RUNNING"
                    finishedState = trial.State;
                else
                    finishedState = storedState;
                end
                error("radia:optuna:TrialState", "Trial %d is already %s.", ...
                    trial.Number, finishedState);
            end
            endTime = now; %#ok<TNOW1> hot-path serial timestamp
            trial.markFinished(state, value, endTime, message);
            elapsed = (endTime-trial.startTimeSerial())*86400;
            objectiveRows=obj.ObjectiveTrialNumberData==trial.Number;
            if any(objectiveRows)
                obj.ObjectiveTrialNumberData(objectiveRows)=[];
                obj.ObjectiveIndexData(objectiveRows)=[];
                obj.ObjectiveValueData(objectiveRows)=[];
            end
            if state == "COMPLETE"
                values = reshape(double(value),[],1);
                count=numel(values);
                newRows=numel(obj.ObjectiveTrialNumberData)+(1:count);
                obj.ObjectiveTrialNumberData(newRows,1)=trial.Number;
                obj.ObjectiveIndexData(newRows,1)=(1:count)';
                obj.ObjectiveValueData(newRows,1)=values;
            end
            obj.ObjectiveTableDirty=true;
            obj.TrialStateData(rows)=state;
            obj.TrialValueData(rows)=value(1);
            obj.TrialEndTimeData(rows)=endTime;
            obj.TrialDurationData(rows)=elapsed;
            obj.TrialErrorData(rows)=message;
            obj.TrialParamsData{rows}=trial.Params;
            obj.TrialIntermediateData{rows}=trial.IntermediateValues;
            obj.TrialTableDirty=true;
            if obj.SamplerHasAfterTrial
                try
                    obj.Sampler.afterTrial(obj, trial);
                catch exception
                    % Optuna commits the requested trial state in a finally
                    % block even when sampler post-processing raises.
                    obj.persist();
                    rethrow(exception)
                end
            end
            if obj.AutoSave && strlength(obj.StoragePath)>0
                obj.save();
            end
            if ~isempty(obj.ProgressFcn)
                obj.reportProgress("TRIAL_FINISHED", trial);
            end
        end

        function updateTrialSnapshot(obj, trial)
            rows = obj.trialRow(trial.Number);
            obj.TrialParamsData{rows}=trial.Params;
            obj.TrialIntermediateData{rows}=trial.IntermediateValues;
            obj.TrialTableDirty=true;
        end

        function persist(obj)
            if obj.AutoSave && strlength(obj.StoragePath) > 0
                obj.save();
            end
        end

        function reportProgress(obj, eventName, trial)
            if isempty(obj.ProgressFcn)
                return;
            end
            snapshot = struct( ...
                "event", string(eventName), ...
                "trial_number", trial.Number, ...
                "trial_state", trial.State, ...
                "trial_value", trial.Value, ...
                "trial_table", obj.TrialTable, ...
                "param_table", obj.ParamTable, ...
                "pareto_front", obj.paretoFront());
            if isscalar(obj.Directions)
                try
                    snapshot.best_trial = obj.bestTrial();
                catch exception
                    if string(exception.identifier) ~= ...
                            "radia:optuna:NoFeasibleTrial"
                        rethrow(exception)
                    end
                    snapshot.best_trial = obj.TrialTable([],:);
                end
            else
                snapshot.best_trial = obj.TrialTable([],:);
            end
            obj.ProgressFcn(snapshot);
        end

        function skip = skipFinishedTrial(obj, trial, skipIfFinished)
            obj.ensureTrialOwnership(trial);
            rows = obj.trialRow(trial.Number);
            if numel(rows) ~= 1
                error("radia:optuna:UnknownTrial", ...
                    "Trial %d does not identify exactly one row in this study.", ...
                    trial.Number);
            end
            storedState = obj.TrialStateData(rows);
            finished = trial.State ~= "RUNNING" || storedState ~= "RUNNING";
            skip = finished && skipIfFinished;
            if finished && ~skip
                if trial.State ~= "RUNNING"
                    finishedState = trial.State;
                else
                    finishedState = storedState;
                end
                error("radia:optuna:TrialState", ...
                    "Trial %d is already %s.", trial.Number, finishedState);
            end
        end

        function [feasible, constrained] = feasibleTrials(obj, trialNumbers)
            % Once any constraint is present, missing constraint rows are
            % unknown rather than implicitly feasible. This protects legacy
            % and partially written studies from selecting an unchecked run.
            trialNumbers = reshape(double(trialNumbers), [], 1);
            constrained = obj.hasConstraintRecords();
            feasible = true(size(trialNumbers));
            if ~constrained
                return
            end
            feasible(:) = false;
            for index = 1:numel(trialNumbers)
                [present,values] = obj.constraintRecord(trialNumbers(index));
                feasible(index) = present && all(values <= 0);
            end
        end

        function result = emptyParetoFront(~)
            result = table('Size',[0,3], ...
                'VariableTypes',{'double','cell','cell'}, ...
                'VariableNames',{'TrialNumber','Values','Params'});
        end

        function path = backupStoragePath(obj)
            path = obj.StoragePath + ".bak";
        end

        function [targetPath, folder] = resolvedStoragePath(obj)
            targetPath = obj.StoragePath;
            folder = string(fileparts(targetPath));
            if strlength(folder) == 0
                folder = string(pwd);
                targetPath = fullfile(folder, targetPath);
            end
        end

        function data = readStudyData(obj, path)
            loaded = builtin("load", path, "StudyData", "-mat");
            if ~isfield(loaded, "StudyData") || ...
                    ~isstruct(loaded.StudyData) || ~isscalar(loaded.StudyData)
                error("radia:optuna:StorageFormat", ...
                    "Storage '%s' does not contain one StudyData struct.", path);
            end
            data = loaded.StudyData;
            hasSchema = isfield(data, "Schema");
            hasVersion = isfield(data, "Version");
            if xor(hasSchema, hasVersion)
                error("radia:optuna:StorageVersion", ...
                    "Storage '%s' has an incomplete schema/version header.", ...
                    path);
            end
            if hasSchema
                schema = string(data.Schema);
                try
                    version = double(data.Version);
                catch
                    error("radia:optuna:StorageVersion", ...
                        "Storage '%s' has a nonnumeric version header.", path);
                end
                if ~isscalar(schema) || ~isscalar(version)
                    error("radia:optuna:StorageVersion", ...
                        "Storage '%s' has nonscalar schema/version headers.", ...
                        path);
                end
                if schema ~= obj.StorageSchema || ~isfinite(version) || ...
                        version ~= floor(version) || ...
                        ~ismember(version, 1:obj.StorageVersion)
                    error("radia:optuna:StorageVersion", ...
                        "Storage '%s' has unsupported schema/version '%s'/%g.", ...
                        path, schema, version);
                end
            end
            required = ["Name", "Directions", "NextTrialNumber", ...
                "TrialTable", "ParamTable", "IntermediateTable", ...
                "UserAttrTable"];
            for field = required
                if ~isfield(data, field)
                    error("radia:optuna:StorageFormat", ...
                        "Storage '%s' is missing StudyData.%s.", path, field);
                end
            end
        end

        function replaceStorageFile(obj, temporaryPath, targetPath, backupPath)
            if isfile(targetPath) && ~isfile(backupPath)
                try
                    obj.readStudyData(targetPath);
                    obj.copyValidStudyFile(targetPath, backupPath);
                catch
                    % Do not seed the backup from a damaged primary.
                end
            end
            obj.atomicMoveFile(temporaryPath, targetPath);
            % Refresh the backup only after the primary atomic move. The
            % backup's own temporary move leaves the previous copy intact
            % if this refresh is interrupted.
            obj.copyValidStudyFile(targetPath, backupPath);
        end

        function copyValidStudyFile(obj, sourcePath, targetPath)
            folder = string(fileparts(targetPath));
            temporaryPath = string(tempname(folder)) + ".mat";
            cleanup = onCleanup(@() obj.removeTemporaryFile(temporaryPath));
            [copied, message, messageId] = copyfile( ...
                sourcePath, temporaryPath, "f");
            if ~copied
                error("radia:optuna:StorageBackup", ...
                    "Unable to create storage backup (%s): %s", ...
                    messageId, message);
            end
            obj.readStudyData(temporaryPath);
            obj.atomicMoveFile(temporaryPath, targetPath);
            clear cleanup
        end

        function atomicMoveFile(~, sourcePath, targetPath)
            % Java NIO requests a same-filesystem atomic replace. movefile
            % remains a compatibility fallback for no-JVM deployments.
            try
                emptyStrings = javaArray('java.lang.String', 0);
                source = javaMethod('get', 'java.nio.file.Paths', ...
                    char(sourcePath), emptyStrings);
                target = javaMethod('get', 'java.nio.file.Paths', ...
                    char(targetPath), emptyStrings);
                options = javaArray('java.nio.file.CopyOption', 2);
                options(1) = javaMethod('valueOf', ...
                    'java.nio.file.StandardCopyOption', 'ATOMIC_MOVE');
                options(2) = javaMethod('valueOf', ...
                    'java.nio.file.StandardCopyOption', 'REPLACE_EXISTING');
                javaMethod('move', 'java.nio.file.Files', ...
                    source, target, options);
            catch javaFailure
                [moved, message, messageId] = movefile( ...
                    sourcePath, targetPath, "f");
                if ~moved
                    storageFailure = MException( ...
                        "radia:optuna:StorageReplace", ...
                        "Unable to replace study storage (%s): %s", ...
                        messageId, message);
                    storageFailure = addCause(storageFailure, javaFailure);
                    throw(storageFailure)
                end
            end
        end

        function removeTemporaryFile(~, path)
            if isfile(path)
                try
                    delete(path);
                catch
                    % Cleanup must not hide the original persistence error.
                end
            end
        end

        function ensureTrialOwnership(obj, trial)
            if ~isequal(trial.Study, obj)
                error("radia:optuna:StudyMismatch", ...
                    "Trial %d belongs to a different study.", trial.Number);
            end
        end

        function result=containsParameterSet(obj,params)
            result=false;
            target=orderfields(params);
            for index=1:numel(obj.TrialNumberData)
                number=obj.TrialNumberData(index);
                if obj.TrialStateData(index)=="WAITING"
                    rows=obj.QueueParamTable.TrialNumber==number;
                    candidate=struct();
                    for row=find(rows)'
                        candidate.(obj.QueueParamTable.Name(row))= ...
                            obj.QueueParamTable.Value{row};
                    end
                else
                    candidate=obj.TrialParamsData{index};
                end
                if isequaln(orderfields(candidate),target)
                    result=true;
                    return
                end
            end
        end

        function callbacks=normalizeCallbacks(~,value)
            if isempty(value)
                callbacks=cell(1,0);
            elseif isa(value,"function_handle")
                callbacks={value};
            elseif iscell(value) && all(cellfun( ...
                    @(item)isa(item,"function_handle"),value))
                callbacks=reshape(value,1,[]);
            else
                error("radia:optuna:Callbacks", ...
                    "Callbacks must be a function handle or cell array of handles.");
            end
        end

        function finishOptimize(obj)
            obj.InOptimize=false;
            obj.StopRequested=false;
        end

        function attrs=attributesForTrial(~,source,trialNumber)
            attrs=struct();
            rows=find(source.TrialNumber==trialNumber)';
            for row=rows
                attrs.(matlab.lang.makeValidName(source.Name(row)))= ...
                    jsondecode(source.ValueJSON(row));
            end
        end

        function appendImportedParameter(obj,number,name,value,distribution)
            if ~radia.optuna.internal.DistributionCodec.isSpec(distribution)
                error("radia:optuna:AddTrialDistribution", ...
                    "Imported distribution for '%s' is invalid.",name);
            end
            distribution=radia.optuna.internal.DistributionCodec.normalize( ...
                distribution);
            numeric=NaN;
            textValue="";
            if isnumeric(value) && isscalar(value)
                numeric=double(value);
            else
                textValue=string(jsonencode(value));
            end
            row=numel(obj.ParamTrialNumberData)+1;
            obj.ParamTrialNumberData(row,1)=number;
            obj.ParamNameData(row,1)=name;
            obj.ParamKindData(row,1)=string(distribution.kind);
            obj.ParamValueNumericData(row,1)=numeric;
            obj.ParamValueTextData(row,1)=textValue;
            obj.ParamDistributionData(row,1)= ...
                radia.optuna.internal.DistributionCodec.encode(distribution);
            obj.ParamTableDirty=true;
        end

        function target=appendImportedAttributes(~,target,number,attrs)
            names=string(fieldnames(attrs));
            for name=reshape(names,1,[])
                target(end+1,:)={number,name, ...
                    string(jsonencode(attrs.(name)))}; %#ok<AGROW>
            end
        end
    end
end
