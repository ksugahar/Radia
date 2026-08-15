classdef Study < handle
    %STUDY Table-backed Optuna-style experiment database.

    properties (SetAccess=private)
        Name (1,1) string
        Directions string
        StoragePath (1,1) string = ""
        Sampler
        Pruner
        AutoSave (1,1) logical = true
        TrialTable table
        ParamTable table
        IntermediateTable table
        UserAttrTable table
        SystemAttrTable table
        ConstraintTable table
        ConstraintCountTable table
        ObjectiveTable table
        SamplerStateTable table
        QueueParamTable table
        ProgressFcn = []
        UserAttrs struct = struct()
        SystemAttrs struct = struct()
        MetricNames string = strings(1,0)
    end

    properties (Access=private)
        NextTrialNumber (1,1) double = 0
        StopRequested (1,1) logical = false
        InOptimize (1,1) logical = false
    end

    properties (Constant, Access=private)
        StorageSchema = "radia.optuna.study"
        StorageVersion = 4
    end

    methods
        function obj = Study(options)
            arguments
                options.Name (1,1) string = "radia-study"
                options.StoragePath (1,1) string = ""
                options.Directions string = "minimize"
                options.Sampler = []
                options.Pruner = []
                options.AutoSave (1,1) logical = true
                options.ProgressFcn = []
            end
            obj.Name = options.Name;
            obj.Directions = options.Directions;
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
                        Seed=0, NStartupTrials=10);
                else
                    obj.Sampler = radia.optuna.MOTPESampler( ...
                        Seed=0, NStartupTrials=20);
                end
            end
            if isempty(obj.Pruner)
                obj.Pruner = radia.optuna.MedianPruner();
            end
            obj.initializeTables();
            if strlength(obj.StoragePath) > 0 && ...
                    (isfile(obj.StoragePath) || isfile(obj.backupStoragePath()))
                obj.loadState();
            end
        end

        function trial = ask(obj)
            waiting=find(obj.TrialTable.State=="WAITING",1);
            if isempty(waiting)
                trial = radia.optuna.Trial(obj, obj.NextTrialNumber);
                obj.NextTrialNumber = obj.NextTrialNumber + 1;
                obj.TrialTable(end+1,:) = {trial.Number, "RUNNING", NaN, ...
                    trial.StartTime, NaT, NaN, {trial.Params}, ...
                    {trial.IntermediateValues}, ""};
            else
                number=obj.TrialTable.TrialNumber(waiting);
                trial=radia.optuna.Trial(obj,number);
                rows=obj.QueueParamTable.TrialNumber==number;
                queuedParams=obj.QueueParamTable(rows,:);
                trial.setFixedParameters(queuedParams.Name,queuedParams.Value);
                userRows=obj.UserAttrTable.TrialNumber==number;
                for index=find(userRows)'
                    trial.setUserAttr(obj.UserAttrTable.Name(index), ...
                        jsondecode(obj.UserAttrTable.ValueJSON(index)));
                end
                obj.TrialTable.State(waiting)="RUNNING";
                obj.TrialTable.StartTime(waiting)=trial.StartTime;
                obj.TrialTable.EndTime(waiting)=NaT;
                obj.TrialTable.Duration_s(waiting)=NaN;
                obj.TrialTable.ErrorMessage(waiting)="";
            end
            if ismethod(obj.Sampler, "beforeTrial")
                obj.Sampler.beforeTrial(obj, trial);
            end
            obj.updateTrialSnapshot(trial);
            obj.persist();
            obj.reportProgress("TRIAL_STARTED", trial);
        end

        function tell(obj, trial, value, options)
            arguments
                obj
                trial (1,1) radia.optuna.Trial
                value double
                options.SkipIfFinished (1,1) logical = false
            end
            if obj.skipFinishedTrial(trial, options.SkipIfFinished)
                return
            end
            value = reshape(double(value), 1, []);
            if numel(value) ~= numel(obj.Directions) || any(~isfinite(value))
                error("radia:optuna:Value", ...
                    "Trial value must contain one finite value per direction.");
            end
            obj.finishTrial(trial, "COMPLETE", value, "");
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
            notATime=datetime(NaT,"TimeZone","local");
            obj.TrialTable(end+1,:)={number,"WAITING",NaN,notATime, ...
                notATime,NaN,{struct()},{table()},""};
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
                    if trial.State == "RUNNING"
                        obj.fail(trial, exception.message);
                    end
                    caught=any(options.Catch==string(exception.identifier)) || ...
                        any(options.Catch=="*");
                    if ~caught
                        rethrow(exception);
                    end
                end
                frozen=obj.freezeTrial(trial.Number);
                for callback=callbacks
                    callback{1}(obj,frozen);
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
                isfinite(obj.TrialTable.Value);
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
            valid = all(isfinite(values),2);
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
            result = obj.bestTrial();
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
            result=obj.paretoFront();
        end

        function result = trials_dataframe(obj)
            result=obj.TrialTable;
        end

        function value = direction(obj)
            if numel(obj.Directions)~=1
                error("radia:optuna:MultiObjectiveDirection", ...
                    "Use Directions for a multi-objective study.");
            end
            value=obj.Directions(1);
        end

        function result = get_trials(obj, states)
            if nargin < 2
                result = obj.trials();
            else
                result = obj.trials(states);
            end
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
            nowTime = datetime("now", "TimeZone", "local");
            ages = seconds(nowTime - obj.TrialTable.StartTime);
            rows = obj.TrialTable.State == "RUNNING" & ...
                isfinite(ages) & ages >= timeoutSeconds;
            if ~any(rows)
                recovered = obj.TrialTable([],:);
                return
            end
            message = options.Message;
            if strlength(message) == 0
                message = "Recovered stale RUNNING trial after timeout.";
            end
            trialNumbers = obj.TrialTable.TrialNumber(rows);
            obj.TrialTable.State(rows) = "FAIL";
            obj.TrialTable.Value(rows) = NaN;
            obj.TrialTable.EndTime(rows) = nowTime;
            obj.TrialTable.Duration_s(rows) = ages(rows);
            obj.TrialTable.ErrorMessage(rows) = message;
            obj.ObjectiveTable(ismember( ...
                obj.ObjectiveTable.TrialNumber, trialNumbers), :) = [];
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

    methods (Hidden=true)
        function stopWhenOptimizing(obj)
            if obj.InOptimize
                obj.StopRequested=true;
            end
        end

        function initializeTables(obj)
            obj.TrialTable = table('Size', [0, 9], ...
                'VariableTypes', {'double','string','double','datetime','datetime', ...
                'double','cell','cell','string'}, ...
                'VariableNames', {'TrialNumber','State','Value','StartTime','EndTime', ...
                'Duration_s','Params','IntermediateValues','ErrorMessage'});
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
            obj.TrialTable.StartTime.TimeZone = "local";
            obj.TrialTable.EndTime.TimeZone = "local";
            obj.IntermediateTable.Timestamp.TimeZone = "local";
            obj.SamplerStateTable.Timestamp.TimeZone = "local";
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
                complete = obj.TrialTable.State == "COMPLETE" & isfinite(obj.TrialTable.Value);
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
            row=obj.TrialTable.TrialNumber==trialNumber;
            if sum(row)~=1
                error("radia:optuna:UnknownTrial", ...
                    "Trial %d does not identify exactly one row.",trialNumber);
            end
            selected=obj.TrialTable(row,:);
            values=NaN;
            objectiveRows=obj.ObjectiveTable.TrialNumber==trialNumber;
            if any(objectiveRows)
                objectives=sortrows(obj.ObjectiveTable(objectiveRows,:), ...
                    "ObjectiveIndex");
                values=reshape(objectives.Value,1,[]);
            elseif isfinite(selected.Value)
                values=selected.Value;
            end
            distributions=struct();
            parameterRows=find(obj.ParamTable.TrialNumber==trialNumber)';
            for index=parameterRows
                key=matlab.lang.makeValidName(obj.ParamTable.Name(index));
                distributions.(key)= ...
                    radia.optuna.internal.DistributionCodec.decode( ...
                    obj.ParamTable.Kind(index), ...
                    obj.ParamTable.Distribution(index));
            end
            userAttrs=obj.attributesForTrial(obj.UserAttrTable,trialNumber);
            systemAttrs=obj.attributesForTrial(obj.SystemAttrTable,trialNumber);
            intermediateRows=obj.IntermediateTable.TrialNumber==trialNumber;
            intermediate=obj.IntermediateTable(intermediateRows, ...
                ["Step","Value","Timestamp"]);
            [constraintPresent,constraints]=obj.constraintRecord(trialNumber);
            frozen=radia.optuna.FrozenTrial(Number=trialNumber, ...
                State=selected.State,Values=values,Params=selected.Params{1}, ...
                Distributions=distributions, ...
                IntermediateValues=intermediate,UserAttrs=userAttrs, ...
                SystemAttrs=systemAttrs,Constraints=constraints, ...
                ConstraintPresent=constraintPresent, ...
                DatetimeStart=selected.StartTime, ...
                DatetimeComplete=selected.EndTime, ...
                ErrorMessage=selected.ErrorMessage);
        end

        function addTrial(obj,frozen)
            arguments
                obj
                frozen (1,1) radia.optuna.FrozenTrial
            end
            if ismember(frozen.State,["RUNNING","WAITING"])
                error("radia:optuna:AddTrialState", ...
                    "addTrial requires a finished trial.");
            end
            if frozen.State=="COMPLETE" && ...
                    (numel(frozen.Values)~=numel(obj.Directions) || ...
                    any(~isfinite(frozen.Values)))
                error("radia:optuna:AddTrialValue", ...
                    "A COMPLETE trial must have one finite value per direction.");
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
            if isnat(startTime), startTime=datetime("now","TimeZone","local"); end
            endTime=frozen.DatetimeComplete;
            if isnat(endTime), endTime=startTime; end
            duration=seconds(endTime-startTime);
            scalarValue=NaN;
            if frozen.State=="COMPLETE", scalarValue=frozen.Values(1); end
            obj.TrialTable(end+1,:)={number,frozen.State,scalarValue, ...
                startTime,endTime,duration,{frozen.Params}, ...
                {frozen.IntermediateValues},frozen.ErrorMessage};
            if frozen.State=="COMPLETE"
                obj.ObjectiveTable=[obj.ObjectiveTable;table( ...
                    repmat(number,numel(frozen.Values),1), ...
                    (1:numel(frozen.Values))',reshape(frozen.Values,[],1), ...
                    'VariableNames',obj.ObjectiveTable.Properties.VariableNames)];
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
            rows = obj.ParamTable.TrialNumber == trial.Number & obj.ParamTable.Name == name;
            obj.ParamTable(rows,:) = [];
            numeric = NaN;
            text = "";
            if isnumeric(value) && isscalar(value)
                numeric = double(value);
            else
                text = string(jsonencode(value));
            end
            obj.ParamTable(end+1,:) = {trial.Number, name, kind, numeric, text, distribution};
            obj.updateTrialSnapshot(trial);
            obj.persist();
        end

        function recordIntermediate(obj, trial, value, step)
            rows = obj.IntermediateTable.TrialNumber == trial.Number & ...
                obj.IntermediateTable.Step == step;
            obj.IntermediateTable(rows,:) = [];
            obj.IntermediateTable(end+1,:) = {trial.Number, step, value, ...
                datetime("now", "TimeZone", "local")};
            obj.updateTrialSnapshot(trial);
            obj.persist();
        end

        function recordUserAttribute(obj, trial, name, value)
            rows = obj.UserAttrTable.TrialNumber == trial.Number & ...
                obj.UserAttrTable.Name == name;
            obj.UserAttrTable(rows,:) = [];
            obj.UserAttrTable(end+1,:) = {trial.Number, name, string(jsonencode(value))};
            obj.persist();
        end

        function recordSystemAttribute(obj,trial,name,value)
            rows=obj.SystemAttrTable.TrialNumber==trial.Number & ...
                obj.SystemAttrTable.Name==name;
            obj.SystemAttrTable(rows,:)=[];
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
            rows = obj.TrialTable.TrialNumber == trial.Number;
            if sum(rows) ~= 1
                error("radia:optuna:UnknownTrial", ...
                    "Trial %d does not identify exactly one row in this study.", ...
                    trial.Number);
            end
            storedState = obj.TrialTable.State(rows);
            if trial.State ~= "RUNNING" || storedState ~= "RUNNING"
                if trial.State ~= "RUNNING"
                    finishedState = trial.State;
                else
                    finishedState = storedState;
                end
                error("radia:optuna:TrialState", "Trial %d is already %s.", ...
                    trial.Number, finishedState);
            end
            endTime = datetime("now", "TimeZone", "local");
            trial.markFinished(state, value, endTime, message);
            elapsed = seconds(endTime - trial.StartTime);
            obj.TrialTable.State(rows) = state;
            obj.TrialTable.Value(rows) = value(1);
            obj.ObjectiveTable(obj.ObjectiveTable.TrialNumber == trial.Number,:) = [];
            if state == "COMPLETE"
                values = reshape(double(value),[],1);
                obj.ObjectiveTable = [obj.ObjectiveTable; table( ...
                    repmat(trial.Number,numel(values),1), (1:numel(values))', values, ...
                    'VariableNames',obj.ObjectiveTable.Properties.VariableNames)];
            end
            obj.TrialTable.EndTime(rows) = endTime;
            obj.TrialTable.Duration_s(rows) = elapsed;
            obj.TrialTable.ErrorMessage(rows) = message;
            obj.updateTrialSnapshot(trial);
            if ismethod(obj.Sampler, "afterTrial")
                try
                    obj.Sampler.afterTrial(obj, trial);
                catch exception
                    % Optuna commits the requested trial state in a finally
                    % block even when sampler post-processing raises.
                    obj.persist();
                    rethrow(exception)
                end
            end
            obj.persist();
            obj.reportProgress("TRIAL_FINISHED", trial);
        end

        function updateTrialSnapshot(obj, trial)
            rows = obj.TrialTable.TrialNumber == trial.Number;
            obj.TrialTable.Params(rows) = {trial.Params};
            obj.TrialTable.IntermediateValues(rows) = {trial.IntermediateValues};
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
            rows = obj.TrialTable.TrialNumber == trial.Number;
            if sum(rows) ~= 1
                error("radia:optuna:UnknownTrial", ...
                    "Trial %d does not identify exactly one row in this study.", ...
                    trial.Number);
            end
            storedState = obj.TrialTable.State(rows);
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
            for index=1:height(obj.TrialTable)
                number=obj.TrialTable.TrialNumber(index);
                if obj.TrialTable.State(index)=="WAITING"
                    rows=obj.QueueParamTable.TrialNumber==number;
                    candidate=struct();
                    for row=find(rows)'
                        candidate.(obj.QueueParamTable.Name(row))= ...
                            obj.QueueParamTable.Value{row};
                    end
                else
                    candidate=obj.TrialTable.Params{index};
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
            required=["kind","low","high","log","step","choices"];
            if ~isstruct(distribution) || ~isscalar(distribution) || ...
                    ~all(isfield(distribution,required))
                error("radia:optuna:AddTrialDistribution", ...
                    "Imported distribution for '%s' is invalid.",name);
            end
            numeric=NaN;
            textValue="";
            if isnumeric(value) && isscalar(value)
                numeric=double(value);
            else
                textValue=string(jsonencode(value));
            end
            obj.ParamTable(end+1,:)={number,name,string(distribution.kind), ...
                numeric,textValue, ...
                radia.optuna.internal.DistributionCodec.encode(distribution)};
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
