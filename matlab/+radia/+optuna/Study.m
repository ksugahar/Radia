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
        ConstraintTable table
        ObjectiveTable table
        SamplerStateTable table
        ProgressFcn = []
        UserAttrs struct = struct()
    end

    properties (Access=private)
        NextTrialNumber (1,1) double = 0
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
                obj.Sampler = radia.optuna.RandomSampler(0);
            end
            if isempty(obj.Pruner)
                obj.Pruner = radia.optuna.MedianPruner();
            end
            obj.initializeTables();
            if strlength(obj.StoragePath) > 0 && isfile(obj.StoragePath)
                obj.loadState();
            end
        end

        function trial = ask(obj)
            trial = radia.optuna.Trial(obj, obj.NextTrialNumber);
            obj.NextTrialNumber = obj.NextTrialNumber + 1;
            if ismethod(obj.Sampler, "beforeTrial")
                obj.Sampler.beforeTrial(obj, trial);
            end
            obj.TrialTable(end+1,:) = {trial.Number, "RUNNING", NaN, ...
                trial.StartTime, NaT, NaN, {trial.Params}, ...
                {trial.IntermediateValues}, ""};
            obj.persist();
            obj.reportProgress("TRIAL_STARTED", trial);
        end

        function tell(obj, trial, value)
            arguments
                obj
                trial (1,1) radia.optuna.Trial
                value double
            end
            value = reshape(double(value), 1, []);
            if numel(value) ~= numel(obj.Directions) || any(~isfinite(value))
                error("radia:optuna:Value", ...
                    "Trial value must contain one finite value per direction.");
            end
            obj.finishTrial(trial, "COMPLETE", value, "");
        end

        function fail(obj, trial, message)
            arguments
                obj
                trial (1,1) radia.optuna.Trial
                message (1,1) string
            end
            obj.finishTrial(trial, "FAIL", NaN, message);
        end

        function results = optimize(obj, objective, nTrials)
            arguments
                obj
                objective (1,1) function_handle
                nTrials (1,1) double {mustBeInteger, mustBeNonnegative}
            end
            for k = 1:nTrials
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
                    rethrow(exception);
                end
            end
            results = obj.TrialTable;
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
                result = table('Size',[0,3], ...
                    'VariableTypes',{'double','cell','cell'}, ...
                    'VariableNames',{'TrialNumber','Values','Params'});
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

        function result = get_trials(obj, states)
            if nargin < 2
                result = obj.trials();
            else
                result = obj.trials(states);
            end
        end

        function values = intermediateValuesAtStep(obj, step)
            rows = obj.IntermediateTable.Step == step;
            trialNumbers = obj.IntermediateTable.TrialNumber(rows);
            completed = obj.TrialTable.State == "COMPLETE";
            completedNumbers = obj.TrialTable.TrialNumber(completed);
            values = obj.IntermediateTable.Value(rows & ...
                ismember(trialNumbers, completedNumbers));
        end

        function save(obj)
            if strlength(obj.StoragePath) == 0
                error("radia:optuna:Storage", "StoragePath is empty.");
            end
            StudyData = struct( ...
                "Name", obj.Name, ...
                "Directions", obj.Directions, ...
                "NextTrialNumber", obj.NextTrialNumber, ...
                "TrialTable", obj.TrialTable, ...
                "ParamTable", obj.ParamTable, ...
                "IntermediateTable", obj.IntermediateTable, ...
                "UserAttrTable", obj.UserAttrTable, ...
                "ConstraintTable", obj.ConstraintTable, ...
                "UserAttrs", obj.UserAttrs, ...
                "ObjectiveTable", obj.ObjectiveTable, ...
                "SamplerStateTable", obj.SamplerStateTable); %#ok<NASGU>
            folder = fileparts(obj.StoragePath);
            if strlength(folder) > 0 && ~isfolder(folder)
                mkdir(folder);
            end
            builtin("save", obj.StoragePath, "StudyData", "-mat");
        end
    end

    methods (Hidden=true)
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
            obj.ConstraintTable = table('Size', [0, 3], ...
                'VariableTypes', {'double','double','double'}, ...
                'VariableNames', {'TrialNumber','ConstraintIndex','Value'});
            obj.ObjectiveTable = table('Size', [0, 3], ...
                'VariableTypes', {'double','double','double'}, ...
                'VariableNames', {'TrialNumber','ObjectiveIndex','Value'});
            obj.SamplerStateTable = table('Size', [0, 7], ...
                'VariableTypes', ...
                {'string','string','double','double','double','cell','datetime'}, ...
                'VariableNames', {'Sampler','Schema','Revision','TrialNumber', ...
                'Generation','State','Timestamp'});
            obj.TrialTable.StartTime.TimeZone = "local";
            obj.TrialTable.EndTime.TimeZone = "local";
            obj.IntermediateTable.Timestamp.TimeZone = "local";
            obj.SamplerStateTable.Timestamp.TimeZone = "local";
        end

        function loadState(obj)
            loaded = builtin("load", obj.StoragePath, "StudyData");
            data = loaded.StudyData;
            obj.Name = string(data.Name);
            obj.Directions = string(data.Directions);
            obj.NextTrialNumber = data.NextTrialNumber;
            obj.TrialTable = data.TrialTable;
            obj.ParamTable = data.ParamTable;
            obj.IntermediateTable = data.IntermediateTable;
            obj.UserAttrTable = data.UserAttrTable;
            if isfield(data, "ConstraintTable")
                obj.ConstraintTable = data.ConstraintTable;
            end
            if isfield(data, "UserAttrs")
                obj.UserAttrs = data.UserAttrs;
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

        function recordConstraints(obj, trial, values)
            values = reshape(double(values), [], 1);
            if any(isnan(values))
                error("radia:optuna:Constraints", ...
                    "Constraint values must not contain NaN.");
            end
            obj.ConstraintTable( ...
                obj.ConstraintTable.TrialNumber == trial.Number, :) = [];
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
            rows = obj.ConstraintTable.TrialNumber == trialNumber;
            selected = sortrows(obj.ConstraintTable(rows,:), "ConstraintIndex");
            values = reshape(selected.Value, 1, []);
        end

        function finishTrial(obj, trial, state, value, message)
            if trial.State ~= "RUNNING"
                error("radia:optuna:TrialState", "Trial %d is already %s.", ...
                    trial.Number, trial.State);
            end
            endTime = datetime("now", "TimeZone", "local");
            trial.markFinished(state, value, endTime, message);
            elapsed = seconds(endTime - trial.StartTime);
            rows = obj.TrialTable.TrialNumber == trial.Number;
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
                obj.Sampler.afterTrial(obj, trial);
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
                snapshot.best_trial = obj.bestTrial();
            else
                snapshot.best_trial = obj.TrialTable([],:);
            end
            obj.ProgressFcn(snapshot);
        end
    end
end
