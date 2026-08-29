classdef OptimizationSession < handle
    %OPTIMIZATIONSESSION Persistent lifecycle for MATLAB/Simulink studies.
    %
    % The session is the unit operated by a script or a Level-2 MATLAB
    % S-Function. It keeps Study as the Optuna state machine and adds the
    % explicit start/pause/resume/cancel/select/apply workflow students need.

    properties (SetAccess=private)
        Parameters
        Options = radia.optuna.OptimizeOptions()
        ObjectiveFcn = []
        Study = []
        State (1,1) string = "configured"
        SelectedTrialNumber (1,1) double = NaN
        StartedAt datetime = NaT
        UpdatedAt datetime = NaT
        StopReason (1,1) string = ""
        CheckpointRevision (1,1) double = 0
        ModelName (1,1) string = ""
    end

    properties (Access=private)
        Searchable
        Fixed
        Names string = strings(1,0)
        Distributions cell = cell(1,0)
        ProgressFcn = []
    end

    properties (Constant, Access=private)
        SessionSchema = "radia.optuna.session"
        SessionVersion = 1
    end

    methods
        function obj = OptimizationSession(objectiveFcn, parameters, options, ...
                sessionOptions)
            arguments
                objectiveFcn (1,1) function_handle
                parameters
                options (1,1) radia.optuna.OptimizeOptions = ...
                    radia.optuna.OptimizeOptions()
                sessionOptions.ModelName (1,1) string = ""
                sessionOptions.ProgressFcn = []
            end

            obj.ObjectiveFcn = objectiveFcn;
            obj.Options = options;
            obj.ModelName = sessionOptions.ModelName;
            obj.ProgressFcn = sessionOptions.ProgressFcn;
            if ~isempty(obj.ProgressFcn) && ...
                    ~isa(obj.ProgressFcn, "function_handle")
                error("radia:optuna:ProgressFcn", ...
                    "ProgressFcn must be a function handle.");
            end
            obj.Parameters = obj.normalizeParameters(parameters);
            obj.Searchable = obj.Parameters([obj.Parameters.Free]);
            obj.Fixed = obj.Parameters(~[obj.Parameters.Free]);
            if isempty(obj.Searchable)
                error("radia:optuna:SessionParameters", ...
                    "Every parameter has Free=false; there is nothing to search.");
            end
            arrayfun(@(p) p.mustBeSearchable(), obj.Searchable);
            for index = 1:numel(obj.Fixed)
                if isempty(obj.Fixed(index).Value)
                    error("radia:optuna:SessionParameters", ...
                        "Fixed parameter '%s' has no Value.", ...
                        obj.Fixed(index).Name);
                end
            end
            obj.Names = [obj.Searchable.Name];
            obj.Distributions = arrayfun( ...
                @(p) p.distribution(), obj.Searchable, UniformOutput=false);

            storageExists = strlength(options.StoragePath) > 0 && ...
                (isfile(options.StoragePath) || ...
                isfile(options.StoragePath + ".bak"));
            if storageExists && ~options.Resume
                error("radia:optuna:Resume", ...
                    "Storage '%s' already exists and Resume=false.", ...
                    options.StoragePath);
            end
            obj.Study = obj.buildStudy();

            if storageExists
                obj.restoreSessionMetadata();
                recovered = obj.Study.recoverStaleRunning(0, ...
                    Message="Recovered orphaned RUNNING trial while " + ...
                    "restoring OptimizationSession.");
                if ~isempty(recovered)
                    obj.State = "paused";
                    obj.StopReason = "recovered-stale-running";
                elseif obj.State == "running"
                    obj.State = "paused";
                    obj.StopReason = "restored-running-session-as-paused";
                end
                if obj.State == "completed" && ~obj.atTrialLimit()
                    obj.State = "paused";
                    obj.StopReason = "trial-budget-extended";
                end
            end
            if obj.atTrialLimit()
                obj.State = "completed";
                obj.StopReason = "max-trials";
            end
            obj.UpdatedAt = datetime("now");
            obj.persist();
        end

        function start(obj)
            if obj.State == "running"
                return
            end
            if ~ismember(obj.State, ["configured", "paused"])
                error("radia:optuna:SessionState", ...
                    "Cannot start a session in state '%s'.", obj.State);
            end
            if obj.atTrialLimit()
                obj.State = "completed";
                obj.StopReason = "max-trials";
                obj.persist();
                return
            end
            obj.State = "running";
            obj.StopReason = "";
            if isnat(obj.StartedAt)
                obj.StartedAt = datetime("now");
            end
            obj.touch();
        end

        function pause(obj)
            if obj.State == "paused"
                return
            end
            obj.requireState("running", "pause");
            obj.State = "paused";
            obj.StopReason = "paused-by-user";
            obj.touch();
        end

        function resume(obj)
            if obj.State == "running"
                return
            end
            obj.requireState("paused", "resume");
            obj.State = "running";
            obj.StopReason = "";
            obj.touch();
        end

        function cancel(obj)
            if obj.State == "cancelled"
                return
            end
            if ~ismember(obj.State, ["configured", "running", "paused"])
                error("radia:optuna:SessionState", ...
                    "Cannot cancel a session in state '%s'.", obj.State);
            end
            obj.State = "cancelled";
            obj.StopReason = "cancelled-by-user";
            obj.touch();
        end

        function frozen = runNext(obj)
            obj.requireState("running", "runNext");
            if obj.atTrialLimit()
                obj.finish("max-trials");
                frozen = radia.optuna.FrozenTrial.empty(1,0);
                return
            end
            if obj.elapsedSeconds() >= obj.Options.MaxTime
                obj.finish("max-time");
                frozen = radia.optuna.FrozenTrial.empty(1,0);
                return
            end

            trial = obj.Study.ask();
            values = obj.valuesForTrial(trial);
            context = radia.optuna.ObjectiveContext(trial);
            failure = [];
            try
                result = obj.callObjective(values, context);
                normalized = obj.normalizeResult(result, context, values);
                if ~isempty(normalized.Constraints)
                    obj.Study.recordConstraints( ...
                        trial, normalized.Constraints);
                end
                frozen = obj.Study.tell(trial, normalized.Objectives);
            catch problem
                if problem.identifier == "radia:optuna:TrialPruned"
                    frozen = obj.Study.tell(trial, State="PRUNED");
                else
                    obj.Study.fail(trial, string(problem.message));
                    frozen = obj.trialByNumber(trial.Number);
                    failure = problem;
                end
            end

            obj.callHooks(frozen);
            if obj.atTrialLimit()
                obj.finish("max-trials");
            elseif obj.elapsedSeconds() >= obj.Options.MaxTime
                obj.finish("max-time");
            else
                obj.touch();
            end
            if ~isempty(failure) && ~obj.Options.CatchObjectiveErrors
                rethrow(failure)
            end
        end

        function trials = run(obj, count)
            arguments
                obj
                count (1,1) double = Inf
            end
            if count <= 0 || (~isinf(count) && count ~= round(count))
                error("radia:optuna:SessionRunCount", ...
                    "count must be a positive integer or Inf.");
            end
            if obj.State == "configured"
                obj.start();
            end
            trials = radia.optuna.FrozenTrial.empty(1,0);
            index = 0;
            while obj.State == "running" && index < count
                frozen = obj.runNext();
                if ~isempty(frozen)
                    index = index + 1;
                    trials(index) = frozen;
                end
            end
        end

        function selectTrial(obj, trialNumber)
            arguments
                obj
                trialNumber (1,1) double {mustBeInteger, mustBeNonnegative}
            end
            frozen = obj.trialByNumber(trialNumber);
            if isempty(frozen) || frozen.State ~= "COMPLETE"
                error("radia:optuna:SessionSelection", ...
                    "Trial %d is not a completed trial.", trialNumber);
            end
            obj.SelectedTrialNumber = trialNumber;
            obj.touch();
        end

        function values = selectedParameters(obj)
            if isnan(obj.SelectedTrialNumber)
                error("radia:optuna:SessionSelection", ...
                    "No trial is selected.");
            end
            frozen = obj.trialByNumber(obj.SelectedTrialNumber);
            values = frozen.Params;
            for index = 1:numel(obj.Fixed)
                name = matlab.lang.makeValidName(obj.Fixed(index).Name);
                values.(name) = obj.Fixed(index).Value;
            end
        end

        function values = applySelectedToModel(obj, options)
            arguments
                obj
                options.ModelName (1,1) string = obj.ModelName
                options.Target (1,1) string {mustBeMember(options.Target, ...
                    ["auto", "model", "base"])} = "auto"
            end
            values = obj.selectedParameters();
            target = options.Target;
            if target == "auto"
                if strlength(options.ModelName) > 0
                    target = "model";
                else
                    target = "base";
                end
            end

            if target == "model"
                if strlength(options.ModelName) == 0
                    error("radia:optuna:SessionModel", ...
                        "ModelName is required when Target='model'.");
                end
                if exist("load_system", "file") ~= 2
                    error("radia:optuna:SimulinkRequired", ...
                        "Applying parameters to a model requires Simulink.");
                end
                if ~bdIsLoaded(options.ModelName)
                    load_system(options.ModelName);
                end
                workspace = get_param(options.ModelName, "ModelWorkspace");
            else
                workspace = [];
            end

            names = string(fieldnames(values));
            for index = 1:numel(names)
                name = names(index);
                value = values.(name);
                if target == "model" && workspace.hasVariable(char(name))
                    workspace.assignin(char(name), value);
                elseif target == "model"
                    assignin("base", char(name), value);
                else
                    assignin("base", char(name), value);
                end
            end
            obj.ModelName = options.ModelName;
            obj.touch();
        end

        function data = snapshot(obj)
            tableValue = obj.Study.TrialTable;
            states = tableValue.State;
            bestValue = NaN;
            bestTrialNumber = NaN;
            if isscalar(obj.Study.Directions)
                try
                    best = obj.Study.bestSolution();
                    if best.available
                        bestValue = best.value;
                        bestTrialNumber = best.trial_number;
                    end
                catch
                    % No feasible completed trial yet.
                end
            end
            running = tableValue.TrialNumber(states == "RUNNING");
            currentTrial = NaN;
            if ~isempty(running)
                currentTrial = running(end);
            end
            pareto = obj.Study.paretoFront();
            data = struct( ...
                "schema", "radia.optuna.session-snapshot.v1", ...
                "state", obj.State, ...
                "stopReason", obj.StopReason, ...
                "attempted", height(tableValue), ...
                "complete", sum(states == "COMPLETE"), ...
                "pruned", sum(states == "PRUNED"), ...
                "failed", sum(states == "FAIL"), ...
                "running", sum(states == "RUNNING"), ...
                "currentTrial", currentTrial, ...
                "bestValue", bestValue, ...
                "bestTrialNumber", bestTrialNumber, ...
                "paretoCount", height(pareto), ...
                "selectedTrialNumber", obj.SelectedTrialNumber, ...
                "elapsedSeconds", obj.elapsedSeconds(), ...
                "checkpointRevision", obj.CheckpointRevision, ...
                "storagePath", obj.Options.StoragePath, ...
                "sampler", obj.Options.Sampler, ...
                "pruner", obj.Options.Pruner, ...
                "seed", obj.Options.Seed);
        end

        function path = save(obj, path)
            arguments
                obj
                path (1,1) string = ""
            end
            if strlength(obj.Options.StoragePath) > 0
                obj.Study.save();
            end
            if strlength(path) == 0
                path = obj.sessionPath();
            end
            if strlength(path) == 0
                error("radia:optuna:SessionStorage", ...
                    "Set Options.StoragePath or pass a session path.");
            end
            folder = fileparts(path);
            if strlength(folder) > 0 && ~isfolder(folder)
                mkdir(folder);
            end
            obj.CheckpointRevision = obj.CheckpointRevision + 1;
            obj.UpdatedAt = datetime("now");
            SessionData = obj.sessionData(); %#ok<NASGU>
            builtin("save", path, "SessionData", "-mat");
        end
    end

    methods (Static)
        function obj = load(path, options)
            arguments
                path (1,1) string
                options.ObjectiveFcn = []
                options.ModelName (1,1) string = ""
            end
            loaded = builtin("load", path, "SessionData", "-mat");
            if ~isfield(loaded, "SessionData")
                error("radia:optuna:SessionStorage", ...
                    "'%s' has no SessionData.", path);
            end
            data = loaded.SessionData;
            if ~isfield(data, "Schema") || ...
                    string(data.Schema) ~= "radia.optuna.session" || ...
                    data.Version ~= 1
                error("radia:optuna:SessionStorage", ...
                    "'%s' is not a supported OptimizationSession.", path);
            end
            objectiveFcn = options.ObjectiveFcn;
            if isempty(objectiveFcn)
                objectiveFcn = data.ObjectiveFcn;
            end
            modelName = options.ModelName;
            if strlength(modelName) == 0 && isfield(data, "ModelName")
                modelName = string(data.ModelName);
            end
            obj = radia.optuna.OptimizationSession( ...
                objectiveFcn, data.Parameters, data.Options, ...
                ModelName=modelName);
        end
    end

    methods (Access=private)
        function parameters = normalizeParameters(~, parameters)
            if isstruct(parameters)
                parameters = ...
                    radia.optuna.OptimizationParameter.fromStruct(parameters);
            end
            if ~isa(parameters, "radia.optuna.OptimizationParameter")
                error("radia:optuna:SessionParameters", ...
                    "Parameters must be OptimizationParameter objects or " + ...
                    "a struct of bounds.");
            end
            parameters = reshape(parameters, 1, []);
        end

        function study = buildStudy(obj)
            samplerName = obj.Options.Sampler;
            if samplerName == "auto"
                spec = struct( ...
                    "fixed_numeric", ...
                        all([obj.Searchable.Type] ~= "categorical"), ...
                    "dimensions", numel(obj.Searchable), ...
                    "has_constraints", ~isempty(obj.Options.ConstraintFcn), ...
                    "constraints_declared", ...
                        ~isempty(obj.Options.ConstraintFcn), ...
                    "has_categorical", ...
                        any([obj.Searchable.Type] == "categorical"), ...
                    "is_conditional", false);
                samplerName = ...
                    radia.optuna.internal.AutoSamplerPolicy.choose( ...
                        spec, numel(obj.Options.Directions), ...
                        obj.Options.MaxTrials);
            end
            sampler = radia.optuna.internal.samplerFromName( ...
                samplerName, obj.Options.Seed);
            settings = {"Sampler", sampler, ...
                "Directions", obj.Options.Directions, ...
                "ProgressFcn", obj.ProgressFcn};
            pruner = radia.optuna.internal.prunerFromName( ...
                obj.Options.Pruner);
            if ~isempty(pruner)
                settings = [settings, {"Pruner", pruner}];
            end
            if strlength(obj.Options.StudyName) > 0
                settings = [settings, {"Name", obj.Options.StudyName}];
            end
            if strlength(obj.Options.StoragePath) > 0
                settings = [settings, ...
                    {"StoragePath", obj.Options.StoragePath}];
            else
                settings = [settings, {"AutoSave", false}];
            end
            study = radia.optuna.Study(settings{:});
        end

        function values = valuesForTrial(obj, trial)
            values = struct();
            for index = 1:numel(obj.Names)
                name = matlab.lang.makeValidName(obj.Names(index));
                values.(name) = ...
                    radia.optuna.internal.suggestFromDistribution( ...
                        trial, obj.Names(index), ...
                        obj.Distributions{index});
            end
            for index = 1:numel(obj.Fixed)
                name = matlab.lang.makeValidName(obj.Fixed(index).Name);
                values.(name) = obj.Fixed(index).Value;
            end
        end

        function result = callObjective(obj, values, context)
            argumentCount = nargin(obj.ObjectiveFcn);
            if argumentCount == 1
                result = obj.ObjectiveFcn(values);
            else
                result = obj.ObjectiveFcn(values, context);
            end
        end

        function normalized = normalizeResult(obj, result, context, values)
            normalized = struct( ...
                "Objectives", double.empty(1,0), ...
                "Constraints", double.empty(1,0));
            if isnumeric(result)
                normalized.Objectives = reshape(double(result), 1, []);
            elseif isstruct(result) && isscalar(result)
                if ~isfield(result, "Objectives")
                    error("radia:optuna:ObjectiveResult", ...
                        "A result struct must contain Objectives.");
                end
                normalized.Objectives = ...
                    reshape(double(result.Objectives), 1, []);
                if isfield(result, "InequalityConstraints")
                    normalized.Constraints = [normalized.Constraints, ...
                        reshape(double( ...
                        result.InequalityConstraints), 1, [])];
                end
                if isfield(result, "EqualityConstraints")
                    normalized.Constraints = [normalized.Constraints, ...
                        abs(reshape(double( ...
                        result.EqualityConstraints), 1, []))];
                end
                if isfield(result, "IntermediateValues")
                    obj.reportIntermediate( ...
                        context, result.IntermediateValues);
                end
                if isfield(result, "UserData")
                    context.UserData = result.UserData;
                end
            else
                error("radia:optuna:ObjectiveResult", ...
                    "Objective must return a numeric vector or a scalar struct.");
            end

            if ~isempty(obj.Options.ConstraintFcn)
                normalized.Constraints = [normalized.Constraints, ...
                    reshape(double(obj.Options.ConstraintFcn(values)), 1, [])];
            end
            if numel(normalized.Objectives) ~= ...
                    numel(obj.Options.Directions) || ...
                    any(~isfinite(normalized.Objectives))
                error("radia:optuna:ObjectiveResult", ...
                    "Objective returned %d values for %d directions, or " + ...
                    "returned a non-finite value.", ...
                    numel(normalized.Objectives), ...
                    numel(obj.Options.Directions));
            end
            if any(~isfinite(normalized.Constraints))
                error("radia:optuna:ObjectiveResult", ...
                    "Constraints must be finite.");
            end
            if ~isempty(context.Trial.IntermediateValues) && ...
                    context.shouldPrune()
                context.prune();
            end
        end

        function reportIntermediate(~, context, intermediate)
            if istable(intermediate)
                if ~all(ismember(["Step", "Value"], ...
                        string(intermediate.Properties.VariableNames)))
                    error("radia:optuna:IntermediateValues", ...
                        "IntermediateValues table needs Step and Value.");
                end
                steps = intermediate.Step;
                values = intermediate.Value;
            elseif isnumeric(intermediate) && size(intermediate,2) == 2
                steps = intermediate(:,1);
                values = intermediate(:,2);
            else
                error("radia:optuna:IntermediateValues", ...
                    "IntermediateValues must be an N-by-2 numeric array " + ...
                    "or a table with Step and Value.");
            end
            for index = 1:numel(steps)
                context.report(values(index), steps(index));
            end
        end

        function callHooks(obj, frozen)
            snapshotValue = obj.snapshot();
            optimValues = struct( ...
                "trial", frozen.Number + 1, ...
                "fval", frozen.Value, ...
                "bestfval", snapshotValue.bestValue, ...
                "bestx", struct(), ...
                "funccount", snapshotValue.attempted, ...
                "elapsedtime", snapshotValue.elapsedSeconds, ...
                "stalled", 0, ...
                "state", frozen.State, ...
                "sampler", obj.Options.Sampler, ...
                "pruned", snapshotValue.pruned, ...
                "failed", snapshotValue.failed);
            hooks = [obj.Options.OutputFcn, obj.Options.PlotFcn];
            for index = 1:numel(hooks)
                stop = hooks{index}(frozen.Params, optimValues, "iter");
                if ~isempty(stop) && islogical(stop) && any(stop)
                    obj.pause();
                    obj.StopReason = "stopped-by-callback";
                    break
                end
            end
        end

        function frozen = trialByNumber(obj, number)
            trials = obj.Study.get_trials();
            if isempty(trials)
                frozen = radia.optuna.FrozenTrial.empty(1,0);
                return
            end
            index = find([trials.Number] == number, 1);
            if isempty(index)
                frozen = radia.optuna.FrozenTrial.empty(1,0);
            else
                frozen = trials(index);
            end
        end

        function result = atTrialLimit(obj)
            result = height(obj.Study.TrialTable) >= obj.Options.MaxTrials;
        end

        function finish(obj, reason)
            obj.State = "completed";
            obj.StopReason = reason;
            obj.touch();
        end

        function requireState(obj, required, operation)
            if obj.State ~= required
                error("radia:optuna:SessionState", ...
                    "Cannot %s a session in state '%s'; required '%s'.", ...
                    operation, obj.State, required);
            end
        end

        function secondsValue = elapsedSeconds(obj)
            if isnat(obj.StartedAt)
                secondsValue = 0;
            else
                secondsValue = seconds(datetime("now") - obj.StartedAt);
            end
        end

        function touch(obj)
            obj.UpdatedAt = datetime("now");
            obj.persist();
        end

        function persist(obj)
            if strlength(obj.Options.StoragePath) == 0
                return
            end
            obj.save();
        end

        function path = sessionPath(obj)
            if strlength(obj.Options.StoragePath) == 0
                path = "";
            else
                path = obj.Options.StoragePath + ".session.mat";
            end
        end

        function data = sessionData(obj)
            data = struct( ...
                "Schema", obj.SessionSchema, ...
                "Version", obj.SessionVersion, ...
                "Parameters", obj.Parameters, ...
                "Options", obj.Options, ...
                "ObjectiveFcn", obj.ObjectiveFcn, ...
                "State", obj.State, ...
                "SelectedTrialNumber", obj.SelectedTrialNumber, ...
                "StartedAt", obj.StartedAt, ...
                "UpdatedAt", obj.UpdatedAt, ...
                "StopReason", obj.StopReason, ...
                "CheckpointRevision", obj.CheckpointRevision, ...
                "ModelName", obj.ModelName);
        end

        function restoreSessionMetadata(obj)
            path = obj.sessionPath();
            if strlength(path) == 0 || ~isfile(path)
                return
            end
            loaded = builtin("load", path, "SessionData", "-mat");
            if ~isfield(loaded, "SessionData")
                error("radia:optuna:SessionStorage", ...
                    "'%s' has no SessionData.", path);
            end
            data = loaded.SessionData;
            if string(data.Schema) ~= obj.SessionSchema || ...
                    data.Version ~= obj.SessionVersion
                error("radia:optuna:SessionStorage", ...
                    "'%s' has an unsupported session schema.", path);
            end
            obj.State = string(data.State);
            obj.SelectedTrialNumber = data.SelectedTrialNumber;
            obj.StartedAt = data.StartedAt;
            obj.UpdatedAt = data.UpdatedAt;
            obj.StopReason = string(data.StopReason);
            obj.CheckpointRevision = data.CheckpointRevision;
            if strlength(obj.ModelName) == 0 && isfield(data, "ModelName")
                obj.ModelName = string(data.ModelName);
            end
        end
    end
end
