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
        SessionVersion = 2
        OptionFields = ["MaxTrials","MaxTime","FunctionTolerance", ...
            "MaxStallTrials","Display","UseParallel","ParallelMode", ...
            "BatchSize","Sampler","Pruner","Seed","StoragePath", ...
            "StudyName","Resume","Directions","CatchObjectiveErrors"]
        ParameterFields = ["Name","Value","Minimum","Maximum","Free", ...
            "Type","Transform","Step","Choices"]
        CallbackFields = ["Slot","Position","Description","HandleType"]
        SessionFields = ["Schema","Version","Parameters","Options", ...
            "UnresolvedCallbacks","State","SelectedTrialNumber", ...
            "StartedAt","UpdatedAt","StopReason","CheckpointRevision", ...
            "ModelName"]
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
            radia.optuna.Study.requireHdf5Path(path);
            obj.CheckpointRevision = obj.CheckpointRevision + 1;
            obj.UpdatedAt = datetime("now");
            radia.optuna.OptimizationSession.writeSessionFile( ...
                path, obj.sessionData());
        end
    end

    methods (Static)
        function obj = load(path, options)
            %LOAD Rebuild a session from its HDF5 checkpoint.
            %   Function handles are never restored from storage. The
            %   caller re-binds every callback slot recorded at save:
            %   ObjectiveFcn always, and OutputFcn, PlotFcn and
            %   ConstraintFcn exactly when the checkpoint recorded them
            %   (same count, in saved order). Missing, extra or miscounted
            %   callbacks raise radia:optuna:SessionCallbackUnresolved.
            arguments
                path (1,1) string
                options.ObjectiveFcn = []
                options.OutputFcn = []
                options.PlotFcn = []
                options.ConstraintFcn = []
                options.ModelName (1,1) string = ""
            end
            data = radia.optuna.OptimizationSession.readSessionFile(path);
            [objectiveFcn, sessionOptions] = radia.optuna. ...
                OptimizationSession.rebindCallbacks(data, path, options);
            modelName = options.ModelName;
            if strlength(modelName) == 0
                modelName = data.ModelName;
            end
            obj = radia.optuna.OptimizationSession( ...
                objectiveFcn, data.Parameters, sessionOptions, ...
                ModelName=modelName);
        end

        function path = sessionPathFor(storagePath)
            %SESSIONPATHFOR Session checkpoint path next to a Study file.
            path = string(storagePath) + ".session.h5";
        end

        function data = readSessionFile(path)
            %READSESSIONFILE Read one typed HDF5 session checkpoint.
            %   Parameters and Options are rebuilt through their public
            %   constructors; UnresolvedCallbacks is descriptive only.
            import radia.optuna.OptimizationSession
            path = string(path);
            fid = radia.optuna.internal.TypedHdf5.openFile(path, ...
                OptimizationSession.SessionSchema, ...
                OptimizationSession.SessionVersion);
            closer = onCleanup(@() H5F.close(fid));
            try
                raw = radia.optuna.internal.TypedHdf5.readTree(fid, "SessionData");
            catch exception
                failure = MException("radia:optuna:SessionStorage", ...
                    "'%s' has invalid typed SessionData: %s", ...
                    path, exception.message);
                throw(addCause(failure, exception))
            end
            OptimizationSession.requireFields(raw, ...
                OptimizationSession.SessionFields, path, "SessionData");
            OptimizationSession.requireFields(raw.Options, ...
                OptimizationSession.OptionFields, path, "Options");
            OptimizationSession.requireFields(raw.Parameters, ...
                OptimizationSession.ParameterFields, path, "Parameters");
            OptimizationSession.requireFields(raw.UnresolvedCallbacks, ...
                OptimizationSession.CallbackFields, path, ...
                "UnresolvedCallbacks");
            if raw.Schema ~= OptimizationSession.SessionSchema || ...
                    raw.Version ~= OptimizationSession.SessionVersion
                error("radia:optuna:SessionStorage", ...
                    "'%s' has an unsupported session schema.", path);
            end
            data = raw;
            data.Parameters = OptimizationSession.parametersFromData( ...
                raw.Parameters);
            optionValues = namedargs2cell(raw.Options);
            data.Options = radia.optuna.OptimizeOptions(optionValues{:});
        end
    end

    methods (Static, Access=private)
        function [objectiveFcn, options] = rebindCallbacks(data, path, supplied)
            % Explicit re-binding at load: storage holds only descriptive
            % callback records; every recorded slot must be supplied by the
            % caller with the recorded count, and nothing else may be.
            callbacks = data.UnresolvedCallbacks;
            objectiveFcn = supplied.ObjectiveFcn;
            if isempty(objectiveFcn) || ~isa(objectiveFcn, "function_handle") || ...
                    ~isscalar(objectiveFcn)
                fail("ObjectiveFcn", sum(callbacks.Slot == "ObjectiveFcn"), ...
                    "a scalar function handle");
            end
            rebound = struct();
            for slot = ["OutputFcn","PlotFcn","ConstraintFcn"]
                value = supplied.(slot);
                if isempty(value)
                    handles = {};
                elseif isa(value, "function_handle") && isscalar(value)
                    handles = {value};
                elseif slot ~= "ConstraintFcn" && iscell(value) && ...
                        all(cellfun(@(h) isa(h, "function_handle") && ...
                        isscalar(h), value))
                    handles = reshape(value, 1, []);
                else
                    fail(slot, NaN, "function handles");
                end
                expected = sum(callbacks.Slot == slot);
                if numel(handles) ~= expected
                    fail(slot, expected, sprintf("%d handle(s)", numel(handles)));
                end
                if slot == "ConstraintFcn"
                    if isempty(handles)
                        rebound.(slot) = [];
                    else
                        rebound.(slot) = handles{1};
                    end
                else
                    rebound.(slot) = handles;
                end
            end
            values = namedargs2cell( ...
                radia.optuna.OptimizationSession.optionData(data.Options));
            options = radia.optuna.OptimizeOptions(values{:}, ...
                OutputFcn=rebound.OutputFcn, PlotFcn=rebound.PlotFcn, ...
                ConstraintFcn=rebound.ConstraintFcn);

            function fail(slot, expected, got)
                rows = callbacks.Slot == slot;
                recorded = strjoin(callbacks.Description(rows), ", ");
                if strlength(recorded) == 0
                    recorded = "<none>";
                end
                error("radia:optuna:SessionCallbackUnresolved", ...
                    "'%s' recorded %g %s callback(s) (%s); load received %s. " + ...
                    "Function handles are never restored from storage: " + ...
                    "re-bind exactly the recorded callbacks by name.", ...
                    path, expected, slot, recorded, got);
            end
        end

        function record = callbackRecord(objectiveFcn, options)
            % Descriptive, non-restorable callback inventory in slot order.
            slots = "ObjectiveFcn";
            positions = 1;
            handles = {objectiveFcn};
            for slot = ["OutputFcn","PlotFcn"]
                values = options.(slot);
                for index = 1:numel(values)
                    slots(end+1,1) = slot; %#ok<AGROW>
                    positions(end+1,1) = index; %#ok<AGROW>
                    handles{end+1,1} = values{index}; %#ok<AGROW>
                end
            end
            if ~isempty(options.ConstraintFcn)
                slots(end+1,1) = "ConstraintFcn";
                positions(end+1,1) = 1;
                handles{end+1,1} = options.ConstraintFcn;
            end
            descriptions = strings(numel(handles),1);
            types = strings(numel(handles),1);
            for index = 1:numel(handles)
                descriptions(index) = string(func2str(handles{index}));
                details = functions(handles{index});
                types(index) = string(details.type);
            end
            record = struct("Slot", slots(:), "Position", positions(:), ...
                "Description", descriptions, "HandleType", types);
        end

        function writeSessionFile(path, data)
            folder = string(fileparts(path));
            if strlength(folder) == 0
                folder = string(pwd);
            end
            temporaryPath = string(tempname(folder)) + ".h5";
            cleanup = onCleanup(@() ...
                radia.optuna.OptimizationSession.removeFile(temporaryPath));
            fid = radia.optuna.internal.TypedHdf5.createFile(temporaryPath, ...
                radia.optuna.OptimizationSession.SessionSchema, ...
                radia.optuna.OptimizationSession.SessionVersion);
            try
                radia.optuna.internal.TypedHdf5.writeTree(fid, "SessionData", data);
                radia.optuna.internal.TypedHdf5.writeAttribute(fid, ...
                    "SessionData", "index_base:SelectedTrialNumber", 0);
                radia.optuna.internal.TypedHdf5.writeAttribute(fid, ...
                    "SessionData", "restorable:UnresolvedCallbacks", 0);
            catch exception
                H5F.close(fid);
                rethrow(exception)
            end
            H5F.close(fid);
            radia.optuna.OptimizationSession.readSessionFile(temporaryPath);
            [moved, message, messageId] = movefile(temporaryPath, path, "f");
            if ~moved
                error("radia:optuna:SessionStorage", ...
                    "Unable to replace session checkpoint (%s): %s", ...
                    messageId, message);
            end
            clear cleanup
        end

        function removeFile(path)
            if isfile(path)
                delete(path);
            end
        end

        function data = parameterData(parameters)
            fields = radia.optuna.OptimizationSession.ParameterFields;
            values = cell(numel(fields), numel(parameters));
            for index = 1:numel(parameters)
                for field = 1:numel(fields)
                    values{field, index} = parameters(index).(fields(field));
                end
            end
            data = reshape(cell2struct(values, cellstr(fields), 1), 1, []);
        end

        function parameters = parametersFromData(data)
            parameters = radia.optuna.OptimizationParameter.empty(1,0);
            for index = 1:numel(data)
                item = data(index);
                parameters(index) = radia.optuna.OptimizationParameter( ...
                    item.Name, Value=item.Value, Minimum=item.Minimum, ...
                    Maximum=item.Maximum, Free=item.Free, Type=item.Type, ...
                    Transform=item.Transform, Step=item.Step, ...
                    Choices=item.Choices);
            end
        end

        function data = optionData(options)
            fields = radia.optuna.OptimizationSession.OptionFields;
            values = cell(numel(fields), 1);
            for index = 1:numel(fields)
                values{index} = options.(fields(index));
            end
            data = cell2struct(values, cellstr(fields), 1);
        end

        function requireFields(value, expected, path, label)
            if ~isstruct(value) || ...
                    ~isequal(string(fieldnames(value)).', expected)
                error("radia:optuna:SessionStorage", ...
                    "'%s' %s fields do not match session schema version %d.", ...
                    path, label, ...
                    radia.optuna.OptimizationSession.SessionVersion);
            end
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
                path = radia.optuna.OptimizationSession.sessionPathFor( ...
                    obj.Options.StoragePath);
            end
        end

        function data = sessionData(obj)
            data = struct( ...
                "Schema", obj.SessionSchema, ...
                "Version", obj.SessionVersion, ...
                "Parameters", obj.parameterData(obj.Parameters), ...
                "Options", obj.optionData(obj.Options), ...
                "UnresolvedCallbacks", ...
                    obj.callbackRecord(obj.ObjectiveFcn, obj.Options), ...
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
            data = radia.optuna.OptimizationSession.readSessionFile(path);
            obj.State = data.State;
            obj.SelectedTrialNumber = data.SelectedTrialNumber;
            obj.StartedAt = data.StartedAt;
            obj.UpdatedAt = data.UpdatedAt;
            obj.StopReason = data.StopReason;
            obj.CheckpointRevision = data.CheckpointRevision;
            if strlength(obj.ModelName) == 0
                obj.ModelName = data.ModelName;
            end
        end
    end
end
