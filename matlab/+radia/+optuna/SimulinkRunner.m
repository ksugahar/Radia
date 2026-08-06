classdef SimulinkRunner < handle
    %SIMULINKRUNNER Run CAE-aware optimization trials inside MATLAB.

    properties (SetAccess=private)
        Model (1,1) string
        ConfigureFcn function_handle
        ScoreFcn function_handle
        ConstraintFcn function_handle
        ValidationFcn function_handle
        ResultFcn function_handle
        FailureClassifierFcn = []
        StopTime (1,1) string = ""
        UseFastRestart (1,1) logical = true
        Context (1,1) struct = struct()
        ModelIdentity (1,1) struct = struct()
    end

    methods
        function obj = SimulinkRunner(model, options)
            arguments
                model (1,1) string
                options.ConfigureFcn = []
                options.ScoreFcn = []
                options.ConstraintFcn = []
                options.ValidationFcn = []
                options.ResultFcn = []
                options.FailureClassifierFcn = []
                options.StopTime (1,1) string = ""
                options.UseFastRestart (1,1) logical = true
                options.Context (1,1) struct = struct()
            end
            if isempty(options.ScoreFcn) || ~isa(options.ScoreFcn, "function_handle")
                error("radia:optuna:SimulinkScore", ...
                    "SimulinkRunner requires a ScoreFcn(simOut, trial) function handle.");
            end
            obj.Model = model;
            obj.ConfigureFcn = obj.optionalFcn(options.ConfigureFcn, ...
                @(simInput, trial) simInput, "ConfigureFcn");
            obj.ScoreFcn = options.ScoreFcn;
            obj.ConstraintFcn = obj.optionalFcn(options.ConstraintFcn, ...
                @(simOut, trial) zeros(1,0), "ConstraintFcn");
            obj.ValidationFcn = obj.optionalFcn(options.ValidationFcn, ...
                @(simOut, trial) struct(), "ValidationFcn");
            obj.ResultFcn = obj.optionalFcn(options.ResultFcn, ...
                @(simOut, trial) struct(), "ResultFcn");
            if ~isempty(options.FailureClassifierFcn) && ...
                    ~isa(options.FailureClassifierFcn, "function_handle")
                error("radia:optuna:FailureClassifier", ...
                    "FailureClassifierFcn must be a function handle.");
            end
            obj.FailureClassifierFcn = options.FailureClassifierFcn;
            obj.StopTime = options.StopTime;
            obj.UseFastRestart = options.UseFastRestart;
            obj.Context = options.Context;
            obj.requireJSON(obj.Context, "Context");
            obj.ModelIdentity = obj.resolveModelIdentity();
        end

        function [score, simOut, record] = evaluate(obj, trial)
            arguments
                obj
                trial (1,1) radia.optuna.Trial
            end
            timing = obj.emptyTiming();
            totalTimer = tic;
            stage = "configuration";
            try
                setupTimer = tic;
                simInput = obj.prepareInput(trial);
                timing.configuration = toc(setupTimer);

                stage = "simulation";
                simulationTimer = tic;
                if obj.UseFastRestart
                    simOut = sim(simInput, "UseFastRestart", "on");
                else
                    simOut = sim(simInput);
                end
                timing.simulation = toc(simulationTimer);

                stage = "postprocess";
                [score, constraints, validation, artifacts, postprocess_s] = ...
                    obj.extractResult(simOut, trial);
                timing.postprocess = postprocess_s;
                timing.total = toc(totalTimer);
                record = obj.completeRecord(trial, score, constraints, ...
                    validation, artifacts, timing, "serial", 1);
                obj.recordCompleteTrial(trial, record, constraints);
            catch exception
                timing.total = toc(totalTimer);
                record = obj.failedRecord(trial, exception, stage, timing, ...
                    "serial", 1);
                obj.safeSetUserAttr(trial, "cae_execution", record);
                obj.safeSetUserAttr(trial, "cae_failure", record.failure);
                rethrow(exception);
            end
        end

        function results = optimize(obj, study, nTrials, options)
            arguments
                obj
                study (1,1) radia.optuna.Study
                nTrials (1,1) double {mustBeInteger, mustBeNonnegative}
                options.ContinueOnError (1,1) logical = true
            end
            obj.registerStudy(study, "serial", 1);
            for index = 1:nTrials
                trial = study.ask();
                try
                    score = obj.evaluate(trial);
                    if trial.State == "RUNNING"
                        study.tell(trial, score);
                        obj.restoreRunnerConstraints(study, trial);
                    end
                catch exception
                    if trial.State == "RUNNING"
                        study.fail(trial, string(exception.message));
                    end
                    if ~options.ContinueOnError
                        rethrow(exception);
                    end
                end
            end
            results = study.TrialTable;
        end

        function results = optimizeParallel(obj, study, nTrials, options)
            arguments
                obj
                study (1,1) radia.optuna.Study
                nTrials (1,1) double {mustBeInteger,mustBeNonnegative}
                options.BatchSize (1,1) double {mustBeInteger,mustBePositive} = 4
                options.ShowProgress (1,1) logical = true
                options.TransferBaseWorkspaceVariables (1,1) logical = true
                options.ContinueOnError (1,1) logical = true
            end
            if isempty(ver("parallel"))
                error("radia:optuna:ParallelUnavailable", ...
                    "Parallel Computing Toolbox is required for Simulink parallel trials.");
            end
            batchSize = min(options.BatchSize, max(1, nTrials));
            obj.registerStudy(study, "parsim", batchSize);
            modelFolder = obj.modelFolder();
            if strlength(modelFolder) > 0
                addpath(modelFolder);
            end
            firstFailure = [];
            completed = 0;
            while completed < nTrials
                count = min(batchSize, nTrials - completed);
                trials = cell(count, 1);
                inputs = repmat(Simulink.SimulationInput(obj.Model), 1, count);
                valid = false(count, 1);
                configuration_s = NaN(count, 1);
                for localIndex = 1:count
                    trial = study.ask();
                    trials{localIndex} = trial;
                    timer = tic;
                    try
                        inputs(localIndex) = obj.prepareInput(trial);
                        configuration_s(localIndex) = toc(timer);
                        valid(localIndex) = true;
                    catch exception
                        configuration_s(localIndex) = toc(timer);
                        timing = obj.emptyTiming();
                        timing.configuration = configuration_s(localIndex);
                        timing.total = configuration_s(localIndex);
                        record = obj.failedRecord(trial, exception, ...
                            "configuration", timing, "parsim", count);
                        obj.safeSetUserAttr(trial, "cae_execution", record);
                        obj.safeSetUserAttr(trial, "cae_failure", record.failure);
                        study.fail(trial, string(exception.message));
                        if isempty(firstFailure)
                            firstFailure = exception;
                        end
                    end
                end

                runIndices = find(valid);
                if ~isempty(runIndices)
                    batchTimer = tic;
                    try
                        outputs = parsim(inputs(runIndices), ...
                            "ShowProgress", options.ShowProgress, ...
                            "UseFastRestart", obj.UseFastRestart, ...
                            "StopOnError", "off", ...
                            "SetupFcn", @()radia.optuna.internal. ...
                                setupSimulinkWorker(modelFolder), ...
                            "TransferBaseWorkspaceVariables", ...
                                options.TransferBaseWorkspaceVariables);
                        batchElapsed = toc(batchTimer);
                        for outputIndex = 1:numel(runIndices)
                            localIndex = runIndices(outputIndex);
                            trial = trials{localIndex};
                            timing = obj.emptyTiming();
                            timing.configuration = configuration_s(localIndex);
                            timing.simulation = obj.simulationSeconds( ...
                                outputs(outputIndex), batchElapsed / numel(runIndices));
                            try
                                if strlength(string(outputs(outputIndex).ErrorMessage)) > 0
                                    error("radia:optuna:SimulinkTrial", "%s", ...
                                        outputs(outputIndex).ErrorMessage);
                                end
                                [score, constraints, validation, artifacts, ...
                                    timing.postprocess] = obj.extractResult( ...
                                        outputs(outputIndex), trial);
                                timing.total = timing.configuration + ...
                                    timing.simulation + timing.postprocess;
                                record = obj.completeRecord(trial, score, ...
                                    constraints, validation, artifacts, timing, ...
                                    "parsim", count);
                                obj.recordCompleteTrial(trial, record, constraints);
                                study.tell(trial, score);
                                obj.restoreRunnerConstraints(study, trial);
                            catch exception
                                timing.total = sum( ...
                                    [timing.configuration, timing.simulation, ...
                                    timing.postprocess], "omitnan");
                                record = obj.failedRecord(trial, exception, ...
                                    "postprocess", timing, "parsim", count);
                                obj.safeSetUserAttr(trial, "cae_execution", record);
                                obj.safeSetUserAttr(trial, "cae_failure", record.failure);
                                if trial.State == "RUNNING"
                                    study.fail(trial, string(exception.message));
                                end
                                if isempty(firstFailure)
                                    firstFailure = exception;
                                end
                            end
                        end
                    catch exception
                        batchElapsed = toc(batchTimer);
                        for localIndex = reshape(runIndices, 1, [])
                            trial = trials{localIndex};
                            if trial.State ~= "RUNNING"
                                continue
                            end
                            timing = obj.emptyTiming();
                            timing.configuration = configuration_s(localIndex);
                            timing.simulation = batchElapsed;
                            timing.total = timing.configuration + batchElapsed;
                            record = obj.failedRecord(trial, exception, ...
                                "simulation", timing, "parsim", count);
                            obj.safeSetUserAttr(trial, "cae_execution", record);
                            obj.safeSetUserAttr(trial, "cae_failure", record.failure);
                            study.fail(trial, string(exception.message));
                        end
                        if isempty(firstFailure)
                            firstFailure = exception;
                        end
                    end
                end
                completed = completed + count;
                if ~options.ContinueOnError && ~isempty(firstFailure)
                    rethrow(firstFailure);
                end
            end
            results = study.TrialTable;
        end
    end

    methods (Access=private)
        function value = optionalFcn(~, candidate, fallback, name)
            if isempty(candidate)
                value = fallback;
            elseif isa(candidate, "function_handle")
                value = candidate;
            else
                error("radia:optuna:SimulinkCallback", ...
                    "%s must be a function handle.", name);
            end
        end

        function simInput = prepareInput(obj, trial)
            simInput = Simulink.SimulationInput(obj.Model);
            if strlength(obj.StopTime) > 0
                simInput = simInput.setModelParameter( ...
                    "StopTime", char(obj.StopTime));
            end
            simInput = obj.ConfigureFcn(simInput, trial);
            if ~isa(simInput, "Simulink.SimulationInput")
                error("radia:optuna:SimulinkInput", ...
                    "ConfigureFcn must return a Simulink.SimulationInput object.");
            end
        end

        function [score, constraints, validation, artifacts, elapsed] = ...
                extractResult(obj, simOut, trial)
            timer = tic;
            score = reshape(double(obj.ScoreFcn(simOut, trial)), 1, []);
            if isempty(score) || any(~isfinite(score))
                error("radia:optuna:SimulinkScore", ...
                    "ScoreFcn must return finite numeric objective values.");
            end
            constraints = reshape(double(obj.ConstraintFcn(simOut, trial)), 1, []);
            if any(isnan(constraints))
                error("radia:optuna:SimulinkConstraints", ...
                    "ConstraintFcn must return numeric c <= 0 values without NaN.");
            end
            validation = obj.invokeOutputFcn(obj.ValidationFcn, simOut, trial, ...
                "ValidationFcn");
            if isempty(validation)
                validation = struct();
            end
            if ~isstruct(validation) || ~isscalar(validation)
                error("radia:optuna:SimulinkValidation", ...
                    "ValidationFcn must return one JSON-serializable struct.");
            end
            artifacts = obj.invokeOutputFcn(obj.ResultFcn, simOut, trial, ...
                "ResultFcn");
            if isempty(artifacts)
                artifacts = struct();
            end
            if ~isstruct(artifacts) || ~isscalar(artifacts)
                error("radia:optuna:SimulinkArtifacts", ...
                    "ResultFcn must return one JSON-serializable artifact manifest struct or no output.");
            end
            obj.requireJSON(validation, "ValidationFcn result");
            obj.requireJSON(artifacts, "ResultFcn result");
            elapsed = toc(timer);
        end

        function value = invokeOutputFcn(~, callback, simOut, trial, name)
            outputCount = nargout(callback);
            if outputCount == 0
                callback(simOut, trial);
                value = struct();
            else
                try
                    value = callback(simOut, trial);
                catch exception
                    if exception.identifier == "MATLAB:maxlhs" || ...
                            contains(exception.message, "Too many output")
                        callback(simOut, trial);
                        value = struct();
                    else
                        rethrow(exception);
                    end
                end
            end
            if isempty(value)
                value = struct();
            end
            if ~(isstruct(value) && isscalar(value))
                error("radia:optuna:SimulinkCallbackResult", ...
                    "%s must return one scalar struct or no output.", name);
            end
        end

        function recordCompleteTrial(~, trial, record, constraints)
            trial.setUserAttr("cae_execution", record);
            if ~isempty(fieldnames(record.validation))
                trial.setUserAttr("cae_validation", record.validation);
            end
            if ~isempty(fieldnames(record.artifacts))
                trial.setUserAttr("cae_artifacts", record.artifacts);
            end
            trial.Study.recordConstraints(trial, constraints);
        end

        function restoreRunnerConstraints(~, study, trial)
            key = matlab.lang.makeValidName("cae_execution");
            if isfield(trial.UserAttrs, key)
                record = trial.UserAttrs.(key);
                if isfield(record, "constraint_values")
                    study.recordConstraints(trial, record.constraint_values);
                end
            end
        end

        function registerStudy(obj, study, mode, batchSize)
            contract = struct( ...
                "schema", "radia.optuna.cae-simulink-runner.v1", ...
                "model", obj.ModelIdentity, ...
                "context", obj.Context, ...
                "execution_mode", string(mode), ...
                "batch_size", double(batchSize), ...
                "use_fast_restart", obj.UseFastRestart, ...
                "constraint_convention", "c <= 0", ...
                "versions", obj.versionRecord());
            study.setUserAttr("cae_execution_contract", contract);
        end

        function record = completeRecord(obj, trial, score, constraints, ...
                validation, artifacts, timing, mode, batchSize)
            record = obj.baseRecord(trial, timing, mode, batchSize);
            record.status = "complete";
            record.objective_values = score;
            record.constraint_values = constraints;
            record.validation = validation;
            record.artifacts = artifacts;
            record.completed_at = obj.timestamp();
            obj.requireJSON(record, "CAE execution record");
        end

        function record = failedRecord(obj, trial, exception, stage, timing, ...
                mode, batchSize)
            record = obj.baseRecord(trial, timing, mode, batchSize);
            record.status = "failed";
            record.objective_values = zeros(1,0);
            record.constraint_values = zeros(1,0);
            record.validation = struct();
            record.artifacts = struct();
            record.failure = obj.classifyFailure(exception, stage, trial);
            record.completed_at = obj.timestamp();
        end

        function record = baseRecord(obj, trial, timing, mode, batchSize)
            record = struct( ...
                "schema", "radia.optuna.cae-trial.v1", ...
                "trial_number", trial.Number, ...
                "model", obj.ModelIdentity, ...
                "context", obj.Context, ...
                "parameters", trial.Params, ...
                "execution_mode", string(mode), ...
                "batch_size", double(batchSize), ...
                "use_fast_restart", obj.UseFastRestart, ...
                "constraint_convention", "c <= 0", ...
                "versions", obj.versionRecord(), ...
                "started_at", obj.datetimeText(trial.StartTime), ...
                "timing_s", timing);
        end

        function failure = classifyFailure(obj, exception, stage, trial)
            if ~isempty(obj.FailureClassifierFcn)
                failure = obj.FailureClassifierFcn(exception, string(stage), trial);
                if ~isstruct(failure) || ~isscalar(failure)
                    error("radia:optuna:FailureClassifier", ...
                        "FailureClassifierFcn must return one scalar struct.");
                end
            else
                failure = obj.defaultFailure(exception, stage);
            end
            if ~isfield(failure, "schema")
                failure.schema = "radia.optuna.cae-failure.v1";
            end
            if ~isfield(failure, "stage")
                failure.stage = string(stage);
            end
            if ~isfield(failure, "retryable")
                failure.retryable = false;
            end
            obj.requireJSON(failure, "failure classification");
        end

        function failure = defaultFailure(~, exception, stage)
            token = lower(string(exception.identifier) + " " + ...
                string(exception.message));
            kind = "unknown";
            retryable = false;
            action = "inspect_exception";
            if contains(token, "license") || contains(token, "checkout")
                kind = "license_or_resource";
                retryable = true;
                action = "wait_for_owned_resource_then_retry_same_parameters";
            elseif contains(token, "timeout") || contains(token, "timed out")
                kind = "timeout";
                retryable = true;
                action = "inspect_progress_then_retry_or_raise_timeout";
            elseif contains(token, "mesh") || contains(token, "element quality")
                kind = "mesh_invalid";
                action = "repair_or_refine_mesh_before_retry";
            elseif contains(token, "converg") || contains(token, "singular") || ...
                    contains(token, "solver")
                kind = "solver_nonconvergence";
                action = "adjust_continuation_mesh_or_solver_policy";
            elseif string(stage) == "configuration"
                kind = "configuration_invalid";
                action = "reject_or_repair_design_parameters";
            elseif string(stage) == "postprocess"
                kind = "observable_invalid";
                action = "repair_score_constraint_or_artifact_contract";
            elseif string(stage) == "simulation"
                kind = "simulation_failed";
                action = "inspect_simulink_diagnostics";
            end
            causeIdentifiers = strings(1, numel(exception.cause));
            for index = 1:numel(exception.cause)
                causeIdentifiers(index) = string(exception.cause{index}.identifier);
            end
            failure = struct( ...
                "schema", "radia.optuna.cae-failure.v1", ...
                "class", kind, ...
                "stage", string(stage), ...
                "retryable", retryable, ...
                "recommended_action", action, ...
                "identifier", string(exception.identifier), ...
                "message", string(exception.message), ...
                "cause_identifiers", causeIdentifiers, ...
                "timestamp", radia.optuna.SimulinkRunner.timestamp());
        end

        function identity = resolveModelIdentity(obj)
            path = "";
            if isfile(obj.Model)
                path = string(java.io.File(char(obj.Model)).getCanonicalPath());
            else
                located = which(char(obj.Model));
                if ~isempty(located)
                    path = string(located);
                end
            end
            identity = struct( ...
                "name", obj.Model, ...
                "path", path, ...
                "sha256", "", ...
                "bytes", NaN, ...
                "modified_at", "");
            if strlength(path) == 0 || ~isfile(path)
                return
            end
            info = dir(path);
            identity.bytes = info.bytes;
            modified = datetime(info.datenum, "ConvertFrom", "datenum", ...
                "TimeZone", "local");
            identity.modified_at = obj.datetimeText(modified);
            try
                identity.sha256 = obj.sha256File(path);
            catch
                identity.sha256 = "";
            end
        end

        function folder = modelFolder(obj)
            if strlength(obj.ModelIdentity.path) > 0
                folder = string(fileparts(obj.ModelIdentity.path));
            else
                folder = string(fileparts(obj.Model));
            end
        end

        function secondsValue = simulationSeconds(~, simOut, fallback)
            secondsValue = double(fallback);
            try
                timing = simOut.SimulationMetadata.TimingInfo;
                candidates = ["TotalElapsedWallTime", "ExecutionElapsedWallTime"];
                for name = candidates
                    if isprop(timing, name)
                        value = timing.(name);
                    elseif isstruct(timing) && isfield(timing, name)
                        value = timing.(name);
                    else
                        continue
                    end
                    if isduration(value)
                        value = seconds(value);
                    end
                    if isnumeric(value) && isscalar(value) && isfinite(value)
                        secondsValue = double(value);
                        return
                    end
                end
            catch
            end
        end

        function safeSetUserAttr(~, trial, name, value)
            try
                trial.setUserAttr(name, value);
            catch
            end
        end

        function requireJSON(~, value, name)
            try
                jsonencode(value);
            catch exception
                error("radia:optuna:JSONContract", ...
                    "%s must be JSON serializable: %s", name, exception.message);
            end
        end

        function versions = versionRecord(~)
            simulinkVersion = "unavailable";
            entry = ver("simulink");
            if ~isempty(entry)
                simulinkVersion = string(entry(1).Version);
            end
            versions = struct( ...
                "matlab", string(version), ...
                "simulink", simulinkVersion);
        end

        function timing = emptyTiming(~)
            timing = struct( ...
                "configuration", NaN, ...
                "simulation", NaN, ...
                "postprocess", NaN, ...
                "total", NaN);
        end

        function value = sha256File(~, path)
            digest = java.security.MessageDigest.getInstance("SHA-256");
            stream = fopen(path, "rb");
            if stream < 0
                error("radia:optuna:ModelIdentity", ...
                    "Cannot open model file for hashing: %s", path);
            end
            cleanup = onCleanup(@()fclose(stream));
            while true
                bytes = fread(stream, 1024 * 1024, "*uint8");
                if isempty(bytes)
                    break
                end
                digest.update(typecast(bytes, "int8"));
            end
            raw = typecast(digest.digest(), "uint8");
            value = lower(string(reshape(dec2hex(raw, 2).', 1, [])));
            clear cleanup
        end
    end

    methods (Static, Access=private)
        function value = timestamp()
            instant = datetime("now", "TimeZone", "UTC", ...
                "Format", "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'");
            value = string(instant);
        end

        function value = datetimeText(input)
            if ischar(input) || isstring(input)
                input = datetime(input);
            end
            input.TimeZone = "UTC";
            input.Format = "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'";
            value = string(input);
        end
    end
end
