classdef OptunaRunner < handle
    %OPTUNARUNNER Optimize explicit Stream Function application trials.
    %   The outer loop is MATLAB Optuna. Each trial writes an isolated
    %   DesignSpec configuration and launches one complete Stream Function
    %   batch analysis. The inner inverse remains the production C++
    %   ACA+/QR/TSVD kernel; this runner does not replace that factorization.

    properties (SetAccess=private)
        BaseSettings (1,1) struct
        SearchSpace (1,1) struct
        ObjectiveKeys (1,:) string
        ConfigureFcn
        ScoreFcn
        ConstraintFcn
        ResultFcn
        OutputRoot (1,1) string
        WorkingDirectory (1,1) string
        Timeout_s (1,1) double
        PythonExecutable (1,1) string
        RunFcn function_handle
    end

    methods
        function obj = OptunaRunner(baseSettings, options)
            arguments
                baseSettings (1,1) struct
                options.SearchSpace (1,1) struct = struct()
                options.ObjectiveKeys string = "rms"
                options.ConfigureFcn = []
                options.ScoreFcn = []
                options.ConstraintFcn = []
                options.ResultFcn = []
                options.OutputRoot (1,1) string = ...
                    "C:\temp\radia_streamfunction_optuna"
                options.WorkingDirectory (1,1) string = ""
                options.Timeout_s (1,1) double {mustBeNonnegative} = 3600
                options.PythonExecutable (1,1) string = "python"
                options.RunFcn = []
            end
            obj.BaseSettings = baseSettings;
            obj.SearchSpace = options.SearchSpace;
            obj.ObjectiveKeys = reshape(string(options.ObjectiveKeys), 1, []);
            obj.ConfigureFcn = obj.optionalFcn( ...
                options.ConfigureFcn, "ConfigureFcn");
            obj.ScoreFcn = obj.optionalFcn(options.ScoreFcn, "ScoreFcn");
            obj.ConstraintFcn = obj.optionalFcn( ...
                options.ConstraintFcn, "ConstraintFcn");
            obj.ResultFcn = obj.optionalFcn(options.ResultFcn, "ResultFcn");
            obj.OutputRoot = options.OutputRoot;
            obj.WorkingDirectory = options.WorkingDirectory;
            obj.Timeout_s = options.Timeout_s;
            obj.PythonExecutable = options.PythonExecutable;
            if isempty(options.RunFcn)
                obj.RunFcn = @(configPath, runRoot, timeout_s, python) ...
                    radia.simulink.runApplication( ...
                        "streamfunction", configPath, RunRoot=runRoot, ...
                        Timeout_s=timeout_s, PythonExecutable=python, ...
                        ThrowOnFailure=false);
            elseif isa(options.RunFcn, "function_handle")
                obj.RunFcn = options.RunFcn;
            else
                error("radia:stream:OptunaRunFcn", ...
                    "RunFcn must be a function handle.");
            end
            if isempty(obj.ScoreFcn) && isempty(obj.ObjectiveKeys)
                error("radia:stream:OptunaObjective", ...
                    "ObjectiveKeys or ScoreFcn must define an objective.");
            end
            obj.validateSearchSpace();
            obj.requireJSON(obj.BaseSettings, "BaseSettings");
        end

        function [score, result] = evaluate(obj, trial)
            arguments
                obj
                trial (1,1) radia.optuna.Trial
            end
            settings = obj.applySearchSpace(trial, obj.BaseSettings);
            if ~isempty(obj.ConfigureFcn)
                settings = obj.ConfigureFcn(trial, settings);
            end
            if ~(isstruct(settings) && isscalar(settings))
                error("radia:stream:OptunaSettings", ...
                    "ConfigureFcn must return one settings struct.");
            end
            obj.requireJSON(settings, "trial settings");

            trialDirectory = obj.createTrialDirectory(trial.Number);
            configPath = fullfile(trialDirectory, "config.json");
            runRoot = fullfile(trialDirectory, "runs");
            primaryKey = "rms";
            if ~isempty(obj.ObjectiveKeys)
                primaryKey = obj.ObjectiveKeys(1);
            end
            radia.simulink.writeApplicationConfig( ...
                "streamfunction", settings, configPath, ...
                PrimaryKey=primaryKey, ...
                WorkingDirectory=obj.WorkingDirectory);

            applicationResult = obj.RunFcn( ...
                configPath, runRoot, obj.Timeout_s, obj.PythonExecutable);
            obj.validateApplicationResult(applicationResult);
            solverResultPath = string(applicationResult.solver_result);
            solverResult = jsondecode(fileread(solverResultPath));

            if isempty(obj.ScoreFcn)
                score = obj.valuesForKeys(solverResult, obj.ObjectiveKeys, ...
                    "objective");
            else
                score = obj.ScoreFcn(solverResult, applicationResult, trial);
            end
            score = obj.validateNumericVector(score, "ScoreFcn", true);

            constraints = zeros(1, 0);
            if ~isempty(obj.ConstraintFcn)
                constraints = obj.ConstraintFcn( ...
                    solverResult, applicationResult, trial);
                constraints = obj.validateNumericVector( ...
                    constraints, "ConstraintFcn", false);
                trial.Study.recordConstraints(trial, constraints);
            end

            result = struct( ...
                "schema", "radia.stream.optuna-trial.v1", ...
                "trial_number", trial.Number, ...
                "settings", settings, ...
                "objective_values", score, ...
                "constraint_values", constraints, ...
                "trial_dir", trialDirectory, ...
                "config_json", string(configPath), ...
                "application_result_json", ...
                    string(applicationResult.result_json), ...
                "solver_result_json", solverResultPath, ...
                "application_run_dir", string(applicationResult.run_dir), ...
                "application_elapsed_s", ...
                    double(applicationResult.elapsed_s));
            summaryPath = fullfile(trialDirectory, "optuna_trial.json");
            obj.writeJSON(summaryPath, result);
            trial.setUserAttr("streamfunction_trial", result);
            trial.setUserAttr("result_json", string(summaryPath));
            if ~isempty(obj.ResultFcn)
                obj.ResultFcn(result, solverResult, applicationResult, trial);
            end
        end

        function results = optimize(obj, study, nTrials, options)
            arguments
                obj
                study (1,1) radia.optuna.Study
                nTrials (1,1) double {mustBeInteger, mustBeNonnegative}
                options.ContinueOnError (1,1) logical = true
            end
            study.setUserAttr("streamfunction_runner", struct( ...
                "schema", "radia.stream.optuna-runner.v1", ...
                "output_root", obj.OutputRoot, ...
                "objective_keys", obj.ObjectiveKeys, ...
                "execution", "one-explicit-batch-analysis-per-trial", ...
                "inner_factorization", "aca_plus_qr_tsvd", ...
                "python_per_simulink_step", false));
            caught = strings(1, 0);
            if options.ContinueOnError
                caught = "*";
            end
            results = study.optimize(@(trial) obj.evaluate(trial), nTrials, ...
                Catch=caught);
        end
    end

    methods (Access=private)
        function settings = applySearchSpace(obj, trial, settings)
            names = string(fieldnames(obj.SearchSpace));
            for index = 1:numel(names)
                name = names(index);
                spec = obj.SearchSpace.(name);
                kind = lower(string(spec.Kind));
                switch kind
                    case "float"
                        value = trial.suggest_float(name, ...
                            double(spec.Low), double(spec.High), ...
                            Log=obj.logicalField(spec, "Log", false), ...
                            Step=obj.numericField(spec, "Step", NaN));
                    case "int"
                        value = trial.suggest_int(name, ...
                            double(spec.Low), double(spec.High), ...
                            Log=obj.logicalField(spec, "Log", false), ...
                            Step=obj.numericField(spec, "Step", 1));
                    case "categorical"
                        value = trial.suggest_categorical(name, spec.Choices);
                    otherwise
                        error("radia:stream:OptunaSearchSpace", ...
                            "SearchSpace.%s.Kind must be float, int, or categorical.", ...
                            name);
                end
                settings.(name) = value;
            end
        end

        function validateSearchSpace(obj)
            names = string(fieldnames(obj.SearchSpace));
            for index = 1:numel(names)
                name = names(index);
                spec = obj.SearchSpace.(name);
                if ~(isstruct(spec) && isscalar(spec) && isfield(spec, "Kind"))
                    error("radia:stream:OptunaSearchSpace", ...
                        "SearchSpace.%s must be one struct with Kind.", name);
                end
                kind = lower(string(spec.Kind));
                if ismember(kind, ["float", "int"])
                    if ~all(isfield(spec, ["Low", "High"])) || ...
                            ~isscalar(spec.Low) || ~isscalar(spec.High) || ...
                            ~isfinite(double(spec.Low)) || ...
                            ~isfinite(double(spec.High)) || ...
                            double(spec.Low) >= double(spec.High)
                        error("radia:stream:OptunaSearchSpace", ...
                            "SearchSpace.%s requires finite Low < High.", name);
                    end
                elseif kind == "categorical"
                    if ~isfield(spec, "Choices") || isempty(spec.Choices)
                        error("radia:stream:OptunaSearchSpace", ...
                            "SearchSpace.%s requires nonempty Choices.", name);
                    end
                else
                    error("radia:stream:OptunaSearchSpace", ...
                        "SearchSpace.%s has unknown Kind '%s'.", name, kind);
                end
            end
        end

        function validateApplicationResult(~, result)
            if ~(isstruct(result) && isscalar(result)) || ...
                    ~isfield(result, "status") || ...
                    string(result.status) ~= "passed"
                message = "Stream Function application trial failed.";
                if isstruct(result) && isfield(result, "error")
                    message = message + " " + string(result.error);
                end
                error("radia:stream:OptunaApplication", "%s", message);
            end
            required = ["solver_result", "result_json", "run_dir", "elapsed_s"];
            missing = required(~isfield(result, required));
            if ~isempty(missing) || ~isfile(string(result.solver_result))
                error("radia:stream:OptunaApplication", ...
                    "Application result is missing its solver artifact: %s", ...
                    strjoin(missing, ", "));
            end
        end

        function values = valuesForKeys(obj, payload, keys, role)
            values = zeros(1, numel(keys));
            for index = 1:numel(keys)
                value = obj.valueForKey(payload, keys(index));
                if ~(isnumeric(value) && isscalar(value) && isfinite(value))
                    error("radia:stream:OptunaObjective", ...
                        "%s key '%s' is not a finite numeric scalar.", ...
                        role, keys(index));
                end
                values(index) = double(value);
            end
        end

        function value = valueForKey(~, payload, key)
            value = payload;
            parts = split(string(key), ".");
            for part = reshape(parts, 1, [])
                if ~(isstruct(value) && isscalar(value) && isfield(value, part))
                    value = [];
                    return
                end
                value = value.(part);
            end
        end

        function value = logicalField(~, spec, name, fallback)
            if isfield(spec, name)
                value = logical(spec.(name));
            else
                value = fallback;
            end
            if ~isscalar(value)
                error("radia:stream:OptunaSearchSpace", ...
                    "%s must be one logical value.", name);
            end
        end

        function value = numericField(~, spec, name, fallback)
            if isfield(spec, name)
                value = double(spec.(name));
            else
                value = fallback;
            end
            if ~isscalar(value) || isnan(value) && ~isnan(fallback)
                error("radia:stream:OptunaSearchSpace", ...
                    "%s must be one numeric value.", name);
            end
        end

        function value = optionalFcn(~, candidate, name)
            if isempty(candidate)
                value = [];
            elseif isa(candidate, "function_handle")
                value = candidate;
            else
                error("radia:stream:OptunaCallback", ...
                    "%s must be a function handle.", name);
            end
        end

        function values = validateNumericVector(~, value, name, requireFinite)
            if ~isnumeric(value) || isempty(value) || ~isvector(value)
                error("radia:stream:OptunaCallback", ...
                    "%s must return a nonempty numeric vector.", name);
            end
            values = reshape(double(value), 1, []);
            if any(isnan(values)) || requireFinite && any(~isfinite(values))
                error("radia:stream:OptunaCallback", ...
                    "%s returned an invalid numeric value.", name);
            end
        end

        function directory = createTrialDirectory(obj, trialNumber)
            if ~isfolder(obj.OutputRoot)
                mkdir(obj.OutputRoot);
            end
            base = fullfile(obj.OutputRoot, ...
                compose("trial_%06d", trialNumber));
            directory = string(base);
            suffix = 1;
            while isfolder(directory)
                suffix = suffix + 1;
                directory = string(base) + compose("_%02d", suffix);
            end
            mkdir(directory);
        end

        function writeJSON(~, path, value)
            file = fopen(path, "w", "n", "UTF-8");
            if file < 0
                error("radia:stream:OptunaArtifact", ...
                    "Cannot write trial artifact: %s", path);
            end
            cleanup = onCleanup(@() fclose(file));
            count = fprintf(file, "%s\n", ...
                jsonencode(value, PrettyPrint=true));
            if count <= 0
                error("radia:stream:OptunaArtifact", ...
                    "Cannot write trial artifact: %s", path);
            end
            clear cleanup
        end

        function requireJSON(~, value, name)
            try
                jsonencode(value);
            catch exception
                error("radia:stream:OptunaJSON", ...
                    "%s must be JSON serializable: %s", ...
                    name, exception.message);
            end
        end
    end
end
