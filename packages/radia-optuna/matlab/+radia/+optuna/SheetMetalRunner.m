classdef SheetMetalRunner < handle
    %SHEETMETALRUNNER Evaluate Optuna trials with native sheet-metal solvers.
    %   A problem factory samples the outer design variables from a Trial and
    %   returns the callback contract consumed by Radia's existing two-level
    %   NGSolve/Cubit sheet-metal optimizer. Each trial receives an isolated
    %   working directory and a compact result.json provenance record.

    properties (SetAccess=private)
        ProblemFcn function_handle
        ScoreFcn function_handle
        ResultFcn function_handle
        Solver (1,1) string
        OutputRoot (1,1) string
    end

    methods
        function obj = SheetMetalRunner(problemFcn, options)
            arguments
                problemFcn (1,1) function_handle
                options.Solver (1,1) string {mustBeMember(options.Solver, ...
                    ["hex-sheet", "hcurl-joule", "hcurl-activation"])} = ...
                    "hex-sheet"
                options.ScoreFcn = []
                options.ResultFcn = []
                options.OutputRoot (1,1) string = ...
                    "C:\temp\radia_sheet_metal_optuna"
            end
            obj.ProblemFcn = problemFcn;
            obj.Solver = options.Solver;
            obj.OutputRoot = options.OutputRoot;
            if isempty(options.ScoreFcn)
                obj.ScoreFcn = @(result, trial) result.state.objective;
            elseif isa(options.ScoreFcn, "function_handle")
                obj.ScoreFcn = options.ScoreFcn;
            else
                error("radia:optuna:SheetMetalScore", ...
                    "ScoreFcn must be a function handle.");
            end
            if isempty(options.ResultFcn)
                obj.ResultFcn = @(result, trial, runDirectory) [];
            elseif isa(options.ResultFcn, "function_handle")
                obj.ResultFcn = options.ResultFcn;
            else
                error("radia:optuna:SheetMetalResult", ...
                    "ResultFcn must be a function handle.");
            end
        end

        function [score, result] = evaluate(obj, trial)
            arguments
                obj
                trial (1,1) radia.optuna.Trial
            end
            problem = obj.ProblemFcn(trial);
            if ~isstruct(problem) || ~isscalar(problem)
                error("radia:optuna:SheetMetalProblem", ...
                    "ProblemFcn must return one problem struct.");
            end
            runDirectory = obj.createTrialDirectory(trial.Number);
            problem = obj.withTrialDirectory(problem, runDirectory);
            result = obj.solve(problem);
            score = obj.ScoreFcn(result, trial);
            if ~isnumeric(score) || isempty(score) || ...
                    ~isvector(score) || any(~isfinite(score), "all")
                error("radia:optuna:SheetMetalScore", ...
                    "ScoreFcn must return finite scalar or vector objectives.");
            end
            score = reshape(double(score), 1, []);

            trial.setUserAttr("sheet_metal_solver", obj.Solver);
            trial.setUserAttr("run_dir", runDirectory);
            obj.invokeResultFcn(result, trial, runDirectory);
            resultPath = fullfile(runDirectory, "result.json");
            trial.setUserAttr("result_json", resultPath);
            obj.writeSummary(resultPath, trial, result, score, runDirectory);
        end

        function results = optimize(obj, study, nTrials)
            arguments
                obj
                study (1,1) radia.optuna.Study
                nTrials (1,1) double {mustBeInteger, mustBeNonnegative}
            end
            results = study.optimize(@(trial) obj.evaluate(trial), nTrials);
        end
    end

    methods (Access=private)
        function problem = withTrialDirectory(~, problem, runDirectory)
            if isfield(problem, "driverOptions")
                driverOptions = problem.driverOptions;
                if ~isstruct(driverOptions) || ~isscalar(driverOptions)
                    error("radia:optuna:SheetMetalProblem", ...
                        "driverOptions must be one struct.");
                end
            else
                driverOptions = struct();
            end
            % Trial isolation is part of the production artifact contract.
            driverOptions.WorkDirectory = runDirectory;
            problem.driverOptions = driverOptions;
        end

        function result = solve(obj, problem)
            common = ["initialState", "deformationFactory", ...
                "evaluateObjective", "rebuildHMatrix", "cubitBackend", ...
                "elementSizes", "driverOptions"];
            switch obj.Solver
                case "hex-sheet"
                    required = [common, "linearizeStep"];
                    obj.requireFields(problem, required);
                    args = namedargs2cell(problem.driverOptions);
                    result = radia.topopt.optimizeHexSheetTopology( ...
                        problem.initialState, problem.linearizeStep, ...
                        problem.deformationFactory, problem.evaluateObjective, ...
                        problem.rebuildHMatrix, problem.cubitBackend, ...
                        problem.elementSizes, args{:});
                case "hcurl-joule"
                    required = [common, "buildStepInputs"];
                    obj.requireFields(problem, required);
                    result = radia.topopt.optimizeHCurlEddyBubbleHexSheet( ...
                        problem.initialState, problem.buildStepInputs, ...
                        problem.deformationFactory, problem.evaluateObjective, ...
                        problem.rebuildHMatrix, problem.cubitBackend, ...
                        problem.elementSizes, problem.driverOptions);
                case "hcurl-activation"
                    required = [common, "buildStepInputs"];
                    obj.requireFields(problem, required);
                    result = ...
                        radia.topopt.optimizeHCurlEddyBubbleActivationHexSheet( ...
                        problem.initialState, problem.buildStepInputs, ...
                        problem.deformationFactory, problem.evaluateObjective, ...
                        problem.rebuildHMatrix, problem.cubitBackend, ...
                        problem.elementSizes, problem.driverOptions);
            end
        end

        function requireFields(~, problem, required)
            missing = required(~isfield(problem, required));
            if ~isempty(missing)
                error("radia:optuna:SheetMetalProblem", ...
                    "Sheet-metal problem is missing: %s", ...
                    strjoin(missing, ", "));
            end
        end

        function invokeResultFcn(obj, result, trial, runDirectory)
            if nargin(obj.ResultFcn) == 2
                obj.ResultFcn(result, trial);
            else
                obj.ResultFcn(result, trial, runDirectory);
            end
        end

        function runDirectory = createTrialDirectory(obj, trialNumber)
            if ~isfolder(obj.OutputRoot)
                mkdir(obj.OutputRoot);
            end
            base = fullfile(obj.OutputRoot, sprintf("trial_%06d", trialNumber));
            runDirectory = string(base);
            suffix = 1;
            while isfolder(runDirectory)
                suffix = suffix + 1;
                runDirectory = string(base) + compose("_%02d", suffix);
            end
            mkdir(runDirectory);
        end

        function writeSummary(obj, path, trial, result, score, runDirectory)
            summary = struct( ...
                "schema", "radia.optuna.sheet-metal-trial.v1", ...
                "trial_number", trial.Number, ...
                "solver", obj.Solver, ...
                "status", "completed", ...
                "objective_values", score, ...
                "params", trial.Params, ...
                "user_attrs", trial.UserAttrs, ...
                "run_dir", runDirectory, ...
                "result_schema", obj.resultField(result, "schema", ""), ...
                "converged", logical(obj.resultField(result, ...
                    "converged", false)), ...
                "inner_iteration_count", double(obj.resultField(result, ...
                    "inner_iteration_count", 0)), ...
                "outer_iteration_count", double(obj.resultField(result, ...
                    "outer_iteration_count", 0)), ...
                "cubit_rebuild_count", double(obj.resultField(result, ...
                    "cubit_rebuild_count", 0)), ...
                "generated_at_utc", string(datetime("now", ...
                    "TimeZone", "UTC", ...
                    "Format", "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'")), ...
                "matlab_version", string(version), ...
                "platform", string(computer));
            file = fopen(path, "w", "n", "UTF-8");
            if file < 0
                error("radia:optuna:SheetMetalArtifact", ...
                    "Cannot write sheet-metal trial result: %s", path);
            end
            cleanup = onCleanup(@() fclose(file));
            count = fprintf(file, "%s\n", ...
                jsonencode(summary, PrettyPrint=true));
            if count <= 0
                error("radia:optuna:SheetMetalArtifact", ...
                    "Cannot write sheet-metal trial result: %s", path);
            end
            clear cleanup
        end

        function value = resultField(~, result, name, fallback)
            if isstruct(result) && isfield(result, name)
                value = result.(name);
            else
                value = fallback;
            end
        end
    end
end
