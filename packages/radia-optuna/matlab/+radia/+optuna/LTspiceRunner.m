classdef LTspiceRunner < handle
    %LTSPICERUNNER Evaluate Optuna-style trials through LTspice.

    properties (SetAccess=private)
        NetlistFile (1,1) string
        ConfigureFcn function_handle
        ScoreFcn function_handle
        ResultFcn function_handle
        Executable (1,1) string = ""
        OutputRoot (1,1) string = "C:\temp\radia_ltspice_optuna"
    end

    methods
        function obj = LTspiceRunner(netlistFile, options)
            arguments
                netlistFile (1,1) string {mustBeFile}
                options.ConfigureFcn = []
                options.ScoreFcn = []
                options.ResultFcn = []
                options.Executable (1,1) string = ""
                options.OutputRoot (1,1) string = "C:\temp\radia_ltspice_optuna"
            end
            if isempty(options.ScoreFcn) || ~isa(options.ScoreFcn, "function_handle")
                error("radia:optuna:LTspiceScore", ...
                    "LTspiceRunner requires ScoreFcn(result, trial).");
            end
            obj.NetlistFile = netlistFile;
            if isempty(options.ConfigureFcn)
                obj.ConfigureFcn = @(trial) struct();
            else
                obj.ConfigureFcn = options.ConfigureFcn;
            end
            obj.ScoreFcn = options.ScoreFcn;
            if isempty(options.ResultFcn)
                obj.ResultFcn = @(result, trial) [];
            else
                obj.ResultFcn = options.ResultFcn;
            end
            obj.Executable = options.Executable;
            obj.OutputRoot = options.OutputRoot;
        end

        function [score, result] = evaluate(obj, trial)
            arguments
                obj
                trial (1,1) radia.optuna.Trial
            end
            parameters = obj.ConfigureFcn(trial);
            if ~isstruct(parameters) || ~isscalar(parameters)
                error("radia:optuna:LTspiceParameters", ...
                    "ConfigureFcn must return one parameter struct.");
            end
            runDirectory = fullfile(obj.OutputRoot, ...
                sprintf("trial_%06d", trial.Number));
            result = radia.ltspice.run(obj.NetlistFile, ...
                Parameters=parameters, Executable=obj.Executable, ...
                OutputDirectory=runDirectory);
            score = obj.ScoreFcn(result, trial);
            if ~isnumeric(score) || isempty(score) || any(~isfinite(score),"all")
                error("radia:optuna:LTspiceScore", ...
                    "ScoreFcn must return finite numeric objective values.");
            end
            obj.ResultFcn(result, trial);
        end

        function results = optimize(obj, study, nTrials)
            arguments
                obj
                study (1,1) radia.optuna.Study
                nTrials (1,1) double {mustBeInteger, mustBeNonnegative}
            end
            results = study.optimize(@(trial) obj.evaluate(trial), nTrials);
        end

        function results = optimizeParallel(obj, study, nTrials, options)
            arguments
                obj
                study (1,1) radia.optuna.Study
                nTrials (1,1) double {mustBeInteger,mustBeNonnegative}
                options.Pool = []
                options.ShowProgress (1,1) logical = true
            end
            if isempty(ver("parallel"))
                error("radia:optuna:ParallelUnavailable", ...
                    "Parallel Computing Toolbox is required for LTspice parallel trials.");
            end
            pool = options.Pool;
            if isempty(pool), pool = gcp(); end
            trials = cell(nTrials,1); futures = parallel.FevalFuture.empty(0,nTrials);
            for k = 1:nTrials
                trial = study.ask(); trials{k} = trial;
                parameters = obj.ConfigureFcn(trial);
                if ~isstruct(parameters) || ~isscalar(parameters)
                    error("radia:optuna:LTspiceParameters", ...
                        "ConfigureFcn must return one parameter struct.");
                end
                runDirectory = fullfile(obj.OutputRoot, ...
                    sprintf("trial_%06d", trial.Number));
                futures(k) = parfeval(pool, ...
                    @radia.optuna.internal.runLTspiceTrial, 1, ...
                    obj.NetlistFile, parameters, obj.Executable, runDirectory);
            end
            cleanup = onCleanup(@() cancel(futures));
            for completed = 1:nTrials
                [index,result] = fetchNext(futures);
                trial = trials{index};
                try
                    score = obj.ScoreFcn(result,trial);
                    study.tell(trial,score);
                    obj.ResultFcn(result,trial);
                catch exception
                    if trial.State == "RUNNING", study.fail(trial,exception.message); end
                    rethrow(exception);
                end
                if options.ShowProgress
                    fprintf("LTspice parallel trials: %d/%d complete\n",completed,nTrials);
                end
            end
            clear cleanup
            results = study.TrialTable;
        end

    end
end
