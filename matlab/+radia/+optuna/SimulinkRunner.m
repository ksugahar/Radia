classdef SimulinkRunner < handle
    %SIMULINKRUNNER Evaluate one Optuna-style Trial through Simulink.

    properties (SetAccess=private)
        Model (1,1) string
        ConfigureFcn function_handle
        ScoreFcn function_handle
        ResultFcn function_handle
        StopTime (1,1) string = ""
        UseFastRestart (1,1) logical = false
    end

    methods
        function obj = SimulinkRunner(model, options)
            arguments
                model (1,1) string
                options.ConfigureFcn = []
                options.ScoreFcn = []
                options.ResultFcn = []
                options.StopTime (1,1) string = ""
                options.UseFastRestart (1,1) logical = false
            end
            if isempty(options.ScoreFcn) || ~isa(options.ScoreFcn, "function_handle")
                error("radia:optuna:SimulinkScore", ...
                    "SimulinkRunner requires a ScoreFcn(simOut, trial) function handle.");
            end
            obj.Model = model;
            if isempty(options.ConfigureFcn)
                obj.ConfigureFcn = @(simInput, trial) simInput;
            else
                obj.ConfigureFcn = options.ConfigureFcn;
            end
            obj.ScoreFcn = options.ScoreFcn;
            if isempty(options.ResultFcn)
                obj.ResultFcn = @(simOut, trial) [];
            else
                obj.ResultFcn = options.ResultFcn;
            end
            obj.StopTime = options.StopTime;
            obj.UseFastRestart = options.UseFastRestart;
        end

        function [score, simOut] = evaluate(obj, trial)
            arguments
                obj
                trial (1,1) radia.optuna.Trial
            end
            simInput = Simulink.SimulationInput(obj.Model);
            if strlength(obj.StopTime) > 0
                simInput = simInput.setModelParameter("StopTime", char(obj.StopTime));
            end
            simInput = obj.ConfigureFcn(simInput, trial);
            if ~isa(simInput, "Simulink.SimulationInput")
                error("radia:optuna:SimulinkInput", ...
                    "ConfigureFcn must return a Simulink.SimulationInput object.");
            end
            if obj.UseFastRestart
                simInput = simInput.setModelParameter("FastRestart", "on");
            end
            simOut = sim(simInput);
            score = obj.ScoreFcn(simOut, trial);
            if ~isnumeric(score) || isempty(score) || any(~isfinite(score),"all")
                error("radia:optuna:SimulinkScore", ...
                    "ScoreFcn must return finite numeric objective values.");
            end
            if ~isempty(obj.ResultFcn)
                obj.ResultFcn(simOut, trial);
            end
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
                options.ShowProgress (1,1) logical = true
                options.TransferBaseWorkspaceVariables (1,1) logical = true
            end
            if isempty(ver("parallel"))
                error("radia:optuna:ParallelUnavailable", ...
                    "Parallel Computing Toolbox is required for Simulink parallel trials.");
            end
            trials = cell(nTrials,1);
            modelFolder = string(fileparts(obj.Model));
            if strlength(modelFolder) > 0
                addpath(modelFolder);
            end
            inputs(1,nTrials) = Simulink.SimulationInput(obj.Model);
            for k = 1:nTrials
                trial = study.ask(); trials{k} = trial;
                simInput = Simulink.SimulationInput(obj.Model);
                if strlength(obj.StopTime) > 0
                    simInput = simInput.setModelParameter("StopTime",char(obj.StopTime));
                end
                simInput = obj.ConfigureFcn(simInput,trial);
                if ~isa(simInput,"Simulink.SimulationInput")
                    error("radia:optuna:SimulinkInput", ...
                        "ConfigureFcn must return a Simulink.SimulationInput object.");
                end
                inputs(k) = simInput;
            end
            outputs = parsim(inputs, ...
                "ShowProgress",options.ShowProgress, ...
                "UseFastRestart",obj.UseFastRestart, ...
                "SetupFcn",@()radia.optuna.internal.setupSimulinkWorker(modelFolder), ...
                "TransferBaseWorkspaceVariables",options.TransferBaseWorkspaceVariables);
            for k = 1:nTrials
                trial = trials{k};
                try
                    if ~isempty(outputs(k).ErrorMessage)
                        error("radia:optuna:SimulinkTrial", "%s", outputs(k).ErrorMessage);
                    end
                    score = obj.ScoreFcn(outputs(k),trial);
                    study.tell(trial,score);
                    obj.ResultFcn(outputs(k),trial);
                catch exception
                    if trial.State == "RUNNING", study.fail(trial,exception.message); end
                    rethrow(exception);
                end
            end
            results = study.TrialTable;
        end
    end
end
