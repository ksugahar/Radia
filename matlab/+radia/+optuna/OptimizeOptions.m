classdef OptimizeOptions
    %OPTIMIZEOPTIONS Solver options, named the way optimoptions names them.
    %
    %   Global Optimization Toolbox solvers are configured through one
    %   options object rather than a long argument list, and the option names
    %   are shared across solvers: MaxTime, FunctionTolerance, Display,
    %   OutputFcn, PlotFcn, UseParallel.  radia.optuna.optimize takes the same
    %   object so an Optuna study is set up the way ga or particleswarm is.
    %
    %   Where an Optuna concept has no toolbox equivalent it keeps its Optuna
    %   name -- Sampler, Pruner, Seed, StoragePath, Directions -- rather than
    %   being bent into a toolbox option that means something else.
    %
    %   Names that differ deliberately:
    %     MaxTrials      ga calls this MaxGenerations and fmincon calls it
    %                    MaxIterations.  A trial is Optuna's unit of work and
    %                    is what the study records, so the option is named
    %                    after it.
    %     MaxStallTrials ga's MaxStallGenerations, counted in trials.
    %
    %   See also radia.optuna.optimoptions, radia.optuna.optimize.

    properties
        MaxTrials (1,1) double {mustBeInteger, mustBePositive} = 100
        MaxTime (1,1) double {mustBePositive} = Inf
        FunctionTolerance (1,1) double {mustBeNonnegative} = 0
        MaxStallTrials (1,1) double {mustBePositive} = Inf
        Display (1,1) string {mustBeMember(Display, ...
            ["off", "none", "final", "iter"])} = "final"
        OutputFcn = []
        PlotFcn = []
        UseParallel (1,1) logical = false
        ParallelMode (1,1) string {mustBeMember(ParallelMode, ...
            ["sequential", "batch", "steady-state"])} = "sequential"
        BatchSize (1,1) double {mustBeInteger, mustBeNonnegative} = 0
        Sampler (1,1) string {mustBeMember(Sampler, ...
            ["auto", "random", "tpe", "cmaes", "motpe", "nsgaii", ...
            "nsgaiii", "gp", "bruteforce", "qmc"])} = "auto"
        Pruner (1,1) string {mustBeMember(Pruner, ...
            ["none", "median", "hyperband", "percentile", "patient", ...
            "successivehalving", "threshold"])} = "none"
        Seed = []
        StoragePath (1,1) string = ""
        StudyName (1,1) string = ""
        Resume (1,1) logical = true
        Directions string = "minimize"
        ConstraintFcn = []
        CatchObjectiveErrors (1,1) logical = true
    end

    methods
        function obj = OptimizeOptions(options)
            arguments
                options.MaxTrials (1,1) double = 100
                options.MaxTime (1,1) double = Inf
                options.FunctionTolerance (1,1) double = 0
                options.MaxStallTrials (1,1) double = Inf
                options.Display (1,1) string = "final"
                options.OutputFcn = []
                options.PlotFcn = []
                options.UseParallel (1,1) logical = false
                options.ParallelMode (1,1) string = "sequential"
                options.BatchSize (1,1) double = 0
                options.Sampler (1,1) string = "auto"
                options.Pruner (1,1) string = "none"
                options.Seed = []
                options.StoragePath (1,1) string = ""
                options.StudyName (1,1) string = ""
                options.Resume (1,1) logical = true
                options.Directions string = "minimize"
                options.ConstraintFcn = []
                options.CatchObjectiveErrors (1,1) logical = true
            end
            names = string(fieldnames(options));
            for index = 1:numel(names)
                obj.(names(index)) = options.(names(index));
            end
            obj.OutputFcn = obj.normalizeCallbacks(obj.OutputFcn, "OutputFcn");
            obj.PlotFcn = obj.normalizeCallbacks(obj.PlotFcn, "PlotFcn");
            if ~isempty(obj.ConstraintFcn) && ...
                    ~isa(obj.ConstraintFcn, "function_handle")
                error("radia:optuna:OptimizeOptions", ...
                    "ConstraintFcn must be a function handle.");
            end
            if isfinite(obj.MaxStallTrials) && ...
                    obj.MaxStallTrials ~= round(obj.MaxStallTrials)
                error("radia:optuna:OptimizeOptions", ...
                    "MaxStallTrials counts trials, so it must be a " + ...
                    "positive integer or Inf; got %g.", obj.MaxStallTrials);
            end
            if any(~ismember(obj.Directions, ["minimize", "maximize"]))
                error("radia:optuna:OptimizeOptions", ...
                    "Directions must contain only minimize or maximize.");
            end
            if obj.UseParallel && obj.ParallelMode == "sequential"
                obj.ParallelMode = "batch";
            elseif obj.ParallelMode ~= "sequential"
                obj.UseParallel = true;
            end
        end

        function tf = isQuiet(obj)
            tf = ismember(obj.Display, ["off", "none"]);
        end
    end

    methods (Static, Access=private)
        function handles = normalizeCallbacks(value, label)
            if isempty(value)
                handles = {};
                return
            end
            if isa(value, "function_handle")
                handles = {value};
                return
            end
            if iscell(value) && all(cellfun( ...
                    @(h) isa(h, "function_handle"), value))
                handles = reshape(value, 1, []);
                return
            end
            error("radia:optuna:OptimizeOptions", ...
                "%s must be a function handle or a cell array of handles.", ...
                label);
        end
    end
end
