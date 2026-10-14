classdef AdjointRunner < handle
    %ADJOINTRUNNER Simulink-facing runner for one checked adjoint study.
    properties
        InitialDesign (:,1) double
        EvaluateFcn
        Options (1,1) struct = struct
        Solver = "mma"
        Metadata (1,1) struct = struct
        Result = []
    end

    methods
        function obj = AdjointRunner(initialDesign, evaluateFcn, options)
            arguments
                initialDesign (:,1) double {mustBeReal,mustBeFinite}
                evaluateFcn (1,1) function_handle
                options.Solver (1,1) string ...
                    {mustBeMember(options.Solver,["mma","sqp"])} = "mma"
                options.OptimizerOptions (1,1) struct = struct
                options.Metadata (1,1) struct = struct
            end
            if ~isa(evaluateFcn,"function_handle")
                error("radia:topopt:AdjointRunner", ...
                    "evaluateFcn must be a function handle.");
            end
            obj.InitialDesign = initialDesign;
            obj.EvaluateFcn = evaluateFcn;
            obj.Solver = options.Solver;
            obj.Options = options.OptimizerOptions;
            obj.Metadata = options.Metadata;
        end

        function setSolver(obj, solver)
            if isnumeric(solver) && isscalar(solver) && isfinite(solver) && ...
                    solver == floor(solver) && solver >= 1 && solver <= 2
                choices = ["mma","sqp"];
                solver = choices(double(solver));
            else
                solver = strip(string(solver));
            end
            if ~isscalar(solver) || ~ismember(solver,["mma","sqp"])
                error("radia:topopt:AdjointRunner", ...
                    "Solver must be mma or sqp.");
            end
            obj.Solver = solver;
        end

        function result = run(obj)
            options = obj.Options;
            obj.setSolver(obj.Solver);
            options.Solver = string(obj.Solver);
            args = namedargs2cell(options);
            result = radia.topopt.optimizeAdjoint( ...
                obj.InitialDesign, obj.EvaluateFcn, args{:});
            obj.Result = result;
        end
    end
end
