classdef OptimizationParameter
    %OPTIMIZATIONPARAMETER One tunable quantity, described the way Simulink
    %Design Optimization describes a model parameter.
    %
    %   sdo.getParameterFromModel returns objects carrying Name, Value,
    %   Minimum, Maximum and Free, and a user narrows the search by editing
    %   Minimum/Maximum before calling sdo.optimize. radia.optuna.optimize
    %   accepts the same shape, so a reader who knows Simulink Design
    %   Optimization can read an Optuna study without learning a second
    %   vocabulary.
    %
    %   Type and Transform carry the part Optuna needs and sdo has no word
    %   for: whether the parameter is continuous, integer or categorical, and
    %   whether it is searched on a log scale.  They map onto
    %   FloatDistribution, IntDistribution and CategoricalDistribution.
    %
    %   Example:
    %       p = radia.optuna.OptimizationParameter("Kp", ...
    %           Minimum=0, Maximum=10);
    %       p(2) = radia.optuna.OptimizationParameter("Ki", ...
    %           Minimum=1e-3, Maximum=1, Transform="log");
    %
    %   See also radia.optuna.optimize, radia.optuna.optimoptions,
    %   radia.optuna.getParameterFromModel.

    properties
        Name (1,1) string
        Value double = []
        Minimum (1,1) double = -Inf
        Maximum (1,1) double = Inf
        Free (1,1) logical = true
        Type (1,1) string {mustBeMember(Type, ...
            ["continuous", "integer", "categorical"])} = "continuous"
        Transform (1,1) string {mustBeMember(Transform, ...
            ["linear", "log"])} = "linear"
        Step double = []
        Choices = []
    end

    methods
        function obj = OptimizationParameter(name, options)
            arguments
                name (1,1) string = "parameter"
                options.Value double = []
                options.Minimum (1,1) double = -Inf
                options.Maximum (1,1) double = Inf
                options.Free (1,1) logical = true
                options.Type (1,1) string = "continuous"
                options.Transform (1,1) string = "linear"
                options.Step double = []
                options.Choices = []
            end
            obj.Name = name;
            obj.Value = options.Value;
            obj.Minimum = options.Minimum;
            obj.Maximum = options.Maximum;
            obj.Free = options.Free;
            obj.Type = options.Type;
            obj.Transform = options.Transform;
            obj.Step = options.Step;
            obj.Choices = options.Choices;
        end

        function distribution = distribution(obj)
            %DISTRIBUTION The Optuna distribution this parameter describes.
            obj.mustBeSearchable();
            switch obj.Type
                case "categorical"
                    distribution = radia.optuna.CategoricalDistribution( ...
                        obj.Choices);
                case "integer"
                    args = {};
                    if ~isempty(obj.Step)
                        args = {"Step", obj.Step};
                    end
                    distribution = radia.optuna.IntDistribution( ...
                        obj.Minimum, obj.Maximum, ...
                        "Log", obj.Transform == "log", args{:});
                otherwise
                    args = {};
                    if ~isempty(obj.Step)
                        args = {"Step", obj.Step};
                    end
                    distribution = radia.optuna.FloatDistribution( ...
                        obj.Minimum, obj.Maximum, ...
                        "Log", obj.Transform == "log", args{:});
            end
        end

        function mustBeSearchable(obj)
            %MUSTBESEARCHABLE Reject a parameter no sampler could search.
            %   A missing bound is not defaulted to something plausible: an
            %   invented range silently changes the study, and the user
            %   cannot tell from the result that it happened.
            if obj.Type == "categorical"
                if isempty(obj.Choices)
                    error("radia:optuna:ParameterChoices", ...
                        "Categorical parameter '%s' has no Choices.", ...
                        obj.Name);
                end
                return
            end
            if ~isfinite(obj.Minimum) || ~isfinite(obj.Maximum)
                error("radia:optuna:ParameterBounds", ...
                    "Parameter '%s' needs finite Minimum and Maximum; " + ...
                    "got [%g, %g]. Set the bounds the study should search.", ...
                    obj.Name, obj.Minimum, obj.Maximum);
            end
            if obj.Minimum > obj.Maximum
                error("radia:optuna:ParameterBounds", ...
                    "Parameter '%s' has Minimum %g above Maximum %g.", ...
                    obj.Name, obj.Minimum, obj.Maximum);
            end
            if obj.Transform == "log" && obj.Minimum <= 0
                error("radia:optuna:ParameterBounds", ...
                    "Parameter '%s' is searched on a log scale, so " + ...
                    "Minimum must be positive; got %g.", ...
                    obj.Name, obj.Minimum);
            end
        end
    end

    methods (Static)
        function parameters = fromStruct(values, options)
            %FROMSTRUCT Build a parameter array from a bounds struct.
            %   radia.optuna.OptimizationParameter.fromStruct( ...
            %       struct("Kp", [0 10], "Ki", [1e-3 1]))
            arguments
                values (1,1) struct
                options.Transform (1,1) string = "linear"
                options.Type (1,1) string = "continuous"
            end
            names = string(fieldnames(values));
            parameters = radia.optuna.OptimizationParameter.empty(1, 0);
            for index = 1:numel(names)
                bounds = values.(names(index));
                if isnumeric(bounds) && numel(bounds) == 2
                    parameters(index) = radia.optuna.OptimizationParameter( ...
                        names(index), Minimum=bounds(1), Maximum=bounds(2), ...
                        Transform=options.Transform, Type=options.Type);
                else
                    parameters(index) = radia.optuna.OptimizationParameter( ...
                        names(index), Type="categorical", Choices=bounds);
                end
            end
        end
    end
end
