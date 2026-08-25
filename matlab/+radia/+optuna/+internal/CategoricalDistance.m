classdef CategoricalDistance
    %CATEGORICALDISTANCE Validate and evaluate TPE categorical distances.

    methods (Static)
        function validate(mapping)
            if isempty(mapping)
                return
            end
            if isa(mapping,"containers.Map")
                values=mapping.values;
                for index=1:numel(values)
                    if ~isa(values{index},"function_handle")
                        error("radia:optuna:TPECategoricalDistance", ...
                            "Every categorical distance must be a function handle.");
                    end
                end
                return
            end
            if isstruct(mapping) && isscalar(mapping)
                names=fieldnames(mapping);
                for index=1:numel(names)
                    if ~isa(mapping.(names{index}),"function_handle")
                        error("radia:optuna:TPECategoricalDistance", ...
                            "Every categorical distance must be a function handle.");
                    end
                end
                return
            end
            error("radia:optuna:TPECategoricalDistance", ...
                "CategoricalDistanceFcn must be a containers.Map or scalar struct of function handles.");
        end

        function distanceFcn=get(mapping,name)
            distanceFcn=[];
            if isempty(mapping)
                return
            end
            if isa(mapping,"containers.Map")
                key=char(string(name));
                if isKey(mapping,key)
                    distanceFcn=mapping(key);
                end
                return
            end
            key=char(string(name));
            if isfield(mapping,key)
                distanceFcn=mapping.(key);
            end
        end
    end
end
