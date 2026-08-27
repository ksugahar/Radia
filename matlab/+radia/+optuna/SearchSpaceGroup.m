classdef SearchSpaceGroup < handle
    %SEARCHSPACEGROUP MATLAB name for Optuna 4.9's _SearchSpaceGroup.

    properties (Access=private)
        Core (1,1) radia.optuna.internal.GroupDecomposedSearchSpace
    end

    methods
        function obj=SearchSpaceGroup(groups)
            obj.Core=radia.optuna.internal.GroupDecomposedSearchSpace();
            if nargin>0
                for index=1:numel(groups)
                    group=groups{index};
                    obj.Core.update(string({group.name}), ...
                        reshape({group.distribution},[],1));
                end
            end
        end

        function add_distributions(obj,distributions)
            if isa(distributions,"containers.Map")
                names=string(distributions.keys);
                values=cell(numel(names),1);
                for index=1:numel(names)
                    values{index}=distributions(char(names(index)));
                end
            elseif isstruct(distributions) && isscalar(distributions)
                names=string(fieldnames(distributions));
                values=cell(numel(names),1);
                for index=1:numel(names)
                    values{index}=distributions.(names(index));
                end
            else
                error("radia:optuna:SearchSpaceGroup", ...
                    "distributions must be a scalar struct or containers.Map.");
            end
            obj.Core.update(names,values);
        end

        function spaces=search_spaces(obj)
            groups=obj.Core.current(ExcludeSingle=false);
            spaces=cell(numel(groups),1);
            for index=1:numel(groups)
                group=groups{index};
                current=containers.Map("KeyType","char","ValueType","any");
                for item=1:numel(group)
                    current(char(group(item).name))=group(item).distribution;
                end
                spaces{index}=current;
            end
        end
    end
end
