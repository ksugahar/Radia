classdef (Abstract) BaseMutation < handle
    %BASEMUTATION Base class for Optuna 5 numerical mutation operators.

    methods (Abstract)
        value = mutation(obj,param,rng,study,searchSpaceBounds)
    end

    methods
        function value=string(obj)
            value=string(class(obj));
            parts=split(value,".");
            value=parts(end);
        end

        function value=char(obj)
            value=char(string(obj));
        end
    end
end
