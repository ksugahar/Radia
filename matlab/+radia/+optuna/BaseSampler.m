classdef (Abstract) BaseSampler < handle
    %BASESAMPLER Common type for Optuna 4.9 sampler implementations.

    methods
        function before_trial(obj,study,trial)
            if ismethod(obj,"beforeTrial")
                obj.beforeTrial(study,trial);
            end
        end

        function after_trial(obj,study,trial,state,values) %#ok<INUSD>
            if ismethod(obj,"afterTrial")
                obj.afterTrial(study,trial);
            end
        end

        function searchSpace=infer_relative_search_space(obj,study,trial)
            if nargin<3
                trial=[];
            end
            if ismethod(obj,"inferRelativeSearchSpace")
                searchSpace=obj.inferRelativeSearchSpace(study,trial);
            else
                searchSpace=struct("name",{},"distribution",{});
            end
        end
    end
end
