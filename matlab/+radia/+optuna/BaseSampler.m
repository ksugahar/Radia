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

        function value=sample_independent(obj,study,trial,param_name,distribution)
            spec=radia.optuna.internal.DistributionCodec.normalize(distribution);
            switch spec.kind
                case "float"
                    value=obj.sampleFloat(study,trial,param_name, ...
                        spec.low,spec.high,struct("Log",spec.log,"Step",spec.step));
                case "integer"
                    if ismethod(obj,"sampleIntegerDetailed")
                        value=obj.sampleIntegerDetailed(study,trial,param_name, ...
                            spec.low,spec.high,spec.step,spec.log);
                    elseif ~spec.log && spec.step==1
                        value=obj.sampleInteger(study,trial,param_name, ...
                            spec.low,spec.high);
                    else
                        value=obj.sampleFloat(study,trial,param_name, ...
                            spec.low,spec.high,struct("Log",spec.log,"Step",spec.step));
                        value=spec.low+round((double(value)-spec.low)/spec.step)*spec.step;
                        value=min(max(value,spec.low),spec.high);
                    end
                case "categorical"
                    value=obj.sampleCategorical(study,trial,param_name,spec.choices);
                otherwise
                    error("radia:optuna:DistributionKind", ...
                        "Unsupported distribution kind '%s'.",spec.kind);
            end
        end

        function params=sample_relative(~,~,trial,search_space)
            % Study.ask has already run each sampler's relative phase. This
            % public upstream adapter returns that trial-visible result.
            params=struct();
            if isempty(search_space) || ~isprop(trial,"relative_params")
                return
            end
            available=trial.relative_params;
            if isa(search_space,"containers.Map")
                names=string(keys(search_space));
            elseif isstruct(search_space) && ...
                    all(isfield(search_space,["name","distribution"]))
                names=reshape(string({search_space.name}),1,[]);
            elseif isstruct(search_space) && isscalar(search_space)
                names=reshape(string(fieldnames(search_space)),1,[]);
            else
                error("radia:optuna:SearchSpace", ...
                    "search_space must be a distribution struct or containers.Map.");
            end
            for name=names
                key=matlab.lang.makeValidName(name);
                if isfield(available,key)
                    params.(key)=available.(key);
                end
            end
        end

        function reseed_rng(~)
            % The upstream BaseSampler default is intentionally a no-op.
        end
    end
end
