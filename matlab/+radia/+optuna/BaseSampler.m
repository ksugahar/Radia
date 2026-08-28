classdef (Abstract) BaseSampler < handle
    %BASESAMPLER Optuna 4.9 sampler interface shared by every sampler.
    %   Optuna drives a sampler through sample_independent, sample_relative,
    %   infer_relative_search_space, before_trial, after_trial and
    %   reseed_rng. radia.optuna's samplers are written around the
    %   suggest-shaped entry points Trial actually calls (sampleFloat,
    %   sampleInteger, sampleCategorical), so this class supplies the
    %   upstream-named interface on top of them rather than duplicating each
    %   sampler's logic.
    %
    %   sample_relative is deliberately NOT defined here. Returning an empty
    %   map would be the right answer for RandomSampler and wrong for
    %   TPESampler, CmaEsSampler, NSGAIISampler and GPSampler, whose relative
    %   proposals are produced inside beforeTrial. A base that answered for
    %   them would be a plausible wrong number, which is worse than a missing
    %   method; the samplers that own a relative search space must implement
    %   it themselves.
    %
    %   reseed_rng is omitted for the same reason. A generic implementation
    %   cannot reach every sampler's stream -- GridSampler keeps its stream
    %   fully private and QMCSampler owns none -- so it would either fail or
    %   silently do nothing on exactly the samplers a caller cares about.

    methods
        function value = sample_independent(obj, study, trial, name, distribution)
            %SAMPLE_INDEPENDENT Sample one parameter from its own distribution.
            %   Dispatches on the distribution kind exactly as Trial does, so
            %   every sampler answers with the same value it would have
            %   produced through suggest_float / suggest_int /
            %   suggest_categorical.
            arguments
                obj
                study (1,1) radia.optuna.Study
                trial
                name (1,1) string
                distribution (1,1) struct
            end
            if ~isfield(distribution, "kind")
                error("radia:optuna:Distribution", ...
                    "A distribution requires a kind field.");
            end
            if radia.optuna.internal.DistributionCodec.isSingle(distribution)
                % Optuna's Trial._suggest never reaches a sampler for a
                % single-valued distribution.
                if distribution.kind == "categorical"
                    value = radia.optuna.internal.DistributionCodec.choiceAt( ...
                        distribution.choices, 1);
                else
                    value = distribution.low;
                end
                return
            end
            switch string(distribution.kind)
                case "float"
                    value = obj.sampleFloat(study, trial, name, ...
                        distribution.low, distribution.high, ...
                        struct("Log", distribution.log, ...
                        "Step", distribution.step));
                case "integer"
                    value = obj.sampleIntegerParameter(study, trial, name, ...
                        distribution);
                case "categorical"
                    value = obj.sampleCategorical(study, trial, name, ...
                        distribution.choices);
                otherwise
                    error("radia:optuna:DistributionKind", ...
                        "Unsupported distribution kind '%s'.", ...
                        string(distribution.kind));
            end
        end

        function searchSpace = infer_relative_search_space(obj, study, trial) %#ok<INUSD>
            %INFER_RELATIVE_SEARCH_SPACE No relative space unless overridden.
            searchSpace = radia.optuna.internal.IntersectionSearchSpace.empty();
        end

        function before_trial(obj, study, trial)
            %BEFORE_TRIAL Upstream-named alias of the beforeTrial hook.
            if ismethod(obj, "beforeTrial")
                obj.beforeTrial(study, trial);
            end
        end

        function after_trial(obj, study, trial, state, values) %#ok<INUSD>
            %AFTER_TRIAL Upstream-named alias of the afterTrial hook.
            %   STATE and VALUES are accepted for signature compatibility;
            %   radia.optuna reads both from the trial itself.
            if ismethod(obj, "afterTrial")
                obj.afterTrial(study, trial);
            end
        end
    end

    methods (Access=private)
        function value = sampleIntegerParameter(obj, study, trial, name, ...
                distribution)
            step = distribution.step;
            if ~distribution.log && step == 1
                value = obj.sampleInteger(study, trial, name, ...
                    distribution.low, distribution.high);
                return
            end
            value = obj.sampleFloat(study, trial, name, distribution.low, ...
                distribution.high, ...
                struct("Log", distribution.log, "Step", step));
            value = distribution.low + ...
                radia.optuna.internal.UpstreamNumerics.roundTiesToEven( ...
                (double(value) - distribution.low) / step) * step;
            value = min(max(value, distribution.low), distribution.high);
        end
    end
end
