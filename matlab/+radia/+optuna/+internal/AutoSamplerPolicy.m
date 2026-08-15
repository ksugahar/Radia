classdef AutoSamplerPolicy
    %AUTOSAMPLERPOLICY Budget/search-space routing aligned with OptunaHub.

    methods (Static)
        function [name,reason]=choose(spec,nObjectives,nTrials)
            isMultiObjective=nObjectives>1;
            if isMultiObjective
                if spec.has_categorical || spec.is_conditional
                    name="motpe";
                    reason="multiobjective_mixed_or_conditional";
                elseif spec.fixed_numeric && nTrials<=250 && nObjectives<4
                    name="gp";
                    reason="fixed_numeric_small_budget_few_objectives";
                elseif spec.fixed_numeric && nTrials>250 && nObjectives>=4
                    name="nsgaiii";
                    reason="fixed_numeric_many_objectives_population_budget";
                elseif spec.fixed_numeric && nTrials>250
                    name="nsgaii";
                    reason="fixed_numeric_population_budget";
                else
                    name="motpe";
                    reason="multiobjective_unknown_space_safe_default";
                end
                return
            end

            minimumCmaTrials=Inf;
            if isfinite(spec.dimensions)
                minimumCmaTrials=max(40,8*spec.dimensions);
            end
            cmaEligible=spec.fixed_numeric && spec.dimensions>=2 && ...
                spec.constraints_declared && ~spec.has_constraints && ...
                nTrials>250 && nTrials>=minimumCmaTrials;
            gpEligible=spec.fixed_numeric && nTrials<=250;
            if gpEligible
                name="gp";
                reason="fixed_numeric_small_budget";
            elseif cmaEligible
                name="cmaes";
                reason="declared_fixed_numeric_correlated_budget";
            elseif spec.has_constraints
                name="tpe";
                reason="constrained_space";
            elseif spec.has_categorical || spec.is_conditional
                name="tpe";
                reason="categorical_or_conditional_space";
            elseif spec.fixed_numeric
                name="tpe";
                reason="fixed_numeric_budget_or_metadata_insufficient_for_cmaes";
            else
                name="tpe";
                reason="unknown_space_safe_default";
            end
        end
    end
end
