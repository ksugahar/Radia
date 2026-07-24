classdef IntersectionSearchSpace
    %INTERSECTIONSEARCHSPACE Infer stable distributions shared by trials.

    methods (Static)
        function searchSpace = calculate(study, options)
            arguments
                study (1,1) radia.optuna.Study
                options.IncludePruned (1,1) logical = true
                options.NumericOnly (1,1) logical = false
            end
            searchSpace = ...
                radia.optuna.internal.IntersectionSearchSpace.empty();
            finished = study.TrialTable.State == "COMPLETE";
            if options.IncludePruned
                finished = finished | study.TrialTable.State == "PRUNED";
            end
            trialNumbers = study.TrialTable.TrialNumber(finished);
            if isempty(trialNumbers)
                return
            end

            params = study.ParamTable;
            commonNames = unique(params.Name( ...
                params.TrialNumber == trialNumbers(1)), "stable");
            for index = 2:numel(trialNumbers)
                names = unique(params.Name( ...
                    params.TrialNumber == trialNumbers(index)), "stable");
                commonNames = intersect(commonNames, names, "stable");
                if isempty(commonNames)
                    return
                end
            end

            for name = reshape(sort(commonNames), 1, [])
                rows = params.Name == name & ...
                    ismember(params.TrialNumber, trialNumbers);
                selected = params(rows,:);
                if height(selected) ~= numel(trialNumbers)
                    continue
                end
                distribution = ...
                    radia.optuna.internal.DistributionCodec.decode( ...
                    selected.Kind(1), selected.Distribution(1));
                compatible = true;
                for row = 2:height(selected)
                    candidate = ...
                        radia.optuna.internal.DistributionCodec.decode( ...
                        selected.Kind(row), selected.Distribution(row));
                    if ~radia.optuna.internal.DistributionCodec.equivalent( ...
                            distribution, candidate)
                        compatible = false;
                        break
                    end
                end
                if ~compatible || ...
                        radia.optuna.internal.DistributionCodec.isSingle( ...
                        distribution)
                    continue
                end
                if options.NumericOnly && distribution.kind == "categorical"
                    continue
                end
                searchSpace(end+1,1) = struct( ...
                    "name", name, ...
                    "distribution", distribution); %#ok<AGROW>
            end
        end

        function searchSpace = empty()
            template = struct( ...
                "name", "", ...
                "distribution", ...
                radia.optuna.internal.DistributionCodec.float( ...
                0, 1, false, NaN));
            searchSpace = reshape(template([]), 0, 1);
        end
    end
end
