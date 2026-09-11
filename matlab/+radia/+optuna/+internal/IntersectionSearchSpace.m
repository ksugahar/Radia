classdef IntersectionSearchSpace
    %INTERSECTIONSEARCHSPACE Infer stable distributions shared by trials.

    methods (Static)
        function searchSpace = calculate(study, options)
            arguments
                study (1,1) radia.optuna.Study
                options.IncludePruned (1,1) logical = true
                options.NumericOnly (1,1) logical = false
                options.ExcludeSingle (1,1) logical = true
            end
            searchSpace = ...
                radia.optuna.internal.IntersectionSearchSpace.empty();
            trials=study.trialData();
            finished = trials.State == "COMPLETE";
            if options.IncludePruned
                finished = finished | trials.State == "PRUNED";
            end
            trialNumbers = trials.TrialNumber(finished);
            if isempty(trialNumbers)
                return
            end

            % Use the Study's columnar arrays directly. Constructing and
            % indexing MATLAB tables here dominated short TPE studies even
            % though the intersection itself contained only a few names.
            params = study.parameterData();
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
                selected = find(rows);
                if numel(selected) ~= numel(trialNumbers)
                    continue
                end
                distribution = ...
                    radia.optuna.internal.DistributionCodec.decode( ...
                    params.Kind(selected(1)), ...
                    params.Distribution(selected(1)));
                compatible = true;
                for row = 2:numel(selected)
                    if params.Kind(selected(row))==params.Kind(selected(1)) && ...
                            params.Distribution(selected(row))== ...
                            params.Distribution(selected(1))
                        continue
                    end
                    candidate = ...
                        radia.optuna.internal.DistributionCodec.decode( ...
                        params.Kind(selected(row)), ...
                        params.Distribution(selected(row)));
                    if ~radia.optuna.internal.DistributionCodec.equivalent( ...
                            distribution, candidate)
                        compatible = false;
                        break
                    end
                end
                if ~compatible || (options.ExcludeSingle && ...
                        radia.optuna.internal.DistributionCodec.isSingle( ...
                        distribution))
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
