classdef TPESampler < radia.optuna.BaseSampler
    %TPESAMPLER Optuna-style tree-structured Parzen estimator sampler.
    %   The default settings and univariate Parzen construction follow the
    %   Optuna 4.9 TPESampler: ten random startup trials, a ten-percent good
    %   set capped at 25 trials, 24 expected-improvement candidates,
    %   observation-specific truncated-normal kernels, history weights,
    %   an explicit prior, and magic clipping.
    %   Constraint-aware splits use c <= 0 as feasible, rank infeasible
    %   trials by total positive violation, and rank missing/nonfinite
    %   constraint records last.

    properties (SetAccess=private)
        Stream
        Seed (1,1) double = 0
        NStartupTrials (1,1) double = 10
        Gamma (1,1) double = 0.1
        MaxGoodTrials (1,1) double = 25
        GammaFcn = []
        WeightsFcn = []
        NumberOfEIChoices (1,1) double = 24
        PriorWeight (1,1) double = 1
        ConsiderMagicClip (1,1) logical = true
        ConsiderEndpoints (1,1) logical = false
        Multivariate (1,1) logical = false
        Group (1,1) logical = false
        WarnIndependentSampling (1,1) logical = false
        ConstantLiar (1,1) logical = false
        ConstraintsFcn = []
        CategoricalDistanceFcn = []
    end

    properties (Access=private)
        AttachedStudy = []
        Restored (1,1) logical = false
        MultiObjectiveSampler
        IndependentSampler
        GroupDecomposition
        EncodingCacheNames string = strings(0,1)
        EncodingCacheDistributions cell = cell(0,1)
        EncodingCacheValues string = strings(0,1)
        NativeHistoryValid (1,1) logical = false
        NativeHistoryCompleteCount (1,1) double = 0
        HistoryDistributionNames string = strings(0,1)
        HistoryDistributions cell = cell(0,1)
        NativeGroupRevision (1,1) double = -1
        NativeGroupMetadata struct = struct()
    end

    properties (Constant, Access=private)
        StateSchema = "radia.optuna.tpe-sampler-state.v1"
        SamplerName = "tpe"
    end

    methods
        function obj = TPESampler(options)
            arguments
                options.Seed double = double.empty(1,0)
                options.NStartupTrials (1,1) double = 10
                options.Gamma (1,1) double = 0.1
                options.MaxGoodTrials (1,1) double = 25
                options.GammaFcn = []
                options.WeightsFcn = []
                options.NumberOfEIChoices (1,1) double = 24
                options.PriorWeight (1,1) double = 1
                options.ConsiderMagicClip (1,1) logical = true
                options.ConsiderEndpoints (1,1) logical = false
                options.Multivariate (1,1) logical = false
                options.Group (1,1) logical = false
                options.WarnIndependentSampling (1,1) logical = false
                options.ConstantLiar (1,1) logical = false
                options.ConstraintsFcn = []
                options.CategoricalDistanceFcn = []
            end
            if options.NStartupTrials < 0 || ...
                    options.NStartupTrials ~= floor(options.NStartupTrials)
                error("radia:optuna:TPEStartup", ...
                    "NStartupTrials must be a nonnegative integer.");
            end
            if ~(options.Gamma > 0 && options.Gamma <= 1)
                error("radia:optuna:TPEGamma", ...
                    "Gamma must be greater than zero and at most one.");
            end
            if options.MaxGoodTrials < 1 || ...
                    options.MaxGoodTrials ~= floor(options.MaxGoodTrials)
                error("radia:optuna:TPEMaxGood", ...
                    "MaxGoodTrials must be a positive integer.");
            end
            if ~isempty(options.GammaFcn) && ...
                    ~isa(options.GammaFcn, "function_handle")
                error("radia:optuna:TPEGamma", ...
                    "GammaFcn must be a function handle.");
            end
            if ~isempty(options.WeightsFcn) && ...
                    ~isa(options.WeightsFcn, "function_handle")
                error("radia:optuna:TPEWeights", ...
                    "WeightsFcn must be a function handle.");
            end
            if options.NumberOfEIChoices < 1 || ...
                    options.NumberOfEIChoices ~= floor(options.NumberOfEIChoices)
                error("radia:optuna:TPECandidates", ...
                    "NumberOfEIChoices must be a positive integer.");
            end
            if options.PriorWeight < 0 || ~isfinite(options.PriorWeight)
                error("radia:optuna:TPEPriorWeight", ...
                    "PriorWeight must be finite and nonnegative.");
            end
            obj.Seed = radia.optuna.internal.resolveSeed(options.Seed);
            obj.Stream = ...
                radia.optuna.internal.NumpyRandomState(obj.Seed);
            obj.IndependentSampler = radia.optuna.RandomSampler(options.Seed);
            obj.NStartupTrials = options.NStartupTrials;
            obj.Gamma = options.Gamma;
            obj.MaxGoodTrials = options.MaxGoodTrials;
            obj.GammaFcn = options.GammaFcn;
            obj.WeightsFcn = options.WeightsFcn;
            obj.NumberOfEIChoices = options.NumberOfEIChoices;
            obj.PriorWeight = options.PriorWeight;
            obj.ConsiderMagicClip = options.ConsiderMagicClip;
            obj.ConsiderEndpoints = options.ConsiderEndpoints;
            obj.Multivariate = options.Multivariate;
            if options.Group && ~options.Multivariate
                error("radia:optuna:TPEGroup", ...
                    "Group can only be enabled when Multivariate is enabled.");
            end
            obj.Group = options.Group;
            obj.WarnIndependentSampling = options.WarnIndependentSampling;
            obj.ConstantLiar = options.ConstantLiar;
            if ~isempty(options.ConstraintsFcn) && ...
                    ~isa(options.ConstraintsFcn, "function_handle")
                error("radia:optuna:ConstraintsFcn", ...
                    "ConstraintsFcn must be a function handle.");
            end
            obj.ConstraintsFcn = options.ConstraintsFcn;
            radia.optuna.internal.CategoricalDistance.validate( ...
                options.CategoricalDistanceFcn);
            obj.CategoricalDistanceFcn=options.CategoricalDistanceFcn;
            obj.GroupDecomposition= ...
                radia.optuna.internal.GroupDecomposedSearchSpace();
            % Optuna's TPESampler owns both the single- and
            % multi-objective TPE paths. Keep MOTPESampler as an internal
            % compatibility implementation, but route an explicitly chosen
            % TPESampler through it for multi-objective studies.
            obj.MultiObjectiveSampler = radia.optuna.MOTPESampler( ...
                Seed=obj.Seed, NStartupTrials=obj.NStartupTrials, ...
                Gamma=obj.Gamma, MaxGoodTrials=obj.MaxGoodTrials, ...
                GammaFcn=obj.GammaFcn, WeightsFcn=obj.WeightsFcn, ...
                NumberOfEIChoices=obj.NumberOfEIChoices, ...
                PriorWeight=obj.PriorWeight, ...
                ConsiderMagicClip=obj.ConsiderMagicClip, ...
                ConsiderEndpoints=obj.ConsiderEndpoints, ...
                ConstraintsFcn=obj.ConstraintsFcn, ...
                CategoricalDistanceFcn=obj.CategoricalDistanceFcn);
        end

        function value = sampleFloat(obj, study, trial, name, low, high, options)
            if numel(study.Directions) > 1
                value = obj.MultiObjectiveSampler.sampleFloat( ...
                    study, trial, name, low, high, options);
                return
            end
            obj.attach(study);
            obj.warnIndependent(study,trial,name);
            obj.validateBounds(low, high, options.Log, options.Step);
            if low == high
                value = low;
                return
            end
            [x, y, trialNumbers, pending] = ...
                obj.numericObservations(study, name);
            valid = isfinite(x) & isfinite(y) & x >= low & x <= high;
            if options.Log
                valid = valid & x > 0;
            end
            x = x(valid);
            y = y(valid);
            trialNumbers = trialNumbers(valid);
            pending = pending(valid);
            finishedCount = sum(~pending);
            if finishedCount == 0 || finishedCount < obj.NStartupTrials
                value = obj.IndependentSampler.sampleFloat( ...
                    study, trial, name, low, high, options);
                return
            end

            [good, bad] = obj.splitObservations( ...
                x, y, study.Directions(1), study, trialNumbers, pending);
            belowWeights=obj.observationWeights(numel(good));
            aboveWeights=obj.observationWeights(numel(bad));
            nativeHandle=obj.Stream.nativeHandle();
            if nativeHandle~=0 && radia.optuna.internal.NativeKernels.has( ...
                    "optuna.tpe.best_numerical_observations")
                value=radia.optuna.internal.NativeKernels.call( ...
                    "optuna.tpe.best_numerical_observations", ...
                    nativeHandle,obj.NumberOfEIChoices,good,bad,low,high, ...
                    options.Log,options.Step,obj.PriorWeight, ...
                    obj.ConsiderMagicClip,obj.ConsiderEndpoints, ...
                    belowWeights,aboveWeights);
                obj.recordState(study,trial.Number);
                return
            end
            estimatorOptions = { ...
                "Log", options.Log, ...
                "Step", options.Step, ...
                "PriorWeight", obj.PriorWeight, ...
                "ConsiderMagicClip", obj.ConsiderMagicClip, ...
                "ConsiderEndpoints", obj.ConsiderEndpoints};
            below = radia.optuna.internal.ParzenEstimator.numerical( ...
                good, low, high, estimatorOptions{:}, ...
                ObservationWeights=belowWeights);
            above = radia.optuna.internal.ParzenEstimator.numerical( ...
                bad, low, high, estimatorOptions{:}, ...
                ObservationWeights=aboveWeights);
            nativeHandle=obj.Stream.nativeHandle();
            if nativeHandle~=0 && radia.optuna.internal.NativeKernels.has( ...
                    "optuna.tpe.best_numerical")
                value=radia.optuna.internal.NativeKernels.call( ...
                    "optuna.tpe.best_numerical",nativeHandle, ...
                    obj.NumberOfEIChoices,below.weights,below.mu,below.sigma, ...
                    above.weights,above.mu,above.sigma,below.internal_low, ...
                    below.internal_high,below.low,below.high,below.log,below.step);
            else
                candidates = ...
                    radia.optuna.internal.ParzenEstimator.sampleNumerical( ...
                    below, obj.Stream, obj.NumberOfEIChoices);
                acquisition = ...
                    radia.optuna.internal.ParzenEstimator.logPdfNumerical( ...
                    below, candidates) - ...
                    radia.optuna.internal.ParzenEstimator.logPdfNumerical( ...
                    above, candidates);
                [~, best] = max(acquisition);
                value = candidates(best);
            end
            obj.recordState(study, trial.Number);
        end

        function value = sampleInteger(obj, study, trial, name, low, high)
            if numel(study.Directions) > 1
                value = obj.MultiObjectiveSampler.sampleInteger( ...
                    study, trial, name, low, high);
                return
            end
            if low ~= floor(low) || high ~= floor(high) || low > high
                error("radia:optuna:Bounds", ...
                    "Integer bounds must be finite integers with low <= high.");
            end
            if low == high
                value = low;
                return
            end
            value = obj.sampleFloat(study, trial, name, low, high, ...
                struct("Log", false, "Step", 1));
            value = min(max(round(value), low), high);
        end

        function value = sampleCategorical(obj, study, trial, name, choices)
            if numel(study.Directions) > 1
                value = obj.MultiObjectiveSampler.sampleCategorical( ...
                    study, trial, name, choices);
                return
            end
            obj.attach(study);
            obj.warnIndependent(study,trial,name);
            if isempty(choices)
                error("radia:optuna:Choices", ...
                    "Categorical choices must not be empty.");
            end
            [tokens, y, trialNumbers, pending] = ...
                obj.categoricalObservations(study, name);
            choiceTokens = obj.choiceTokens(choices);
            observed = zeros(numel(tokens), 1);
            valid = isfinite(y);
            for index = 1:numel(tokens)
                match = find(choiceTokens == tokens(index), 1);
                if isempty(match)
                    valid(index) = false;
                else
                    observed(index) = match;
                end
            end
            observed = observed(valid);
            y = y(valid);
            trialNumbers = trialNumbers(valid);
            pending = pending(valid);
            count = numel(choiceTokens);
            finishedCount = sum(~pending);
            if finishedCount == 0 || finishedCount < obj.NStartupTrials
                value = obj.IndependentSampler.sampleCategorical( ...
                    study, trial, name, choices);
                return
            end

            [good, bad] = obj.splitObservations( ...
                observed, y, study.Directions(1), study, ...
                trialNumbers, pending);
            below = radia.optuna.internal.ParzenEstimator.categorical( ...
                good, count, PriorWeight=obj.PriorWeight, ...
                ObservationWeights=obj.observationWeights(numel(good)), ...
                DistanceFcn=radia.optuna.internal.CategoricalDistance. ...
                get(obj.CategoricalDistanceFcn,name),Choices=choices);
            above = radia.optuna.internal.ParzenEstimator.categorical( ...
                bad, count, PriorWeight=obj.PriorWeight, ...
                ObservationWeights=obj.observationWeights(numel(bad)), ...
                DistanceFcn=radia.optuna.internal.CategoricalDistance. ...
                get(obj.CategoricalDistanceFcn,name),Choices=choices);
            candidates = ...
                radia.optuna.internal.ParzenEstimator.sampleCategorical( ...
                below, obj.Stream, obj.NumberOfEIChoices);
            acquisition = ...
                radia.optuna.internal.ParzenEstimator.logPdfCategorical( ...
                below, candidates) - ...
                radia.optuna.internal.ParzenEstimator.logPdfCategorical( ...
                above, candidates);
            [~, best] = max(acquisition);
            value = obj.choiceAt(choices, candidates(best));
            obj.recordState(study, trial.Number);
        end

        function values = sampleJoint(obj, study, trial, names, lows, highs, options)
            %SAMPLEJOINT Multivariate TPE with a shared mixture component.
            if numel(study.Directions) > 1
                values = zeros(1, numel(names));
                for index = 1:numel(names)
                    values(index) = obj.MultiObjectiveSampler.sampleFloat( ...
                        study, trial, names(index), lows(index), ...
                        highs(index), struct("Log", options.Log(index), ...
                        "Step", NaN));
                end
                return
            end
            obj.attach(study);
            trials=study.trialData();
            finished = trials.State == "COMPLETE" | ...
                trials.State == "PRUNED";
            if sum(finished) < obj.NStartupTrials
                values = zeros(1,numel(names));
                for index = 1:numel(names)
                    values(index) = obj.IndependentSampler.sampleFloat( ...
                        study,trial,names(index),lows(index),highs(index), ...
                        struct("Log",options.Log(index),"Step",NaN));
                end
                return
            end
            template = struct( ...
                "name", "", ...
                "distribution", ...
                    radia.optuna.internal.DistributionCodec.float( ...
                    0, 1, false, NaN));
            searchSpace = repmat(template, 1, numel(names));
            for index = 1:numel(names)
                searchSpace(index).name = names(index);
                searchSpace(index).distribution = ...
                    radia.optuna.internal.DistributionCodec.float( ...
                    lows(index), highs(index), options.Log(index), NaN);
            end
            sampled = obj.sampleRelativeSpace(study, searchSpace);
            values = zeros(1, numel(sampled));
            for index = 1:numel(sampled)
                values(index) = double(sampled{index});
            end
            obj.recordState(study, trial.Number);
        end

        function searchSpace = inferRelativeSearchSpace(obj, study, trial) %#ok<INUSD>
            if ~obj.Multivariate
                searchSpace = obj.emptySearchSpace();
                return
            end
            if obj.Group
                groups=obj.groupSearchSpaces(study);
                searchSpace=obj.emptySearchSpace();
                for index=1:numel(groups)
                    searchSpace=[searchSpace;groups{index}(:)]; %#ok<AGROW>
                end
            else
                searchSpace = obj.intersectionSearchSpace(study);
            end
        end

        function searchSpace = infer_relative_search_space(obj, study, trial)
            if nargin < 3
                trial = [];
            end
            searchSpace = obj.inferRelativeSearchSpace(study, trial);
        end

        function beforeTrial(obj, study, trial)
            if numel(study.Directions) > 1
                obj.MultiObjectiveSampler.beforeTrial(study, trial);
                return
            end
            if strlength(study.StoragePath)>0
                obj.IndependentSampler.beforeTrial(study, trial);
            end
            obj.attach(study);
            if ~obj.Multivariate || numel(study.Directions) ~= 1
                return
            end
            trials=study.trialData();
            finished = trials.State == "COMPLETE" | ...
                trials.State == "PRUNED";
            if sum(finished) < obj.NStartupTrials
                return
            end
            if obj.Group
                groups=obj.groupSearchSpaces(study);
                if obj.canUseNativeGroupedHistory(study) && ...
                        obj.sampleNativeGroupedHistory(trial,groups)
                    obj.recordState(study, trial.Number);
                    return
                end
                allSpace=obj.emptySearchSpace();
                for index=1:numel(groups)
                    allSpace=[allSpace;groups{index}(:)]; %#ok<AGROW>
                end
                [allValues,allObjectives,allNumbers,allPending]= ...
                    obj.relativeObservations(study,allSpace,true);
                globalValid=isfinite(allObjectives);
                globalData=struct('observations',zeros(sum(globalValid),0), ...
                    'objectives',allObjectives(globalValid), ...
                    'trialNumbers',allNumbers(globalValid), ...
                    'pending',allPending(globalValid));
                globalGoodTrialNumbers=obj.globalGoodTrialNumbers( ...
                    study,globalData);
                allNames=string({allSpace.name});
                for index=1:numel(groups)
                    searchSpace=groups{index};
                    if isempty(searchSpace), continue; end
                    groupNames=string({searchSpace.name});
                    columns=zeros(size(groupNames));
                    for columnIndex=1:numel(groupNames)
                        match=find(allNames==groupNames(columnIndex),1);
                        if ~isempty(match)
                            columns(columnIndex)=match;
                        end
                    end
                    present=columns>0;
                    if ~all(present)
                        error("radia:optuna:GroupSearchSpace", ...
                            "Grouped search-space cache is inconsistent.");
                    end
                    valid=isfinite(allObjectives) & ...
                        all(isfinite(allValues(:,columns)),2);
                    precomputed=struct( ...
                        'observations',allValues(valid,columns), ...
                        'objectives',allObjectives(valid), ...
                        'trialNumbers',allNumbers(valid), ...
                        'pending',allPending(valid));
                    values=obj.sampleRelativeSpace( ...
                        study,searchSpace,EnforceStartup=false, ...
                        GlobalSplit=true, ...
                        GlobalGoodTrialNumbers=globalGoodTrialNumbers, ...
                        Precomputed=precomputed);
                    trial.setRelativeParameters(searchSpace,values,"");
                end
            else
                searchSpace = obj.inferRelativeSearchSpace(study, trial);
                if isempty(searchSpace), return; end
                values = obj.sampleRelativeSpace(study, searchSpace);
                trial.setRelativeParameters(searchSpace,values,"");
            end
            obj.recordState(study, trial.Number);
        end

        function afterTrial(obj, study, trial)
            if numel(study.Directions) > 1
                obj.MultiObjectiveSampler.afterTrial(study, trial);
                return
            end
            if trial.State == "COMPLETE" && ~isempty(obj.ConstraintsFcn)
                study.recordConstraints(trial, obj.ConstraintsFcn(trial));
            end
            if obj.Group && obj.Multivariate
                obj.updateNativeHistory(study,trial);
            end
        end
    end

    methods (Access=private)
        function groups=groupSearchSpaces(obj,study)
            if obj.NativeHistoryValid
                groups=obj.GroupDecomposition.current(ExcludeSingle=true);
            else
                groups=obj.GroupDecomposition.calculate( ...
                    study,IncludePruned=true,ExcludeSingle=true);
            end
        end

        function goodNumbers=globalGoodTrialNumbers(obj,study,precomputed)
            if nargin<3 || isempty(fieldnames(precomputed))
                [~,objectives,numbers,pending]=obj.relativeObservations( ...
                    study,obj.emptySearchSpace());
            else
                objectives=precomputed.objectives;
                numbers=precomputed.trialNumbers;
                pending=precomputed.pending;
            end
            nGood=obj.goodTrialCount(sum(~pending));
            order=obj.rankObservations(objectives, ...
                study.Directions(1),study,numbers,pending);
            goodNumbers=numbers(order(1:nGood));
        end

        function warnIndependent(obj,study,trial,name)
            if ~obj.Multivariate || ~obj.WarnIndependentSampling || obj.Group
                return
            end
            trials=study.trialData();
            parameters=study.parameterData();
            eligible=ismember(trials.State,["COMPLETE","PRUNED"]);
            if sum(eligible)<obj.NStartupTrials
                return
            end
            previous=parameters.Name==string(name) & ...
                ismember(parameters.TrialNumber,trials.TrialNumber(eligible));
            if any(previous)
                warning("radia:optuna:TPEIndependentSampling", ...
                    "Parameter '%s' in trial %d is sampled independently because the dynamic search space is not grouped.", ...
                    name,trial.Number);
            end
        end

        function [good, bad] = splitObservations(obj, values, objectives, ...
                direction, study, trialNumbers, pending)
            % Match Optuna's _split_complete_trials: n_below may equal the
            % number of observations.  In that case the above density is
            % represented by its prior component only.
            finishedCount = sum(~pending);
            nGood = obj.goodTrialCount(finishedCount);
            order = obj.rankObservations( ...
                objectives, direction, study, trialNumbers, pending);
            isGood = false(numel(values), 1);
            isGood(order(1:nGood)) = true;
            % Keep chronological order so Optuna's history weights attach to
            % the same observations after the objective-based split.
            good = values(isGood);
            bad = values(~isGood);
        end

        function count = goodTrialCount(obj, finishedCount)
            if isempty(obj.GammaFcn)
                count = min(obj.MaxGoodTrials, ...
                    ceil(obj.Gamma * finishedCount));
                count = max(1, min(finishedCount, count));
                return
            end
            count = obj.GammaFcn(finishedCount);
            if ~(isnumeric(count) && isreal(count) && isscalar(count) && ...
                    isfinite(count) && count == floor(count) && ...
                    count >= 0 && count <= finishedCount)
                error("radia:optuna:TPEGamma", ...
                    "GammaFcn(%d) must return an integer from 0 through %d.", ...
                    finishedCount, finishedCount);
            end
            count = double(count);
        end

        function weights = observationWeights(obj, count)
            if isempty(obj.WeightsFcn)
                weights = zeros(0,1);
                return
            end
            weights = obj.WeightsFcn(count);
            if ~(isnumeric(weights) && isreal(weights))
                error("radia:optuna:TPEWeights", ...
                    "WeightsFcn(%d) must return a real numeric vector.", count);
            end
            weights = reshape(double(weights),[],1);
            if numel(weights) ~= count
                error("radia:optuna:TPEWeights", ...
                    "WeightsFcn(%d) returned %d weights; expected %d.", ...
                    count, numel(weights), count);
            end
            if any(~isfinite(weights)) || any(weights < 0) || ...
                    (count > 0 && sum(weights) <= 0)
                error("radia:optuna:TPEWeights", ...
                    "WeightsFcn(%d) must return finite nonnegative weights " + ...
                    "with positive total mass.", count);
            end
        end

        function [x, y, trialNumbers, pending] = ...
                numericObservations(obj, study, name)
            p = study.parameterData();
            t = study.trialData();
            rows = p.Name == string(name) & isfinite(p.ValueNumeric);
            indices = find(rows);
            trialRows=p.TrialNumber(indices)+1;
            known=trialRows>=1 & trialRows<=numel(t.TrialNumber);
            known(known)=t.TrialNumber(trialRows(known))== ...
                p.TrialNumber(indices(known));
            indices=indices(known);
            trialRows=trialRows(known);
            states=t.State(trialRows);
            complete=states=="COMPLETE" & isfinite(t.Value(trialRows));
            running=states=="RUNNING" & obj.ConstantLiar;
            keep=complete | running;
            indices=indices(keep);
            trialRows=trialRows(keep);
            pending=running(keep);
            x=p.ValueNumeric(indices);
            y=t.Value(trialRows);
            if any(pending)
                y(pending)=obj.liarObjective(study);
            end
            trialNumbers=p.TrialNumber(indices);
        end

        function [tokens, y, trialNumbers, pending] = ...
                categoricalObservations(obj, study, name)
            p = study.parameterData();
            t = study.trialData();
            rows = p.Name == string(name) & p.Kind == "categorical";
            indices = find(rows);
            trialRows=p.TrialNumber(indices)+1;
            known=trialRows>=1 & trialRows<=numel(t.TrialNumber);
            known(known)=t.TrialNumber(trialRows(known))== ...
                p.TrialNumber(indices(known));
            indices=indices(known);
            trialRows=trialRows(known);
            states=t.State(trialRows);
            complete=states=="COMPLETE" & isfinite(t.Value(trialRows));
            running=states=="RUNNING" & obj.ConstantLiar;
            keep=complete | running;
            indices=indices(keep);
            trialRows=trialRows(keep);
            pending=running(keep);
            tokens=p.ValueText(indices);
            numeric=isfinite(p.ValueNumeric(indices));
            for index=reshape(find(numeric),1,[])
                tokens(index)=obj.token(p.ValueNumeric(indices(index)));
            end
            y=t.Value(trialRows);
            if any(pending)
                y(pending)=obj.liarObjective(study);
            end
            trialNumbers=p.TrialNumber(indices);
        end

        function value = uniform(obj, low, high, logScale)
            u = rand(obj.Stream, 1, 1);
            if logScale
                value = exp(log(low) + u * (log(high) - log(low)));
            else
                value = low + u * (high - low);
            end
        end

        function indices = sampleComponents(obj, weights, count)
            cumulative = cumsum(weights(:));
            cumulative(end) = 1;
            u = rand(obj.Stream, count, 1);
            indices = 1 + sum(cumulative.' < u, 2);
        end

        function searchSpace = intersectionSearchSpace(~, study)
            searchSpace = ...
                radia.optuna.internal.IntersectionSearchSpace.calculate( ...
                study, IncludePruned=true);
        end

        function searchSpace = emptySearchSpace(~)
            template = struct( ...
                "name", "", ...
                "distribution", ...
                    radia.optuna.internal.DistributionCodec.float( ...
                    0, 1, false, NaN));
            searchSpace = reshape(template([]), 0, 1);
        end

        function values = sampleRelativeSpace(obj, study, searchSpace, options)
            arguments
                obj
                study
                searchSpace
                options.EnforceStartup (1,1) logical = true
                options.GlobalSplit (1,1) logical = false
                options.GlobalGoodTrialNumbers double = NaN
                options.Precomputed struct = struct()
            end
            if isempty(fieldnames(options.Precomputed))
                [observations, objectives, trialNumbers, pending] = ...
                    obj.relativeObservations(study, searchSpace);
            else
                observations=options.Precomputed.observations;
                objectives=options.Precomputed.objectives;
                trialNumbers=options.Precomputed.trialNumbers;
                pending=options.Precomputed.pending;
            end
            finishedCount = sum(~pending);
            if finishedCount == 0 || ...
                    (options.EnforceStartup && ...
                    finishedCount < obj.NStartupTrials) || ...
                    isempty(objectives)
                values = cell(1, numel(searchSpace));
                for index = 1:numel(searchSpace)
                    values{index} = obj.randomRelativeValue( ...
                        searchSpace(index).distribution);
                end
                return
            end

            if options.GlobalSplit
                goodNumbers=options.GlobalGoodTrialNumbers;
                if isscalar(goodNumbers) && isnan(goodNumbers)
                    goodNumbers=obj.globalGoodTrialNumbers(study);
                end
                maximum=max([trialNumbers;goodNumbers(:);0]);
                goodLookup=false(maximum+1,1);
                goodLookup(goodNumbers+1)=true;
                isGood=goodLookup(trialNumbers+1);
                good=observations(isGood,:);
                bad=observations(~isGood,:);
            else
                [good, bad] = obj.splitJointObservations( ...
                    observations, objectives, study.Directions(1), study, ...
                    trialNumbers, pending);
            end
            dimension = numel(searchSpace);
            belowWeights = obj.observationWeights(size(good,1));
            aboveWeights = obj.observationWeights(size(bad,1));
            nativeHandle=obj.Stream.nativeHandle();
            nativeFastPath=nativeHandle~=0 && ...
                isempty(obj.CategoricalDistanceFcn) && ...
                radia.optuna.internal.NativeKernels.has( ...
                "optuna.tpe.best_joint_observations");
            if nativeFastPath
                categorical=false(1,dimension);
                lows=zeros(1,dimension);
                highs=ones(1,dimension);
                logScale=false(1,dimension);
                steps=NaN(1,dimension);
                choiceCounts=zeros(1,dimension);
                for index=1:dimension
                    distribution=searchSpace(index).distribution;
                    categorical(index)=distribution.kind=="categorical";
                    if categorical(index)
                        choiceCounts(index)=numel(distribution.choices);
                    else
                        lows(index)=distribution.low;
                        highs(index)=distribution.high;
                        logScale(index)=distribution.log;
                        steps(index)=distribution.step;
                    end
                end
                bestCandidates=radia.optuna.internal.NativeKernels.call( ...
                    "optuna.tpe.best_joint_observations",nativeHandle, ...
                    obj.NumberOfEIChoices,categorical,lows,highs,logScale, ...
                    steps,choiceCounts,good,bad,obj.PriorWeight, ...
                    obj.ConsiderMagicClip,belowWeights,aboveWeights);
                values=cell(1,dimension);
                for index=1:dimension
                    distribution=searchSpace(index).distribution;
                    if categorical(index)
                        values{index}=radia.optuna.internal. ...
                            DistributionCodec.choiceAt( ...
                            distribution.choices,bestCandidates(index));
                    else
                        values{index}=bestCandidates(index);
                    end
                end
                return
            end
            below = cell(1, dimension);
            above = cell(1, dimension);
            for index = 1:dimension
                distribution = searchSpace(index).distribution;
                if distribution.kind == "categorical"
                    choiceCount = numel(distribution.choices);
                    distanceFcn=radia.optuna.internal.CategoricalDistance. ...
                        get(obj.CategoricalDistanceFcn, ...
                        searchSpace(index).name);
                    below{index} = ...
                        radia.optuna.internal.ParzenEstimator.categorical( ...
                        good(:,index), choiceCount, ...
                        PriorWeight=obj.PriorWeight, ...
                        ObservationWeights=belowWeights, ...
                        DistanceFcn=distanceFcn, ...
                        Choices=distribution.choices);
                    above{index} = ...
                        radia.optuna.internal.ParzenEstimator.categorical( ...
                        bad(:,index), choiceCount, ...
                        PriorWeight=obj.PriorWeight, ...
                        ObservationWeights=aboveWeights, ...
                        DistanceFcn=distanceFcn, ...
                        Choices=distribution.choices);
                else
                    estimatorOptions = { ...
                        "Log", distribution.log, ...
                        "Step", distribution.step, ...
                        "PriorWeight", obj.PriorWeight, ...
                        "ConsiderMagicClip", obj.ConsiderMagicClip, ...
                        "ConsiderEndpoints", obj.ConsiderEndpoints, ...
                        "MultivariateDimension", dimension, ...
                        "ObservationWeights", belowWeights};
                    below{index} = ...
                        radia.optuna.internal.ParzenEstimator.numerical( ...
                        good(:,index), distribution.low, distribution.high, ...
                        estimatorOptions{:});
                    above{index} = ...
                        radia.optuna.internal.ParzenEstimator.numerical( ...
                        bad(:,index), distribution.low, distribution.high, ...
                        "Log", distribution.log, ...
                        "Step", distribution.step, ...
                        "PriorWeight", obj.PriorWeight, ...
                        "ConsiderMagicClip", obj.ConsiderMagicClip, ...
                        "ConsiderEndpoints", obj.ConsiderEndpoints, ...
                        "MultivariateDimension", dimension, ...
                        "ObservationWeights", aboveWeights);
                end
            end

            nativeHandle=obj.Stream.nativeHandle();
            if nativeHandle~=0 && radia.optuna.internal.NativeKernels.has( ...
                    "optuna.tpe.best_joint")
                categorical=false(1,dimension);
                for index=1:dimension
                    categorical(index)= ...
                        searchSpace(index).distribution.kind=="categorical";
                end
                bestCandidates=radia.optuna.internal.NativeKernels.call( ...
                    "optuna.tpe.best_joint",nativeHandle, ...
                    obj.NumberOfEIChoices,categorical,below,above);
                values=cell(1,dimension);
                for index=1:dimension
                    distribution=searchSpace(index).distribution;
                    if categorical(index)
                        values{index}=radia.optuna.internal. ...
                            DistributionCodec.choiceAt( ...
                            distribution.choices,bestCandidates(index));
                    else
                        values{index}=bestCandidates(index);
                    end
                end
                return
            end

            count = obj.NumberOfEIChoices;
            components = obj.sampleComponents(below{1}.weights, count);
            candidates = zeros(count, dimension);
            % Optuna's _MixtureOfProductDistribution samples every
            % categorical dimension first, then draws all numerical
            % dimensions in search-space order through one vectorized
            % truncnorm call. Preserve that RandomState consumption order.
            for index = 1:dimension
                distribution = searchSpace(index).distribution;
                if distribution.kind == "categorical"
                    candidates(:,index) = ...
                        radia.optuna.internal.ParzenEstimator. ...
                        sampleCategoricalComponents( ...
                        below{index}, obj.Stream, components);
                end
            end
            for index = 1:dimension
                distribution = searchSpace(index).distribution;
                if distribution.kind ~= "categorical"
                    candidates(:,index) = ...
                        radia.optuna.internal.ParzenEstimator. ...
                        sampleNumericalComponents( ...
                        below{index}, obj.Stream, components);
                end
            end
            % Optuna models a mixture of product distributions: one
            % component index is shared by every dimension. Therefore the
            % per-component log densities must be added before the mixture
            % log-sum-exp. Summing already-mixed marginal densities changes
            % the acquisition ranking and breaks seeded parity.
            belowComponentLogPdf = zeros(count, numel(below{1}.weights));
            aboveComponentLogPdf = zeros(count, numel(above{1}.weights));
            for index = 1:dimension
                distribution = searchSpace(index).distribution;
                if distribution.kind == "categorical"
                    belowPart = ...
                        radia.optuna.internal.ParzenEstimator. ...
                        componentLogPdfCategorical( ...
                        below{index}, candidates(:,index));
                    abovePart = ...
                        radia.optuna.internal.ParzenEstimator. ...
                        componentLogPdfCategorical( ...
                        above{index}, candidates(:,index));
                else
                    belowPart = ...
                        radia.optuna.internal.ParzenEstimator. ...
                        componentLogPdfNumerical( ...
                        below{index}, candidates(:,index));
                    abovePart = ...
                        radia.optuna.internal.ParzenEstimator. ...
                        componentLogPdfNumerical( ...
                        above{index}, candidates(:,index));
                end
                belowComponentLogPdf = belowComponentLogPdf + belowPart;
                aboveComponentLogPdf = aboveComponentLogPdf + abovePart;
            end
            acquisition = obj.mixtureLogPdf( ...
                belowComponentLogPdf, below{1}.weights) - ...
                obj.mixtureLogPdf( ...
                aboveComponentLogPdf, above{1}.weights);
            [~, best] = max(acquisition);
            values = cell(1, dimension);
            for index = 1:dimension
                distribution = searchSpace(index).distribution;
                if distribution.kind == "categorical"
                    values{index} = ...
                        radia.optuna.internal.DistributionCodec.choiceAt( ...
                        distribution.choices, candidates(best,index));
                else
                    values{index} = candidates(best,index);
                end
            end
        end

        function values = mixtureLogPdf(~, componentLogPdf, weights)
            weighted = componentLogPdf + ...
                log(reshape(weights, 1, []));
            maximum = max(weighted, [], 2);
            maximum(isinf(maximum) & maximum < 0) = 0;
            values = log(sum(exp(weighted - maximum), 2)) + maximum;
        end

        function [x, y, observationTrialNumbers, pending] = ...
                relativeObservations(obj, study, searchSpace, keepIncomplete)
            if nargin<4
                keepIncomplete=false;
            end
            trials=study.trialData();
            states = trials.State;
            usable = states == "COMPLETE" | states == "PRUNED" | ...
                (states == "RUNNING" & obj.ConstantLiar);
            trialRows=find(usable);
            trialNumbers = trials.TrialNumber(trialRows);
            finished = trials.Value(usable);
            finiteFinished = finished(isfinite(finished));
            if isempty(finiteFinished)
                liar = 0;
            elseif study.Directions(1) == "minimize"
                liar = max(finiteFinished);
            else
                liar = min(finiteFinished);
            end
            count=numel(trialNumbers);
            x=NaN(count,numel(searchSpace));
            valid=true(count,1);
            objectiveValid=true(count,1);
            parameters=study.parameterData();
            maximumTrial=max([trials.TrialNumber;0]);
            positionByNumber=zeros(maximumTrial+1,1);
            positionByNumber(trialNumbers+1)=1:count;
            for index=1:numel(searchSpace)
                distribution=searchSpace(index).distribution;
                candidateRows=find(parameters.Name==searchSpace(index).name);
                candidateNumbers=parameters.TrialNumber(candidateRows);
                inRange=candidateNumbers>=0 & candidateNumbers<=maximumTrial;
                positions=zeros(size(candidateNumbers));
                positions(inRange)=positionByNumber(candidateNumbers(inRange)+1);
                found=positions>0;
                candidateRows=candidateRows(found);
                positions=positions(found);
                rowForTrial=zeros(count,1);
                rowForTrial(positions)=candidateRows;
                present=rowForTrial>0;
                presentPositions=find(present);
                rows=rowForTrial(present);
                compatible=false(count,1);
                expectedEncoding=obj.distributionEncoding( ...
                    searchSpace(index).name,distribution);
                exact=parameters.Kind(rows)==distribution.kind & ...
                    parameters.Distribution(rows)==expectedEncoding;
                compatible(presentPositions(exact))=true;
                for fallbackIndex=reshape(find(~exact),1,[])
                    row=rows(fallbackIndex);
                    if parameters.Kind(row)~=distribution.kind
                        continue
                    end
                    stored=radia.optuna.internal.DistributionCodec.decode( ...
                        parameters.Kind(row),parameters.Distribution(row));
                    compatible(presentPositions(fallbackIndex))= ...
                        radia.optuna.internal.DistributionCodec.equivalent( ...
                        distribution,stored);
                end
                if distribution.kind=="categorical"
                    tokens=parameters.ValueText(rows);
                    numeric=isfinite(parameters.ValueNumeric(rows));
                    for numericIndex=reshape(find(numeric),1,[])
                        tokens(numericIndex)=obj.token( ...
                            parameters.ValueNumeric(rows(numericIndex)));
                    end
                    [validValues,indices]=ismember(tokens, ...
                        radia.optuna.internal.DistributionCodec. ...
                        choiceTokens(distribution.choices));
                else
                    indices=parameters.ValueNumeric(rows);
                    validValues=isfinite(indices) & ...
                        indices>=distribution.low & ...
                        indices<=distribution.high & ...
                        (~distribution.log | indices>0);
                end
                internalValues=NaN(count,1);
                internalValues(presentPositions)=indices;
                mappedValid=false(count,1);
                mappedValid(presentPositions)=validValues;
                columnValid=present & compatible & mappedValid;
                internalValues(~columnValid)=NaN;
                valid=valid & columnValid;
                x(:,index)=internalValues;
            end
            y=trials.Value(trialRows);
            selectedStates=trials.State(trialRows);
            pruned=find(selectedStates=="PRUNED");
            for index=reshape(pruned,1,[])
                intermediate=study.IntermediateTable.TrialNumber== ...
                    trialNumbers(index);
                if any(intermediate)
                    y(index)=study.IntermediateTable.Value( ...
                        find(intermediate,1,"last"));
                else
                    valid(index)=false;
                    objectiveValid(index)=false;
                end
            end
            pending=selectedStates=="RUNNING";
            y(pending)=liar;
            if keepIncomplete
                y(~objectiveValid)=NaN;
                observationTrialNumbers=trialNumbers;
            else
                x=x(valid,:);
                y=y(valid);
                observationTrialNumbers=trialNumbers(valid);
                pending=pending(valid);
            end
        end

        function value = randomRelativeValue(obj, distribution)
            if distribution.kind == "categorical"
                [~,index] = max(rand(obj.Stream, ...
                    numel(distribution.choices),1));
                value = ...
                    radia.optuna.internal.DistributionCodec.choiceAt( ...
                    distribution.choices, index);
                return
            end
            value = obj.uniform( ...
                distribution.low, distribution.high, distribution.log);
            value = obj.quantize(value, distribution.low, ...
                distribution.high, distribution.step);
            if distribution.kind == "integer"
                value = round(value);
            end
        end

        function [good, bad] = splitJointObservations(obj, values, ...
                objectives, direction, study, trialNumbers, pending)
            finishedCount = sum(~pending);
            nGood = obj.goodTrialCount(finishedCount);
            order = obj.rankObservations( ...
                objectives, direction, study, trialNumbers, pending);
            mask = false(size(values,1),1);
            mask(order(1:nGood)) = true;
            good = values(mask,:);
            bad = values(~mask,:);
        end

        function order = rankObservations(obj, objectives, direction, ...
                study, trialNumbers, pending)
            % Feasible completed/pruned trials rank before infeasible ones;
            % infeasible trials rank by total positive violation. Missing or
            % nonfinite constraint data is deliberately worst (Inf). Pending
            % constant-liar observations can shape only the above density.
            violations = obj.constraintViolations(study, trialNumbers);
            objectiveRank = reshape(double(objectives), [], 1);
            if direction == "maximize"
                objectiveRank = -objectiveRank;
            end
            sequence = (1:numel(objectiveRank))';
            rankKeys = [double(reshape(pending, [], 1)), ...
                double(violations > 0), violations, objectiveRank, sequence];
            [~, order] = sortrows(rankKeys, 1:size(rankKeys, 2));
        end

        function violations = constraintViolations(obj, study, trialNumbers)
            trialNumbers = reshape(double(trialNumbers), [], 1);
            violations = zeros(size(trialNumbers));
            constraintsEnabled = ~isempty(obj.ConstraintsFcn) || ...
                study.hasConstraintRecords();
            if ~constraintsEnabled
                return
            end
            for index = 1:numel(trialNumbers)
                [present,constraintValues] = ...
                    study.constraintRecord(trialNumbers(index));
                if ~present
                    violations(index) = Inf;
                    continue
                end
                violations(index) = sum(max(constraintValues, 0));
            end
        end

        function value = liarObjective(~, study)
            trials=study.trialData();
            complete = trials.State == "COMPLETE" & ...
                isfinite(trials.Value);
            finished = trials.Value(complete);
            if isempty(finished)
                value = 0;
            elseif study.Directions(1) == "minimize"
                value = max(finished);
            else
                value = min(finished);
            end
        end

        function encoding=distributionEncoding(obj,name,distribution)
            candidates=find(obj.EncodingCacheNames==string(name));
            for index=reshape(candidates,1,[])
                if radia.optuna.internal.DistributionCodec.equivalent( ...
                        obj.EncodingCacheDistributions{index},distribution)
                    encoding=obj.EncodingCacheValues(index);
                    return
                end
            end
            encoding=radia.optuna.internal.DistributionCodec.encode(distribution);
            row=numel(obj.EncodingCacheNames)+1;
            obj.EncodingCacheNames(row,1)=string(name);
            obj.EncodingCacheDistributions{row,1}=distribution;
            obj.EncodingCacheValues(row,1)=encoding;
        end

        function result=canUseNativeGroupedHistory(obj,study)
            result=obj.NativeHistoryValid && obj.Group && obj.Multivariate && ...
                isscalar(study.Directions) && isempty(obj.GammaFcn) && ...
                isempty(obj.WeightsFcn) && isempty(obj.ConstraintsFcn) && ...
                isempty(obj.CategoricalDistanceFcn) && ~obj.ConstantLiar && ...
                strlength(study.StoragePath)==0 && ...
                obj.Stream.nativeHandle()~=0 && ...
                radia.optuna.internal.NativeKernels.has( ...
                    "optuna.tpe.best_grouped_history");
            if ~result
                return
            end
            trials=study.trialData();
            result=sum(trials.State=="COMPLETE")== ...
                obj.NativeHistoryCompleteCount && ...
                ~any(trials.State=="PRUNED") && ...
                sum(trials.State=="RUNNING")==1;
            if ~result
                obj.NativeHistoryValid=false;
            end
        end

        function updateNativeHistory(obj,~,trial)
            if ~obj.NativeHistoryValid || ~obj.Group || ~obj.Multivariate
                return
            end
            if trial.State=="PRUNED"
                obj.NativeHistoryValid=false;
                return
            end
            if trial.State~="COMPLETE"
                return
            end
            if numel(trial.Values)~=1 || ~isfinite(trial.Value) || ...
                    obj.Stream.nativeHandle()==0 || ...
                    ~radia.optuna.internal.NativeKernels.has( ...
                        "optuna.tpe.history.append_complete")
                obj.NativeHistoryValid=false;
                return
            end
            [names,values,distributions]=trial.parameterRecords();
            distributionIds=zeros(numel(names),1,"int32");
            for index=1:numel(names)
                distributionIds(index)=obj.historyDistributionId( ...
                    names(index),distributions{index});
            end
            radia.optuna.internal.NativeKernels.call( ...
                "optuna.tpe.history.append_complete", ...
                obj.Stream.nativeHandle(),trial.Number,trial.Value, ...
                distributionIds,values);
            obj.NativeHistoryCompleteCount= ...
                obj.NativeHistoryCompleteCount+1;
            obj.GroupDecomposition.update(names,distributions);
        end

        function identifier=historyDistributionId(obj,name,distribution)
            candidates=find(obj.HistoryDistributionNames==string(name));
            for index=reshape(candidates,1,[])
                if radia.optuna.internal.DistributionCodec.equivalent( ...
                        obj.HistoryDistributions{index},distribution)
                    identifier=int32(index);
                    return
                end
            end
            row=numel(obj.HistoryDistributionNames)+1;
            obj.HistoryDistributionNames(row,1)=string(name);
            obj.HistoryDistributions{row,1}=distribution;
            identifier=int32(row);
        end

        function used=sampleNativeGroupedHistory(obj,trial,groups)
            used=false;
            if isempty(groups)
                return
            end
            metadata=obj.nativeGroupedMetadata(groups);
            best=radia.optuna.internal.NativeKernels.call( ...
                "optuna.tpe.best_grouped_history", ...
                obj.Stream.nativeHandle(),obj.NumberOfEIChoices, ...
                metadata.groupOffsets,metadata.distributionIds, ...
                metadata.categorical,metadata.lows,metadata.highs, ...
                metadata.logScale,metadata.steps,metadata.choiceCounts, ...
                trial.Study.Directions(1)=="minimize",obj.Gamma, ...
                obj.MaxGoodTrials,obj.PriorWeight,obj.ConsiderMagicClip);
            position=0;
            for groupIndex=1:numel(groups)
                searchSpace=groups{groupIndex};
                values=cell(1,numel(searchSpace));
                for index=1:numel(searchSpace)
                    position=position+1;
                    distribution=searchSpace(index).distribution;
                    if distribution.kind=="categorical"
                        values{index}=radia.optuna.internal. ...
                            DistributionCodec.choiceAt( ...
                                distribution.choices,best(position));
                    else
                        values{index}=best(position);
                    end
                end
                trial.setRelativeParameters(searchSpace,values,"");
            end
            used=true;
        end

        function metadata=nativeGroupedMetadata(obj,groups)
            revision=obj.GroupDecomposition.revision();
            if obj.NativeGroupRevision==revision && ...
                    ~isempty(fieldnames(obj.NativeGroupMetadata))
                metadata=obj.NativeGroupMetadata;
                return
            end
            groupCount=numel(groups);
            dimensions=reshape(cellfun(@numel,groups),1,[]);
            total=sum(dimensions);
            metadata=struct( ...
                "groupOffsets",int32([0,cumsum(dimensions)]), ...
                "distributionIds",zeros(1,total,"int32"), ...
                "categorical",false(1,total), ...
                "lows",zeros(1,total), ...
                "highs",ones(1,total), ...
                "logScale",false(1,total), ...
                "steps",NaN(1,total), ...
                "choiceCounts",zeros(1,total,"int32"));
            position=0;
            for groupIndex=1:groupCount
                searchSpace=groups{groupIndex};
                for index=1:numel(searchSpace)
                    position=position+1;
                    distribution=searchSpace(index).distribution;
                    metadata.distributionIds(position)= ...
                        obj.historyDistributionId( ...
                            searchSpace(index).name,distribution);
                    metadata.categorical(position)= ...
                        distribution.kind=="categorical";
                    if metadata.categorical(position)
                        metadata.choiceCounts(position)= ...
                            numel(distribution.choices);
                    else
                        metadata.lows(position)=distribution.low;
                        metadata.highs(position)=distribution.high;
                        metadata.logScale(position)=distribution.log;
                        metadata.steps(position)=distribution.step;
                    end
                end
            end
            obj.NativeGroupRevision=revision;
            obj.NativeGroupMetadata=metadata;
        end

        function value = randomNumerical(obj, low, high, logScale, step)
            if isfinite(step) && ~logScale
                count = floor((high - low) / step + 1e-12) + 1;
                index = floor(rand(obj.Stream, 1, 1) * count);
                value = low + index * step;
                value = min(value, high);
                return
            end
            value = obj.uniform(low, high, logScale);
            value = obj.quantize(value, low, high, step);
        end

        function attach(obj, study)
            changed = isempty(obj.AttachedStudy) || ...
                ~isequal(obj.AttachedStudy, study);
            if changed
                obj.AttachedStudy = study;
                obj.Stream = ...
                    radia.optuna.internal.NumpyRandomState(obj.Seed);
                obj.GroupDecomposition= ...
                    radia.optuna.internal.GroupDecomposedSearchSpace();
                obj.NativeHistoryCompleteCount=0;
                obj.HistoryDistributionNames=strings(0,1);
                obj.HistoryDistributions=cell(0,1);
                obj.NativeGroupRevision=-1;
                obj.NativeGroupMetadata=struct();
                trials=study.trialData();
                terminal=trials.State=="COMPLETE" | ...
                    trials.State=="PRUNED";
                obj.NativeHistoryValid=strlength(study.StoragePath)==0 && ...
                    ~any(terminal) && sum(trials.State=="RUNNING")<=1;
                obj.Restored = false;
            end
            if obj.Restored
                return
            end
            state = study.samplerState(obj.SamplerName, obj.StateSchema);
            if ~isempty(state)
                obj.restoreState(state);
            end
            obj.Restored = true;
        end

        function restoreState(obj, state)
            required = ["schema", "seed", "random_state"];
            if ~isstruct(state) || ~isscalar(state) || ...
                    any(~isfield(state, required)) || ...
                    string(state.schema) ~= obj.StateSchema
                error("radia:optuna:TPEState", ...
                    "Stored TPE sampler state is invalid or unsupported.");
            end
            if ~isnumeric(state.seed) || ~isscalar(state.seed) || ...
                    ~isfinite(double(state.seed)) || ...
                    double(state.seed) ~= obj.Seed
                error("radia:optuna:TPEStateSeed", ...
                    "Stored TPE sampler seed (%g) does not match the " + ...
                    "configured seed (%g).", double(state.seed), obj.Seed);
            end
            try
                obj.Stream.State = state.random_state;
            catch exception
                error("radia:optuna:TPEState", ...
                    "Stored TPE random state is invalid: %s", ...
                    exception.message);
            end
        end

        function recordState(obj, study, trialNumber)
            % In-memory studies have nowhere to resume from. Avoid copying
            % 624 MT words into a table on every proposal when no storage
            % contract exists; persisted studies retain the exact snapshot.
            if strlength(study.StoragePath)==0
                return
            end
            state = struct( ...
                "schema", obj.StateSchema, ...
                "seed", obj.Seed, ...
                "random_state", obj.Stream.State);
            trials=study.trialData();
            generation = sum(trials.State == "COMPLETE");
            study.recordSamplerState(obj.SamplerName, obj.StateSchema, ...
                trialNumber, generation, state);
        end

        function value = quantize(~, value, low, high, step)
            if isfinite(step)
                value = low + round((value - low) / step) * step;
            end
            value = min(max(value, low), high);
        end

        function validateBounds(~, low, high, logScale, step)
            if ~(isfinite(low) && isfinite(high) && low <= high)
                error("radia:optuna:Bounds", ...
                    "Float bounds must satisfy low <= high.");
            end
            if logScale && low <= 0
                error("radia:optuna:LogBounds", ...
                    "Log-uniform bounds must be positive.");
            end
            if isfinite(step) && step <= 0
                error("radia:optuna:Step", "Step must be positive.");
            end
            if logScale && isfinite(step) && low - step / 2 <= 0
                error("radia:optuna:LogBounds", ...
                    "Expanded log-distribution support must be positive.");
            end
        end

        function tokens = choiceTokens(obj, choices)
            tokens = strings(numel(choices), 1);
            for index = 1:numel(choices)
                tokens(index) = obj.token(obj.choiceAt(choices, index));
            end
        end

        function token = token(~, value)
            if isnumeric(value) && isscalar(value)
                token = "numeric:" + string(value, "%.17g");
            else
                token = string(jsonencode(value));
            end
        end

        function value = choiceAt(~, choices, index)
            if iscell(choices)
                value = choices{index};
            else
                value = choices(index);
            end
        end
    end
end
