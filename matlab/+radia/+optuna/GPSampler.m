classdef GPSampler < handle
    %GPSAMPLER Matérn-5/2 ARD Bayesian optimization for CAE objectives.
    %   Supports mixed stable search spaces, single-objective expected
    %   improvement, multi-objective expected hypervolume improvement,
    %   soft constraints c<=0, and pending-trial repulsion. The hot loop is
    %   MATLAB-native and does not require Statistics or Optimization Toolbox.

    properties (SetAccess=private)
        Seed (1,1) double = 0
        NStartupTrials (1,1) double = 10
        DeterministicObjective (1,1) logical = false
        CandidateCount (1,1) double = 2048
        LocalSearchCount (1,1) double = 10
        MonteCarloSamples (1,1) double = 128
        ConstraintsFcn = []
        Stream
    end

    properties (Access=private)
        AttachedStudy = []
        Restored (1,1) logical = false
        ObjectiveTheta cell = cell(0,1)
        ConstraintTheta cell = cell(0,1)
    end

    properties (Constant, Access=private)
        StateSchema = "radia.optuna.gp-sampler-state.v1"
        SamplerName = "gp"
    end

    methods
        function obj=GPSampler(options)
            arguments
                options.Seed (1,1) double = 0
                options.NStartupTrials (1,1) double ...
                    {mustBeInteger,mustBeNonnegative} = 10
                options.DeterministicObjective (1,1) logical = false
                options.CandidateCount (1,1) double ...
                    {mustBeInteger,mustBePositive} = 2048
                options.LocalSearchCount (1,1) double ...
                    {mustBeInteger,mustBeNonnegative} = 10
                options.MonteCarloSamples (1,1) double ...
                    {mustBeInteger,mustBePositive} = 128
                options.ConstraintsFcn = []
            end
            if options.CandidateCount<16
                error("radia:optuna:GPCandidates", ...
                    "CandidateCount must be at least 16.");
            end
            if ~isempty(options.ConstraintsFcn) && ...
                    ~isa(options.ConstraintsFcn,"function_handle")
                error("radia:optuna:GPConstraints", ...
                    "ConstraintsFcn must be a function handle.");
            end
            obj.Seed=double(options.Seed);
            obj.NStartupTrials=double(options.NStartupTrials);
            obj.DeterministicObjective=options.DeterministicObjective;
            obj.CandidateCount=double(options.CandidateCount);
            obj.LocalSearchCount=double(options.LocalSearchCount);
            obj.MonteCarloSamples=double(options.MonteCarloSamples);
            obj.ConstraintsFcn=options.ConstraintsFcn;
            obj.Stream=RandStream("mt19937ar","Seed",obj.Seed);
        end

        function searchSpace=inferRelativeSearchSpace(~,study,trial) %#ok<INUSD>
            searchSpace=radia.optuna.internal.IntersectionSearchSpace. ...
                calculate(study,IncludePruned=false);
        end

        function searchSpace=infer_relative_search_space(obj,study,trial)
            if nargin<3, trial=[]; end
            searchSpace=obj.inferRelativeSearchSpace(study,trial);
        end

        function beforeTrial(obj,study,trial)
            obj.attach(study);
            completed=sum(study.TrialTable.State=="COMPLETE");
            if completed<obj.NStartupTrials
                trial.setSystemAttr("gp_sampling_mode","startup_random");
                obj.recordState(study,trial.Number);
                return
            end
            searchSpace=obj.inferRelativeSearchSpace(study,trial);
            if isempty(searchSpace)
                trial.setSystemAttr("gp_sampling_mode", ...
                    "independent_dynamic_space");
                obj.recordState(study,trial.Number);
                return
            end
            [observations,objectives,numbers,pending,constraints, ...
                constraintPresent]=obj.observations(study,searchSpace,trial.Number);
            finished=~pending;
            if sum(finished)<obj.NStartupTrials || ...
                    size(objectives,1)<obj.NStartupTrials
                trial.setSystemAttr("gp_sampling_mode","startup_random");
                obj.recordState(study,trial.Number);
                return
            end
            values=obj.sampleRelative(study,searchSpace,observations, ...
                objectives,numbers,pending,constraints,constraintPresent);
            trial.setRelativeParameters(searchSpace,values,"gp");
            if isscalar(study.Directions)
                acquisition="expected_improvement";
            else
                acquisition="expected_hypervolume_improvement";
            end
            trial.setSystemAttr("gp_sampling_mode","matern52_ard");
            trial.setSystemAttr("gp_acquisition",acquisition);
            trial.setSystemAttr("gp_pending_count",sum(pending));
            obj.recordState(study,trial.Number);
        end

        function value=sampleFloat(obj,study,trial,name,low,high,options) %#ok<INUSD>
            obj.attach(study);
            distribution=radia.optuna.internal.DistributionCodec.float( ...
                low,high,options.Log,options.Step);
            value=obj.randomValue(distribution);
            obj.recordState(study,trial.Number);
        end

        function value=sampleInteger(obj,study,trial,name,low,high) %#ok<INUSD>
            obj.attach(study);
            distribution=radia.optuna.internal.DistributionCodec.integer( ...
                low,high,false,1);
            value=obj.randomValue(distribution);
            obj.recordState(study,trial.Number);
        end

        function value=sampleCategorical(obj,study,trial,name,choices) %#ok<INUSD>
            obj.attach(study);
            distribution=radia.optuna.internal.DistributionCodec. ...
                categorical(choices);
            value=obj.randomValue(distribution);
            obj.recordState(study,trial.Number);
        end

        function values=sampleJoint(obj,study,trial,names,lows,highs,options)
            values=zeros(1,numel(names));
            for index=1:numel(names)
                distribution=radia.optuna.internal.DistributionCodec.float( ...
                    lows(index),highs(index),options.Log(index),NaN);
                values(index)=obj.randomValue(distribution);
            end
            obj.recordState(study,trial.Number);
        end

        function afterTrial(obj,study,trial)
            if trial.State=="COMPLETE" && ~isempty(obj.ConstraintsFcn)
                study.recordConstraints(trial,obj.ConstraintsFcn(trial));
            end
            obj.recordState(study,trial.Number);
        end
    end

    methods (Access=private)
        function values=sampleRelative(obj,study,searchSpace,observations, ...
                objectives,trialNumbers,pending,constraints,constraintPresent)
            finished=~pending;
            x=observations(finished,:);
            y=objectives(finished,:);
            finishedNumbers=trialNumbers(finished);
            categorical=obj.categoricalMask(searchSpace);
            candidates=obj.candidatePool(x,categorical,searchSpace);
            objectiveModels=cell(1,size(y,2));
            if numel(obj.ObjectiveTheta)~=size(y,2)
                obj.ObjectiveTheta=cell(size(y,2),1);
            end
            for objective=1:size(y,2)
                [objectiveModels{objective},obj.ObjectiveTheta{objective}]= ...
                    radia.optuna.internal.GaussianProcess.fit( ...
                    x,y(:,objective),categorical, ...
                    obj.DeterministicObjective,obj.ObjectiveTheta{objective});
            end
            [probability,constraintModels]=obj.probabilityFeasible( ...
                candidates,x,constraints(finished,:), ...
                constraintPresent(finished),categorical);
            obj.ConstraintTheta=cell(numel(constraintModels),1);
            for index=1:numel(constraintModels)
                obj.ConstraintTheta{index}=constraintModels{index}.theta;
            end
            feasible=obj.feasibleMask(study,finishedNumbers);
            if size(y,2)==1
                acquisition=obj.expectedImprovement( ...
                    objectiveModels{1},candidates,y(:,1), ...
                    study.Directions(1),feasible,probability);
                running=observations(pending,:);
                acquisition=acquisition.*obj.pendingPenalty( ...
                    candidates,running,categorical);
                [refined,refinedAcquisition]=obj.refineSingleAcquisition( ...
                    candidates,acquisition,objectiveModels{1},y(:,1), ...
                    study.Directions(1),feasible,constraintModels, ...
                    running,categorical,searchSpace);
                candidates=[candidates;refined];
                acquisition=[acquisition;refinedAcquisition];
            else
                acquisition=obj.expectedHypervolumeImprovement( ...
                    objectiveModels,candidates,y,study.Directions, ...
                    feasible,probability);
                running=observations(pending,:);
                acquisition=acquisition.*obj.pendingPenalty( ...
                    candidates,running,categorical);
            end
            if all(~isfinite(acquisition)) || all(acquisition<=0)
                uncertainty=zeros(size(candidates,1),1);
                for objective=1:numel(objectiveModels)
                    [~,standardDeviation]=radia.optuna.internal. ...
                        GaussianProcess.predict( ...
                        objectiveModels{objective},candidates);
                    uncertainty=uncertainty+standardDeviation;
                end
                acquisition=probability.*uncertainty.* ...
                    obj.pendingPenalty(candidates,running,categorical);
            end
            [~,best]=max(acquisition);
            values=obj.decodePoint(candidates(best,:),searchSpace);
        end

        function candidates=candidatePool( ...
                obj,observations,categorical,searchSpace)
            count=obj.CandidateCount;
            dimensions=size(observations,2);
            % Match Optuna's 2048-point QMC preliminary search where the
            % native Sobol table covers the requested dimensionality.
            if dimensions<=32
                qmcSeed=floor(rand(obj.Stream)*(2^31-1));
                qmc=radia.optuna.QMCSampler(QMCType="sobol", ...
                    Scramble=true,Seed=qmcSeed);
                candidates=qmc.unitPoints(dimensions,count);
            else
                % Explicit bounded fallback above the native direction
                % table: randomized Latin hypercube, not iid sampling.
                candidates=zeros(count,dimensions);
                for dimension=1:dimensions
                    order=randperm(obj.Stream,count)';
                    candidates(:,dimension)=((order-1)+ ...
                        rand(obj.Stream,count,1))/count;
                end
            end
            candidates=obj.quantizeCategorical( ...
                candidates,categorical,searchSpace);
            if ~isempty(observations)
                candidates=[candidates;observations];
            end
            candidates=unique(candidates,"rows","stable");
        end

        function [refined,acquisition]=refineSingleAcquisition( ...
                obj,candidates,initialAcquisition,model,values,direction, ...
                feasible,constraintModels,pending,categorical,searchSpace)
            refined=zeros(0,size(candidates,2));
            acquisition=zeros(0,1);
            continuous=find(~categorical);
            if obj.LocalSearchCount==0 || isempty(continuous)
                return
            end
            finite=initialAcquisition;
            finite(~isfinite(finite))=-Inf;
            [~,order]=sort(finite,"descend");
            count=min(obj.LocalSearchCount,numel(order));
            anchors=candidates(order(1:count),:);
            anchorAcquisition=initialAcquisition(order(1:count));
            steps=[0.2 0.08 0.032 0.0128 0.00512];
            for step=steps
                probes=zeros(count*(1+2*numel(continuous)),size(anchors,2));
                owners=zeros(size(probes,1),1);
                cursor=1;
                for anchorIndex=1:count
                    probes(cursor,:)=anchors(anchorIndex,:);
                    owners(cursor)=anchorIndex;
                    cursor=cursor+1;
                    for dimension=reshape(continuous,1,[])
                        lower=anchors(anchorIndex,:);
                        upper=anchors(anchorIndex,:);
                        lower(dimension)=max(0,lower(dimension)-step);
                        upper(dimension)=min(1,upper(dimension)+step);
                        probes(cursor,:)=lower;
                        owners(cursor)=anchorIndex;
                        probes(cursor+1,:)=upper;
                        owners(cursor+1)=anchorIndex;
                        cursor=cursor+2;
                    end
                end
                probes=obj.quantizeCategorical(probes,categorical,searchSpace);
                probabilities=obj.probabilityFromModels( ...
                    constraintModels,probes);
                probeAcquisition=obj.expectedImprovement( ...
                    model,probes,values,direction,feasible,probabilities).* ...
                    obj.pendingPenalty(probes,pending,categorical);
                for anchorIndex=1:count
                    local=find(owners==anchorIndex);
                    [bestValue,bestOffset]=max(probeAcquisition(local));
                    if bestValue>anchorAcquisition(anchorIndex)
                        anchors(anchorIndex,:)=probes(local(bestOffset),:);
                        anchorAcquisition(anchorIndex)=bestValue;
                    end
                end
            end
            refined=anchors;
            acquisition=anchorAcquisition;
        end

        function [probability,models]=probabilityFeasible( ...
                obj,candidates,x,constraints,present,categorical)
            probability=ones(size(candidates,1),1);
            models=cell(0,1);
            if isempty(constraints) || ~any(present)
                return
            end
            counts=sum(isfinite(constraints),2);
            expected=max(counts(present));
            if any(counts(present)~=expected)
                error("radia:optuna:ConstraintShape", ...
                    "Trials with different numbers of constraints cannot be compared.");
            end
            models=cell(expected,1);
            if numel(obj.ConstraintTheta)~=expected
                obj.ConstraintTheta=cell(expected,1);
            end
            for constraint=1:expected
                rows=present & isfinite(constraints(:,constraint));
                [models{constraint},theta]=radia.optuna.internal. ...
                    GaussianProcess.fit(x(rows,:),constraints(rows,constraint), ...
                    categorical,obj.DeterministicObjective, ...
                    obj.ConstraintTheta{constraint});
                models{constraint}.theta=theta;
                [meanValue,stdValue]=radia.optuna.internal. ...
                    GaussianProcess.predict(models{constraint},candidates);
                z=(0-meanValue)./max(stdValue,1e-12);
                probability=probability.*(0.5*erfc(-z/sqrt(2)));
            end
            probability=max(probability,realmin("double"));
        end

        function probability=probabilityFromModels(~,models,candidates)
            probability=ones(size(candidates,1),1);
            for constraint=1:numel(models)
                [meanValue,stdValue]=radia.optuna.internal. ...
                    GaussianProcess.predict(models{constraint},candidates);
                z=(0-meanValue)./max(stdValue,1e-12);
                probability=probability.*(0.5*erfc(-z/sqrt(2)));
            end
            probability=max(probability,realmin("double"));
        end

        function acquisition=expectedImprovement(~,model,candidates,values, ...
                direction,feasible,probability)
            [meanValue,stdValue]=radia.optuna.internal. ...
                GaussianProcess.predict(model,candidates);
            if string(direction)=="minimize"
                meanValue=-meanValue;
                signedValues=-values;
            else
                signedValues=values;
            end
            if any(feasible)
                best=max(signedValues(feasible));
                improvement=meanValue-best;
                z=improvement./max(stdValue,1e-12);
                expected=improvement.*(0.5*erfc(-z/sqrt(2)))+ ...
                    stdValue.*exp(-0.5*z.^2)/sqrt(2*pi);
                expected(stdValue<1e-12)=max(improvement(stdValue<1e-12),0);
            else
                expected=ones(size(meanValue));
            end
            acquisition=max(expected,0).*probability;
        end

        function acquisition=expectedHypervolumeImprovement(obj,models, ...
                candidates,values,directions,feasible,probability)
            if ~any(feasible)
                acquisition=probability;
                return
            end
            signs=ones(1,numel(directions));
            signs(string(directions)=="maximize")=-1;
            front=values(feasible,:).*signs;
            worst=max(front,[],1);
            span=max(front,[],1)-min(front,[],1);
            reference=worst+max(0.1*max(span,abs(worst)),1e-9);
            means=zeros(size(candidates,1),numel(models));
            deviations=zeros(size(means));
            for objective=1:numel(models)
                [means(:,objective),deviations(:,objective)]= ...
                    radia.optuna.internal.GaussianProcess.predict( ...
                    models{objective},candidates);
            end
            if numel(models)==2
                uniforms=min(max(rand(obj.Stream,obj.MonteCarloSamples, ...
                    numel(models)),eps),1-eps);
                normals=sqrt(2)*erfinv(2*uniforms-1);
                first=means(:,1)+deviations(:,1).*normals(:,1)';
                second=means(:,2)+deviations(:,2).*normals(:,2)';
                samples=[first(:),second(:)].*signs;
                improvement=radia.optuna.internal.ParetoSupport. ...
                    hypervolumeImprovement2D( ...
                    samples,front,reference);
                acquisition=mean(reshape(improvement, ...
                    size(candidates,1),obj.MonteCarloSamples),2).*probability;
                return
            end
            % For three or more objectives, integrate the probability of
            % dominating objective-space Sobol nodes.  This is the same
            % EHVI integral, evaluated in a vectorized bounded domain, and
            % avoids one recursive hypervolume solve per posterior draw.
            signedMeans=means.*signs;
            lower=min([front;signedMeans-4*deviations],[],1);
            span=max(reference-lower,1e-12);
            qmcSeed=floor(rand(obj.Stream)*(2^31-1));
            qmc=radia.optuna.QMCSampler(QMCType="sobol", ...
                Scramble=true,Seed=qmcSeed);
            integration=lower+qmc.unitPoints( ...
                numel(models),obj.MonteCarloSamples).*span;
            dominated=false(size(integration,1),1);
            for point=1:size(front,1)
                dominated=dominated | all( ...
                    integration>=front(point,:),2);
            end
            integration=integration(~dominated,:);
            if isempty(integration)
                acquisition=zeros(size(candidates,1),1);
            else
                dominanceProbability=ones( ...
                    size(candidates,1),size(integration,1));
                for objective=1:numel(models)
                    z=(integration(:,objective)'- ...
                        signedMeans(:,objective))./ ...
                        max(deviations(:,objective),1e-12);
                    dominanceProbability=dominanceProbability.* ...
                        (0.5*erfc(-z/sqrt(2)));
                end
                acquisition=prod(span)*mean(dominanceProbability,2);
            end
            acquisition=acquisition.*probability;
        end

        function penalty=pendingPenalty(~,candidates,pending,categorical)
            penalty=ones(size(candidates,1),1);
            if isempty(pending), return, end
            for candidate=1:size(candidates,1)
                difference=candidates(candidate,:)-pending;
                difference(:,categorical)= ...
                    difference(:,categorical)~=0;
                distance=min(sum(difference.^2,2));
                % Kriging-believer in Optuna excludes the immediate basin
                % around each pending point.  This bounded surrogate keeps
                % the same anti-duplicate behavior without refitting the GP.
                exclusionRadiusSquared=0.01;
                if distance<=exclusionRadiusSquared
                    penalty(candidate)=1e-12;
                else
                    penalty(candidate)=max(1-exp( ...
                        -(distance-exclusionRadiusSquared)/0.04),1e-12);
                end
            end
        end

        function feasible=feasibleMask(~,study,trialNumbers)
            feasible=true(numel(trialNumbers),1);
            if ~study.hasConstraintRecords(), return, end
            feasible(:)=false;
            expected=NaN;
            for index=1:numel(trialNumbers)
                [present,values]=study.constraintRecord(trialNumbers(index));
                if ~present, continue, end
                if isnan(expected), expected=numel(values); end
                if numel(values)~=expected
                    error("radia:optuna:ConstraintShape", ...
                        "Trials with different numbers of constraints cannot be compared.");
                end
                feasible(index)=all(values<=0);
            end
        end

        function [x,y,numbers,pending,constraints,present]= ...
                observations(obj,study,searchSpace,excludedNumber)
            states=study.TrialTable.State;
            rows=find(states=="COMPLETE" | states=="RUNNING");
            rows=rows(study.TrialTable.TrialNumber(rows)~=excludedNumber);
            numbers=study.TrialTable.TrialNumber(rows);
            pending=states(rows)=="RUNNING";
            x=zeros(numel(rows),numel(searchSpace));
            y=NaN(numel(rows),numel(study.Directions));
            constraintCells=cell(numel(rows),1);
            present=false(numel(rows),1);
            maximumConstraints=0;
            for index=1:numel(rows)
                x(index,:)=obj.encodeTrial(study,numbers(index),searchSpace);
                if ~pending(index)
                    for objective=1:numel(study.Directions)
                        mask=study.ObjectiveTable.TrialNumber==numbers(index) & ...
                            study.ObjectiveTable.ObjectiveIndex==objective;
                        if sum(mask)==1
                            y(index,objective)=study.ObjectiveTable.Value(mask);
                        end
                    end
                    [present(index),constraintCells{index}]= ...
                        study.constraintRecord(numbers(index));
                    maximumConstraints=max(maximumConstraints, ...
                        numel(constraintCells{index}));
                end
            end
            constraints=NaN(numel(rows),maximumConstraints);
            for index=1:numel(rows)
                if present(index) && ~isempty(constraintCells{index})
                    constraints(index,1:numel(constraintCells{index}))= ...
                        constraintCells{index};
                end
            end
            completeRows=~pending & all(isfinite(y),2);
            keep=pending | completeRows;
            x=x(keep,:); y=y(keep,:); numbers=numbers(keep);
            pending=pending(keep); constraints=constraints(keep,:);
            present=present(keep);
        end

        function encoded=encodeTrial(obj,study,trialNumber,searchSpace)
            encoded=zeros(1,numel(searchSpace));
            for dimension=1:numel(searchSpace)
                row=study.ParamTable.TrialNumber==trialNumber & ...
                    study.ParamTable.Name==searchSpace(dimension).name;
                if sum(row)~=1
                    error("radia:optuna:GPObservations", ...
                        "A GP intersection parameter is missing from a trial.");
                end
                if isfinite(study.ParamTable.ValueNumeric(row))
                    value=study.ParamTable.ValueNumeric(row);
                else
                    value=jsondecode(study.ParamTable.ValueText(row));
                end
                encoded(dimension)=obj.encodeValue( ...
                    value,searchSpace(dimension).distribution);
            end
        end

        function encoded=encodeValue(~,value,distribution)
            if distribution.kind=="categorical"
                tokens=radia.optuna.internal.DistributionCodec. ...
                    choiceTokens(distribution.choices);
                token=radia.optuna.internal.DistributionCodec.choiceToken(value);
                index=find(tokens==token,1);
                if isempty(index)
                    error("radia:optuna:GPObservations", ...
                        "Observed categorical value is outside its distribution.");
                end
                encoded=(index-1)/max(1,numel(tokens)-1);
            elseif distribution.log
                encoded=(log(double(value))-log(distribution.low))/ ...
                    (log(distribution.high)-log(distribution.low));
            else
                encoded=(double(value)-distribution.low)/ ...
                    (distribution.high-distribution.low);
            end
            encoded=min(max(encoded,0),1);
        end

        function values=decodePoint(~,point,searchSpace)
            values=cell(1,numel(searchSpace));
            for dimension=1:numel(searchSpace)
                distribution=searchSpace(dimension).distribution;
                unit=min(max(point(dimension),0),1);
                if distribution.kind=="categorical"
                    count=numel(distribution.choices);
                    index=min(max(round(unit*max(1,count-1))+1,1),count);
                    values{dimension}=radia.optuna.internal. ...
                        DistributionCodec.choiceAt(distribution.choices,index);
                else
                    if distribution.log
                        value=exp(log(distribution.low)+unit* ...
                            (log(distribution.high)-log(distribution.low)));
                    else
                        value=distribution.low+unit* ...
                            (distribution.high-distribution.low);
                    end
                    if isfinite(distribution.step)
                        value=distribution.low+round( ...
                            (value-distribution.low)/distribution.step)* ...
                            distribution.step;
                    end
                    value=min(max(value,distribution.low),distribution.high);
                    if distribution.kind=="integer", value=round(value); end
                    values{dimension}=value;
                end
            end
        end

        function mask=categoricalMask(~,searchSpace)
            mask=false(1,numel(searchSpace));
            for index=1:numel(searchSpace)
                mask(index)=searchSpace(index).distribution.kind=="categorical";
            end
        end

        function candidates=quantizeCategorical( ...
                ~,candidates,categorical,searchSpace)
            for dimension=find(categorical)
                levelCount=numel(searchSpace(dimension).distribution.choices);
                levels=linspace(0,1,levelCount)';
                indices=min(floor(candidates(:,dimension)*numel(levels))+1, ...
                    numel(levels));
                candidates(:,dimension)=levels(indices);
            end
        end

        function value=randomValue(obj,distribution)
            if distribution.kind=="categorical"
                index=randi(obj.Stream,numel(distribution.choices));
                value=radia.optuna.internal.DistributionCodec. ...
                    choiceAt(distribution.choices,index);
                return
            end
            unit=rand(obj.Stream);
            if distribution.log
                value=exp(log(distribution.low)+unit* ...
                    (log(distribution.high)-log(distribution.low)));
            else
                value=distribution.low+unit* ...
                    (distribution.high-distribution.low);
            end
            if isfinite(distribution.step)
                value=distribution.low+round( ...
                    (value-distribution.low)/distribution.step)*distribution.step;
            end
            value=min(max(value,distribution.low),distribution.high);
            if distribution.kind=="integer", value=round(value); end
        end

        function attach(obj,study)
            changed=isempty(obj.AttachedStudy) || ~isequal(obj.AttachedStudy,study);
            if changed
                obj.AttachedStudy=study;
                obj.Stream=RandStream("mt19937ar","Seed",obj.Seed);
                obj.ObjectiveTheta=cell(0,1);
                obj.ConstraintTheta=cell(0,1);
                obj.Restored=false;
            end
            if obj.Restored, return, end
            state=study.samplerState(obj.SamplerName,obj.StateSchema);
            if ~isempty(state), obj.restoreState(state); end
            obj.Restored=true;
        end

        function restoreState(obj,state)
            required=["schema","seed","random_state", ...
                "objective_theta","constraint_theta"];
            if ~isstruct(state) || ~isscalar(state) || ...
                    any(~isfield(state,required)) || ...
                    string(state.schema)~=obj.StateSchema || ...
                    double(state.seed)~=obj.Seed
                error("radia:optuna:GPState", ...
                    "Stored GP sampler state is invalid or incompatible.");
            end
            try
                obj.Stream.State=state.random_state;
            catch exception
                error("radia:optuna:GPState", ...
                    "Stored GP random state is invalid: %s",exception.message);
            end
            obj.ObjectiveTheta=state.objective_theta;
            obj.ConstraintTheta=state.constraint_theta;
        end

        function recordState(obj,study,trialNumber)
            state=struct("schema",obj.StateSchema,"seed",obj.Seed, ...
                "random_state",obj.Stream.State, ...
                "objective_theta",{obj.ObjectiveTheta}, ...
                "constraint_theta",{obj.ConstraintTheta});
            generation=sum(study.TrialTable.State=="COMPLETE");
            study.recordSamplerState(obj.SamplerName,obj.StateSchema, ...
                trialNumber,generation,state);
        end
    end
end
