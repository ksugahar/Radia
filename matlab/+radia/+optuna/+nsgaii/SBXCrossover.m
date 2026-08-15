classdef SBXCrossover < radia.optuna.nsgaii.BaseCrossover
    %SBXCROSSOVER Simulated binary crossover matching Optuna 4.9.

    properties (SetAccess=private)
        Eta (1,1) double = NaN
        UniformCrossoverProbability (1,1) double = 0.5
        UseChildGeneProbability (1,1) double = 0.5
    end

    methods
        function obj = SBXCrossover(options)
            arguments
                options.Eta (1,1) double = NaN
                options.UniformCrossoverProbability (1,1) double = 0.5
                options.UseChildGeneProbability (1,1) double = 0.5
            end
            if ~isnan(options.Eta) && options.Eta < 0
                error("radia:optuna:NSGAIICrossover", ...
                    "Eta must be nonnegative or NaN.");
            end
            validateProbabilities(options.UniformCrossoverProbability, ...
                options.UseChildGeneProbability);
            obj.NParents = 2;
            obj.Name = "SBXCrossover";
            obj.Eta = options.Eta;
            obj.UniformCrossoverProbability = ...
                options.UniformCrossoverProbability;
            obj.UseChildGeneProbability = options.UseChildGeneProbability;
        end

        function child = crossover(obj,parents,stream,study,bounds)
            lower = bounds(:,1).';
            upper = bounds(:,2).';
            parentsMin = min(parents,[],1);
            parentsMax = max(parents,[],1);
            eta = obj.Eta;
            if isnan(eta)
                if numel(study.Directions)>1
                    eta = 20;
                else
                    eta = 2;
                end
            end
            difference = max(parentsMax-parentsMin,1e-10);
            beta1 = 1+2*(parentsMin-lower)./difference;
            beta2 = 1+2*(upper-parentsMax)./difference;
            alpha1 = 2-beta1.^(-(eta+1));
            alpha2 = 2-beta2.^(-(eta+1));
            uniforms = rand(stream,1,size(parents,2));
            betaQ1 = (uniforms.*alpha1).^(1/(eta+1));
            mask = uniforms>1./alpha1;
            betaQ1(mask) = (1./(2-uniforms(mask).*alpha1(mask))).^ ...
                (1/(eta+1));
            betaQ2 = (uniforms.*alpha2).^(1/(eta+1));
            mask = uniforms>1./alpha2;
            betaQ2(mask) = (1./(2-uniforms(mask).*alpha2(mask))).^ ...
                (1/(eta+1));
            firstChild = 0.5*((parentsMin+parentsMax)- ...
                betaQ1.*difference);
            secondChild = 0.5*((parentsMin+parentsMax)+ ...
                betaQ2.*difference);
            child = obj.selectChildGenes( ...
                firstChild,secondChild,parents,stream);
        end

        function config = configuration(obj)
            config = configuration@radia.optuna.nsgaii.BaseCrossover(obj);
            config.eta = obj.Eta;
            config.uniform_crossover_probability = ...
                obj.UniformCrossoverProbability;
            config.use_child_gene_probability = ...
                obj.UseChildGeneProbability;
        end
    end

    methods (Access=protected)
        function child = selectChildGenes(obj,firstChild,secondChild, ...
                parents,stream)
            dimension = size(parents,2);
            candidate1 = zeros(1,dimension);
            candidate2 = zeros(1,dimension);
            for index = 1:dimension
                useChild = rand(stream)<obj.UseChildGeneProbability;
                swap = rand(stream)<obj.UniformCrossoverProbability;
                if useChild
                    pair = [firstChild(index),secondChild(index)];
                else
                    pair = parents(:,index).';
                end
                if swap
                    pair = pair([2,1]);
                end
                candidate1(index) = pair(1);
                candidate2(index) = pair(2);
            end
            if rand(stream)<0.5
                child = candidate1;
            else
                child = candidate2;
            end
        end
    end
end

function validateProbabilities(uniformProbability,useChildProbability)
if ~isfinite(uniformProbability) || uniformProbability<0 || ...
        uniformProbability>1
    error("radia:optuna:NSGAIIProbability", ...
        "UniformCrossoverProbability must be in [0,1].");
end
if ~isfinite(useChildProbability) || useChildProbability<=0 || ...
        useChildProbability>1
    error("radia:optuna:NSGAIIProbability", ...
        "UseChildGeneProbability must be in (0,1].");
end
end
