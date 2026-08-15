classdef VSBXCrossover < radia.optuna.nsgaii.SBXCrossover
    %VSBXCROSSOVER Modified simulated binary crossover from Optuna 4.9.

    methods
        function obj = VSBXCrossover(options)
            arguments
                options.Eta (1,1) double = NaN
                options.UniformCrossoverProbability (1,1) double = 0.5
                options.UseChildGeneProbability (1,1) double = 0.5
            end
            obj@radia.optuna.nsgaii.SBXCrossover( ...
                Eta=options.Eta, ...
                UniformCrossoverProbability= ...
                options.UniformCrossoverProbability, ...
                UseChildGeneProbability=options.UseChildGeneProbability);
            obj.Name = "VSBXCrossover";
        end

        function child = crossover(obj,parents,stream,study,~)
            eta = obj.Eta;
            if isnan(eta)
                if numel(study.Directions)>1
                    eta = 20;
                else
                    eta = 2;
                end
            end
            dimension = size(parents,2);
            uniforms = rand(stream,1,dimension);
            beta1 = (1./max(2*uniforms,1e-10)).^(1/(eta+1));
            beta2 = (1./max(2*(1-uniforms),1e-10)).^(1/(eta+1));
            if rand(stream)<=0.5
                firstChild = 0.5*((1+beta1).*parents(1,:)+ ...
                    (1-beta2).*parents(2,:));
            else
                firstChild = 0.5*((1-beta1).*parents(1,:)+ ...
                    (1+beta2).*parents(2,:));
            end
            if rand(stream)<=0.5
                secondChild = 0.5*((3-beta1).*parents(1,:)- ...
                    (1-beta2).*parents(2,:));
            else
                secondChild = 0.5*(-(1-beta1).*parents(1,:)+ ...
                    (3-beta2).*parents(2,:));
            end
            child = obj.selectChildGenes( ...
                firstChild,secondChild,parents,stream);
        end
    end
end
