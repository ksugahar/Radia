classdef UNDXCrossover < radia.optuna.nsgaii.BaseCrossover
    %UNDXCROSSOVER Three-parent unimodal normal distribution crossover.

    properties (SetAccess=private)
        SigmaXi (1,1) double = 0.5
        SigmaEta (1,1) double = NaN
    end

    methods
        function obj = UNDXCrossover(options)
            arguments
                options.SigmaXi (1,1) double = 0.5
                options.SigmaEta (1,1) double = NaN
            end
            obj.NParents = 3;
            obj.Name = "UNDXCrossover";
            obj.SigmaXi = options.SigmaXi;
            obj.SigmaEta = options.SigmaEta;
        end

        function child = crossover(obj,parents,stream,study,bounds) %#ok<INUSD>
            dimension = size(bounds,1);
            midpoint = (parents(1,:)+parents(2,:))/2;
            primary = parents(1,:)-parents(2,:);
            sigmaEta = obj.SigmaEta;
            if isnan(sigmaEta)
                sigmaEta = 0.35/sqrt(dimension);
            end
            etas = sigmaEta^2*randn(stream,1,dimension);
            xi = obj.SigmaXi^2*randn(stream);
            direction = parents(2,:)-parents(1,:);
            direction = direction/max(norm(direction,2),1e-10);
            basis = eye(dimension);
            if any(direction~=0)
                basis(1,:) = direction;
            end
            [orthogonal,~] = qr(basis.');
            perpendicularBasis = orthogonal.';
            perpendicularBasis = perpendicularBasis(2:end,:);
            child = midpoint+xi*primary;
            if dimension>1
                thirdVector = parents(3,:)-parents(1,:);
                distance = norm(thirdVector- ...
                    dot(thirdVector,direction)*direction,2);
                subsearch = zeros(1,dimension);
                for index = 1:dimension-1
                    subsearch = subsearch+ ...
                        etas(index)*perpendicularBasis(index,:);
                end
                child = child+distance*subsearch;
            end
        end

        function config = configuration(obj)
            config = configuration@radia.optuna.nsgaii.BaseCrossover(obj);
            config.sigma_xi = obj.SigmaXi;
            config.sigma_eta = obj.SigmaEta;
        end
    end
end
