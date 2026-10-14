classdef SeparableCMAEvolutionStrategy < handle
    %SEPARABLECMAEVOLUTIONSTRATEGY Optuna/cmaes 0.12 SepCMA state.

    properties (SetAccess=private)
        Dimension (1,1) double
        PopulationSize (1,1) double
        Generation (1,1) double = 0
        Mean (1,:) double
        Sigma (1,1) double
        Covariance (1,:) double
        PSigma (1,:) double
        PC (1,:) double
        Bounds double
        Stream
    end

    properties (Access=private)
        Mu (1,1) double
        MuEff (1,1) double
        Weights (:,1) double
        CC (1,1) double
        C1 (1,1) double
        CMu (1,1) double
        CSigma (1,1) double
        DSigma (1,1) double
        ChiN (1,1) double
        MaxResampling (1,1) double = 100
    end

    methods
        function obj=SeparableCMAEvolutionStrategy(mean,sigma,options)
            arguments
                mean (1,:) double
                sigma (1,1) double
                options.Bounds double = zeros(0,2)
                options.PopulationSize (1,1) double = 0
                options.Seed (1,1) double = 0
                options.Covariance double = zeros(0,0)
                options.MaxResampling (1,1) double = 100
            end
            dimension=numel(mean);
            if dimension<=1 || any(~isfinite(mean))
                error("radia:optuna:CMAMean", ...
                    "Separable CMA-ES requires at least two finite dimensions.");
            end
            if ~(isfinite(sigma)&&sigma>0)
                error("radia:optuna:CMASigma", ...
                    "CMA-ES sigma must be positive and finite.");
            end
            population=options.PopulationSize;
            if population==0,population=4+floor(3*log(dimension));end
            if population<2 || population~=floor(population)
                error("radia:optuna:CMAPopulation", ...
                    "CMA-ES population size must be an integer of at least two.");
            end
            bounds=options.Bounds;
            if isempty(bounds),bounds=[-inf(dimension,1),inf(dimension,1)];end
            if ~isequal(size(bounds),[dimension,2]) || ...
                    any(mean(:)<bounds(:,1)) || any(mean(:)>bounds(:,2))
                error("radia:optuna:CMABounds","Invalid CMA-ES bounds.");
            end
            covariance=options.Covariance;
            if isempty(covariance),covariance=ones(1,dimension);end
            if ismatrix(covariance) && isequal(size(covariance), ...
                    [dimension,dimension])
                covariance=diag(covariance).';
            end
            covariance=reshape(double(covariance),1,[]);
            if numel(covariance)~=dimension || any(~isfinite(covariance))
                error("radia:optuna:CMACovariance", ...
                    "Invalid separable CMA-ES covariance.");
            end

            obj.Dimension=dimension;
            obj.PopulationSize=population;
            obj.Mean=double(mean);
            obj.Sigma=double(sigma);
            obj.Covariance=covariance;
            obj.PSigma=zeros(1,dimension);
            obj.PC=zeros(1,dimension);
            obj.Bounds=double(bounds);
            obj.Stream=radia.optuna.internal.NumpyRandomState(double(options.Seed));
            obj.MaxResampling=options.MaxResampling;
            obj.initializeCoefficients();
        end

        function point=ask(obj)
            for attempt=1:obj.MaxResampling
                point=obj.samplePoint();
                if all(point(:)>=obj.Bounds(:,1)) && ...
                        all(point(:)<=obj.Bounds(:,2)),return,end
            end
            point=obj.samplePoint();
            point=min(max(point,obj.Bounds(:,1).'),obj.Bounds(:,2).');
        end

        function reseed(obj,seed)
            obj.Stream=radia.optuna.internal.NumpyRandomState(double(seed));
        end

        function tell(obj,points,fitness)
            points=double(points);
            fitness=reshape(double(fitness),[],1);
            if ~isequal(size(points),[obj.PopulationSize,obj.Dimension]) || ...
                    numel(fitness)~=obj.PopulationSize
                error("radia:optuna:CMATell", ...
                    "CMA-ES tell requires one full population.");
            end
            [~,order]=sort(fitness,"ascend");
            points=points(order,:);
            obj.Generation=obj.Generation+1;
            scales=sqrt(max(obj.Covariance,1e-8));
            y=(points-obj.Mean)/obj.Sigma;
            yWeighted=obj.Weights.'*y(1:obj.Mu,:);
            obj.Mean=obj.Mean+obj.Sigma*yWeighted;
            obj.PSigma=(1-obj.CSigma)*obj.PSigma+sqrt( ...
                obj.CSigma*(2-obj.CSigma)*obj.MuEff)*(yWeighted./scales);
            normPath=norm(obj.PSigma);
            obj.Sigma=min(1e32,obj.Sigma*exp(obj.CSigma/obj.DSigma* ...
                (normPath/obj.ChiN-1)));
            left=normPath/sqrt(1-(1-obj.CSigma)^(2*(obj.Generation+1)));
            right=(1.4+2/(obj.Dimension+1))*obj.ChiN;
            hSigma=double(left<right);
            obj.PC=(1-obj.CC)*obj.PC+hSigma*sqrt( ...
                obj.CC*(2-obj.CC)*obj.MuEff)*yWeighted;
            delta=(1-hSigma)*obj.CC*(2-obj.CC);
            rankMu=sum((obj.Weights.*ones(1,obj.Dimension)).* ...
                (y(1:obj.Mu,:).^2),1);
            obj.Covariance=(1+obj.C1*delta-obj.C1- ...
                obj.CMu*sum(obj.Weights))*obj.Covariance+ ...
                obj.C1*(obj.PC.^2)+obj.CMu*rankMu;
        end

        function state=snapshot(obj)
            state=struct("schema","radia.optuna.separable-cma-state.v1", ...
                "dimension",obj.Dimension,"population_size",obj.PopulationSize, ...
                "generation",obj.Generation,"mean",obj.Mean, ...
                "sigma",obj.Sigma,"covariance",obj.Covariance, ...
                "p_sigma",obj.PSigma,"p_c",obj.PC,"bounds",obj.Bounds, ...
                "max_resampling",obj.MaxResampling, ...
                "random_state",obj.Stream.State);
        end

        function result=shouldStop(obj)
            scales=sqrt(max(obj.Covariance,1e-8));
            result=any(~isfinite(scales)) || ...
                obj.Sigma*max(scales)>1e4 || ...
                max(scales)/max(min(scales),eps)>1e14;
        end
    end

    methods (Static)
        function obj=fromSnapshot(state)
            if ~isstruct(state) || string(state.schema)~= ...
                    "radia.optuna.separable-cma-state.v1"
                error("radia:optuna:CMAState", ...
                    "Invalid separable CMA-ES state.");
            end
            obj=radia.optuna.internal.SeparableCMAEvolutionStrategy( ...
                reshape(double(state.mean),1,[]),double(state.sigma), ...
                Bounds=double(state.bounds), ...
                PopulationSize=double(state.population_size), ...
                Covariance=double(state.covariance), ...
                MaxResampling=double(state.max_resampling));
            obj.Generation=double(state.generation);
            obj.PSigma=reshape(double(state.p_sigma),1,[]);
            obj.PC=reshape(double(state.p_c),1,[]);
            obj.Stream.State=state.random_state;
        end
    end

    methods (Access=private)
        function initializeCoefficients(obj)
            obj.Mu=floor(obj.PopulationSize/2);
            raw=log(obj.Mu+1)-log((1:obj.Mu).');
            obj.Weights=raw/sum(raw);
            obj.MuEff=1/sum(obj.Weights.^2);
            obj.C1=2/((obj.Dimension+1.3)^2+obj.MuEff);
            full=2/obj.MuEff/((obj.Dimension+sqrt(2))^2)+ ...
                (1-1/obj.MuEff)*min(1,(2*obj.MuEff-1)/ ...
                ((obj.Dimension+2)^2+obj.MuEff));
            obj.CMu=(obj.Dimension+2)/3*full;
            obj.CSigma=(obj.MuEff+2)/(obj.Dimension+obj.MuEff+3);
            obj.DSigma=1+2*max(0,sqrt((obj.MuEff-1)/ ...
                (obj.Dimension+1))-1)+obj.CSigma;
            obj.CC=4/(obj.Dimension+4);
            obj.ChiN=sqrt(obj.Dimension)*(1-1/(4*obj.Dimension)+ ...
                1/(21*obj.Dimension^2));
        end

        function point=samplePoint(obj)
            z=randn(obj.Stream,obj.Dimension,1);
            point=obj.Mean+obj.Sigma*(sqrt(max(obj.Covariance,1e-8)).'.*z).';
        end
    end
end
