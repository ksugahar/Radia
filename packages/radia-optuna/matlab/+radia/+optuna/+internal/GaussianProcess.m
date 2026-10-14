classdef GaussianProcess
    %GAUSSIANPROCESS Toolbox-free Matérn-5/2 ARD Gaussian process.

    methods (Static)
        function [model,theta]=fit(x,y,categorical,deterministic,initialTheta)
            arguments
                x double
                y double
                categorical logical
                deterministic (1,1) logical = false
                initialTheta double = zeros(0,1)
            end
            x=double(x);
            y=reshape(double(y),[],1);
            categorical=reshape(logical(categorical),1,[]);
            if isempty(x) || size(x,1)~=numel(y) || ...
                    size(x,2)~=numel(categorical) || ...
                    any(~isfinite(x),"all") || any(~isfinite(y))
                error("radia:optuna:GPObservations", ...
                    "GP observations must be finite and aligned.");
            end
            center=mean(y);
            scale=std(y,1);
            if ~isfinite(scale) || scale<1e-12
                scale=1;
            end
            standardized=(y-center)/scale;
            dimensions=size(x,2);
            parameterCount=dimensions+1+double(~deterministic);
            if numel(initialTheta)~=parameterCount || ...
                    any(~isfinite(initialTheta))
                % Optuna initializes inverse squared length scales and
                % kernel scale to one. theta stores their natural logs.
                initialTheta=zeros(dimensions+1,1);
                if ~deterministic
                    initialTheta(end+1,1)=0;
                end
            else
                initialTheta=reshape(double(initialTheta),[],1);
            end
            objective=@(candidate)radia.optuna.internal. ...
                GaussianProcess.negativeLogLikelihood(candidate,x, ...
                standardized,categorical,deterministic);
            optimizerOptions=optimset("Display","off","MaxIter",300, ...
                "MaxFunEvals",2000,"TolX",1e-4,"TolFun",1e-4);
            theta=fminsearch(objective,initialTheta,optimizerOptions);
            theta=max(min(reshape(theta,[],1),6),-7);
            [kernel,noiseVariance]=radia.optuna.internal. ...
                GaussianProcess.trainingKernel( ...
                theta,x,categorical,deterministic);
            [factor,jitter]=radia.optuna.internal. ...
                GaussianProcess.stableCholesky(kernel+ ...
                noiseVariance*eye(size(kernel)));
            alpha=factor'\(factor\standardized);
            model=struct("x",x,"categorical",categorical, ...
                "theta",theta,"factor",factor,"alpha",alpha, ...
                "center",center,"scale",scale,"jitter",jitter, ...
                "deterministic",deterministic);
        end

        function [meanValue,stdValue]=predict(model,query)
            query=double(query);
            if size(query,2)~=size(model.x,2) || ...
                    any(~isfinite(query),"all")
                error("radia:optuna:GPPrediction", ...
                    "GP query points must be finite and dimension-compatible.");
            end
            dimensions=size(model.x,2);
            inverseSquaredLengthScales=exp(model.theta(1:dimensions));
            lengthScales=1./sqrt(inverseSquaredLengthScales);
            signalVariance=exp(model.theta(dimensions+1));
            cross=radia.optuna.internal.GaussianProcess.kernel( ...
                model.x,query,lengthScales,signalVariance, ...
                model.categorical);
            standardizedMean=cross'*model.alpha;
            projected=model.factor\cross;
            latentVariance=max(signalVariance-sum(projected.^2,1)',1e-14);
            meanValue=model.center+model.scale*standardizedMean;
            stdValue=model.scale*sqrt(latentVariance);
        end
    end

    methods (Static, Access=private)
        function value=negativeLogLikelihood(theta,x,y,categorical,deterministic)
            theta=max(min(reshape(double(theta),[],1),6),-7);
            try
                [kernel,noiseVariance]=radia.optuna.internal. ...
                    GaussianProcess.trainingKernel( ...
                    theta,x,categorical,deterministic);
                factor=radia.optuna.internal.GaussianProcess. ...
                    stableCholesky(kernel+noiseVariance*eye(size(kernel)));
                alpha=factor'\(factor\y);
                value=0.5*(y'*alpha)+sum(log(diag(factor)))+ ...
                    0.5*numel(y)*log(2*pi)+ ...
                    radia.optuna.internal.GaussianProcess.negativeLogPrior( ...
                    theta,deterministic);
                if ~isfinite(value), value=realmax("double")/4; end
            catch
                value=realmax("double")/4;
            end
        end

        function [matrix,noiseVariance]=trainingKernel( ...
                theta,x,categorical,deterministic)
            dimensions=size(x,2);
            inverseSquaredLengthScales=exp(theta(1:dimensions));
            lengthScales=1./sqrt(inverseSquaredLengthScales);
            signalVariance=exp(theta(dimensions+1));
            if deterministic
                noiseVariance=1e-6;
            else
                noiseVariance=1e-6+exp(theta(dimensions+2));
            end
            matrix=radia.optuna.internal.GaussianProcess.kernel( ...
                x,x,lengthScales,signalVariance,categorical);
        end

        function value=negativeLogPrior(theta,deterministic)
            % Match Optuna's default GP hyperparameter prior.  q is the
            % inverse squared length scale and k is the kernel scale:
            %   log p(q) = -0.1/q - 0.1*q
            %   k ~ Gamma(shape=2, rate=1)
            %   noise ~ Gamma(shape=1.1, rate=30)
            dimensions=numel(theta)-1-double(~deterministic);
            inverseSquared=exp(theta(1:dimensions));
            kernelScale=exp(theta(dimensions+1));
            value=sum(0.1./inverseSquared+0.1*inverseSquared)+ ...
                kernelScale-log(max(kernelScale,realmin("double")));
            if ~deterministic
                noiseVariance=1e-6+exp(theta(dimensions+2));
                value=value+30*noiseVariance- ...
                    0.1*log(max(noiseVariance,realmin("double")));
            end
        end

        function matrix=kernel(left,right,lengthScales,signalVariance,categorical)
            distance=zeros(size(left,1),size(right,1));
            for dimension=1:size(left,2)
                difference=left(:,dimension)-right(:,dimension)';
                if categorical(dimension)
                    distance=distance+(difference~=0)/ ...
                        max(lengthScales(dimension)^2,eps);
                else
                    distance=distance+(difference/ ...
                        max(lengthScales(dimension),eps)).^2;
                end
            end
            radius=sqrt(5*max(distance,0));
            matrix=signalVariance*(1+radius+radius.^2/3).*exp(-radius);
        end

        function [factor,jitter]=stableCholesky(matrix)
            jitter=1e-12*max(1,mean(diag(matrix)));
            for attempt=1:10
                [factor,status]=chol(matrix+jitter*eye(size(matrix)),"lower");
                if status==0
                    return
                end
                jitter=jitter*10;
            end
            error("radia:optuna:GPCholesky", ...
                "GP covariance remained singular after jitter escalation.");
        end
    end
end
