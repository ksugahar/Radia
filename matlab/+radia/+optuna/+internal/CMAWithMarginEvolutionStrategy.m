classdef CMAWithMarginEvolutionStrategy < handle
    %CMAWITHMARGINEVOLUTIONSTRATEGY Optuna/cmaes 0.12 CMAwM state.

    properties (SetAccess=private)
        Base
        Dimension (1,1) double
        PopulationSize (1,1) double
        Generation (1,1) double
        Sigma (1,1) double
    end

    properties (Access=private)
        MaxResampling (1,1) double = 100
        DiscreteIndices double = zeros(1,0)
        ContinuousIndices double = zeros(1,0)
        ContinuousBounds double = zeros(0,2)
        ZSpace double = zeros(0,0)
        ZLimits double = zeros(0,0)
        LowerLimits double = zeros(1,0)
        UpperLimits double = zeros(1,0)
        ScaleA double = zeros(1,0)
        Margin (1,1) double = NaN
    end

    methods
        function obj=CMAWithMarginEvolutionStrategy(mean,sigma,options)
            arguments
                mean (1,:) double
                sigma (1,1) double
                options.Bounds double
                options.Steps (1,:) double
                options.PopulationSize (1,1) double = 0
                options.Seed (1,1) double = 0
                options.Covariance double = zeros(0,0)
                options.MaxResampling (1,1) double = 100
                options.Margin (1,1) double = NaN
            end
            obj.Base=radia.optuna.internal.CMAEvolutionStrategy(mean,sigma, ...
                Bounds=options.Bounds,PopulationSize=options.PopulationSize, ...
                Seed=options.Seed,Covariance=options.Covariance, ...
                MaxResampling=options.MaxResampling);
            obj.Dimension=obj.Base.Dimension;
            obj.PopulationSize=obj.Base.PopulationSize;
            obj.Generation=obj.Base.Generation;
            obj.Sigma=obj.Base.Sigma;
            obj.MaxResampling=options.MaxResampling;
            steps=reshape(double(options.Steps),1,[]);
            if numel(steps)~=obj.Dimension || any(isnan(steps))
                error("radia:optuna:CMAMarginSteps", ...
                    "Margin steps must match the CMA-ES dimension.");
            end
            obj.DiscreteIndices=find(steps>0);
            obj.ContinuousIndices=find(steps<=0);
            obj.ContinuousBounds=options.Bounds(obj.ContinuousIndices,:);
            obj.ScaleA=ones(1,obj.Dimension);
            if isempty(obj.DiscreteIndices),return,end

            counts=zeros(1,numel(obj.DiscreteIndices));
            values=cell(1,numel(obj.DiscreteIndices));
            for row=1:numel(obj.DiscreteIndices)
                index=obj.DiscreteIndices(row);
                values{row}=options.Bounds(index,1):steps(index): ...
                    (options.Bounds(index,2)+steps(index)/2);
                counts(row)=numel(values{row});
            end
            width=max(counts);
            obj.ZSpace=NaN(numel(values),width);
            for row=1:numel(values)
                obj.ZSpace(row,1:counts(row))=values{row};
            end
            obj.ZLimits=(obj.ZSpace(:,2:end)+obj.ZSpace(:,1:end-1))/2;
            for row=1:size(obj.ZSpace,1)
                obj.ZSpace(row,isnan(obj.ZSpace(row,:)))= ...
                    max(obj.ZSpace(row,:),[],"omitnan");
                obj.ZLimits(row,isnan(obj.ZLimits(row,:)))= ...
                    max(obj.ZLimits(row,:),[],"omitnan");
            end
            obj.updateLimits(obj.Base.Mean(obj.DiscreteIndices));
            if isfinite(options.Margin)
                obj.Margin=options.Margin;
            else
                obj.Margin=1/(obj.Dimension*obj.PopulationSize);
            end
            if obj.Margin<=0
                error("radia:optuna:CMAMargin", ...
                    "CMA-ES margin must be positive.");
            end
        end

        function [encoded,raw]=ask(obj)
            for attempt=1:obj.MaxResampling
                raw=obj.Base.sampleRaw();
                if obj.continuousFeasible(raw(obj.ContinuousIndices))
                    encoded=obj.encode(raw);
                    return
                end
            end
            raw=obj.Base.sampleRaw();
            raw(obj.ContinuousIndices)=min(max( ...
                raw(obj.ContinuousIndices), ...
                obj.ContinuousBounds(:,1).'),obj.ContinuousBounds(:,2).');
            encoded=obj.encode(raw);
        end

        function reseed(obj,seed)
            obj.Base.reseed(seed);
        end

        function tell(obj,points,fitness)
            obj.Base.tell(points,fitness);
            obj.Generation=obj.Base.Generation;
            obj.Sigma=obj.Base.Sigma;
            if isempty(obj.DiscreteIndices),return,end
            mean=obj.Base.Mean;
            covariance=obj.Base.Covariance;
            discreteMean=mean(obj.DiscreteIndices);
            obj.updateLimits(discreteMean);
            diagonal=sqrt(diag(covariance)).';
            scale=obj.Base.Sigma*obj.ScaleA(obj.DiscreteIndices).* ...
                diagonal(obj.DiscreteIndices);
            low=obj.normalCdf(obj.LowerLimits,discreteMean,scale);
            upper=1-obj.normalCdf(obj.UpperLimits,discreteMean,scale);
            middle=1-(low+upper);
            edge=max(low,upper)>0.5;
            side=~edge;
            if any(edge)
                modify=min(low,upper)<obj.Margin;
                direction=sign(discreteMean-obj.UpperLimits);
                distance=obj.Base.Sigma*obj.ScaleA(obj.DiscreteIndices).* ...
                    sqrt(obj.chiSquareOneQuantile(1-2*obj.Margin).* ...
                    diag(covariance(obj.DiscreteIndices, ...
                    obj.DiscreteIndices)).');
                discreteMean=discreteMean+modify.*edge.*( ...
                    obj.UpperLimits+direction.*distance-discreteMean);
            end
            low=max(low,obj.Margin/2);
            upper=max(upper,obj.Margin/2);
            denominator=low+middle+upper-3*obj.Margin/2;
            modifiedLow=low+(1-low-upper-middle).* ...
                (low-obj.Margin/2)./denominator;
            modifiedUpper=upper+(1-low-upper-middle).* ...
                (upper-obj.Margin/2)./denominator;
            modifiedLow=min(max(modifiedLow,1e-10),0.5-1e-10);
            modifiedUpper=min(max(modifiedUpper,1e-10),0.5-1e-10);
            chiLow=sqrt(obj.chiSquareOneQuantile(1-2*modifiedLow));
            chiUpper=sqrt(obj.chiSquareOneQuantile(1-2*modifiedUpper));
            diagonal=diagonal(obj.DiscreteIndices);
            obj.ScaleA(obj.DiscreteIndices)=obj.ScaleA(obj.DiscreteIndices)+ ...
                side.*((obj.UpperLimits-obj.LowerLimits)./ ...
                ((chiLow+chiUpper)*obj.Base.Sigma.*diagonal)- ...
                obj.ScaleA(obj.DiscreteIndices));
            discreteMean=discreteMean+side.*((obj.LowerLimits.*chiUpper+ ...
                obj.UpperLimits.*chiLow)./(chiLow+chiUpper)-discreteMean);
            mean(obj.DiscreteIndices)=discreteMean;
            obj.Base.setMean(mean);
            obj.Sigma=obj.Base.Sigma;
        end

        function state=snapshot(obj)
            state=struct("schema","radia.optuna.cma-with-margin-state.v1", ...
                "base",obj.Base.snapshot(),"max_resampling",obj.MaxResampling, ...
                "discrete_indices",obj.DiscreteIndices, ...
                "continuous_indices",obj.ContinuousIndices, ...
                "continuous_bounds",obj.ContinuousBounds, ...
                "z_space",obj.ZSpace,"z_limits",obj.ZLimits, ...
                "lower_limits",obj.LowerLimits,"upper_limits",obj.UpperLimits, ...
                "scale_a",obj.ScaleA,"margin",obj.Margin);
        end

        function result=shouldStop(obj)
            result=obj.Base.shouldStop();
        end
    end

    methods (Static)
        function obj=fromSnapshot(state)
            if ~isstruct(state) || string(state.schema)~= ...
                    "radia.optuna.cma-with-margin-state.v1"
                error("radia:optuna:CMAState", ...
                    "Invalid CMA-ES-with-margin state.");
            end
            obj=radia.optuna.internal.CMAWithMarginEvolutionStrategy.empty();
            obj.Base=radia.optuna.internal.CMAEvolutionStrategy. ...
                fromSnapshot(state.base);
            obj.Dimension=obj.Base.Dimension;
            obj.PopulationSize=obj.Base.PopulationSize;
            obj.Generation=obj.Base.Generation;
            obj.Sigma=obj.Base.Sigma;
            obj.MaxResampling=double(state.max_resampling);
            obj.DiscreteIndices=reshape(double(state.discrete_indices),1,[]);
            obj.ContinuousIndices=reshape(double(state.continuous_indices),1,[]);
            obj.ContinuousBounds=double(state.continuous_bounds);
            obj.ZSpace=double(state.z_space);
            obj.ZLimits=double(state.z_limits);
            obj.LowerLimits=reshape(double(state.lower_limits),1,[]);
            obj.UpperLimits=reshape(double(state.upper_limits),1,[]);
            obj.ScaleA=reshape(double(state.scale_a),1,[]);
            obj.Margin=double(state.margin);
        end
    end

    methods (Access=private)
        function encoded=encode(obj,raw)
            encoded=raw;
            if isempty(obj.DiscreteIndices),return,end
            mean=obj.Base.Mean;
            values=(raw(obj.DiscreteIndices)-mean(obj.DiscreteIndices)).* ...
                obj.ScaleA(obj.DiscreteIndices)+mean(obj.DiscreteIndices);
            for row=1:numel(values)
                position=sum(obj.ZLimits(row,:)<values(row))+1;
                encoded(obj.DiscreteIndices(row))=obj.ZSpace(row,position);
            end
        end

        function result=continuousFeasible(obj,values)
            result=all(values(:)>=obj.ContinuousBounds(:,1)) && ...
                all(values(:)<=obj.ContinuousBounds(:,2));
        end

        function updateLimits(obj,values)
            count=numel(values);
            obj.LowerLimits=zeros(1,count);
            obj.UpperLimits=zeros(1,count);
            last=size(obj.ZLimits,2);
            for row=1:count
                position=sum(obj.ZLimits(row,:)<values(row));
                low=min(max(position,1),last);
                high=min(max(position+1,1),last);
                obj.LowerLimits(row)=obj.ZLimits(row,low);
                obj.UpperLimits(row)=obj.ZLimits(row,high);
            end
        end
    end

    methods (Static,Access=private)
        function value=normalCdf(x,mean,scale)
            value=0.5*erfc(-(x-mean)./(scale*sqrt(2)));
        end

        function value=chiSquareOneQuantile(probability)
            % chi2inv(p,1) without a Statistics Toolbox dependency.
            value=2*erfinv(probability).^2;
        end
    end
end
