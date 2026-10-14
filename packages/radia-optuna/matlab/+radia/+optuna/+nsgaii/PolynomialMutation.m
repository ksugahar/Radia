classdef PolynomialMutation < radia.optuna.nsgaii.BaseMutation
    %POLYNOMIALMUTATION Optuna 5 polynomial mutation for numerical values.

    properties (SetAccess=private)
        Eta (1,1) double = 20.0
    end

    methods
        function obj=PolynomialMutation(options)
            arguments
                options.Eta (1,1) double = 20.0
            end
            if ~isfinite(options.Eta) || options.Eta<0
                error("radia:optuna:PolynomialMutationEta", ...
                    "Eta must be a nonnegative finite scalar.");
            end
            obj.Eta=double(options.Eta);
        end

        function value=mutation(obj,param,rng,study,searchSpaceBounds) %#ok<INUSD>
            arguments
                obj
                param (1,1) double
                rng (1,1) radia.optuna.internal.NumpyRandomState
                study
                searchSpaceBounds (1,2) double
            end
            u=rand(rng);
            lower=searchSpaceBounds(1);
            upper=searchSpaceBounds(2);
            width=upper-lower;
            if width<=0
                value=param;
                return
            end

            delta1=(param-lower)/width;
            delta2=(upper-param)/width;
            mutationPower=1/(obj.Eta+1);
            if u<=0.5
                xy=1-delta1;
                transformed=2*u+(1-2*u)*xy^(obj.Eta+1);
                deltaQ=transformed^mutationPower-1;
            else
                xy=1-delta2;
                transformed=2*(1-u)+2*(u-0.5)*xy^(obj.Eta+1);
                deltaQ=1-transformed^mutationPower;
            end
            value=param+deltaQ*width;
        end
    end
end
