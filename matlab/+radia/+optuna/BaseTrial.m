classdef (Abstract) BaseTrial < handle
    %BASETRIAL Common Optuna 4.9 trial interface.

    methods
        function value=datetime_start(obj)
            if isprop(obj,"DatetimeStart")
                value=obj.DatetimeStart; %#ok<MCNPN>
            else
                value=obj.StartTime; %#ok<MCNPN>
            end
        end

        function value=number(obj)
            value=obj.Number; %#ok<MCNPN>
        end

        function value=params(obj)
            value=obj.Params; %#ok<MCNPN>
        end

        function value=distributions(obj)
            value=obj.Distributions; %#ok<MCNPN>
        end

        function value=user_attrs(obj)
            value=obj.UserAttrs; %#ok<MCNPN>
        end

        function value=system_attrs(obj)
            value=obj.SystemAttrs; %#ok<MCNPN>
        end

        function report(~,~,~)
            error("radia:optuna:AbstractTrial", ...
                "The concrete trial must implement report.");
        end

        function set_user_attr(obj,name,value)
            obj.setUserAttr(name,value); %#ok<MCNPN>
        end

        function set_system_attr(obj,name,value)
            obj.setSystemAttr(name,value); %#ok<MCNPN>
        end

        function value=should_prune(obj)
            value=obj.shouldPrune(); %#ok<MCNPN>
        end

        function value=suggest_float(obj,name,low,high,varargin)
            value=obj.suggestFloat(name,low,high,varargin{:}); %#ok<MCNPN>
        end

        function value=suggest_int(obj,name,low,high,varargin)
            value=obj.suggestInteger(name,low,high,varargin{:}); %#ok<MCNPN>
        end

        function value=suggest_categorical(obj,name,choices)
            value=obj.suggestCategorical(name,choices); %#ok<MCNPN>
        end

        function value=suggest_uniform(obj,name,low,high)
            value=obj.suggestFloat(name,low,high); %#ok<MCNPN>
        end

        function value=suggest_loguniform(obj,name,low,high)
            value=obj.suggestFloat(name,low,high,Log=true); %#ok<MCNPN>
        end

        function value=suggest_discrete_uniform(obj,name,low,high,q)
            value=obj.suggestFloat(name,low,high,Step=q); %#ok<MCNPN>
        end
    end
end
