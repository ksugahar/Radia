classdef RawRead
    %RAWREAD MATLAB-native LTspice RAW reader facade.
    properties (SetAccess=private), Data (1,1) struct, end
    methods
        function obj=RawRead(path), obj.Data=radia.ltspice.readRaw(path); end
        function names=getTraceNames(obj), names=obj.Data.names; end
        function wave=getTrace(obj,name)
            j=find(obj.Data.names==string(name),1); if isempty(j), error("radia:ltspice:TraceNotFound","Trace not found: %s",name); end
            wave=obj.Data.values(:,j);
        end
        function t=getTime(obj), t=obj.getTrace("time"); end
        function axis=getAxis(obj), axis=obj.Data.values(:,1); end
        function answer=isComplex(obj), answer=~isreal(obj.Data.values); end
        function count=getStepCount(obj), count=obj.Data.step_count; end
        function wave=getStep(obj,name,step)
            r=obj.Data.step_ranges(step,:); all=obj.getTrace(name); wave=all(r(1):r(2));
        end
    end
end
