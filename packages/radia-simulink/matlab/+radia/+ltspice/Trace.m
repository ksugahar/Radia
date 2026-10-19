classdef Trace
    %TRACE One RAW trace with PyLTSpice-compatible step slicing.
    properties (SetAccess=private), Name (1,1) string; Data (:,1); StepRanges (:,2) double; end
    methods
        function obj=Trace(name,data,stepRanges), obj.Name=string(name); obj.Data=data(:); obj.StepRanges=stepRanges; end
        function wave=getWave(obj,step)
            if nargin<2, wave=obj.Data; return, end
            index=double(step)+1;
            if index<1||index>size(obj.StepRanges,1), error("radia:ltspice:StepOutOfRange","Step %d is out of range.",step); end
            r=obj.StepRanges(index,:); wave=obj.Data(r(1):r(2));
        end
        function wave=get_wave(obj,varargin),wave=obj.getWave(varargin{:});end
        function n=length(obj),n=numel(obj.Data);end
    end
end
