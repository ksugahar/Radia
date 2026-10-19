classdef RunTask < handle
    %RUNTASK Asynchronous or completed LTspice task.
    properties (SetAccess=private), Result struct=struct(); Future=[]; Callback=[]; CallbackInvoked (1,1) logical=false; end
    properties (Dependent), RawFile; LogFile; Status; end
    methods
        function obj=RunTask(resultOrFuture,callback)
            if nargin<2,callback=[];end,obj.Callback=callback;if isa(resultOrFuture,'parallel.FevalFuture'),obj.Future=resultOrFuture;else,obj.Result=resultOrFuture;end
        end
        function x=get.RawFile(obj),obj.collect(false);if isempty(fieldnames(obj.Result)),x="";else,x=string(obj.Result.raw_file);end,end
        function x=get.LogFile(obj),obj.collect(false);if isempty(fieldnames(obj.Result)),x="";else,x=string(obj.Result.log_file);end,end
        function x=get.Status(obj),if isempty(obj.Future)||strcmpi(obj.Future.State,'finished'),x="completed";else,x=lower(string(obj.Future.State));end,end
        function answer=isAlive(obj),answer=~isempty(obj.Future)&&~strcmpi(obj.Future.State,'finished');end
        function answer=wait(obj,timeout)
            if nargin<2,timeout=Inf;end,t=tic;while obj.isAlive()&&toc(t)<timeout,pause(0.02);end,answer=~obj.isAlive();if answer,obj.collect(true);end
        end
        function x=waitResults(obj,timeout),if nargin<2,timeout=Inf;end,if ~obj.wait(timeout),error("radia:ltspice:Timeout","LTspice task timed out.");end,x={obj.RawFile,obj.LogFile};end
        function cancel(obj),if obj.isAlive(),cancel(obj.Future);end,end
        function x=wait_results(obj),x=obj.waitResults();end
        function x=is_alive(obj),x=obj.isAlive();end
    end
    methods (Access=private)
        function collect(obj,invokeCallback)
            if ~isempty(obj.Future)&&strcmpi(obj.Future.State,'finished')&&isempty(fieldnames(obj.Result)),obj.Result=fetchOutputs(obj.Future);end
            if invokeCallback&&~obj.CallbackInvoked&&~isempty(obj.Callback)&&~isempty(fieldnames(obj.Result)),feval(obj.Callback,string(obj.Result.raw_file),string(obj.Result.log_file));obj.CallbackInvoked=true;end
        end
    end
end
