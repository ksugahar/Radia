classdef RawRead
    %RAWREAD MATLAB-native PyLTSpice-compatible RAW reader facade.
    properties (SetAccess=private), Data (1,1) struct; Path (1,1) string; end
    methods
        function obj=RawRead(path)
            obj.Path=string(path); obj.Data=radia.ltspice.readRaw(path);
            [folder,name]=fileparts(path);logPath=fullfile(folder,name+".log");if isfile(logPath),obj.Data.step_parameters=localStepParameters(logPath,obj.Data.step_count);end
        end
        function names=getTraceNames(obj), names=[obj.Data.names,obj.Data.alias_names(:)']; end
        function wave=getTrace(obj,reference)
            if ~isnumeric(reference)&&~any(strcmpi(obj.Data.names,string(reference))),wave=obj.evaluateAlias(string(reference));return,end
            j=obj.traceIndex(reference); wave=obj.Data.values(:,j);
        end
        function trace=getTraceObject(obj,reference)
            if ~isnumeric(reference)&&~any(strcmpi(obj.Data.names,string(reference))),trace=radia.ltspice.Trace(string(reference),obj.getTrace(reference),obj.Data.step_ranges);return,end
            j=obj.traceIndex(reference); trace=radia.ltspice.Trace(obj.Data.names(j),obj.Data.values(:,j),obj.Data.step_ranges);
        end
        function wave=getWave(obj,reference,step)
            arguments, obj; reference; step (1,1) double {mustBeInteger,mustBeNonnegative}=0; end
            trace=obj.getTraceObject(reference); wave=trace.getWave(step);
        end
        function t=getTime(obj,step), if nargin<2,step=0;end, t=obj.getWave("time",step); end
        function axis=getAxis(obj,step), if nargin<2,step=0;end, axis=obj.getWave(1,step); end
        function answer=isComplex(obj), answer=~isreal(obj.Data.values); end
        function count=getStepCount(obj), count=obj.Data.step_count; end
        function wave=getStep(obj,name,step)
            arguments, obj; name; step (1,1) double {mustBeInteger,mustBePositive}; end
            wave=obj.getWave(name,step-1);
        end
        function n=getLen(obj,step), if nargin<2,step=0;end, n=numel(obj.getAxis(step)); end
        function n=getNrPlots(obj),if isfield(obj.Data,'plots'),n=numel(obj.Data.plots);else,n=1;end,end
        function name=getPlotName(obj), name=obj.Data.plot_name; end
        function names=getPlotNames(obj), names=obj.getPlotName(); end
        function properties=getRawProperties(obj), properties=obj.Data.raw_properties; end
        function value=getRawProperty(obj,name)
            if nargin<2, value=obj.getRawProperties(); return, end
            key=matlab.lang.makeValidName(char(name)); if ~isfield(obj.Data.raw_properties,key), error("radia:ltspice:RawPropertyNotFound","RAW property not found: %s",name); end
            value=obj.Data.raw_properties.(key);
        end
        function steps=getSteps(obj,conditions)
            arguments, obj; conditions (1,1) struct=struct(); end
            steps=0:obj.Data.step_count-1;
            if isempty(fieldnames(conditions)), return, end
            if ~isfield(obj.Data,'step_parameters'), steps=[]; return, end
            keep=true(size(steps)); names=fieldnames(conditions);
            for k=1:numel(steps)
                p=obj.Data.step_parameters(k);
                for j=1:numel(names), keep(k)=keep(k)&&isfield(p,names{j})&&isequal(p.(names{j}),conditions.(names{j})); end
            end
            steps=steps(keep);
        end
        function out=export(obj,columns,step)
            if nargin<2||isempty(columns), columns=cellstr(obj.getTraceNames()); end
            if nargin<3, step=-1; end
            out=struct();
            for k=1:numel(columns)
                key=matlab.lang.makeValidName(char(columns{k}));
                if isequal(step,-1), out.(key)=obj.getTrace(columns{k});
                else, out.(key)=obj.getWave(columns{k},step); end
            end
        end
        function tableData=toTable(obj,columns,step)
            if nargin<2,columns=[];end, if nargin<3,step=-1;end
            values=obj.export(columns,step); tableData=struct2table(values);
        end
        function toCsv(obj,filename,columns,step,options)
            arguments, obj; filename (1,1) string; columns=[]; step=-1; options.Separator (1,1) string=","; end
            writetable(obj.toTable(columns,step),filename,'Delimiter',char(options.Separator));
        end
        function toExcel(obj,filename,columns,step)
            if nargin<3,columns=[];end, if nargin<4,step=-1;end, writetable(obj.toTable(columns,step),filename);
        end
        % PyLTSpice spelling aliases.
        function x=get_trace_names(obj),x=obj.getTraceNames();end
        function x=get_trace(obj,r),x=obj.getTraceObject(r);end
        function x=get_wave(obj,r,s),if nargin<3,s=0;end,x=obj.getWave(r,s);end
        function x=get_axis(obj,s),if nargin<2,s=0;end,x=obj.getAxis(s);end
        function x=get_time_axis(obj,s),if nargin<2,s=0;end,x=obj.getTime(s);end
        function x=get_len(obj,s),if nargin<2,s=0;end,x=obj.getLen(s);end
        function x=get_nr_plots(obj),x=obj.getNrPlots();end
        function x=get_plot_name(obj),x=obj.getPlotName();end
        function x=get_plot_names(obj),x=obj.getPlotNames();end
        function x=get_raw_properties(obj),x=obj.getRawProperties();end
        function x=get_raw_property(obj,n),if nargin<2,x=obj.getRawProperty();else,x=obj.getRawProperty(n);end,end
        function x=get_steps(obj,varargin),if isempty(varargin),x=obj.getSteps();else,x=obj.getSteps(varargin{1});end,end
        function x=aliases(obj),x=obj.Data.aliases;end
        function x=backannotations(obj),x=obj.Data.backannotations;end
        function to_csv(obj,varargin),obj.toCsv(varargin{:});end
        function to_excel(obj,varargin),obj.toExcel(varargin{:});end
        function x=to_dataframe(obj,varargin),x=obj.toTable(varargin{:});end
    end
    methods (Access=private)
        function j=traceIndex(obj,reference)
            if isnumeric(reference)
                if ~isscalar(reference)||~isfinite(reference)||fix(reference)~=reference
                    error("radia:ltspice:TraceNotFound", ...
                        "Trace index must be a finite integer scalar.");
                end
                % PyLTSpice integer trace references are zero-based.
                j=double(reference)+1;
            else
                j=find(strcmpi(obj.Data.names,string(reference)),1);
            end
            if isempty(j)||j<1||j>numel(obj.Data.names), error("radia:ltspice:TraceNotFound","Trace not found: %s",string(reference)); end
        end
        function wave=evaluateAlias(obj,name)
            index=find(strcmpi(obj.Data.alias_names,string(name)),1);
            if isempty(index)
                error("radia:ltspice:TraceNotFound","Trace not found: %s",name);
            end
            formula=char(obj.Data.alias_formulas(index));
            aliasValues={};
            [pairs,pairMatches]=regexp(formula, ...
                'V\(\s*([^,()]+?)\s*,\s*([^()]+?)\s*\)', ...
                'tokens','match');
            for k=1:numel(pairs)
                aliasValues{end+1}=obj.nodeWave(pairs{k}{1}) ...
                    -obj.nodeWave(pairs{k}{2}); %#ok<AGROW>
                formula=strrep(formula,pairMatches{k}, ...
                    sprintf('x%d',numel(aliasValues)));
            end
            singles=unique(string(regexp(formula,'[VIP]\([^()]+\)', ...
                'match')),'stable');
            for k=1:numel(singles)
                aliasValues{end+1}=obj.getTrace(singles(k)); %#ok<AGROW>
                formula=strrep(formula,char(singles(k)), ...
                    sprintf('x%d',numel(aliasValues)));
            end
            formula=regexprep(formula,'(?i)(mho|ohm)','');
            try
                wave=localSafeAliasFormula(formula,aliasValues);
            catch cause
                error("radia:ltspice:AliasEvaluation", ...
                    "Cannot evaluate alias %s: %s",name,cause.message);
            end
            if isscalar(wave)
                wave=repmat(wave,size(obj.Data.values,1),1);
            else
                wave=wave(:);
            end
            if numel(wave)~=size(obj.Data.values,1)
                error("radia:ltspice:AliasEvaluation", ...
                    "Alias %s produced %d values; expected %d.",name, ...
                    numel(wave),size(obj.Data.values,1));
            end
        end
        function wave=nodeWave(obj,node),node=strtrim(string(node));if node=="0",wave=zeros(size(obj.Data.values,1),1);else,wave=obj.getTrace("V("+node+")");end,end
    end
end
function parameters=localStepParameters(path,count)
text=string(fileread(path));tokens=regexp(text,'(?im)^\.step\s+([^\r\n]+)','tokens');parameters=repmat(struct(),1,count);
for k=1:min(count,numel(tokens)),pairs=regexp(tokens{k}{1},'([A-Za-z_]\w*)\s*=\s*([^\s]+)','tokens');for j=1:numel(pairs),key=matlab.lang.makeValidName(pairs{j}{1});number=str2double(pairs{j}{2});if isfinite(number),parameters(k).(key)=number;else,parameters(k).(key)=string(pairs{j}{2});end,end,end
end

function value=localSafeAliasFormula(formula,aliasValues)
% Evaluate arithmetic aliases without executing arbitrary RAW-file text.
compact=regexprep(char(formula),'\s+','');
tokenPattern='x\d+|pi|e|(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?|[()+\-*/^]';
tokens=regexp(compact,tokenPattern,'match');
if isempty(tokens)||~strcmp(strjoin(tokens,''),compact)
    error("radia:ltspice:AliasSyntax", ...
        "Alias contains an unsupported token.");
end
position=1;
value=parseExpression();
if position<=numel(tokens)
    error("radia:ltspice:AliasSyntax","Unexpected token: %s",tokens{position});
end

    function result=parseExpression()
        result=parseTerm();
        while position<=numel(tokens)&&any(strcmp(tokens{position},{'+','-'}))
            operation=tokens{position};position=position+1;rhs=parseTerm();
            if operation=="+",result=result+rhs;else,result=result-rhs;end
        end
    end

    function result=parseTerm()
        result=parsePower();
        while position<=numel(tokens)&&any(strcmp(tokens{position},{'*','/'}))
            operation=tokens{position};position=position+1;rhs=parsePower();
            if operation=="*",result=result.*rhs;else,result=result./rhs;end
        end
    end

    function result=parsePower()
        result=parseUnary();
        if position<=numel(tokens)&&strcmp(tokens{position},'^')
            position=position+1;result=result.^parsePower();
        end
    end

    function result=parseUnary()
        if position<=numel(tokens)&&strcmp(tokens{position},'+')
            position=position+1;result=parseUnary();return
        end
        if position<=numel(tokens)&&strcmp(tokens{position},'-')
            position=position+1;result=-parseUnary();return
        end
        result=parsePrimary();
    end

    function result=parsePrimary()
        if position>numel(tokens)
            error("radia:ltspice:AliasSyntax","Unexpected end of alias.");
        end
        token=tokens{position};position=position+1;
        if strcmp(token,'(')
            result=parseExpression();
            if position>numel(tokens)||~strcmp(tokens{position},')')
                error("radia:ltspice:AliasSyntax","Missing closing parenthesis.");
            end
            position=position+1;return
        end
        if startsWith(token,'x')
            index=str2double(token(2:end));
            if ~isfinite(index)||index<1||index>numel(aliasValues)
                error("radia:ltspice:AliasSyntax","Invalid trace reference.");
            end
            result=aliasValues{index};return
        end
        if strcmp(token,'pi'),result=pi;return,end
        if strcmp(token,'e'),result=exp(1);return,end
        result=str2double(token);
        if ~isfinite(result)
            error("radia:ltspice:AliasSyntax","Invalid numeric literal: %s",token);
        end
    end
end
