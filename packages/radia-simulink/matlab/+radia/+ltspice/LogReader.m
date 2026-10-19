classdef LogReader
    %LOGREADER Parse LTspice log text and .meas scalar results.
    properties (SetAccess=private), Path (1,1) string; Text (1,1) string; Measures table; StepParameters table; end
    methods
        function obj=LogReader(path)
            arguments, path (1,1) string {mustBeFile}, end
            obj.Path=path; obj.Text=string(fileread(path)); obj.Measures=obj.parseMeasures();obj.StepParameters=obj.parseSteps();
        end
        function value=getMeasure(obj,name)
            row=strcmpi(obj.Measures.Name,string(name)); if ~any(row), error("radia:ltspice:MeasureNotFound","Measure not found: %s",name); end
            value=obj.Measures.Value(find(row,1));
        end
        function names=getMeasureNames(obj),names=obj.Measures.Name;end
        function value=getMeasureValue(obj,name,step)
            rows=find(strcmpi(obj.Measures.Name,string(name)));if isempty(rows),error("radia:ltspice:MeasureNotFound","Measure not found: %s",name);end
            if nargin<3||isempty(step),step=0;end,index=double(step)+1;if index>numel(rows),value=[];else,value=obj.Measures.Value(rows(index));end
        end
        function values=getMeasureValuesAtSteps(obj,name,steps)
            if nargin<3||isempty(steps),steps=0:sum(strcmpi(obj.Measures.Name,string(name)))-1;end
            values=arrayfun(@(s)obj.getMeasureValue(name,s),steps);
        end
        function value=avgMeasureValue(obj,name,steps),if nargin<3,steps=[];end,value=mean(obj.getMeasureValuesAtSteps(name,steps));end
        function value=minMeasureValue(obj,name,steps),if nargin<3,steps=[];end,value=min(obj.getMeasureValuesAtSteps(name,steps));end
        function value=maxMeasureValue(obj,name,steps),if nargin<3,steps=[];end,value=max(obj.getMeasureValuesAtSteps(name,steps));end
        function answer=hasSteps(obj),answer=height(obj.StepParameters)>0;end
        function names=getStepVars(obj),names=string(obj.StepParameters.Properties.VariableNames);end
        function steps=stepsWithConditions(obj,conditions)
            arguments,obj;conditions (1,1) struct;end,steps=(0:height(obj.StepParameters)-1)';names=fieldnames(conditions);keep=true(size(steps));
            for k=1:numel(names),key=matlab.lang.makeValidName(names{k});if ~ismember(key,obj.StepParameters.Properties.VariableNames),keep(:)=false;break,end,column=obj.StepParameters.(key);target=conditions.(names{k});if isnumeric(column),keep=keep&(column==target);else,keep=keep&strcmpi(string(column),string(target));end,end,steps=steps(keep)';
        end
        function steps=stepsWithParameterEqualTo(obj,name,value),s=struct();s.(matlab.lang.makeValidName(name))=value;steps=obj.stepsWithConditions(s);end
        function exportData(obj,path),writetable(obj.Measures,path);end
        function obtainAmplitudeAndPhaseFromComplexValues(obj)
            names=obj.Measures.Properties.VariableNames;for k=1:numel(names),column=obj.Measures.(names{k});if isnumeric(column)&&~isreal(column),obj.Measures.(names{k}+"_mag")=abs(column);obj.Measures.(names{k}+"_ph")=rad2deg(angle(column));end,end
        end
        function chart=plotHistogram(obj,name,steps,bins,normalized,varargin)
            if nargin<3,steps=[];end,if nargin<4,bins=50;end,if nargin<5,normalized=true;end,values=obj.getMeasureValuesAtSteps(name,steps);if normalized,mode='pdf';else,mode='count';end,chart=histogram(values,bins,'Normalization',mode);grid on
        end
        function x=get_measure_names(obj),x=obj.getMeasureNames();end
        function x=get_measure_value(obj,varargin),x=obj.getMeasureValue(varargin{:});end
        function x=get_measure_values_at_steps(obj,varargin),x=obj.getMeasureValuesAtSteps(varargin{:});end
        function x=avg_measure_value(obj,varargin),x=obj.avgMeasureValue(varargin{:});end
        function x=min_measure_value(obj,varargin),x=obj.minMeasureValue(varargin{:});end
        function x=max_measure_value(obj,varargin),x=obj.maxMeasureValue(varargin{:});end
        function x=has_steps(obj),x=obj.hasSteps();end
        function x=get_step_vars(obj),x=obj.getStepVars();end
        function x=steps_with_conditions(obj,varargin),x=obj.stepsWithConditions(varargin{:});end
        function x=steps_with_parameter_equal_to(obj,varargin),x=obj.stepsWithParameterEqualTo(varargin{:});end
        function export_data(obj,varargin),obj.exportData(varargin{:});end
        function obtain_amplitude_and_phase_from_complex_values(obj),obj.obtainAmplitudeAndPhaseFromComplexValues();end
        function split_complex_values_on_datasets(obj),obj.obtainAmplitudeAndPhaseFromComplexValues();end
        function x=plot_histogram(obj,varargin),x=obj.plotHistogram(varargin{:});end
    end
    methods (Access=private)
        function out=parseMeasures(obj)
            rows=regexp(char(obj.Text),'(?m)^\s*([A-Za-z_]\w*)\s*(?::[^\r\n=]*)?=\s*([-+0-9.eE]+)','tokens');
            names=strings(0,1); values=zeros(0,1);
            for k=1:numel(rows), v=str2double(rows{k}{2}); if isfinite(v), names(end+1,1)=string(rows{k}{1}); values(end+1,1)=v; end, end %#ok<AGROW>
            out=table(names,values,'VariableNames',{'Name','Value'});
        end
        function out=parseSteps(obj)
            tokens=regexp(obj.Text,'(?im)^\.step\s+([^\r\n]+)','tokens');if isempty(tokens),out=table();return,end
            rows=cell(numel(tokens),1);allNames=strings(0,1);
            for k=1:numel(tokens),pairs=regexp(tokens{k}{1},'([A-Za-z_]\w*)\s*=\s*([^\s]+)','tokens');row=struct();for j=1:numel(pairs),key=matlab.lang.makeValidName(pairs{j}{1});number=str2double(pairs{j}{2});if isfinite(number),row.(key)=number;else,row.(key)=string(pairs{j}{2});end,allNames(end+1,1)=string(key);end,rows{k}=row;end %#ok<AGROW>
            allNames=unique(allNames,'stable');out=table();for j=1:numel(allNames),name=allNames(j);values=cell(numel(rows),1);for k=1:numel(rows),if isfield(rows{k},name),values{k}=rows{k}.(name);else,values{k}=missing;end,end,if all(cellfun(@(x)isnumeric(x)&&isscalar(x),values)),out.(name)=cell2mat(values);else,out.(name)=string(values);end,end
        end
    end
end
