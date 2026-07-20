classdef LogReader
    %LOGREADER Parse LTspice log text and .meas scalar results.
    properties (SetAccess=private), Path (1,1) string; Text (1,1) string; Measures table, end
    methods
        function obj=LogReader(path)
            arguments, path (1,1) string {mustBeFile}, end
            obj.Path=path; obj.Text=string(fileread(path)); obj.Measures=obj.parseMeasures();
        end
        function value=getMeasure(obj,name)
            row=strcmpi(obj.Measures.Name,string(name)); if ~any(row), error("radia:ltspice:MeasureNotFound","Measure not found: %s",name); end
            value=obj.Measures.Value(find(row,1));
        end
    end
    methods (Access=private)
        function out=parseMeasures(obj)
            rows=regexp(char(obj.Text),'(?m)^\s*([A-Za-z_]\w*)\s*(?::[^\r\n=]*)?=\s*([-+0-9.eE]+)','tokens');
            names=strings(0,1); values=zeros(0,1);
            for k=1:numel(rows), v=str2double(rows{k}{2}); if isfinite(v), names(end+1,1)=string(rows{k}{1}); values(end+1,1)=v; end, end %#ok<AGROW>
            out=table(names,values,'VariableNames',{'Name','Value'});
        end
    end
end
