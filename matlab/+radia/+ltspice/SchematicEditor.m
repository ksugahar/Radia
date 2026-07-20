classdef SchematicEditor < handle
    %SCHEMATICEDITOR Edit component values and directives in LTspice ASC text.
    properties(SetAccess=private), SourcePath (1,1) string; Text (1,1) string; end
    methods
        function obj=SchematicEditor(path)
            arguments, path (1,1) string {mustBeFile}, end
            [~,~,extension]=fileparts(path);
            if lower(string(extension))~=".asc", error("radia:ltspice:SchematicRequired","SchematicEditor requires .asc."); end
            obj.SourcePath=path; obj.Text=string(fileread(path));
        end
        function setComponentValue(obj,reference,value)
            value=localValue(value); lines=splitlines(obj.Text); names=find(startsWith(strtrim(lines),"SYMATTR InstName "));
            target=[];
            for k=names', if strtrim(extractAfter(strtrim(lines(k)),"SYMATTR InstName "))==reference, target=k; break, end, end
            if isempty(target), error("radia:ltspice:ComponentNotFound","Component not found: %s",reference); end
            nextSymbol=find(startsWith(strtrim(lines(target+1:end)),"SYMBOL "),1); if isempty(nextSymbol), stop=numel(lines); else, stop=target+nextSymbol-1; end
            valueRow=target+find(startsWith(strtrim(lines(target+1:stop)),"SYMATTR Value "),1);
            if isempty(valueRow), lines=[lines(1:target);"SYMATTR Value "+value;lines(target+1:end)]; else, lines(valueRow)="SYMATTR Value "+value; end
            obj.Text=join(lines,newline);
        end
        function addDirective(obj,directive,x,y)
            arguments, obj; directive (1,1) string; x (1,1) double=64; y (1,1) double=400; end
            if ~startsWith(directive,"."), error("radia:ltspice:Directive","Directive must start with '.'."); end
            obj.Text=obj.Text+newline+sprintf("TEXT %d %d Left 2 !%s",round(x),round(y),directive);
        end
        function saveAs(obj,path)
            folder=fileparts(path); if strlength(folder)>0&&~isfolder(folder), mkdir(folder); end
            f=fopen(path,'w'); if f<0,error("radia:ltspice:Write","Cannot write %s",path);end
            c=onCleanup(@()fclose(f)); fprintf(f,'%s',obj.Text); clear c
        end
    end
end
function text=localValue(value)
if isnumeric(value)&&isscalar(value)&&isfinite(value),text=string(sprintf('%.17g',value));
elseif (isstring(value)&&isscalar(value))||(ischar(value)&&isrow(value)),text=string(value);
else,error("radia:ltspice:Value","Value must be finite scalar or text.");end
end
