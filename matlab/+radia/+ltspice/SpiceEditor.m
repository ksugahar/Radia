classdef SpiceEditor < handle
    %SPICEEDITOR Edit LTspice netlists using a MATLAB-native API.
    properties (SetAccess=private)
        SourcePath (1,1) string
        Text (1,1) string
    end
    methods
        function obj=SpiceEditor(path)
            arguments, path (1,1) string {mustBeFile}, end
            obj.SourcePath=path; obj.Text=string(fileread(path));
        end
        function setParameter(obj,name,value)
            arguments, obj; name (1,1) string; value, end
            v=radia.ltspice.SpiceEditor.formatValue(value);
            p="(?im)(^\s*\.param\s+"+regexptranslate('escape',char(name))+"\s*=\s*)([^\s;]+)";
            if isempty(regexp(obj.Text,p,'once')), error("radia:ltspice:ParameterNotFound","Parameter not found: %s",name); end
            obj.Text=regexprep(obj.Text,p,"$1"+v);
        end
        function setComponentValue(obj,reference,value)
            arguments, obj; reference (1,1) string; value, end
            v=radia.ltspice.SpiceEditor.formatValue(value);
            p="(?im)(^\s*"+regexptranslate('escape',char(reference))+"\s+\S+(?:\s+\S+)+?\s+)(\S+)(\s*(?:;.*)?$)";
            if isempty(regexp(obj.Text,p,'once')), error("radia:ltspice:ComponentNotFound","Component not found: %s",reference); end
            obj.Text=regexprep(obj.Text,p,"$1"+v+"$3",'once');
        end
        function addInstruction(obj,instruction)
            arguments, obj; instruction (1,1) string, end
            obj.Text=regexprep(obj.Text,"(?im)^\s*\.end\s*$",instruction+newline+".end",'once');
        end
        function saveAs(obj,path)
            arguments, obj; path (1,1) string, end
            folder=fileparts(path); if strlength(folder)>0 && ~isfolder(folder), mkdir(folder); end
            f=fopen(path,'w'); if f<0, error("radia:ltspice:Write","Cannot write %s",path); end
            c=onCleanup(@()fclose(f)); fprintf(f,'%s',obj.Text); clear c
        end
    end
    methods (Static,Access=private)
        function v=formatValue(value)
            if isnumeric(value)&&isscalar(value)&&isfinite(value), v=string(sprintf('%.17g',value));
            elseif (isstring(value)&&isscalar(value))||(ischar(value)&&isrow(value)), v=string(value);
            else, error("radia:ltspice:Value","Value must be finite scalar or text."); end
        end
    end
end
