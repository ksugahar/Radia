classdef SchematicEditor < handle
    %SCHEMATICEDITOR Edit component values and directives in LTspice ASC text.
    properties(SetAccess=private), SourcePath (1,1) string; Text (1,1) string; LibraryPaths (1,:) string=strings(1,0); end
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
        function value=getComponentValue(obj,reference),[lines,target,stop]=obj.componentBlock(reference);row=target+find(startsWith(strtrim(lines(target+1:stop)),"SYMATTR Value "),1);if isempty(row),value="";else,value=strtrim(extractAfter(strtrim(lines(row)),"SYMATTR Value "));end,end
        function value=getComponentAttribute(obj,reference,attribute)
            [lines,target,stop]=obj.componentBlock(reference);needle="SYMATTR "+string(attribute)+" ";row=target+find(startsWith(strtrim(lines(target+1:stop)),needle),1);
            if isempty(row),value="";else,value=strtrim(extractAfter(strtrim(lines(row)),needle));end
        end
        function setComponentAttribute(obj,reference,attribute,value)
            [lines,target,stop]=obj.componentBlock(reference);needle="SYMATTR "+string(attribute)+" ";relative=find(startsWith(strtrim(lines(target+1:stop)),needle),1);
            if isempty(relative),lines=[lines(1:target);needle+string(value);lines(target+1:end)];else,lines(target+relative)=needle+string(value);end,obj.Text=join(lines,newline);
        end
        function [position,rotation]=getComponentPosition(obj,reference)
            [lines,target]=obj.componentBlock(reference);symbol=find(startsWith(strtrim(lines(1:target)),"SYMBOL "),1,'last');parts=split(strtrim(lines(symbol)));position=[str2double(parts(3)),str2double(parts(4))];rotation=parts(5);
        end
        function setComponentPosition(obj,reference,position,rotation)
            if nargin<4,rotation="R0";end,[lines,target]=obj.componentBlock(reference);symbol=find(startsWith(strtrim(lines(1:target)),"SYMBOL "),1,'last');parts=split(strtrim(lines(symbol)));parts(3)=string(round(position(1)));parts(4)=string(round(position(2)));parts(5)=string(rotation);lines(symbol)=join(parts," ");obj.Text=join(lines,newline);
        end
        function references=getComponents(obj,prefixes)
            if nargin<2,prefixes="*";end,lines=splitlines(obj.Text);rows=startsWith(strtrim(lines),"SYMATTR InstName ");references=strtrim(extractAfter(strtrim(lines(rows)),"SYMATTR InstName "));if prefixes~="*",references=references(contains(upper(prefixes),upper(extractBefore(references,2))));end
        end
        function component=getComponent(obj,reference),[position,rotation]=obj.getComponentPosition(reference);component=struct('reference',string(reference),'value',obj.getComponentValue(reference),'position',position,'rotation',rotation);end
        function value=getComponentFloatValue(obj,reference),value=localEngineeringNumber(obj.getComponentValue(reference));end
        function nodes=getComponentNodes(obj,reference)
            temporary=fullfile("C:\temp","radia_asc_nodes_"+string(char(java.util.UUID.randomUUID()))+".asc");cleanup=onCleanup(@()localDelete(temporary));obj.saveAs(temporary);converted=radia.ltspice.schematicToNetlist(temporary,OutputDirectory=fileparts(temporary));net=radia.ltspice.SpiceEditor(converted.netlist);nodes=net.getComponentNodes(reference);clear cleanup
        end
        function parameters=getComponentParameters(obj,reference),parameters=struct();line=obj.getComponentAttribute(reference,"SpiceLine");pairs=regexp(line,'([A-Za-z_]\w*)\s*=\s*([^\s]+)','tokens');for k=1:numel(pairs),parameters.(matlab.lang.makeValidName(pairs{k}{1}))=string(pairs{k}{2});end,end
        function setComponentParameters(obj,reference,values),names=fieldnames(values);line="";for k=1:numel(names),line=line+names{k}+"="+string(values.(names{k}))+" ";end,obj.setComponentAttribute(reference,"SpiceLine",strtrim(line));end
        function removeComponent(obj,reference),[lines,target,stop]=obj.componentBlock(reference);symbol=find(startsWith(strtrim(lines(1:target)),"SYMBOL "),1,'last');lines(symbol:stop)=[];obj.Text=join(lines,newline);end
        function addComponent(obj,symbol,reference,position,rotation,value)
            arguments,obj;symbol (1,1) string;reference (1,1) string;position (1,2) double;rotation (1,1) string="R0";value (1,1) string="";end
            if any(obj.getComponents()==reference),error("radia:ltspice:DuplicateComponent","Component already exists: %s",reference);end
            addition=sprintf("SYMBOL %s %d %d %s\nSYMATTR InstName %s",symbol,round(position(1)),round(position(2)),rotation,reference);if strlength(value)>0,addition=addition+newline+"SYMATTR Value "+value;end,obj.Text=obj.Text+newline+addition;
        end
        function addLibraryPaths(obj,varargin),obj.LibraryPaths=[obj.LibraryPaths,string(varargin)];end
        function permission=beginUpdate(~),permission="allow";end
        function endUpdate(varargin),end
        function copyFrom(obj,editor),obj.Text=editor.Text;end
        function updated(~),end
        function prepareForSimulator(~,varargin),end
        function removed=removeXInstruction(obj,pattern),lines=splitlines(obj.Text);mask=startsWith(strtrim(lines),"TEXT ")&contains(lines,"!")&~cellfun(@isempty,regexp(cellstr(lines),char(pattern),'once'));removed=any(mask);lines(mask)=[];obj.Text=join(lines,newline);end
        function setElementModel(obj,reference,model),obj.setComponentValue(reference,model);end
        function addWire(obj,p1,p2),obj.Text=obj.Text+newline+sprintf("WIRE %d %d %d %d",round(p1(1)),round(p1(2)),round(p2(1)),round(p2(2)));end
        function scale(obj,offsetX,offsetY,scaleX,scaleY)
            lines=splitlines(obj.Text);pattern='^(WIRE|FLAG|SYMBOL|TEXT)\s+(.+)$';
            for k=1:numel(lines),t=regexp(char(strtrim(lines(k))),pattern,'tokens','once');if isempty(t),continue,end,parts=split(string(t{2}));switch string(t{1}),case "WIRE",idx=1:4;otherwise,idx=1:2;end,for q=idx,number=str2double(parts(q));if mod(q,2)==1,number=offsetX+scaleX*number;else,number=offsetY+scaleY*number;end,parts(q)=string(round(number));end,lines(k)=string(t{1})+" "+join(parts," ");end,obj.Text=join(lines,newline);
        end
        function addDirective(obj,directive,x,y)
            arguments, obj; directive (1,1) string; x (1,1) double=64; y (1,1) double=400; end
            if ~startsWith(directive,"."), error("radia:ltspice:Directive","Directive must start with '.'."); end
            obj.Text=obj.Text+newline+sprintf("TEXT %d %d Left 2 !%s",round(x),round(y),directive);
        end
        function names=getAllParameterNames(obj),t=regexp(obj.Text,'(?im)^TEXT\s+\d+\s+\d+\s+\S+\s+\d+\s+!\.param\s+([A-Za-z_]\w*)\s*=','tokens');names=sort(string(cellfun(@(x)x{1},t,'UniformOutput',false)));end
        function value=getParameter(obj,name),p="(?im)^TEXT\s+\d+\s+\d+\s+\S+\s+\d+\s+!\.param\s+"+regexptranslate('escape',char(name))+"\s*=\s*([^\s;]+)";t=regexp(obj.Text,p,'tokens','once');if isempty(t),error("radia:ltspice:ParameterNotFound","Parameter not found: %s",name);end,value=string(t{1});end
        function setParameter(obj,name,value)
            text=localValue(value);p="(?im)(^TEXT\s+\d+\s+\d+\s+\S+\s+\d+\s+!\.param\s+"+regexptranslate('escape',char(name))+"\s*=\s*)([^\s;]+)";
            if isempty(regexp(obj.Text,p,'once')),obj.addDirective(".param "+name+"="+text);else,obj.Text=regexprep(obj.Text,p,"$1"+text,'once');end
        end
        function removed=removeInstruction(obj,instruction),p="(?im)^TEXT\s+\d+\s+\d+\s+\S+\s+\d+\s+!"+regexptranslate('escape',char(instruction))+"\s*(?:\r?\n)?";removed=~isempty(regexp(obj.Text,p,'once'));if removed,obj.Text=regexprep(obj.Text,p,"",'once');end,end
        function saveAs(obj,path)
            folder=fileparts(path); if strlength(folder)>0&&~isfolder(folder), mkdir(folder); end
            f=fopen(path,'w'); if f<0,error("radia:ltspice:Write","Cannot write %s",path);end
            c=onCleanup(@()fclose(f)); fprintf(f,'%s',obj.Text); clear c
        end
        function x=get_component_value(obj,varargin),x=obj.getComponentValue(varargin{:});end
        function x=get_component(obj,varargin),x=obj.getComponent(varargin{:});end
        function x=get_component_floatvalue(obj,varargin),x=obj.getComponentFloatValue(varargin{:});end
        function x=get_component_nodes(obj,varargin),x=obj.getComponentNodes(varargin{:});end
        function x=get_component_parameters(obj,varargin),x=obj.getComponentParameters(varargin{:});end
        function x=get_component_attribute(obj,varargin),x=obj.getComponentAttribute(varargin{:});end
        function set_component_attribute(obj,varargin),obj.setComponentAttribute(varargin{:});end
        function [p,r]=get_component_position(obj,varargin),[p,r]=obj.getComponentPosition(varargin{:});end
        function set_component_position(obj,varargin),obj.setComponentPosition(varargin{:});end
        function add_component(obj,varargin),obj.addComponent(varargin{:});end
        function x=get_components(obj,varargin),x=obj.getComponents(varargin{:});end
        function set_component_value(obj,varargin),obj.setComponentValue(varargin{:});end
        function set_component_values(obj,values),n=fieldnames(values);for k=1:numel(n),obj.setComponentValue(n{k},values.(n{k}));end,end
        function set_component_parameters(obj,varargin),obj.setComponentParameters(varargin{:});end
        function x=get_parameter(obj,varargin),x=obj.getParameter(varargin{:});end
        function x=get_all_parameter_names(obj),x=obj.getAllParameterNames();end
        function set_parameter(obj,varargin),obj.setParameter(varargin{:});end
        function set_parameters(obj,values),n=fieldnames(values);for k=1:numel(n),obj.setParameter(n{k},values.(n{k}));end,end
        function add_instruction(obj,instruction),obj.addDirective(instruction);end
        function add_instructions(obj,varargin),for k=1:numel(varargin),obj.addDirective(string(varargin{k}));end,end
        function add_library_paths(obj,varargin),obj.addLibraryPaths(varargin{:});end
        function set_custom_library_paths(obj,varargin),obj.LibraryPaths=string(varargin);end
        function x=begin_update(obj),x=obj.beginUpdate();end
        function end_update(obj,varargin),obj.endUpdate(varargin{:});end
        function copy_from(obj,varargin),obj.copyFrom(varargin{:});end
        function prepare_for_simulator(obj,varargin),obj.prepareForSimulator(varargin{:});end
        function x=remove_Xinstruction(obj,varargin),x=obj.removeXInstruction(varargin{:});end
        function set_element_model(obj,varargin),obj.setElementModel(varargin{:});end
        function x=remove_instruction(obj,varargin),x=obj.removeInstruction(varargin{:});end
        function remove_component(obj,varargin),obj.removeComponent(varargin{:});end
        function save_as(obj,varargin),obj.saveAs(varargin{:});end
        function save_netlist(obj,varargin),obj.saveAs(varargin{:});end
        function write_netlist(obj,varargin),obj.saveAs(varargin{:});end
        function x=is_read_only(~),x=false;end
        function x=reset_netlist(obj),obj.Text=string(fileread(obj.SourcePath));x=true;end
    end
    methods (Access=private)
        function [lines,target,stop]=componentBlock(obj,reference)
            lines=splitlines(obj.Text);names=find(startsWith(strtrim(lines),"SYMATTR InstName "));target=[];for k=names',if strtrim(extractAfter(strtrim(lines(k)),"SYMATTR InstName "))==reference,target=k;break,end,end
            if isempty(target),error("radia:ltspice:ComponentNotFound","Component not found: %s",reference);end,next=find(startsWith(strtrim(lines(target+1:end)),"SYMBOL "),1);if isempty(next),stop=numel(lines);else,stop=target+next-1;end
        end
    end
end
function text=localValue(value)
if isnumeric(value)&&isscalar(value)&&isfinite(value)
    text=string(sprintf('%.17g',value));
elseif (isstring(value)&&isscalar(value))||(ischar(value)&&isrow(value))
    text=string(value);
else
    error("radia:ltspice:Value","Value must be finite scalar or text.");
end
end
function value=localEngineeringNumber(text),s=lower(char(text));t=regexp(s,'^([-+0-9.eE]+)(meg|[tgkmunpf]?)','tokens','once');if isempty(t),value=NaN;return,end,value=str2double(t{1});factor=struct('t',1e12,'g',1e9,'meg',1e6,'k',1e3,'m',1e-3,'u',1e-6,'n',1e-9,'p',1e-12,'f',1e-15);if ~isempty(t{2}),value=value*factor.(t{2});end,end
function localDelete(path),if isfile(path),delete(path);end,[folder,name]=fileparts(path);net=fullfile(folder,name+".net");if isfile(net),delete(net);end,end
