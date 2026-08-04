classdef SpiceEditor < handle
    %SPICEEDITOR Edit LTspice netlists using a MATLAB-native API.
    properties (SetAccess=private)
        SourcePath (1,1) string
        Text (1,1) string
        LibraryPaths (1,:) string=strings(1,0)
        CircuitName (1,1) string=""
    end
    methods
        function obj=SpiceEditor(path)
            arguments, path (1,1) string {mustBeFile}, end
            obj.SourcePath=path; obj.Text=string(fileread(path));[~,obj.CircuitName]=fileparts(path);
        end
        function setParameter(obj,name,value)
            arguments, obj; name (1,1) string; value, end
            v=radia.ltspice.SpiceEditor.formatValue(value);
            p="(?im)(^\s*\.param\s+"+regexptranslate('escape',char(name))+"\s*=\s*)([^\s;]+)";
            if isempty(regexp(obj.Text,p,'once')), error("radia:ltspice:ParameterNotFound","Parameter not found: %s",name); end
            obj.Text=regexprep(obj.Text,p,"$1"+v);
        end
        function value=getParameter(obj,name)
            p="(?im)^\s*\.param\s+"+regexptranslate('escape',char(name))+"\s*=\s*([^\s;]+)"; t=regexp(obj.Text,p,'tokens','once');
            if isempty(t),error("radia:ltspice:ParameterNotFound","Parameter not found: %s",name);end,value=string(t{1});
        end
        function names=getAllParameterNames(obj)
            t=regexp(obj.Text,'(?im)^\s*\.param\s+([A-Za-z_]\w*)\s*=','tokens'); names=string(cellfun(@(x)x{1},t,'UniformOutput',false));
        end
        function references=getComponents(obj,prefixes)
            if nargin<2,prefixes="*";end
            lines=splitlines(obj.Text); references=strings(0,1);
            for k=1:numel(lines),s=strtrim(lines(k));if strlength(s)==0||startsWith(s,[".";"*";";"]),continue,end,r=extractBefore(s," ");if prefixes=="*"||contains(upper(prefixes),upper(extractBefore(r,2))),references(end+1,1)=r;end,end %#ok<AGROW>
        end
        function value=getComponentValue(obj,reference)
            [parts,valueIndex]=obj.componentTokens(reference);value=parts(valueIndex);
        end
        function component=getComponent(obj,reference),component=struct('reference',string(reference),'nodes',obj.getComponentNodes(reference),'value',obj.getComponentValue(reference),'line',obj.componentLine(reference));end
        function value=getComponentAttribute(obj,reference,attribute)
            line=obj.componentLine(reference);t=regexp(line,"(?i)(?:^|\s)"+regexptranslate('escape',char(attribute))+"\s*=\s*([^\s]+)",'tokens','once');if isempty(t),value="";else,value=string(t{1});end
        end
        function setComponentAttribute(obj,reference,attribute,value)
            line=obj.componentLine(reference);p="(?i)(\s"+regexptranslate('escape',char(attribute))+"\s*=\s*)([^\s]+)";if isempty(regexp(line,p,'once')),replacement=line+" "+attribute+"="+string(value);else,replacement=regexprep(line,p,"$1"+string(value),'once');end,obj.Text=replace(obj.Text,line,replacement);
        end
        function parameters=getComponentParameters(obj,reference)
            line=obj.componentLine(reference);pairs=regexp(line,'\s([A-Za-z_]\w*)\s*=\s*([^\s]+)','tokens');parameters=struct();for k=1:numel(pairs),parameters.(matlab.lang.makeValidName(pairs{k}{1}))=string(pairs{k}{2});end
        end
        function setComponentParameters(obj,reference,values),names=fieldnames(values);for k=1:numel(names),obj.setComponentAttribute(reference,names{k},values.(names{k}));end,end
        function value=getComponentFloatValue(obj,reference),value=localEngineeringNumber(obj.getComponentValue(reference));end
        function nodes=getComponentNodes(obj,reference)
            [parts,valueIndex]=obj.componentTokens(reference);if valueIndex<=2,nodes=strings(0,1);else,nodes=parts(2:valueIndex-1);end
        end
        function removeComponent(obj,reference)
            p="(?im)^\s*"+regexptranslate('escape',char(reference))+"\s+[^\r\n]*(?:\r?\n)?";
            if isempty(regexp(obj.Text,p,'once')),error("radia:ltspice:ComponentNotFound","Component not found: %s",reference);end,obj.Text=regexprep(obj.Text,p,"",'once');
        end
        function addComponent(obj,reference,nodes,value,parameters)
            arguments,obj;reference (1,1) string;nodes (:,1) string;value (1,1) string;parameters (1,1) struct=struct();end
            if any(obj.getComponents()==reference),error("radia:ltspice:DuplicateComponent","Component already exists: %s",reference);end,line=reference+" "+join(nodes," ")+" "+value;names=fieldnames(parameters);for k=1:numel(names),line=line+" "+names{k}+"="+string(parameters.(names{k}));end,obj.Text=regexprep(obj.Text,"(?im)^\s*\.end\s*$",line+newline+".end",'once');
        end
        function setElementModel(obj,reference,model),obj.setComponentValue(reference,model);end
        function nodes=getAllNodes(obj),references=obj.getComponents();nodes=strings(0,1);for k=1:numel(references),nodes=[nodes;obj.getComponentNodes(references(k))];end,nodes=unique(nodes,'stable');end %#ok<AGROW>
        function names=getSubcircuitNames(obj),t=regexp(obj.Text,'(?im)^\s*\.subckt\s+(\S+)','tokens');names=string(cellfun(@(x)x{1},t,'UniformOutput',false));end
        function sections=getControlSections(obj),sections=string(regexp(obj.Text,'(?ims)^\s*\.control\s*$.*?^\s*\.endc\s*$','match'))';end
        function addControlSection(obj,instruction),obj.Text=regexprep(obj.Text,"(?im)^\s*\.end\s*$",".control"+newline+instruction+newline+".endc"+newline+".end",'once');end
        function removed=removeControlSection(obj,index),if nargin<2,index=1;end,sections=obj.getControlSections();if index<1||index>numel(sections),removed=false;else,obj.Text=replace(obj.Text,sections(index),"");removed=true;end,end
        function addLibrarySearchPaths(obj,varargin),obj.LibraryPaths=[obj.LibraryPaths,string(varargin)];end
        function path=findLibrary(obj,name),candidates=[fullfile(fileparts(obj.SourcePath),name),fullfile(obj.LibraryPaths,name)];index=find(isfile(candidates),1);if isempty(index),path="";else,path=candidates(index);end,end
        function circuit=getSubcircuitNamed(obj,name)
            p="(?ims)^\s*\.subckt\s+"+regexptranslate('escape',char(name))+"(?:\s+[^\r\n]*)?\r?\n.*?^\s*\.ends(?:\s+"+regexptranslate('escape',char(name))+")?\s*$";block=string(regexp(obj.Text,p,'match','once'));if strlength(block)==0,circuit=[];return,end,path=fullfile("C:\temp","radia_subckt_"+matlab.lang.makeValidName(char(name))+".cir");f=fopen(path,'w');c=onCleanup(@()fclose(f));fprintf(f,'%s\n.end',block);clear c;circuit=radia.ltspice.SpiceCircuit(path);
        end
        function circuit=getSubcircuit(obj,instanceName),line=obj.componentLine(instanceName);parts=split(strtrim(line));nonParameter=find(~contains(parts,"="));model=parts(nonParameter(end));circuit=obj.getSubcircuitNamed(model);if isempty(circuit),error("radia:ltspice:SubcircuitNotFound","Subcircuit not found: %s",model);end,end
        function circuits=modifiedSubcircuits(~),circuits={};end
        function value=name(obj),value=obj.CircuitName;end
        function setname(obj,value),obj.CircuitName=string(value);end
        function permission=beginUpdate(~),permission="allow";end
        function endUpdate(varargin),end
        function cloned=clone(obj,varargin),path=fullfile("C:\temp","radia_spice_clone_"+string(char(java.util.UUID.randomUUID()))+".cir");obj.saveAs(path);cloned=radia.ltspice.SpiceCircuit(path);end
        function lineClass=classForInstruction(~,instruction,varargin),lineClass=upper(extractBefore(strtrim(string(instruction))+" "," "));end
        function removed=removeXInstruction(obj,pattern),lines=splitlines(obj.Text);mask=startsWith(upper(strtrim(lines)),"X")&~cellfun(@isempty,regexp(cellstr(lines),char(pattern),'once'));removed=any(mask);lines(mask)=[];obj.Text=join(lines,newline);end
        function prepareForSimulator(~,varargin),end
        function writeLines(obj,stream),fprintf(stream,'%s',obj.Text);end
        function setComponentValue(obj,reference,value)
            arguments, obj; reference (1,1) string; value, end
            v=radia.ltspice.SpiceEditor.formatValue(value);
            line=obj.componentLine(reference);[parts,index]=obj.componentTokens(reference);parts(index)=v;obj.Text=replace(obj.Text,line,join(parts," "));
        end
        function addInstruction(obj,instruction)
            arguments, obj; instruction (1,1) string, end
            obj.Text=regexprep(obj.Text,"(?im)^\s*\.end\s*$",instruction+newline+".end",'once');
        end
        function removed=removeInstruction(obj,instruction)
            p="(?im)^\s*"+regexptranslate('escape',char(instruction))+"\s*(?:\r?\n)?";removed=~isempty(regexp(obj.Text,p,'once'));if removed,obj.Text=regexprep(obj.Text,p,"",'once');end
        end
        function addInstructions(obj,varargin),for k=1:numel(varargin),obj.addInstruction(string(varargin{k}));end,end
        function setParameters(obj,values),names=fieldnames(values);for k=1:numel(names),obj.setParameter(names{k},values.(names{k}));end,end
        function setComponentValues(obj,values),names=fieldnames(values);for k=1:numel(names),obj.setComponentValue(names{k},values.(names{k}));end,end
        function resetNetlist(obj),obj.Text=string(fileread(obj.SourcePath));end
        function answer=isReadOnly(~),answer=false;end
        function saveAs(obj,path)
            arguments, obj; path (1,1) string, end
            folder=fileparts(path); if strlength(folder)>0 && ~isfolder(folder), mkdir(folder); end
            f=fopen(path,'w'); if f<0, error("radia:ltspice:Write","Cannot write %s",path); end
            c=onCleanup(@()fclose(f)); fprintf(f,'%s',obj.Text); clear c
        end
        function writeNetlist(obj,path),obj.saveAs(path);end
        function saveNetlist(obj,path),obj.saveAs(path);end
        % PyLTSpice spelling aliases.
        function set_parameter(obj,varargin),obj.setParameter(varargin{:});end
        function x=get_parameter(obj,varargin),x=obj.getParameter(varargin{:});end
        function x=get_all_parameter_names(obj),x=obj.getAllParameterNames();end
        function set_parameters(obj,varargin),obj.setParameters(varargin{:});end
        function set_component_value(obj,varargin),obj.setComponentValue(varargin{:});end
        function set_component_values(obj,varargin),obj.setComponentValues(varargin{:});end
        function x=get_component_value(obj,varargin),x=obj.getComponentValue(varargin{:});end
        function x=get_component_floatvalue(obj,varargin),x=obj.getComponentFloatValue(varargin{:});end
        function x=get_component_nodes(obj,varargin),x=obj.getComponentNodes(varargin{:});end
        function x=get_component(obj,varargin),x=obj.getComponent(varargin{:});end
        function x=get_component_attribute(obj,varargin),x=obj.getComponentAttribute(varargin{:});end
        function set_component_attribute(obj,varargin),obj.setComponentAttribute(varargin{:});end
        function x=get_component_parameters(obj,varargin),x=obj.getComponentParameters(varargin{:});end
        function set_component_parameters(obj,varargin),obj.setComponentParameters(varargin{:});end
        function x=get_components(obj,varargin),x=obj.getComponents(varargin{:});end
        function add_component(obj,varargin),obj.addComponent(varargin{:});end
        function remove_component(obj,varargin),obj.removeComponent(varargin{:});end
        function set_element_model(obj,varargin),obj.setElementModel(varargin{:});end
        function x=get_all_nodes(obj),x=obj.getAllNodes();end
        function x=get_subcircuit_names(obj),x=obj.getSubcircuitNames();end
        function x=get_subcircuit(obj,varargin),x=obj.getSubcircuit(varargin{:});end
        function x=get_subcircuit_named(obj,varargin),x=obj.getSubcircuitNamed(varargin{:});end
        function x=modified_subcircuits(obj),x=obj.modifiedSubcircuits();end
        function x=get_control_sections(obj),x=obj.getControlSections();end
        function add_control_section(obj,varargin),obj.addControlSection(varargin{:});end
        function x=remove_control_section(obj,varargin),x=obj.removeControlSection(varargin{:});end
        function add_library_search_paths(obj,varargin),obj.addLibrarySearchPaths(varargin{:});end
        function set_custom_library_paths(obj,varargin),obj.LibraryPaths=string(varargin);end
        function x=find_library(obj,varargin),x=obj.findLibrary(varargin{:});end
        function x=find_subckt_in_included_libs(obj,name),x=obj.getSubcircuitNamed(name);end
        function x=find_subckt_in_lib(~,library,name),editor=radia.ltspice.SpiceEditor(string(library));x=editor.getSubcircuitNamed(name);end
        function x=begin_update(obj),x=obj.beginUpdate();end
        function end_update(obj,varargin),obj.endUpdate(varargin{:});end
        function x=class_for_instruction(obj,varargin),x=obj.classForInstruction(varargin{:});end
        function x=remove_Xinstruction(obj,varargin),x=obj.removeXInstruction(varargin{:});end
        function prepare_for_simulator(obj,varargin),obj.prepareForSimulator(varargin{:});end
        function write_lines(obj,varargin),obj.writeLines(varargin{:});end
        function add_instruction(obj,varargin),obj.addInstruction(varargin{:});end
        function add_instructions(obj,varargin),obj.addInstructions(varargin{:});end
        function x=remove_instruction(obj,varargin),x=obj.removeInstruction(varargin{:});end
        function x=reset_netlist(obj,varargin),obj.resetNetlist(varargin{:});x=true;end
        function x=is_read_only(obj),x=obj.isReadOnly();end
        function save_as(obj,varargin),obj.saveAs(varargin{:});end
        function save_netlist(obj,varargin),obj.saveAs(varargin{:});end
        function write_netlist(obj,varargin),obj.saveAs(varargin{:});end
    end
    methods (Access=private)
        function line=componentLine(obj,reference)
            p="(?im)^\s*"+regexptranslate('escape',char(reference))+"\s+[^\r\n]*";line=string(regexp(obj.Text,p,'match','once'));
            if strlength(line)==0,error("radia:ltspice:ComponentNotFound","Component not found: %s",reference);end
        end
        function [parts,valueIndex]=componentTokens(obj,reference)
            parts=split(strtrim(obj.componentLine(reference)));prefix=upper(extractBefore(parts(1),2));fixed=struct('R',4,'C',4,'L',4,'D',4,'V',4,'I',4,'B',4,'E',6,'G',6,'F',5,'H',5,'Q',5,'J',5,'M',6,'K',4,'T',6);
            if isfield(fixed,prefix),valueIndex=min(fixed.(prefix),numel(parts));elseif prefix=="X",nonParameter=find(~contains(parts,"="));valueIndex=nonParameter(end);else,nonParameter=find(~contains(parts,"="));valueIndex=nonParameter(end);end
        end
    end
    methods (Static,Access=private)
        function v=formatValue(value)
            if isnumeric(value)&&isscalar(value)&&isfinite(value)
                v=string(sprintf('%.17g',value));
            elseif (isstring(value)&&isscalar(value))||(ischar(value)&&isrow(value))
                v=string(value);
            else
                error("radia:ltspice:Value","Value must be finite scalar or text.");
            end
        end
    end
end
function value=localEngineeringNumber(text)
s=lower(char(text)); tokens=regexp(s,'^([-+0-9.eE]+)(meg|[tgkmunpf]?)','tokens','once');if isempty(tokens),value=NaN;return,end
value=str2double(tokens{1}); factors=struct('t',1e12,'g',1e9,'meg',1e6,'k',1e3,'m',1e-3,'u',1e-6,'n',1e-9,'p',1e-12,'f',1e-15);suffix=tokens{2};if ~isempty(suffix),value=value*factors.(suffix);end
end
