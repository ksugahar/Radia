classdef SpiceEditor < handle
    %SPICEEDITOR Edit LTspice netlists using a MATLAB-native API.
    properties (SetAccess=private)
        SourcePath (1,1) string
        Text (1,1) string
        LibraryPaths (1,:) string=strings(1,0)
        CircuitName (1,1) string=""
        LineEnding (1,1) string=newline
    end
    methods
        function obj=SpiceEditor(path)
            arguments, path (1,1) string {mustBeFile}, end
            obj.SourcePath=path; raw=fileread(path); obj.Text=string(raw); obj.LineEnding=obj.detectLineEnding(raw);[~,obj.CircuitName]=fileparts(path);
        end
        function setParameter(obj,name,value)
            arguments, obj; name (1,1) string; value, end
            v=radia.ltspice.SpiceEditor.formatValue(value);
            [line,index,lines]=obj.parameterLine(name); p="(?i)(^\s*\.param\s+"+regexptranslate('escape',char(name))+"\s*=\s*)([^\s;]+)";
            lines(index)=regexprep(line,p,"$1"+v,'once'); obj.Text=join(lines,obj.LineEnding);
        end
        function value=getParameter(obj,name)
            line=obj.parameterLine(name); p="(?i)^\s*\.param\s+"+regexptranslate('escape',char(name))+"\s*=\s*([^\s;]+)"; t=regexp(line,p,'tokens','once'); value=string(t{1});
        end
        function names=getAllParameterNames(obj)
            lines=splitlines(obj.Text); [topMask,~]=obj.lineScopes(lines); names=strings(0,1);
            for k=2:numel(lines),if ~topMask(k),continue,end,t=regexp(lines(k),'(?i)^\s*\.param\s+([A-Za-z_]\w*)\s*=','tokens','once');if ~isempty(t),names(end+1,1)=string(t{1});end,end %#ok<AGROW>
            folded=lower(names); duplicate=unique(folded(arrayfun(@(x)sum(folded==x)>1,folded)),'stable'); if ~isempty(duplicate),error("radia:ltspice:AmbiguousParameter","Parameter is defined more than once: %s",join(duplicate,", "));end
        end
        function references=getComponents(obj,prefixes)
            if nargin<2,prefixes="*";end
            lines=splitlines(obj.Text); [topMask,~]=obj.lineScopes(lines); references=strings(0,1);
            for k=2:numel(lines),if ~topMask(k),continue,end,s=strtrim(lines(k));if strlength(s)==0||startsWith(s,[".";"*";";";"+"]),continue,end,r=extractBefore(s+" "," ");designator=upper(extractBetween(r,1,1));if isempty(regexp(char(designator),'^[A-Z]$','once')),continue,end,if prefixes=="*"||any(startsWith(upper(string(prefixes)),designator)),references(end+1,1)=r;end,end %#ok<AGROW>
        end
        function value=getComponentValue(obj,reference)
            [parts,valueIndex]=obj.componentTokens(reference);value=parts(valueIndex);
        end
        function component=getComponent(obj,reference),component=struct('reference',string(reference),'nodes',obj.getComponentNodes(reference),'value',obj.getComponentValue(reference),'line',obj.componentLine(reference));end
        function value=getComponentAttribute(obj,reference,attribute)
            line=obj.componentLine(reference);t=regexp(line,"(?i)(?:^|\s)"+regexptranslate('escape',char(attribute))+"\s*=\s*([^\s]+)",'tokens','once');if isempty(t),value="";else,value=string(t{1});end
        end
        function setComponentAttribute(obj,reference,attribute,value)
            [line,index,lines]=obj.componentLine(reference);p="(?i)(\s"+regexptranslate('escape',char(attribute))+"\s*=\s*)([^\s]+)";if isempty(regexp(line,p,'once')),replacement=line+" "+attribute+"="+string(value);else,replacement=regexprep(line,p,"$1"+string(value),'once');end,lines(index)=replacement;obj.Text=join(lines,obj.LineEnding);
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
            [~,index,lines]=obj.componentLine(reference);stop=index;while stop<numel(lines)&&startsWith(strtrim(lines(stop+1)),"+"),stop=stop+1;end,lines(index:stop)=[];obj.Text=join(lines,obj.LineEnding);
        end
        function addComponent(obj,reference,nodes,value,parameters)
            arguments,obj;reference (1,1) string;nodes (:,1) string;value (1,1) string;parameters (1,1) struct=struct();end
            if any(strcmpi(obj.getComponents(),reference)),error("radia:ltspice:DuplicateComponent","Component already exists: %s",reference);end,line=reference+" "+join(nodes," ")+" "+value;names=fieldnames(parameters);for k=1:numel(names),line=line+" "+names{k}+"="+string(parameters.(names{k}));end,if isempty(regexp(obj.Text,'(?im)^\s*\.end\s*$','once')),error("radia:ltspice:MissingEnd","Cannot add a component because the netlist has no .end directive.");end,obj.Text=regexprep(obj.Text,"(?im)^\s*\.end\s*$",line+obj.LineEnding+".end",'once');
        end
        function setElementModel(obj,reference,model),obj.setComponentValue(reference,model);end
        function nodes=getAllNodes(obj),references=obj.getComponents();nodes=strings(0,1);for k=1:numel(references),nodes=[nodes;obj.getComponentNodes(references(k))];end,nodes=unique(nodes,'stable');end %#ok<AGROW>
        function names=getSubcircuitNames(obj)
            lines=splitlines(obj.Text);definitions=obj.subcircuitDefinitions(lines);
            if isempty(definitions),names=strings(0,1);return,end
            allNames=[definitions.name];depths=[definitions.depth];names=allNames(depths==1)';
        end
        function sections=getControlSections(obj),sections=string(regexp(obj.Text,'(?ims)^\s*\.control\s*$.*?^\s*\.endc\s*$','match'))';end
        function addControlSection(obj,instruction),if isempty(regexp(obj.Text,'(?im)^\s*\.end\s*$','once')),error("radia:ltspice:MissingEnd","Cannot add a control section because the netlist has no .end directive.");end,obj.Text=regexprep(obj.Text,"(?im)^\s*\.end\s*$",".control"+obj.LineEnding+instruction+obj.LineEnding+".endc"+obj.LineEnding+".end",'once');end
        function removed=removeControlSection(obj,index),if nargin<2,index=1;end,sections=obj.getControlSections();if index<1||index>numel(sections),removed=false;else,obj.Text=replace(obj.Text,sections(index),"");removed=true;end,end
        function addLibrarySearchPaths(obj,varargin),obj.LibraryPaths=[obj.LibraryPaths,string(varargin)];end
        function path=findLibrary(obj,name),candidates=[fullfile(fileparts(obj.SourcePath),name),fullfile(obj.LibraryPaths,name)];index=find(isfile(candidates),1);if isempty(index),path="";else,path=candidates(index);end,end
        function circuit=getSubcircuitNamed(obj,name)
            lines=splitlines(obj.Text);definitions=obj.subcircuitDefinitions(lines);matches=find(strcmpi([definitions.name],name));if isempty(matches),circuit=[];return,end,if numel(matches)>1,error("radia:ltspice:AmbiguousSubcircuit","Subcircuit is defined more than once: %s",name);end,definition=definitions(matches);block=join(lines(definition.first:definition.last),obj.LineEnding);path=string(tempname("C:\temp"))+"_"+matlab.lang.makeValidName(char(name))+".cir";f=fopen(path,'w');c=onCleanup(@()fclose(f));fprintf(f,'%s%s.end',block,obj.LineEnding);clear c;circuit=radia.ltspice.SpiceCircuit(path);
        end
        function circuit=getSubcircuit(obj,instanceName),line=obj.componentLine(instanceName);parts=split(strtrim(line));parameter=find(contains(parts,"="),1);marker=find(strcmpi(parts,"params:"),1);cutoff=[parameter,marker];cutoff=cutoff(cutoff>0);if isempty(cutoff),index=numel(parts);else,index=min(cutoff)-1;end,model=parts(index);circuit=obj.getSubcircuitNamed(model);if isempty(circuit),error("radia:ltspice:SubcircuitNotFound","Subcircuit not found: %s",model);end,end
        function circuits=modifiedSubcircuits(~),circuits={};end
        function value=name(obj),value=obj.CircuitName;end
        function setname(obj,value),obj.CircuitName=string(value);end
        function permission=beginUpdate(~),permission="allow";end
        function endUpdate(varargin),end
        function cloned=clone(obj,varargin),path=fullfile("C:\temp","radia_spice_clone_"+string(char(java.util.UUID.randomUUID()))+".cir");obj.saveAs(path);cloned=radia.ltspice.SpiceCircuit(path);end
        function lineClass=classForInstruction(~,instruction,varargin),lineClass=upper(extractBefore(strtrim(string(instruction))+" "," "));end
        function removed=removeXInstruction(obj,pattern),lines=splitlines(obj.Text);mask=startsWith(upper(strtrim(lines)),"X")&~cellfun(@isempty,regexp(cellstr(lines),char(pattern),'once'));removed=any(mask);lines(mask)=[];obj.Text=join(lines,obj.LineEnding);end
        function prepareForSimulator(~,varargin),end
        function writeLines(obj,stream),fprintf(stream,'%s',obj.Text);end
        function setComponentValue(obj,reference,value)
            arguments, obj; reference (1,1) string; value, end
            v=radia.ltspice.SpiceEditor.formatValue(value);
            [~,lineIndex,lines]=obj.componentLine(reference);[parts,index]=obj.componentTokens(reference);parts(index)=v;lines(lineIndex)=join(parts," ");obj.Text=join(lines,obj.LineEnding);
        end
        function addInstruction(obj,instruction)
            arguments, obj; instruction (1,1) string, end
            if isempty(regexp(obj.Text,'(?im)^\s*\.end\s*$','once')),error("radia:ltspice:MissingEnd","Cannot add an instruction because the netlist has no .end directive.");end,obj.Text=regexprep(obj.Text,"(?im)^\s*\.end\s*$",instruction+obj.LineEnding+".end",'once');
        end
        function removed=removeInstruction(obj,instruction)
            p="(?im)^\s*"+regexptranslate('escape',char(instruction))+"\s*(?:\r?\n)?";removed=~isempty(regexp(obj.Text,p,'once'));if removed,obj.Text=regexprep(obj.Text,p,"",'once');end
        end
        function addInstructions(obj,varargin),for k=1:numel(varargin),obj.addInstruction(string(varargin{k}));end,end
        function setParameters(obj,values),names=fieldnames(values);for k=1:numel(names),obj.setParameter(names{k},values.(names{k}));end,end
        function setComponentValues(obj,values),names=fieldnames(values);for k=1:numel(names),obj.setComponentValue(names{k},values.(names{k}));end,end
        function resetNetlist(obj),raw=fileread(obj.SourcePath);obj.Text=string(raw);obj.LineEnding=obj.detectLineEnding(raw);end
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
        function [line,index,lines]=componentLine(obj,reference)
            lines=splitlines(obj.Text);[topMask,scopeNames]=obj.lineScopes(lines);p="(?i)^\s*"+regexptranslate('escape',char(reference))+"(?:\s|$)";matches=find(~cellfun(@isempty,regexp(cellstr(lines),char(p),'once')));matches(matches==1)=[];topMatches=matches(topMask(matches));
            if numel(topMatches)>1,error("radia:ltspice:AmbiguousComponent","Top-level component reference occurs more than once: %s",reference);end
            if isempty(topMatches)
                nested=matches(~topMask(matches));if ~isempty(nested),scopes=unique(scopeNames(nested));scopes(scopes=="")=[];error("radia:ltspice:ComponentInSubcircuit","Component %s exists only inside subcircuit(s): %s. Edit that subcircuit explicitly.",reference,join(scopes,", "));end
                error("radia:ltspice:ComponentNotFound","Top-level component not found: %s",reference);
            end
            index=topMatches(1);line=lines(index);
        end
        function [line,index,lines]=parameterLine(obj,name)
            lines=splitlines(obj.Text);[topMask,scopeNames]=obj.lineScopes(lines);p="(?i)^\s*\.param\s+"+regexptranslate('escape',char(name))+"\s*=";matches=find(~cellfun(@isempty,regexp(cellstr(lines),char(p),'once')));topMatches=matches(topMask(matches));
            if numel(topMatches)>1,error("radia:ltspice:AmbiguousParameter","Top-level parameter is defined more than once: %s",name);end
            if isempty(topMatches)
                nested=matches(~topMask(matches));if ~isempty(nested),scopes=unique(scopeNames(nested));scopes(scopes=="")=[];error("radia:ltspice:ParameterInSubcircuit","Parameter %s exists only inside subcircuit(s): %s. Edit that subcircuit explicitly.",name,join(scopes,", "));end
                error("radia:ltspice:ParameterNotFound","Top-level parameter not found: %s",name);
            end
            index=topMatches(1);line=lines(index);
        end
        function [topMask,scopeNames]=lineScopes(obj,lines)
            topMask=true(numel(lines),1);scopeNames=strings(numel(lines),1);definitions=obj.subcircuitDefinitions(lines);
            for index=find([definitions.depth]==1)
                first=definitions(index).first;last=definitions(index).last;topMask(first:last)=false;scopeNames(first:last)=definitions(index).name;
            end
        end
        function definitions=subcircuitDefinitions(~,lines)
            definitions=struct('name',{},'first',{},'last',{},'depth',{});stack=zeros(0,1);
            for k=2:numel(lines)
                s=strtrim(lines(k));opening=regexp(s,'(?i)^\.subckt\s+(\S+)','tokens','once');
                if ~isempty(opening)
                    index=numel(definitions)+1;definitions(index)=struct('name',string(opening{1}),'first',k,'last',0,'depth',numel(stack)+1);stack(end+1,1)=index; %#ok<AGROW>
                elseif ~isempty(regexp(s,'(?i)^\.ends(?:\s|$)','once'))&&~isempty(stack)
                    definitions(stack(end)).last=k;stack(end)=[];
                end
            end
            if ~isempty(stack),names=join([definitions(stack).name],", ");error("radia:ltspice:UnterminatedSubcircuit","Subcircuit definition has no matching .ends: %s",names);end
        end
        function [parts,valueIndex]=componentTokens(obj,reference)
            parts=split(strtrim(obj.componentLine(reference)));prefix=upper(extractBetween(parts(1),1,1));fixed=struct('R',4,'C',4,'L',4,'D',4,'V',4,'I',4,'B',4,'E',6,'G',6,'F',5,'H',5,'Q',5,'J',5,'M',6,'K',4,'T',6);
            if isfield(fixed,prefix),valueIndex=min(fixed.(prefix),numel(parts));elseif prefix=="X",parameter=find(contains(parts,"="),1);marker=find(strcmpi(parts,"params:"),1);cutoff=[parameter,marker];cutoff=cutoff(cutoff>0);if isempty(cutoff),valueIndex=numel(parts);else,valueIndex=min(cutoff)-1;end;else,nonParameter=find(~contains(parts,"="));valueIndex=nonParameter(end);end
        end
    end
    methods (Static,Access=private)
        function ending=detectLineEnding(raw)
            if contains(raw,sprintf('\r\n')),ending=string(sprintf('\r\n'));elseif contains(raw,sprintf('\r')),ending=string(sprintf('\r'));else,ending=string(newline);end
        end
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
