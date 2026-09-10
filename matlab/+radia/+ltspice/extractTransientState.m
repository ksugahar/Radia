function state=extractTransientState(raw,options)
%EXTRACTTRANSIENTSTATE Extract capacitor-node voltages and inductor currents.
arguments, raw; options.Step (1,1) double {mustBeInteger,mustBePositive}=1; options.NetlistFile (1,1) string=""; end
if isa(raw,"radia.ltspice.RawRead"),data=raw.Data;elseif isstruct(raw)&&isfield(raw,"waveform"),data=raw.waveform;elseif isstruct(raw),data=raw;else,error("radia:ltspice:StateInput","raw must be RawRead or RAW/run struct.");end
range=data.step_ranges(options.Step,:); names=data.names; values=data.values(range(2),:);
nodeNames=strings(0,1); nodeValues=zeros(0,1); inductorNames=strings(0,1); inductorValues=zeros(0,1);
for k=2:numel(names)
 name=names(k); voltage=regexp(char(name),'^V\((.+)\)$','tokens','once','ignorecase'); current=regexp(char(name),'^I\(((?:[^():]+:)*L[^)]+)\)$','tokens','once','ignorecase');
 if ~isempty(voltage),nodeNames(end+1,1)=string(voltage{1});nodeValues(end+1,1)=real(values(k));
 elseif ~isempty(current),inductorNames(end+1,1)=string(current{1});inductorValues(end+1,1)=real(values(k));end
end
inductorLike=~cellfun(@isempty,regexp(cellstr(names),'^I\((?:L|.*:L)','once','ignorecase'));
captured="I("+inductorNames+")";missing=names(inductorLike&~ismember(upper(names),upper(captured)));
if ~isempty(missing),error("radia:ltspice:StateTraceUnsupported","Inductor state trace(s) could not be represented: %s.",join(missing,", "));end
if strlength(options.NetlistFile)>0
 instances=statefulTopLevelInstances(options.NetlistFile);
 for k=1:numel(instances)
  instance=instances(k);
  escaped=regexptranslate('escape',char(instance));
  voltagePattern=['^V\(' escaped ':[^)]+\)$'];
  inductorPattern=['^I\(' escaped ':(?:[^():]+:)*L[^)]*\)$'];
  carriesState=any(~cellfun(@isempty,regexp(cellstr(names),voltagePattern,'once','ignorecase'))) || ...
   any(~cellfun(@isempty,regexp(cellstr(names),inductorPattern,'once','ignorecase')));
  if ~carriesState
   error("radia:ltspice:MissingSubcircuitState", ...
    "No internal state traces were saved for subcircuit instance %s. Add explicit .save traces before interval handoff.",instance);
  end
 end
end
state=struct("schema","radia.ltspice.transient_state.v1","time_s",real(values(1)), ...
 "node_names",nodeNames,"node_voltages_V",nodeValues, ...
 "inductor_names",inductorNames,"inductor_currents_A",inductorValues);
end

function instances=statefulTopLevelInstances(netlistFile)
manifest=radia.ltspice.collectDependencies(netlistFile);
definitions=containers.Map('KeyType','char','ValueType','any');top=strings(0,2);
for file=manifest.local_files(:).'
 lines=splitlines(string(fileread(file)));current="";
 for row=1:numel(lines)
  line=strtrim(regexprep(lines(row),';.*$',''));if line==""||startsWith(line,"*"),continue,end
  tokens=split(line);tokens(tokens=="")=[];head=lower(tokens(1));
  if head==".subckt"&&numel(tokens)>=2
   current=lower(tokens(2));if ~isKey(definitions,char(current)),definitions(char(current))=struct("direct",false,"unsupported",false,"children",strings(0,1));end
  elseif head==".ends",current="";
  elseif startsWith(head,["l","c"])&&current~=""
   item=definitions(char(current));item.direct=true;definitions(char(current))=item;
  elseif startsWith(head,["d","q","m","j","z","t","b"])&&current~=""
   item=definitions(char(current));item.unsupported=true;definitions(char(current))=item;
  elseif startsWith(head,"x")
   child=subcircuitName(tokens);
   if current=="",top(end+1,:)=[tokens(1),child];else,item=definitions(char(current));item.children(end+1,1)=child;definitions(char(current))=item;end
  end
 end
end
stateful=containers.Map('KeyType','char','ValueType','double');keysList=string(keys(definitions));
for key=keysList
 item=definitions(char(key));stateful(char(key))=double(item.direct)+2*double(item.unsupported);
end
changed=true;
while changed
 changed=false;
 for key=keysList
  if stateful(char(key))~=0,continue,end
  children=definitions(char(key)).children;
  known=children(arrayfun(@(x)isKey(stateful,char(x)),children));
  childStates=arrayfun(@(x)stateful(char(x)),known);
  if any(childStates>=2),stateful(char(key))=2;changed=true;
  elseif any(childStates==1),stateful(char(key))=1;changed=true;end
 end
end
keep=false(size(top,1),1);
for k=1:size(top,1)
 type=top(k,2);
 if isKey(stateful,char(type))&&stateful(char(type))>=2
  error("radia:ltspice:UnsupportedSubcircuitState", ...
   "Subcircuit instance %s contains device-internal state that cannot be completely restored by node-voltage/inductor-current handoff.",top(k,1));
 end
 keep(k)=~isKey(stateful,char(type))||stateful(char(type))==1;
end
instances=top(keep,1);
end

function name=subcircuitName(tokens)
parameter=find(contains(tokens,"="),1);marker=find(strcmpi(tokens,"params:"),1);
cutoff=[parameter,marker];cutoff=cutoff(cutoff>0);
if isempty(cutoff),index=numel(tokens);else,index=min(cutoff)-1;end
if index<2,name="";else,name=lower(tokens(index));end
end
