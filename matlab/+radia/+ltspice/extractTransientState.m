function state=extractTransientState(raw,options)
%EXTRACTTRANSIENTSTATE Extract capacitor-node voltages and inductor currents.
arguments, raw; options.Step (1,1) double {mustBeInteger,mustBePositive}=1; end
if isa(raw,"radia.ltspice.RawRead"),data=raw.Data;elseif isstruct(raw)&&isfield(raw,"waveform"),data=raw.waveform;elseif isstruct(raw),data=raw;else,error("radia:ltspice:StateInput","raw must be RawRead or RAW/run struct.");end
range=data.step_ranges(options.Step,:); names=data.names; values=data.values(range(2),:);
nodeNames=strings(0,1); nodeValues=zeros(0,1); inductorNames=strings(0,1); inductorValues=zeros(0,1);
for k=2:numel(names)
 name=names(k); voltage=regexp(char(name),'^V\((.+)\)$','tokens','once','ignorecase'); current=regexp(char(name),'^I\((L[^)]+)\)$','tokens','once','ignorecase');
 if ~isempty(voltage),nodeNames(end+1,1)=string(voltage{1});nodeValues(end+1,1)=real(values(k));
 elseif ~isempty(current),inductorNames(end+1,1)=string(current{1});inductorValues(end+1,1)=real(values(k));end
end
state=struct("schema","radia.ltspice.transient_state.v1","time_s",real(values(1)), ...
 "node_names",nodeNames,"node_voltages_V",nodeValues, ...
 "inductor_names",inductorNames,"inductor_currents_A",inductorValues);
end
