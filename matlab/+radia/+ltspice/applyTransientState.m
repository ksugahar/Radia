function info=applyTransientState(netlistFile,state,destination,options)
%APPLYTRANSIENTSTATE Write .ic state and a reset-time transient directive.
arguments
 netlistFile (1,1) string {mustBeFile}; state (1,1) struct; destination (1,1) string
 options.Duration_s (1,1) double {mustBePositive}; options.MaxStep_s (1,1) double {mustBePositive}=inf
end
if ~isfield(state,"schema")||state.schema~="radia.ltspice.transient_state.v1",error("radia:ltspice:StateContract","Invalid transient state.");end
required=["node_names","node_voltages_V","inductor_names","inductor_currents_A"];
if ~all(isfield(state,required)),error("radia:ltspice:StateContract","Transient state fields are incomplete.");end
text=string(fileread(netlistFile)); directives=strings(0,1);
for k=1:numel(state.node_names),directives(end+1)=".ic V("+state.node_names(k)+")="+sprintf("%.17g",state.node_voltages_V(k));end
for k=1:numel(state.inductor_names),directives(end+1)=".ic I("+state.inductor_names(k)+")="+sprintf("%.17g",state.inductor_currents_A(k));end
if isfinite(options.MaxStep_s),tran=sprintf(".tran 0 %.17g 0 %.17g uic",options.Duration_s,options.MaxStep_s);else,tran=sprintf(".tran 0 %.17g uic",options.Duration_s);end
hadTran=~isempty(regexp(text,'(?im)^\s*\.tran(?=\s|$)','once'));
text=regexprep(text,'(?im)^\s*\.tran(?=\s|$)[^\r\n]*',tran,'once'); if ~hadTran,text=regexprep(text,'(?im)^\s*\.end\s*$',tran+newline+".end",'once');end
text=regexprep(text,'(?im)^\s*\.ic(?=\s|$)[^\r\n]*\r?\n?'," ");
if isempty(directives),ending=".end";else,ending=join(directives,newline)+newline+".end";end
text=regexprep(text,'(?im)^\s*\.end\s*$',ending,'once');
folder=fileparts(destination);if strlength(folder)>0&&~isfolder(folder),mkdir(folder);end
f=fopen(destination,'w');if f<0,error("radia:ltspice:Write","Cannot write %s",destination);end;c=onCleanup(@()fclose(f));fprintf(f,'%s',text);clear c
info=struct("schema","radia.ltspice.state_netlist.v1","path",destination,"state_time_s",state.time_s,"duration_s",options.Duration_s);
end
