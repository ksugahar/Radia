function result=validate_interval_state_handoff(options)
%VALIDATE_INTERVAL_STATE_HANDOFF Quantify reset-time handoff error growth.
arguments
 options.OutputDirectory (1,1) string="C:\temp\radia_ltspice_interval_validation"
 options.ResultPath (1,1) string=""
end
startedAt=datetime("now","TimeZone","local");stopwatch=tic;
repoRoot=fileparts(fileparts(fileparts(mfilename("fullpath"))));originalPath=path;pathCleanup=onCleanup(@()path(originalPath));addpath(fullfile(repoRoot,"matlab"));
root=options.OutputDirectory;if isfolder(root),rmdir(root,'s');end;mkdir(root);
executable=radia.ltspice.findExecutable();versionInfo=System.Diagnostics.FileVersionInfo.GetVersionInfo(char(executable));
host=string(getenv("COMPUTERNAME"));if strlength(host)==0,host=string(System.Net.Dns.GetHostName());end
netlist=fullfile(root,"ring.cir");writeText(netlist,[ ...
 "* linear RLC ring, Q approximately 32";"V1 in 0 PULSE(0 1 0 1u 1u 5m 10m)"; ...
 "L1 in out 1m";"C1 out 0 1u";"R1 out 0 1";".tran 0 2m 0 1u";".end"]);
continuous=radia.ltspice.run(netlist,OutputDirectory=fullfile(root,"continuous"));
trace=find(continuous.waveform.names=="V(out)",1);reference=continuous.waveform.values(end,trace);
scale=max(abs(continuous.waveform.values(:,trace)));counts=[2,4,8,16];errors=zeros(size(counts));
for k=1:numel(counts)
 n=counts(k);split=radia.ltspice.runIntervals(netlist,repmat(2e-3/n,n,1),OutputDirectory=fullfile(root,"split_"+n),MaxStep_s=1e-6);
 errors(k)=abs(split.runs{end}.waveform.values(end,trace)-reference)/max(scale,eps);
end
handoffs=counts-1;slope=handoffs(:)\errors(:);correlation=corrcoef(handoffs(:),errors(:));
gates=struct("finite",all(isfinite(errors)),"bounded",all(errors<1e-2), ...
 "approximately_linear",correlation(1,2)>0.98,"nondecreasing",all(diff(errors)>0));
result=struct("schema","radia.validation.ltspice.interval_handoff.v1", ...
 "circuit","linear RLC ring","duration_s",2e-3,"split_counts",counts, ...
 "handoff_counts",handoffs,"relative_endpoint_errors",errors, ...
 "relative_error_per_handoff_fit",slope,"fit_model","least squares constrained through origin", ...
 "host",host,"started_at",string(startedAt,"yyyy-MM-dd'T'HH:mm:ssXXX"), ...
 "elapsed_s",toc(stopwatch),"matlab_version",string(version),"ltspice_version",string(versionInfo.FileVersion), ...
 "ltspice_executable",executable,"gates",gates);
if strlength(options.ResultPath)==0,options.ResultPath=fullfile(fileparts(mfilename("fullpath")),"results_interval_state_handoff.json");end
writeText(options.ResultPath,jsonencode(result,PrettyPrint=true));
if ~all(structfun(@(x)x,gates)),error("radia:validation:LTspiceIntervalHandoff","LTspice interval handoff validation failed.");end
fprintf("LTSPICE_INTERVAL_HANDOFF_OK slope=%.6g max_error=%.6g\n",slope,max(errors));
clear pathCleanup
end

function writeText(path,text)
folder=fileparts(path);if ~isfolder(folder),mkdir(folder);end
file=fopen(path,'w');assert(file>=0);cleanup=onCleanup(@()fclose(file));fprintf(file,'%s',join(string(text),newline));clear cleanup
end
