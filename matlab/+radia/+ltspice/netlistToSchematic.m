function result=netlistToSchematic(netlistFile,options)
%NETLISTTOSCHEMATIC Wrap radia-spice-lab's canonical Python CIR->ASC converter.
arguments
 netlistFile (1,1) string {mustBeFile}
 options.OutputFile (1,1) string=""
 options.PythonExecutable (1,1) string="python"
 options.AsySearchDirs (:,1) string=strings(0,1)
 options.ValidateRoundTrip (1,1) logical=true
 options.LTspiceExecutable (1,1) string=""
end
[sourceFolder,stem,extension]=fileparts(netlistFile);
if ~ismember(lower(string(extension)),[".cir",".net",".sp",".spi"]),error("radia:ltspice:NetlistRequired","Input must be a SPICE netlist.");end
if strlength(options.OutputFile)==0,outputFile=fullfile(sourceFolder,stem+".asc");else,outputFile=options.OutputFile;end
folder=fileparts(outputFile);if strlength(folder)>0&&~isfolder(folder),mkdir(folder);end
parts=[quote(options.PythonExecutable),"-m","ltspice_converter.cli",quote(netlistFile),"-o",quote(outputFile),"--to","asc"];
for directory=options.AsySearchDirs(:)',parts(end+1:end+2)=["--asy-dir",quote(directory)];end
[status,commandOutput]=system(join(parts," "));
if status~=0||~isfile(outputFile)
 error("radia:ltspice:ConverterFailed","radia-spice-lab failed (status %d).\n%s",status,commandOutput);
end
validation=struct.empty;
if options.ValidateRoundTrip
 validationFolder=string(tempname("C:\temp"));mkdir(validationFolder);
 validation=radia.ltspice.schematicToNetlist(outputFile,Executable=options.LTspiceExecutable,OutputDirectory=validationFolder);
 pythonCode="from pathlib import Path; from spice_circuit_lab import topology_equivalent; import json,sys; ok,details=topology_equivalent(Path(sys.argv[1]).read_text(encoding='utf-8',errors='replace'),Path(sys.argv[2]).read_text(encoding='utf-8',errors='replace')); print(json.dumps(dict(equivalent=ok,details=details)))";
 [compareStatus,compareOutput]=system(join([quote(options.PythonExecutable),"-c",quote(pythonCode),quote(netlistFile),quote(validation.netlist)]," "));
 if compareStatus~=0,error("radia:ltspice:TopologyCheckFailed","Python topology check failed.\n%s",compareOutput);end
 comparison=jsondecode(compareOutput);
 validation.topology=comparison;
 if ~comparison.equivalent,error("radia:ltspice:TopologyMismatch","Generated ASC does not preserve netlist connectivity.");end
end
result=struct("schema","radia.ltspice.netlist_to_schematic.v1","source_netlist",netlistFile, ...
 "schematic_file",string(outputFile),"converter","radia-spice-lab", ...
 "converter_output",string(commandOutput),"validation",validation);
end

function value=quote(value)
value='"'+replace(string(value),'"','\"')+'"';
end
