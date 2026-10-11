function result=prepareLTspice(schematicFile,options)
%PREPARELTSPICE Export KiCad SPICE and create an editable LTspice schematic.
arguments
 schematicFile (1,1) string {mustBeFile}
 options.OutputDirectory (1,1) string=""
 options.KiCadExecutable (1,1) string=""
 options.LTspiceExecutable (1,1) string=""
 options.PythonExecutable (1,1) string="python"
 options.AsySearchDirs (:,1) string=strings(0,1)
 options.ValidateRoundTrip (1,1) logical=true
end
[sourceFolder,stem]=fileparts(schematicFile);
if strlength(options.OutputDirectory)==0,outputDirectory=sourceFolder;else,outputDirectory=options.OutputDirectory;end
if ~isfolder(outputDirectory),mkdir(outputDirectory);end
netlist=radia.kicad.exportSpiceNetlist(schematicFile,OutputFile=fullfile(outputDirectory,stem+".cir"),Executable=options.KiCadExecutable);
schematic=radia.ltspice.netlistToSchematic(netlist.netlist_file,OutputFile=fullfile(outputDirectory,stem+".asc"), ...
 PythonExecutable=options.PythonExecutable,AsySearchDirs=options.AsySearchDirs, ...
 ValidateRoundTrip=options.ValidateRoundTrip,LTspiceExecutable=options.LTspiceExecutable);
result=struct("schema","radia.kicad.ltspice_preparation.v1","source_schematic",string(schematicFile), ...
 "netlist_file",netlist.netlist_file,"ltspice_schematic",schematic.schematic_file, ...
 "kicad_export",netlist,"ltspice_conversion",schematic);
end
