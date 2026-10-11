function result=exportSpiceNetlist(schematicFile,options)
%EXPORTSPICENETLIST Export a KiCad schematic as a SPICE circuit netlist.
arguments
 schematicFile (1,1) string {mustBeFile}
 options.OutputFile (1,1) string=""
 options.Executable (1,1) string=""
 options.ModelOnly (1,1) logical=false
end
[sourceFolder,stem,extension]=fileparts(schematicFile);
if lower(string(extension))~=".kicad_sch",error("radia:kicad:SchematicRequired","Input must be a .kicad_sch schematic.");end
if strlength(options.OutputFile)==0,outputFile=fullfile(sourceFolder,stem+".cir");else,outputFile=options.OutputFile;end
folder=fileparts(outputFile);if strlength(folder)>0&&~isfolder(folder),mkdir(folder);end
executable=radia.kicad.findExecutable(Executable=options.Executable);
if options.ModelOnly,format="spicemodel";else,format="spice";end
command=join([quote(executable),"sch","export","netlist","--format",format,"--output",quote(outputFile),quote(schematicFile)]," ");
[status,commandOutput]=system(command);
if status~=0||~isfile(outputFile),error("radia:kicad:ExportFailed","KiCad SPICE export failed (status %d).\n%s",status,commandOutput);end
result=struct("schema","radia.kicad.spice_netlist.v1","source_schematic",string(schematicFile), ...
 "netlist_file",string(outputFile),"format",format,"kicad_cli",executable,"command",command, ...
 "command_output",string(commandOutput),"dependencies",radia.ltspice.collectDependencies(outputFile));
end
function value=quote(value),value='"'+replace(string(value),'"','\"')+'"';end
