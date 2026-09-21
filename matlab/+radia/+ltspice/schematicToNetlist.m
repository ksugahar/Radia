function info=schematicToNetlist(schematicFile,options)
%SCHEMATICTONETLIST Convert an LTspice .asc schematic using LTspice itself.
arguments
 schematicFile (1,1) string {mustBeFile}
 options.Executable (1,1) string=""
 options.OutputDirectory (1,1) string=""
 options.Timeout_s (1,1) double {mustBePositive}=60
end
[sourceFolder,stem,extension]=fileparts(schematicFile);
if lower(string(extension))~=".asc", error("radia:ltspice:SchematicRequired","schematicFile must be .asc."); end
output=options.OutputDirectory; if strlength(output)==0, output=string(tempname("C:\temp")); end
if ~isfolder(output), mkdir(output); end
dependencies=radia.ltspice.collectSchematicDependencies(schematicFile);
sourceFolder=string(java.io.File(char(sourceFolder)).getCanonicalPath());
for dependency=dependencies.local_files(:)'
 relative=extractAfter(dependency,strlength(string(sourceFolder))+1); if startsWith(relative,["\","/"]),relative=extractAfter(relative,1);end
 destination=fullfile(output,relative); destinationFolder=fileparts(destination); if ~isfolder(destinationFolder),mkdir(destinationFolder);end
 copyfile(dependency,destination);
end
staged=fullfile(output,stem+extension);
executable=radia.ltspice.findExecutable(Executable=options.Executable);
if ispc
 script="$p=Start-Process -FilePath '"+replace(executable,"'","''")+"' -ArgumentList @('-netlist','"+replace(string(staged),"'","''")+"') -Wait -PassThru -WindowStyle Hidden; exit $p.ExitCode";
 encoded=matlab.net.base64encode(unicode2native(char(script),'UTF-16LE'));
 p=System.Diagnostics.Process(); p.StartInfo.FileName='pwsh'; p.StartInfo.Arguments='-NoLogo -NoProfile -NonInteractive -EncodedCommand '+string(encoded); p.StartInfo.UseShellExecute=false; p.StartInfo.CreateNoWindow=true;
 if ~p.Start(), error("radia:ltspice:ProcessStart","Could not start LTspice schematic conversion."); end
 if ~p.WaitForExit(round(options.Timeout_s*1000))
  terminated=radia.ltspice.internal.terminateProcessTree(p);
  if ~terminated,error("radia:ltspice:TimeoutCleanup","LTspice schematic conversion timed out and its process tree could not be confirmed terminated.");end
  error("radia:ltspice:Timeout","LTspice schematic conversion exceeded Timeout_s=%g.",options.Timeout_s);
 end
 status=double(p.ExitCode);
else
 [status,~]=system(sprintf('"%s" -netlist "%s"',executable,staged));
end
netlist=fullfile(output,stem+".net");
if status~=0||~isfile(netlist), error("radia:ltspice:SchematicConversion","LTspice did not generate a netlist from %s.",schematicFile); end
info=struct("schema","radia.ltspice.schematic_netlist.v1","source",schematicFile,"staged_schematic",string(staged),"netlist",string(netlist),"executable",executable,"dependencies",dependencies);
end
