function manifest=collectSchematicDependencies(schematicFile)
%COLLECTSCHEMATICDEPENDENCIES Collect custom symbols, hierarchy and model files.
arguments, schematicFile (1,1) string {mustBeFile}, end
root=canonical(schematicFile); rootFolder=string(fileparts(root)); files=strings(0,1); unresolved=strings(0,1); visited=strings(0,1);
walk(root);
manifest=struct("schema","radia.ltspice.schematic_dependencies.v1","root",root, ...
 "local_files",unique(files,"stable"),"unresolved",unique(unresolved,"stable"));
 function walk(path)
  path=canonical(path); if any(visited==path),return,end; visited(end+1)=path; files(end+1)=path;
  folder=string(fileparts(path)); rows=splitlines(string(fileread(path)));
  for i=1:numel(rows)
   row=strtrim(rows(i)); candidates=strings(0,1);
   symbol=regexp(char(row),'^SYMBOL\s+(\S+)','tokens','once','ignorecase');
   if ~isempty(symbol)
    token=string(symbol{1}); if ~endsWith(lower(token),".asy"),token=token+".asy";end
    candidates(end+1)=token;
   end
   attribute=regexp(char(row),'^SYMATTR\s+(?:ModelFile|SpiceModel)\s+(.+)$','tokens','once','ignorecase');
   if ~isempty(attribute)
    value=strtrim(string(attribute{1})); if startsWith(value,'"')&&endsWith(value,'"'),value=extractBetween(value,2,strlength(value)-1);end
    candidates(end+1)=value;
   end
   directive=regexp(char(row),'!\s*\.(?:include|inc|lib)\s+(?:"([^"]+)"|(\S+))','tokens','once','ignorecase');
   if ~isempty(directive),candidates(end+1)=string(directive{find(~cellfun(@isempty,directive),1)});end
   for candidate=candidates
    located=fullfile(folder,candidate);
    if isfile(located),walk(located);
    elseif any(contains(candidate,["/","\"]))||any(endsWith(lower(candidate),[".asc",".asy",".lib",".sub",".inc",".cir"])),unresolved(end+1)=candidate;end
   end
  end
 end
 function path=canonical(path),path=string(java.io.File(char(path)).getCanonicalPath());end
end
