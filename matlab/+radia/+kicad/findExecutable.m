function executable=findExecutable(options)
%FINDEXECUTABLE Locate the current KiCad command-line executable.
arguments,options.Executable (1,1) string="";end
if strlength(options.Executable)>0
 if ~isfile(options.Executable),error("radia:kicad:ExecutableNotFound","kicad-cli does not exist: %s",options.Executable);end
 executable=options.Executable;return
end
parts=split(string(getenv("PATH")),pathsep);
for folder=parts(:)',candidate=fullfile(folder,"kicad-cli.exe");if isfile(candidate),executable=string(candidate);return,end,end
candidates=strings(0,1);
if ispc
 root="C:\Program Files\KiCad";
 if isfolder(root)
  entries=dir(root);entries=entries([entries.isdir]&~startsWith({entries.name},'.'));
  candidates=strings(numel(entries),1);
  for k=1:numel(entries),candidates(k)=fullfile(entries(k).folder,entries(k).name,"bin","kicad-cli.exe");end
 end
end
match=candidates(isfile(candidates));
if isempty(match),error("radia:kicad:ExecutableNotFound","Current KiCad kicad-cli was not found. Install KiCad or pass Executable=... explicitly.");end
executable=match(1);
end
