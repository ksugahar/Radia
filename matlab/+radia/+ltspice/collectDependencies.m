function manifest=collectDependencies(netlistFile)
%COLLECTDEPENDENCIES Recursively collect local .include/.inc/.lib dependencies.
arguments, netlistFile (1,1) string {mustBeFile}, end
root=string(netlistFile); visited=strings(0,1); local=strings(0,1); external=strings(0,1); unresolved=strings(0,1);
walk(root);
manifest=struct("schema","radia.ltspice.dependencies.v1","root",root,"local_files",unique(local,"stable"),"absolute_external",unique(external,"stable"),"unresolved_library_names",unique(unresolved,"stable"));
    function walk(path)
        path=string(java.io.File(char(path)).getCanonicalPath()); if any(visited==path),return,end; visited(end+1)=path; local(end+1)=path;
        text=string(fileread(path)); rows=splitlines(text); base=string(fileparts(path));
        for i=1:numel(rows)
            token=regexp(char(rows(i)),'^\s*\.(include|inc|lib)\s+(.+?)\s*(?:;.*)?$','tokens','once','ignorecase');
            if isempty(token),continue,end; value=strtrim(string(token{2}));
            if (startsWith(value,'"')&&endsWith(value,'"'))||(startsWith(value,"'")&&endsWith(value,"'")),value=extractBetween(value,2,strlength(value)-1);end
            candidate=value; if ~isfile(candidate),candidate=fullfile(base,value);end
            if isfile(candidate)
                canonical=string(java.io.File(char(candidate)).getCanonicalPath());
                if startsWith(lower(canonical),lower(base)), walk(canonical); else, external(end+1)=canonical; walk(canonical); end
            elseif contains(value,["/","\"])||startsWith(value,".")
                error("radia:ltspice:MissingDependency","Missing LTspice dependency %s referenced by %s.",value,path);
            else
                unresolved(end+1)=value;
            end
        end
    end
end
