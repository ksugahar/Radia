function inventory = inspectVolBoundaries(meshFile)
%INSPECTVOLBOUNDARIES Read ordered boundary names from a Netgen .vol file.

arguments
    meshFile (1,1) string {mustBeFile}
end
if ~endsWith(meshFile,".vol",IgnoreCase=true)
    error("radia:simulink:VolExtension","MeshFile must be a Netgen .vol file.");
end
lines=splitlines(string(fileread(meshFile))); trimmed=strip(lines);
section=find(strcmpi(trimmed,"bcnames"),1,"first");
if isempty(section) || section==numel(lines)
    error("radia:simulink:VolBoundaries","The .vol file has no bcnames section.");
end
count=str2double(trimmed(section+1));
if ~isscalar(count) || ~isfinite(count) || count<1 || fix(count)~=count || ...
        section+1+count>numel(lines)
    error("radia:simulink:VolBoundaries","The .vol boundary section is invalid.");
end
ids=zeros(count,1,"uint32"); names=strings(count,1);
for k=1:count
    tokens=regexp(char(lines(section+1+k)),"^\s*(\d+)\s+(.+?)\s*$", ...
        "tokens","once");
    if isempty(tokens)
        error("radia:simulink:VolBoundaries","Invalid boundary row %d.",k);
    end
    ids(k)=uint32(str2double(tokens{1})); names(k)=string(tokens{2});
end
if any(ids<1) || numel(unique(ids))~=count || numel(unique(names))~=count || ...
        numel(unique(lower(names)))~=count
    error("radia:simulink:VolBoundaries", ...
        "Boundary ids and case-insensitive names must be unique.");
end
[ids,order]=sort(ids); names=names(order);
inventory=struct("schema","radia.simulink.vol-boundaries.v1", ...
    "mesh_file",meshFile,"boundary_ids",ids,"boundary_names",names, ...
    "boundary_count",uint16(count));
end
