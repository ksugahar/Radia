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
    id=str2double(tokens{1});
    if ~isfinite(id) || id<1 || fix(id)~=id || id>double(intmax("uint32"))
        error("radia:simulink:VolBoundaries", ...
            "Boundary id at row %d is outside the uint32 contract.",k);
    end
    ids(k)=uint32(id); names(k)=string(tokens{2});
end
if any(ids<1) || numel(unique(ids))~=count || any(strlength(names)==0)
    error("radia:simulink:VolBoundaries", ...
        "Boundary ids must be unique and names must be nonempty.");
end
[ids,order]=sort(ids); names=names(order);
uniqueIds=zeros(0,1,"uint32");uniqueNames=strings(0,1);idGroups=cell(0,1);
for k=1:count
    index=find(strcmpi(uniqueNames,names(k)),1);
    if isempty(index)
        uniqueIds(end+1,1)=ids(k); %#ok<AGROW>
        uniqueNames(end+1,1)=names(k); %#ok<AGROW>
        idGroups{end+1,1}=ids(k); %#ok<AGROW>
    else
        if uniqueNames(index)~=names(k)
            error("radia:simulink:VolBoundaries", ...
                "Boundary names that differ only by case are ambiguous.");
        end
        idGroups{index}(end+1,1)=ids(k);
    end
end
inventory=struct("schema","radia.simulink.vol-boundaries.v1", ...
    "mesh_file",meshFile,"boundary_ids",uniqueIds, ...
    "boundary_names",uniqueNames,"boundary_id_groups",{idGroups}, ...
    "raw_boundary_ids",ids,"raw_boundary_names",names, ...
    "boundary_count",uint16(numel(uniqueNames)), ...
    "raw_boundary_id_count",uint16(count));
end
