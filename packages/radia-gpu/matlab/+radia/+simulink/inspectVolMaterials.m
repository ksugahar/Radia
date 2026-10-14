function inventory = inspectVolMaterials(meshFile)
%INSPECTVOLMATERIALS Read ordered material regions and identity from a .vol.

arguments
    meshFile (1,1) string {mustBeFile}
end
if ~endsWith(meshFile,".vol",IgnoreCase=true)
    error("radia:simulink:VolExtension","MeshFile must be a Netgen .vol file.");
end

lines = splitlines(string(fileread(meshFile)));
trimmed = strip(lines);
section = find(strcmpi(trimmed,"materials"),1,"first");
if isempty(section) || section == numel(lines)
    error("radia:simulink:VolMaterials","The .vol file has no materials section.");
end
count = str2double(trimmed(section+1));
if ~isscalar(count) || ~isfinite(count) || count < 1 || fix(count) ~= count
    error("radia:simulink:VolMaterials","The .vol materials count is invalid.");
end
if section + 1 + count > numel(lines)
    error("radia:simulink:VolMaterials","The .vol materials section is truncated.");
end

ids = zeros(count,1,"uint32");
names = strings(count,1);
for k = 1:count
    tokens = regexp(char(lines(section+1+k)),"^\s*(\d+)\s+(.+?)\s*$", ...
        "tokens","once");
    if isempty(tokens)
        error("radia:simulink:VolMaterials", ...
            "Invalid material row %d in %s.",k,meshFile);
    end
    id = str2double(tokens{1});
    if ~isfinite(id) || id < 1 || fix(id) ~= id || id > intmax("uint32")
        error("radia:simulink:VolMaterials","Invalid material id at row %d.",k);
    end
    ids(k) = uint32(id);
    names(k) = string(tokens{2});
end
if numel(unique(ids)) ~= count || numel(unique(names)) ~= count
    error("radia:simulink:VolMaterials","Material ids and names must be unique.");
end
if numel(unique(lower(names))) ~= count
    error("radia:simulink:VolMaterialCaseCollision", ...
        "Material names that differ only by case are not allowed.");
end
[ids,order] = sort(ids);
names = names(order);

file = fopen(meshFile,"rb");
if file < 0, error("radia:simulink:VolRead","Cannot open %s.",meshFile); end
cleanup = onCleanup(@() fclose(file));
bytes = fread(file,inf,"*uint8");
digest = java.security.MessageDigest.getInstance("SHA-256");
digest.update(typecast(bytes,"int8"));
hashBytes = typecast(int8(digest.digest()),"uint8");
sha256 = lower(string(reshape(dec2hex(hashBytes,2).',1,[])));

inventory = struct("schema","radia.simulink.vol-materials.v1", ...
    "mesh_file",string(meshFile),"mesh_sha256",sha256, ...
    "region_ids",ids,"region_names",names,"region_count",uint16(count));
end
