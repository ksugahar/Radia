function contract = compileWindingDictionary(windings, materialContract, options)
%COMPILEWINDINGDICTIONARY Compile winding/terminal setup to a fixed ABI.

arguments
    windings
    materialContract (1,1) struct
    options.MaxWindings (1,1) double {mustBeInteger,mustBePositive} = 16
    options.MaxRegionsPerWinding (1,1) double {mustBeInteger,mustBePositive} = 16
    options.MaxTerminals (1,1) double {mustBeInteger,mustBePositive} = 32
    options.AllowSharedRegions (1,1) logical = false
end
if ~isfield(materialContract,"schema") || ...
        string(materialContract.schema) ~= "radia.simulink.material-dictionary.v1" || ...
        ~isfield(materialContract,"region_names") || ~isfield(materialContract,"mesh")
    error("radia:simulink:WindingMaterialContract", ...
        "materialContract must come from compileMaterialDictionary.");
end
[names,values] = mappingEntries(windings);
if isempty(names)
    error("radia:simulink:WindingDictionaryEmpty", ...
        "The winding dictionary must not be empty.");
end
assertNames(names,"winding");
[names,order] = sort(names); values = values(order);
if numel(names) > options.MaxWindings
    error("radia:simulink:WindingCapacity", ...
        "Winding count %d exceeds MaxWindings=%d.",numel(names),options.MaxWindings);
end

specs = cell(numel(names),1);
terminalNames = strings(0,1);
allRegions = strings(0,1);
for k = 1:numel(names)
    specs{k} = normalizeSpec(values{k});
    if numel(specs{k}.regions) > options.MaxRegionsPerWinding
        error("radia:simulink:WindingRegionCapacity", ...
            "Winding '%s' exceeds MaxRegionsPerWinding=%d.", ...
            names(k),options.MaxRegionsPerWinding);
    end
    missing = setdiff(specs{k}.regions,materialContract.region_names,"stable");
    if ~isempty(missing)
        error("radia:simulink:WindingUnknownRegion", ...
            "Winding '%s' refers to unknown .vol regions: %s.", ...
            names(k),join(missing,","));
    end
    terminalNames = [terminalNames;specs{k}.positive_terminal;specs{k}.negative_terminal]; %#ok<AGROW>
    allRegions = [allRegions;specs{k}.regions]; %#ok<AGROW>
end
if ~options.AllowSharedRegions && numel(unique(allRegions)) ~= numel(allRegions)
    error("radia:simulink:WindingSharedRegion", ...
        "A .vol region may belong to only one winding unless AllowSharedRegions=true.");
end
terminalNames = sort(unique(terminalNames));
assertNames(terminalNames,"terminal");
if numel(terminalNames) > options.MaxTerminals
    error("radia:simulink:WindingTerminalCapacity", ...
        "Terminal count %d exceeds MaxTerminals=%d.", ...
        numel(terminalNames),options.MaxTerminals);
end

maxWindings = options.MaxWindings;
maxRegions = options.MaxRegionsPerWinding;
runtime = struct("schema_version",uint16(1), ...
    "winding_count",uint16(numel(names)), ...
    "winding_active",false(maxWindings,1), ...
    "turn_count",zeros(maxWindings,1,"uint32"), ...
    "polarity",zeros(maxWindings,1,"int8"), ...
    "parallel_paths",zeros(maxWindings,1,"uint16"), ...
    "effective_turns",zeros(maxWindings,1), ...
    "resistance_ohm",zeros(maxWindings,1), ...
    "positive_terminal_index",zeros(maxWindings,1,"uint16"), ...
    "negative_terminal_index",zeros(maxWindings,1,"uint16"), ...
    "region_count",zeros(maxWindings,1,"uint16"), ...
    "region_id",zeros(maxRegions,maxWindings,"uint32"), ...
    "region_polarity",zeros(maxRegions,maxWindings,"int8"), ...
    "terminal_count",uint16(numel(terminalNames)));

meshRegionIds = materialContract.mesh.region_ids;
meshRegionNames = materialContract.mesh.region_names;
if isempty(meshRegionIds)
    meshRegionIds = uint32((1:numel(materialContract.region_names)).');
    meshRegionNames = materialContract.region_names;
end
for k = 1:numel(names)
    spec = specs{k};
    runtime.winding_active(k) = true;
    runtime.turn_count(k) = spec.turns;
    runtime.polarity(k) = spec.polarity;
    runtime.parallel_paths(k) = spec.parallel_paths;
    runtime.effective_turns(k) = double(spec.polarity)*double(spec.turns)/double(spec.parallel_paths);
    runtime.resistance_ohm(k) = spec.resistance_ohm;
    runtime.positive_terminal_index(k) = uint16(find(terminalNames==spec.positive_terminal,1));
    runtime.negative_terminal_index(k) = uint16(find(terminalNames==spec.negative_terminal,1));
    runtime.region_count(k) = uint16(numel(spec.regions));
    for j = 1:numel(spec.regions)
        index = find(meshRegionNames==spec.regions(j),1);
        runtime.region_id(j,k) = meshRegionIds(index);
        runtime.region_polarity(j,k) = spec.region_polarity(j);
    end
end

contract = struct("schema","radia.simulink.winding-dictionary.v1", ...
    "mesh_sha256",materialContract.mesh.mesh_sha256, ...
    "winding_names",names,"terminal_names",terminalNames, ...
    "runtime",runtime,"runtime_policy",struct( ...
        "fixed_width",true,"strings_per_step",false, ...
        "dictionary_lookup_per_step",false,"python_per_step",false), ...
    "sign_convention",struct( ...
        "positive_voltage","positive terminal minus negative terminal", ...
        "positive_current","enters positive terminal", ...
        "positive_torque","increasing rotor angle"));
end

function spec = normalizeSpec(value)
if iscell(value) && isscalar(value),value=value{1};end
if ~isstruct(value) || ~isscalar(value)
    error("radia:simulink:WindingValue", ...
        "Every winding dictionary value must be a scalar struct.");
end
required = ["regions","turns","polarity","parallel_paths", ...
    "resistance_ohm","positive_terminal","negative_terminal"];
for field = required
    if ~isfield(value,field)
        error("radia:simulink:WindingMissingField", ...
            "Winding value is missing '%s'.",field);
    end
end
description = "";
if isfield(value,"description"),description=string(value.description);end
regionPolarity=zeros(0,1);
if isfield(value,"region_polarity"),regionPolarity=double(value.region_polarity);end
spec = radia.simulink.makeWindingSpec(Regions=string(value.regions), ...
    RegionPolarity=regionPolarity, ...
    Turns=double(value.turns),Polarity=double(value.polarity), ...
    ParallelPaths=double(value.parallel_paths), ...
    Resistance_ohm=double(value.resistance_ohm), ...
    PositiveTerminal=string(value.positive_terminal), ...
    NegativeTerminal=string(value.negative_terminal),Description=description);
end

function [names,values] = mappingEntries(mapping)
if isa(mapping,"dictionary")
    names=string(keys(mapping)); names=names(:); values=cell(numel(names),1);
    for k=1:numel(names),values{k}=mapping(names(k));end
elseif isstruct(mapping) && isscalar(mapping)
    names=string(fieldnames(mapping)); values=cell(numel(names),1);
    for k=1:numel(names),values{k}=mapping.(names(k));end
else
    error("radia:simulink:WindingMappingType", ...
        "windings must be a MATLAB dictionary or scalar struct.");
end
end

function assertNames(names,label)
if any(strlength(names)==0) || numel(unique(names))~=numel(names) || ...
        numel(unique(lower(names)))~=numel(names)
    error("radia:simulink:WindingNames", ...
        "%s names must be nonempty, unique, and case-distinct.",label);
end
end
