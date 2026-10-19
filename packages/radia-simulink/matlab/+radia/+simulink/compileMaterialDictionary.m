function contract = compileMaterialDictionary(materials, options)
%COMPILEMATERIALDICTIONARY Compile named materials to a fixed-width ABI.
%   The MATLAB dictionary and region names are setup metadata. The returned
%   runtime struct contains only fixed-size numeric/logical values suitable
%   for a Simulink.Bus, MEX S-Function, or generated code.

arguments
    materials
    options.RegionMaterials = []
    options.MeshFile (1,1) string = ""
    options.ExpectedMeshSHA256 (1,1) string = ""
    options.MaxMaterials (1,1) double {mustBeInteger,mustBePositive} = 32
    options.MaxRegions (1,1) double {mustBeInteger,mustBePositive} = 128
    options.MaxBHPoints (1,1) double {mustBeInteger,mustBePositive} = 256
    options.MaxHysteresisParameters (1,1) double {mustBeInteger,mustBePositive} = 16
end

[materialNames,materialValues] = mappingEntries(materials,"materials");
materialCount = numel(materialNames);
if materialCount == 0
    error("radia:simulink:MaterialDictionaryEmpty", ...
        "The material dictionary must not be empty.");
end
if materialCount > options.MaxMaterials
    error("radia:simulink:MaterialCapacity", ...
        "Material count %d exceeds MaxMaterials=%d.",materialCount,options.MaxMaterials);
end
assertNames(materialNames,"material");
[materialNames,order] = sort(materialNames);
materialValues = materialValues(order);

maxMaterials = options.MaxMaterials;
maxRegions = options.MaxRegions;
maxBH = options.MaxBHPoints;
maxHysteresis = options.MaxHysteresisParameters;
runtime = struct( ...
    "schema_version",uint16(1), ...
    "material_count",uint16(materialCount), ...
    "material_active",false(maxMaterials,1), ...
    "kind_flags",zeros(maxMaterials,1,"uint16"), ...
    "mu_r",ones(maxMaterials,1), ...
    "conductivity_S_per_m",zeros(maxMaterials,1), ...
    "relative_permittivity",ones(maxMaterials,1), ...
    "remanence_T",zeros(maxMaterials,3), ...
    "density_kg_per_m3",zeros(maxMaterials,1), ...
    "specific_heat_J_per_kgK",zeros(maxMaterials,1), ...
    "thermal_conductivity_W_per_mK",zeros(maxMaterials,1), ...
    "bh_count",zeros(maxMaterials,1,"uint16"), ...
    "bh_B_T",zeros(maxBH,maxMaterials), ...
    "bh_H_A_per_m",zeros(maxBH,maxMaterials), ...
    "hysteresis_model_code",zeros(maxMaterials,1,"uint16"), ...
    "hysteresis_parameter_count",zeros(maxMaterials,1,"uint16"), ...
    "hysteresis_parameters",zeros(maxHysteresis,maxMaterials), ...
    "region_count",uint16(0), ...
    "region_active",false(maxRegions,1), ...
    "region_id",zeros(maxRegions,1,"uint32"), ...
    "region_material_index",zeros(maxRegions,1,"uint16"));

for k = 1:materialCount
    spec = normalizeSpec(materialValues{k});
    B = spec.bh_B_T(:); H = spec.bh_H_A_per_m(:);
    if numel(B) > maxBH
        error("radia:simulink:MaterialBHCapacity", ...
            "Material '%s' has %d B-H points; MaxBHPoints=%d.", ...
            materialNames(k),numel(B),maxBH);
    end
    parameters = spec.hysteresis_parameters(:);
    if numel(parameters) > maxHysteresis
        error("radia:simulink:MaterialHysteresisCapacity", ...
            "Material '%s' has %d hysteresis parameters; capacity=%d.", ...
            materialNames(k),numel(parameters),maxHysteresis);
    end
    flags = uint16(0);
    if isempty(B) && abs(spec.mu_r-1) > eps, flags = bitor(flags,uint16(1)); end
    if ~isempty(B), flags = bitor(flags,uint16(2)); end
    if spec.conductivity_S_per_m > 0, flags = bitor(flags,uint16(4)); end
    if any(spec.remanence_T ~= 0), flags = bitor(flags,uint16(8)); end
    if abs(spec.relative_permittivity-1) > eps, flags = bitor(flags,uint16(16)); end
    if spec.density_kg_per_m3 > 0, flags = bitor(flags,uint16(32)); end
    if spec.hysteresis_model == "energy", flags = bitor(flags,uint16(64)); end

    runtime.material_active(k) = true;
    runtime.kind_flags(k) = flags;
    runtime.mu_r(k) = spec.mu_r;
    runtime.conductivity_S_per_m(k) = spec.conductivity_S_per_m;
    runtime.relative_permittivity(k) = spec.relative_permittivity;
    runtime.remanence_T(k,:) = spec.remanence_T;
    runtime.density_kg_per_m3(k) = spec.density_kg_per_m3;
    runtime.specific_heat_J_per_kgK(k) = spec.specific_heat_J_per_kgK;
    runtime.thermal_conductivity_W_per_mK(k) = spec.thermal_conductivity_W_per_mK;
    runtime.bh_count(k) = uint16(numel(B));
    runtime.bh_B_T(1:numel(B),k) = B;
    runtime.bh_H_A_per_m(1:numel(H),k) = H;
    runtime.hysteresis_model_code(k) = uint16(spec.hysteresis_model == "energy");
    runtime.hysteresis_parameter_count(k) = uint16(numel(parameters));
    runtime.hysteresis_parameters(1:numel(parameters),k) = parameters;
end

if strlength(options.MeshFile) > 0
    mesh = radia.simulink.inspectVolMaterials(options.MeshFile);
    regionNames = mesh.region_names;
    regionIds = mesh.region_ids;
    if strlength(options.ExpectedMeshSHA256) > 0 && ...
            ~strcmpi(mesh.mesh_sha256,options.ExpectedMeshSHA256)
        error("radia:simulink:MeshDigest", ...
            "The .vol digest differs from ExpectedMeshSHA256.");
    end
else
    mesh = struct("schema","radia.simulink.vol-materials.v1", ...
        "mesh_file","","mesh_sha256","", ...
        "region_ids",zeros(0,1,"uint32"),"region_names",strings(0,1), ...
        "region_count",uint16(0));
    if isEmptyMapping(options.RegionMaterials)
        regionNames = materialNames;
    else
        [regionNames,~] = mappingEntries(options.RegionMaterials,"RegionMaterials");
        regionNames = sort(regionNames);
    end
    regionIds = uint32((1:numel(regionNames)).');
end
regionCount = numel(regionNames);
if regionCount > maxRegions
    error("radia:simulink:RegionCapacity", ...
        "Region count %d exceeds MaxRegions=%d.",regionCount,maxRegions);
end
assertNames(regionNames,"region");

if isEmptyMapping(options.RegionMaterials)
    targetNames = regionNames;
else
    [mapNames,mapValues] = mappingEntries(options.RegionMaterials,"RegionMaterials");
    assertNames(mapNames,"region assignment");
    missing = setdiff(regionNames,mapNames,"stable");
    extra = setdiff(mapNames,regionNames,"stable");
    if ~isempty(missing) || ~isempty(extra)
        error("radia:simulink:RegionCoverage", ...
            "RegionMaterials must exactly cover .vol regions. Missing=[%s], extra=[%s].", ...
            displayList(missing),displayList(extra));
    end
    targetNames = strings(regionCount,1);
    for k = 1:regionCount
        index = find(mapNames == regionNames(k),1);
        value = unwrap(mapValues{index});
        if ~(ischar(value) || (isstring(value) && isscalar(value)))
            error("radia:simulink:RegionMaterialValue", ...
                "RegionMaterials values must be material-name strings.");
        end
        targetNames(k) = string(value);
    end
end

indices = zeros(regionCount,1,"uint16");
for k = 1:regionCount
    index = find(materialNames == targetNames(k),1);
    if isempty(index)
        error("radia:simulink:UnknownRegionMaterial", ...
            "Region '%s' refers to unknown material '%s'.",regionNames(k),targetNames(k));
    end
    indices(k) = uint16(index);
end
runtime.region_count = uint16(regionCount);
runtime.region_active(1:regionCount) = true;
runtime.region_id(1:regionCount) = regionIds;
runtime.region_material_index(1:regionCount) = indices;

contract = struct("schema","radia.simulink.material-dictionary.v1", ...
    "material_names",materialNames,"region_names",regionNames, ...
    "region_material_names",targetNames,"mesh",mesh,"runtime",runtime, ...
    "runtime_policy",struct("fixed_width",true,"strings_per_step",false, ...
        "dictionary_lookup_per_step",false,"python_per_step",false));
end

function spec = normalizeSpec(value)
value = unwrap(value);
if ~isstruct(value) || ~isscalar(value)
    error("radia:simulink:MaterialValue", ...
        "Every material dictionary value must be a scalar struct.");
end
allowed = ["mu_r","conductivity_S_per_m","relative_permittivity", ...
    "remanence_T","density_kg_per_m3","specific_heat_J_per_kgK", ...
    "thermal_conductivity_W_per_mK","bh_B_T","bh_H_A_per_m", ...
    "hysteresis_model","hysteresis_parameters","description"];
unknown = setdiff(string(fieldnames(value)),allowed);
if ~isempty(unknown)
    error("radia:simulink:MaterialUnknownField", ...
        "Unknown material field(s): %s",join(unknown,","));
end
spec = radia.simulink.makeMaterialSpec( ...
    MuR=fieldOr(value,"mu_r",1), ...
    Conductivity_S_per_m=fieldOr(value,"conductivity_S_per_m",0), ...
    RelativePermittivity=fieldOr(value,"relative_permittivity",1), ...
    Remanence_T=fieldOr(value,"remanence_T",[0 0 0]), ...
    Density_kg_per_m3=fieldOr(value,"density_kg_per_m3",0), ...
    SpecificHeat_J_per_kgK=fieldOr(value,"specific_heat_J_per_kgK",0), ...
    ThermalConductivity_W_per_mK=fieldOr(value,"thermal_conductivity_W_per_mK",0), ...
    BH_B_T=fieldOr(value,"bh_B_T",zeros(0,1)), ...
    BH_H_A_per_m=fieldOr(value,"bh_H_A_per_m",zeros(0,1)), ...
    HysteresisModel=string(fieldOr(value,"hysteresis_model","none")), ...
    HysteresisParameters=fieldOr(value,"hysteresis_parameters",zeros(0,1)), ...
    Description=string(fieldOr(value,"description","")));
end

function value = fieldOr(data,name,defaultValue)
if isfield(data,name), value = data.(name); else, value = defaultValue; end
end

function value = unwrap(value)
if iscell(value) && isscalar(value), value = value{1}; end
end

function [names,values] = mappingEntries(mapping,label)
if isa(mapping,"dictionary")
    names = string(keys(mapping));
    names = names(:);
    values = cell(numel(names),1);
    for k = 1:numel(names), values{k} = unwrap(mapping(names(k))); end
elseif isstruct(mapping) && isscalar(mapping)
    names = string(fieldnames(mapping));
    values = cell(numel(names),1);
    for k = 1:numel(names), values{k} = mapping.(names(k)); end
else
    error("radia:simulink:MappingType", ...
        "%s must be a MATLAB dictionary or scalar struct.",label);
end
end

function tf = isEmptyMapping(mapping)
tf = isempty(mapping) || (isa(mapping,"dictionary") && isempty(keys(mapping))) || ...
    (isstruct(mapping) && isscalar(mapping) && isempty(fieldnames(mapping)));
end

function assertNames(names,label)
if any(strlength(names)==0) || numel(unique(names)) ~= numel(names)
    error("radia:simulink:DictionaryNames","%s names must be nonempty and unique.",label);
end
if numel(unique(lower(names))) ~= numel(names)
    error("radia:simulink:DictionaryCaseCollision", ...
        "%s names that differ only by case are not allowed.",label);
end
end

function text = displayList(values)
if isempty(values), text = ""; else, text = join(values,","); end
end
