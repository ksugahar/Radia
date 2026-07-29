function contract = registerMaterialDictionary(materialsVariable, regionMapVariable, meshFile, options)
%REGISTERMATERIALDICTIONARY Compile workspace dictionaries for Simulink.

arguments
    materialsVariable (1,1) string
    regionMapVariable (1,1) string = ""
    meshFile (1,1) string = ""
    options.ExpectedMeshSHA256 (1,1) string = ""
    options.MaxMaterials (1,1) double {mustBeInteger,mustBePositive} = 32
    options.MaxRegions (1,1) double {mustBeInteger,mustBePositive} = 128
    options.MaxBHPoints (1,1) double {mustBeInteger,mustBePositive} = 256
    options.MaxHysteresisParameters (1,1) double {mustBeInteger,mustBePositive} = 16
    options.ValidateMesh (1,1) logical = true
    options.MeshContract (1,1) string = ""
    options.ReportDirectory (1,1) string = ""
end
if ~isvarname(materialsVariable)
    error("radia:simulink:MaterialVariable", ...
        "MaterialsVariable must name one base-workspace variable.");
end
materials = evalin("base",materialsVariable);
regionMaterials = [];
if strlength(regionMapVariable) > 0
    if ~isvarname(regionMapVariable)
        error("radia:simulink:RegionMapVariable", ...
            "RegionMapVariable must name one base-workspace variable.");
    end
    regionMaterials = evalin("base",regionMapVariable);
end
if strlength(meshFile) == 0
    error("radia:simulink:MaterialMesh", ...
        "Material Dictionary requires a Netgen .vol mesh.");
end
if options.ValidateMesh
    radia.simulink.validateVolFiles(meshFile, ...
        ReportDirectory=options.ReportDirectory,Contract=options.MeshContract);
end
contract = radia.simulink.compileMaterialDictionary(materials, ...
    RegionMaterials=regionMaterials,MeshFile=meshFile, ...
    ExpectedMeshSHA256=options.ExpectedMeshSHA256, ...
    MaxMaterials=options.MaxMaterials,MaxRegions=options.MaxRegions, ...
    MaxBHPoints=options.MaxBHPoints, ...
    MaxHysteresisParameters=options.MaxHysteresisParameters);
assignin("base","radia_material_contract",contract);
assignin("base","radia_material_bus",contract.runtime);
radia.simulink.makeMaterialBusObject(contract.runtime,Name="RadiaMaterialBus");
end
