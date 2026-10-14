function contract = registerWindingDictionary(windingsVariable, materialContractVariable, options)
%REGISTERWINDINGDICTIONARY Compile workspace winding data for Simulink.

arguments
    windingsVariable (1,1) string
    materialContractVariable (1,1) string = "radia_material_contract"
    options.MaxWindings (1,1) double {mustBeInteger,mustBePositive} = 16
    options.MaxRegionsPerWinding (1,1) double {mustBeInteger,mustBePositive} = 16
    options.MaxTerminals (1,1) double {mustBeInteger,mustBePositive} = 32
end
if ~isvarname(windingsVariable) || ~isvarname(materialContractVariable)
    error("radia:simulink:WindingVariable", ...
        "Winding and material-contract inputs must name base-workspace variables.");
end
windings=evalin("base",windingsVariable);
materials=evalin("base",materialContractVariable);
contract=radia.simulink.compileWindingDictionary(windings,materials, ...
    MaxWindings=options.MaxWindings, ...
    MaxRegionsPerWinding=options.MaxRegionsPerWinding, ...
    MaxTerminals=options.MaxTerminals);
assignin("base","radia_winding_contract",contract);
assignin("base","radia_winding_bus",contract.runtime);
radia.simulink.makeWindingBusObject(contract.runtime,Name="RadiaWindingBus");
radia.simulink.makeElectromechanicalBusObjects(MaxWindings=options.MaxWindings);
end
