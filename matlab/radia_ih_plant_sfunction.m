function radia_ih_plant_sfunction(block)
%RADIA_IH_PLANT_SFUNCTION Simulink loader for the packaged Radia block.
%   Simulink's Level-2 S-function parameter accepts a plain function name;
%   delegate to the namespaced implementation used by the MATLAB API.

radia.simulink.ihPlantSFunction(block);
end
