classdef SpiceCircuit < radia.ltspice.SpiceEditor
    %SPICECIRCUIT Editable SPICE circuit compatibility facade.
    methods
        function obj = SpiceCircuit(path)
            obj@radia.ltspice.SpiceEditor(path);
        end
    end
end
