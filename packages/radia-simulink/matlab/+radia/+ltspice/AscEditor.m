classdef AscEditor < radia.ltspice.SchematicEditor
    %ASCEDITOR PyLTSpice-compatible name for the MATLAB schematic editor.
    methods
        function obj = AscEditor(path)
            obj@radia.ltspice.SchematicEditor(path);
        end
    end
end
