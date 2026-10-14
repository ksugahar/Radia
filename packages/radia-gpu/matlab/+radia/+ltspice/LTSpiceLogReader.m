classdef LTSpiceLogReader < radia.ltspice.LogReader
    %LTSPICELOGREADER PyLTSpice-compatible log-reader name.
    methods
        function obj = LTSpiceLogReader(path)
            obj@radia.ltspice.LogReader(path);
        end
    end
end
