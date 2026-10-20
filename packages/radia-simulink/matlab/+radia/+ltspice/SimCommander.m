classdef SimCommander < radia.ltspice.SpiceEditor
    %SIMCOMMANDER Deprecated PyLTSpice convenience API retained for parity.
    properties (SetAccess=private), Runner (1,1) radia.ltspice.SimRunner; end
    methods
        function obj=SimCommander(path,options)
            arguments,path (1,1) string {mustBeFile};options.OutputFolder (1,1) string="C:\temp\radia_ltspice_runs";end
            obj@radia.ltspice.SpiceEditor(path);obj.Runner=radia.ltspice.SimRunner(OutputFolder=options.OutputFolder);
        end
        function task=run(obj,varargin),temporary=fullfile(obj.Runner.OutputFolder,"simcommander.cir");obj.saveAs(temporary);task=obj.Runner.run(temporary,varargin{:});end
        function answer=wait_completion(obj,varargin),answer=obj.Runner.waitCompletion(varargin{:});end
    end
end
