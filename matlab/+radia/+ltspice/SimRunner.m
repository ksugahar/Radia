classdef SimRunner < handle
    %SIMRUNNER MATLAB-native LTspice run manager.
    properties, OutputFolder (1,1) string="C:\temp\radia_ltspice_runs"; Executable (1,1) string=""; end
    methods
        function obj=SimRunner(options)
            arguments, options.OutputFolder (1,1) string="C:\temp\radia_ltspice_runs"; options.Executable (1,1) string=""; end
            obj.OutputFolder=options.OutputFolder; obj.Executable=options.Executable;
        end
        function result=runNow(obj,netlist,options)
            arguments, obj; netlist (1,1) string {mustBeFile}; options.Parameters (1,1) struct=struct(); options.RunName (1,1) string="run", end
            result=radia.ltspice.run(netlist,Parameters=options.Parameters,Executable=obj.Executable,OutputDirectory=fullfile(obj.OutputFolder,options.RunName));
            result.raw=radia.ltspice.RawRead(result.raw_file); result.log_reader=radia.ltspice.LogReader(result.log_file);
        end
        function results=runMany(obj,netlist,parameterSets,options)
            arguments, obj; netlist (1,1) string {mustBeFile}; parameterSets (1,:) cell; options.UseParallel (1,1) logical=false, end
            results=cell(size(parameterSets));
            if options.UseParallel
                if isempty(ver('parallel')), error("radia:ltspice:ParallelToolbox","Parallel Computing Toolbox is required."); end
                root=obj.OutputFolder; exe=obj.Executable;
                parfor k=1:numel(parameterSets)
                    results{k}=radia.ltspice.run(netlist,Parameters=parameterSets{k},Executable=exe,OutputDirectory=fullfile(root,sprintf("run_%06d",k)));
                end
            else
                for k=1:numel(parameterSets), results{k}=obj.runNow(netlist,Parameters=parameterSets{k},RunName=sprintf("run_%06d",k)); end
            end
        end
    end
end
