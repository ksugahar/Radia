function result=buildLTspiceBlock(modelName,schematicFile,options)
%BUILDLTSPICEBLOCK Export KiCad and add its LTspice circuit to Simulink.
arguments
 modelName (1,1) string
 schematicFile (1,1) string {mustBeFile}
 options.OutputDirectory (1,1) string="C:\temp\radia_kicad_ltspice"
 options.KiCadExecutable (1,1) string=""
 options.InputNames (:,1) string="control"
 options.OutputTraces (:,1) string="V(out)"
 options.SampleTime_s (1,1) double {mustBePositive}=1e-3
 options.MaxStep_s (1,1) double {mustBePositive}=inf
 options.Timeout_s (1,1) double {mustBePositive}=300
 options.LTspiceExecutable (1,1) string=""
 options.Save (1,1) logical=true
end
prepared=radia.kicad.prepareLTspice(schematicFile,OutputDirectory=options.OutputDirectory, ...
 KiCadExecutable=options.KiCadExecutable,LTspiceExecutable=options.LTspiceExecutable);
block=radia.simulink.buildLTspiceBlock(modelName,Netlist=prepared.netlist_file, ...
 InputNames=options.InputNames,OutputTraces=options.OutputTraces,SampleTime_s=options.SampleTime_s, ...
 MaxStep_s=options.MaxStep_s,Timeout_s=options.Timeout_s,Executable=options.LTspiceExecutable,Save=options.Save);
result=prepared;result.schema="radia.kicad.simulink_ltspice.v1";result.block_path=block;
end
