function configureFieldStudyBlock(blockPath)
%CONFIGUREFIELDSTUDYBLOCK Configure the public study-contract library block.

arguments
    blockPath (1,1) string
end
mask=Simulink.Mask.get(blockPath);if isempty(mask),mask=Simulink.Mask.create(blockPath);end
mask.Type="Radia Field Study";
mask.Description= ...
    "Compile electrostatic, current-flow, steady-heat, or harmonic-eddy setup " + ...
    "to RadiaStudyBus. The Python/NGSolve solve runs once per explicit trigger, not per step.";
addParameter(mask,"study_contract_variable","Field study contract","'radia_field_study_contract'");
addParameter(mask,"request_file","Request JSON","''");
mask.Initialization= ...
    "radia.simulink.registerFieldStudy(string(study_contract_variable)," + ...
    "RequestFile=string(request_file));";
mask.Display="disp('Field Study');port_label('output',1,'study');";
set_param(blockPath,"UserData",struct("schema","radia.simulink.field-study.v1", ...
    "physics","electrostatic,current_flow,steady_heat,harmonic_eddy", ...
    "runtime_bus","RadiaStudyBus","mesh_format","Netgen .vol", ...
    "dictionary_lookup_per_step",false,"python_per_step",false, ...
    "batch_python_per_trigger",true),"UserDataPersistent","on");
end

function addParameter(mask,name,prompt,value)
if isempty(mask.getParameter(name))
    mask.addParameter(Type="edit",Name=name,Prompt=prompt,Value=value,Evaluate="on");
end
end
