function contract = registerFieldStudy(contractVariable, options)
%REGISTERFIELDSTUDY Register setup metadata and fixed runtime Bus in base workspace.

arguments
    contractVariable (1,1) string = "radia_field_study_contract"
    options.RequestFile (1,1) string = ""
end
contract=evalin("base",contractVariable);
if ~isfield(contract,"schema") || string(contract.schema)~="radia.simulink.field-study.v1"
    error("radia:simulink:FieldStudyContract", ...
        "The selected variable must come from compileFieldStudy.");
end
radia.simulink.makeFieldStudyBusObject(contract.runtime,Name="RadiaStudyBus");
assignin("base","radia_field_study_bus",contract.runtime);
if strlength(options.RequestFile)>0
    radia.simulink.writeFieldStudyRequest(contract,options.RequestFile);
end
end
