function requestPath = writeFieldStudyRequest(contract, requestPath)
%WRITEFIELDSTUDYREQUEST Write the closed-world request consumed by Field Study.

arguments
    contract (1,1) struct
    requestPath (1,1) string
end
if ~isfield(contract,"schema") || string(contract.schema)~="radia.simulink.field-study.v1"
    error("radia:simulink:FieldStudyContract", ...
        "contract must come from compileFieldStudy.");
end
parent=fileparts(requestPath);if strlength(parent)>0 && ~isfolder(parent),mkdir(parent);end
request=jsonReadyRequest(contract.request);
text=jsonencode(request,PrettyPrint=true);
file=fopen(requestPath,"w","n","UTF-8");
if file<0,error("radia:simulink:FieldStudyWrite","Cannot open %s.",requestPath);end
cleanup=onCleanup(@()fclose(file));count=fprintf(file,"%s\n",text);
if count<=0,error("radia:simulink:FieldStudyWrite","Could not write %s.",requestPath);end
clear cleanup
requestPath=string(requestPath);
end

function request=jsonReadyRequest(request)
if ~isfield(request,"branches"),return;end
request.branches=num2cell(request.branches(:));
currents=request.branch_current_a;
request.branch_current_a=mat2cell(currents,ones(size(currents,1),1),size(currents,2));
end
