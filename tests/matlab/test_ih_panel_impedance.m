function tests = test_ih_panel_impedance
tests = functiontests(localfunctions);
end

function testWritesComplexValuesWithoutAveraging(testCase)
folder = fullfile("C:\temp", "radia-zs-" + string(java.util.UUID.randomUUID));
mkdir(folder);
cleanup = onCleanup(@() rmdir(folder,"s"));
template = fullfile(folder,"template.json");
output = fullfile(folder,"zs.json");
data = struct("schema","radia.panel_surface_impedance.v1", ...
    "layout","BND-element-order","unit","ohm","frequency_hz",1000, ...
    "surface_sha256","test-identity","centroids_m",[0 0 0;1 0 0]);
file = fopen(template,"w"); fwrite(file,jsonencode(data)); fclose(file);
radia.simulink.writeIHPanelImpedance(template,output,[.01+.02i;.03+.04i]);
result = jsondecode(fileread(output));
verifyEqual(testCase,result.real_ohm,[.01;.03]);
verifyEqual(testCase,result.imag_ohm,[.02;.04]);
verifyEqual(testCase,string(result.surface_sha256),data.surface_sha256);
verifyError(testCase,@() radia.simulink.writeIHPanelImpedance(template,output,-ones(2,1)), ...
    "radia:simulink:IHPanelZsValues");
verifyError(testCase,@() radia.simulink.writeIHPanelImpedance(template,output,1), ...
    "radia:simulink:IHPanelZsValues");
end
