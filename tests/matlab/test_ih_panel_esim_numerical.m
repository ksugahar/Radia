function tests = test_ih_panel_esim_numerical
% The strong Python reference invokes the existing solve_panel_esim helper.
tests=functiontests(localfunctions);
end

function testNamedAssemblyMatchesIndependentPythonReference(testCase)
exerciseAssembly(testCase,0);
end

function testGenusOneStrongAssemblyMatchesIndependentPythonReference(testCase)
exerciseAssembly(testCase,1);
end

function exerciseAssembly(testCase,genus)
work=string(tempname()); mkdir(work);
cleanup=onCleanup(@() rmdir(work,"s"));
repo=string(fileparts(fileparts(fileparts(mfilename("fullpath")))));
python=radia.internal.resolveAssemblyPython("");
script=fullfile(repo,"validation_test","induction_heating","make_simulink_esim_fixture.py");
[status,output]=system('"'+python+'" "'+script+'" "'+work+'" --genus '+string(genus));
verifyEqual(testCase,status,0,output);
reference=jsondecode(fileread(fullfile(work,"reference.json")));
verifyEqual(testCase,reference.wp_genus,genus);
if genus==1, verifyTrue(testCase,reference.wp_loop_dof); end
options=struct("esim_bh_file",fullfile(work,"bh.txt"),"esim_reference_current_A",2, ...
    "coupling_mode","strong","python_executable",python);
config=radia.simulink.assembleIHPanelESIMFromGeometry(fullfile(work,"workpiece.vol"), ...
    fullfile(work,"coil.vol"),fullfile(work,"native.json"),options);
zs=jsondecode(fileread(string(config.surface_impedance.artifact_file)));
expected=complex(reference.esim_per_panel_Z_s_real,reference.esim_per_panel_Z_s_imag);
actual=complex(zs.real_ohm,zs.imag_ohm);
% Both routes solve identical discrete equations; relative vector norm gate.
verifyLessThanOrEqual(testCase,norm(actual-expected)/norm(expected),1e-9);
verifyEqual(testCase,config.unit_current.electromagnetic_power_W*4,reference.P_wp_W,"RelTol",1e-9);
verifyEqual(testCase,config.surface_impedance.zs_semantics,'frozen_at_reference_current');
metrics=struct("schema","radia.simulink.panel_esim.acceptance.v1", ...
    "workpiece_genus",genus,"coupling_mode","strong","backend","intree-dense", ...
    "evaluator","table","reference_current_A",2,"reference_phase_rad",0,"frequency_hz",7000, ...
    "panel_count",numel(expected),"relative_zs_l2_error",norm(actual-expected)/norm(expected), ...
    "relative_heat_error",abs(config.unit_current.electromagnetic_power_W*4-reference.P_wp_W)/reference.P_wp_W, ...
    "comparison_tolerance",1e-9,"constitutive_residual",config.surface_impedance.final_residual, ...
    "final_certification",config.surface_impedance.final_certification, ...
    "body_residual",reference.body_residual,"body_power_balance_relative_error",reference.body_power_balance_relative_error);
fid=fopen(fullfile(pwd,compose("panel-esim-acceptance-genus%d.json",genus)),"w");
fprintf(fid,"%s",jsonencode(metrics)); fclose(fid);
if genus==1, return; end
% Geometry update fingerprints the operating state and performs no repeat solve.
model="r09_esim_"+string(java.util.UUID.randomUUID()); model=replace(model,"-","_");
new_system(model);
modelCleanup=onCleanup(@() close_system(model,0));
block=radia.simulink.addIHGeometryUpdateBlock(model);
set_param(block,"wp_vol",char(fullfile(work,"workpiece.vol")), ...
    "coil_file",char(fullfile(work,"coil.vol")),"config_file",char(fullfile(work,"lifecycle.json")), ...
    "zs_mode","per-panel-esim","esim_bh_file",char(fullfile(work,"bh.txt")), ...
    "esim_reference_current_A","2","coupling_mode","strong","python_executable",char(python));
first=radia.simulink.updateIHGeometry(model);
second=radia.simulink.updateIHGeometry(model);
verifyTrue(testCase,first.rebuilt); verifyFalse(testCase,second.rebuilt);
verifyEqual(testCase,first.revision,second.revision);
set_param(block,"esim_reference_current_A","2.05");
third=radia.simulink.updateIHGeometry(model);
verifyTrue(testCase,third.rebuilt); verifyEqual(testCase,third.revision,first.revision+1);
current=jsondecode(fileread(fullfile(work,"lifecycle.json")));
writelines("corrupt",string(current.surface_impedance.artifact_file));
verifyError(testCase,@() radia.simulink.updateIHGeometry(model),"radia:simulink:IHPanelZsProvenance");
end
