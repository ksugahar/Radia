function tests = test_ih_panel_esim_contract
% Frozen constitutive semantics and stale-certificate rejection; no fixture skips.
tests = functiontests(localfunctions);
end

function testDriveBandUsesPeakAmplitude(testCase)
config.surface_impedance = struct("mode","per-panel-esim", ...
    "reference_current_A",2,"drive_relative_band",.1);
radia.simulink.verifyIHESIMDrive(config,2);
radia.simulink.verifyIHESIMDrive(config,-2);
radia.simulink.verifyIHESIMDrive(config,2.2);
verifyError(testCase,@() radia.simulink.verifyIHESIMDrive(config,2.3),"radia:simulink:IHESIMDriveBand");
verifyError(testCase,@() radia.simulink.verifyIHESIMDrive(config,0),"radia:simulink:IHESIMDriveBand");
end

function testProvenanceRejectsStaleArtifactFrequencyAndMaterial(testCase)
work=string(tempname()); mkdir(work);
cleanup=onCleanup(@() rmdir(work,"s"));
bh=fullfile(work,"bh.txt"); writelines("0 0",bh);
identity=radia.simulink.fileFingerprint(bh);
options=struct("zs_mode","per-panel-esim","esim_bh_file",bh, ...
    "esim_reference_current_A","2","frequency_hz","7000", ...
    "esim_drive_relative_band","0.1","esim_panel_evaluator","table","esim_tolerance","0.001","esim_max_iter","15");
state=struct("mode","per-panel-esim","surface_sha256","synthetic-surface","bh_sha256",identity.sha256,"bh_source_file",bh, ...
    "reference_current_A",2,"reference_phase_rad",0,"frequency_hz",7000, ...
    "drive_relative_band",.1,"iterations",3,"final_residual",1e-5,"tolerance",.001, ...
    "final_certification","direct-all-panels","evaluator","table","zs_semantics","frozen_at_reference_current");
artifact=fullfile(work,"panel.json");
data=struct("schema","radia.panel_surface_impedance.v1","layout","BND-element-order", ...
    "unit","ohm","frequency_hz",7000,"surface_sha256","synthetic-surface","provenance",state);
writelines(jsonencode(data),artifact);
identity=radia.simulink.fileFingerprint(artifact);
state.artifact_file=artifact; state.artifact_sha256=identity.sha256;
config.surface_impedance=state;
file=fullfile(work,"config.json"); writelines(jsonencode(config),file);
radia.simulink.verifyIHPanelImpedanceProvenance(file,"","",options);
for invalid = [-1,NaN]
    malformed=config; malformed.surface_impedance.final_residual=invalid; writelines(jsonencode(malformed),file);
    verifyError(testCase,@() radia.simulink.verifyIHPanelImpedanceProvenance(file,"","",options),"radia:simulink:IHPanelZsProvenance");
end
for invalid = [1.5,16,Inf]
    malformed=config; malformed.surface_impedance.iterations=invalid; writelines(jsonencode(malformed),file);
    verifyError(testCase,@() radia.simulink.verifyIHPanelImpedanceProvenance(file,"","",options),"radia:simulink:IHPanelZsProvenance");
end
malformed=config; malformed.surface_impedance=rmfield(malformed.surface_impedance,"mode"); writelines(jsonencode(malformed),file);
verifyError(testCase,@() radia.simulink.verifyIHPanelImpedanceProvenance(file,"","",options),"radia:simulink:IHPanelZsProvenance");
writelines(jsonencode(config),file);
changed=options; changed.frequency_hz="8000";
verifyError(testCase,@() radia.simulink.verifyIHPanelImpedanceProvenance(file,"","",changed),"radia:simulink:IHPanelZsProvenance");
writelines("changed",bh);
verifyError(testCase,@() radia.simulink.verifyIHPanelImpedanceProvenance(file,"","",options),"radia:simulink:IHPanelZsProvenance");
writelines("0 0",bh);
writelines("changed",artifact);
verifyError(testCase,@() radia.simulink.verifyIHPanelImpedanceProvenance(file,"","",options),"radia:simulink:IHPanelZsProvenance");
end
