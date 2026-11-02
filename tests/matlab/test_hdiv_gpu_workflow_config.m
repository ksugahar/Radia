function tests=test_hdiv_gpu_workflow_config
% Guards: fast pre-allocation config/setup contracts; no MATLAB GPU context is requested.
tests=functiontests(localfunctions);
end
function setupOnce(t)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'packages','radia-simulink','matlab'));
end
function testExplicitGPUOnly(t)
c=fixture();c.Backend='cpu';
verifyError(t,@()radia.simulink.validateHDivGPUWorkflowConfig(c),'radia:simulink:HDivBackend');
end
function testDeclaredPortMatchesSetup(t)
c=fixture();c.NumDoF=5;
verifyError(t,@()radia.simulink.validateHDivGPUWorkflowConfig(c),'radia:simulink:HDivDimensions');
end
function testNoSilentOldOptions(t)
c=fixture();c.HostSetup=struct();
verifyError(t,@()radia.simulink.validateHDivGPUWorkflowConfig(c),'radia:simulink:HDivConfig');
c=fixture();c.Preconditioner='mass-cholesky';
verifyError(t,@()radia.simulink.validateHDivGPUWorkflowConfig(c),'radia:simulink:HDivRoute');
end
function testSetupSchemaIsChecked(t)
c=fixture();c.Setup.schema='radia.hdiv-host-setup.v1';
verifyError(t,@()radia.simulink.validateHDivGPUWorkflowConfig(c),'radia:gpu:DeviceSetup');
c=fixture();c.Setup.charges.cell_dofs=[0.5;1];
verifyError(t,@()radia.simulink.validateHDivGPUWorkflowConfig(c),'radia:gpu:DeviceSetup');
end
function testIntegerTablesRestoredToInt32(t)
out=radia.simulink.validateHDivGPUWorkflowConfig(fixture());
verifyClass(t,out.Setup.charges.cell_dofs,'int32');verifyClass(t,out.Setup.evaluation.element_points,'int32');
verifyClass(t,out.Setup.charges.compression,'char');verifyEqual(t,out.Setup.dofs,4);
end
function c=fixture()
% Structural payload only: the native MEX recomputes and checks every dimension.
charges=struct('cell_vertices',zeros(12,1),'face_vertices',zeros(0,1),'tet_rule',ones(8,1), ...
    'tri_rule',ones(4,1),'tet_low_rule',zeros(0,1),'tri_low_rule',zeros(0,1),'cell_tables',zeros(6,1), ...
    'face_tables',zeros(0,1),'charge_host',[0;0],'charge_kind',[0;0],'charge_exponents',[0;0;0;1;0;0], ...
    'cell_class',0,'cell_dofs',[3;1;0;2],'cell_exponents',[0;0;0;1;0;0],'face_class',zeros(0,1), ...
    'face_dofs',zeros(0,1),'face_exponents',zeros(0,1),'far_factor',1e30,'compression',"aca",'dofs',4, ...
    'aca',struct('eps',1e-8,'eta',2,'leaf_size',32,'max_rank',64,'rank_cap',"fail"));
evaluation=struct('rule',ones(12,1),'reference_shapes',zeros(36,1),'vertices',zeros(12,1), ...
    'element_class',0,'element_dofs',[3;1;0;2],'element_points',0);
setup=struct('schema',"radia.hdiv-device-setup.v1",'backend',"gpu",'family',"affine-tet",'order',1,'dofs',4, ...
    'mesh_sha256',repmat('a',1,64),'label_contract_sha256',repmat('b',1,64),'charges',charges, ...
    'evaluation',evaluation,'rhs_operator',zeros(4,3),'mean_operator',zeros(3,4));
c=struct('Backend','gpu','Setup',setup,'AppliedField',[1 0 0],'BHTable',[0 0;1 1],'NumDoF',4, ...
    'Preconditioner','reference-diagonal','SampleTime',.01,'RunDirectory',tempdir);
end
