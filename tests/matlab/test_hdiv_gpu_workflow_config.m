function tests=test_hdiv_gpu_workflow_config
% Fast pre-allocation contracts; no MATLAB GPU context is requested.
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
function testFixedPortBudget(t)
c=fixture();c.NumDoF=4097;
verifyError(t,@()radia.simulink.validateHDivGPUWorkflowConfig(c),'radia:simulink:HDivDimensions');
c=fixture();c.NumDoF=4;
verifyError(t,@()radia.simulink.validateHDivGPUWorkflowConfig(c),'radia:simulink:HDivDimensions');
end
function testNoSilentOldOptions(t)
c=fixture();c.FinalResidualOnly=true;
verifyError(t,@()radia.simulink.validateHDivGPUWorkflowConfig(c),'radia:simulink:HDivConfig');
c=fixture();c.Preconditioner='mass-cholesky';
verifyError(t,@()radia.simulink.validateHDivGPUWorkflowConfig(c),'radia:simulink:HDivRoute');
end
function testPreliminaryConfigDoesNotPublishFEAcceptance(t)
c=fixture();out=radia.simulink.validateHDivGPUWorkflowConfig(c);
verifyEqual(t,out,c); % Complete host FE payload checks remain in PreparedHDivNonlinear.
end
function c=fixture()
c=struct('Backend','gpu','MeshPath','','Order',2,'AppliedField',[1 0 0], ...
    'BHTable',[0 0;1 1],'NumDoF',3,'Assembly','ngsolve-host-normalized', ...
    'HostSetup',struct('mass',speye(3)),'Preconditioner','reference-diagonal', ...
    'SampleTime',.01,'RunDirectory',tempdir);
end
