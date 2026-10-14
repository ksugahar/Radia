function model = install_radia_ih(options)
%INSTALL_RADIA_IH Add the release folder to MATLAB and open radia_ih.slx.
arguments
    options.Open (1,1) logical = true
end
matlabRoot = fileparts(mfilename("fullpath"));
addpath(matlabRoot);
radia.setup(RequireMex=true, Force=true);
radia.simulink.requireIHNativeRuntime();
model = radia.simulink.openIH(Open=options.Open);
end
