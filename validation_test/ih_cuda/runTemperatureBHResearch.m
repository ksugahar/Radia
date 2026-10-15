% Original public research case. No measured material data or Python solver.
here=fileparts(mfilename('fullpath'));
if isfile(fullfile(here,'makeTemperatureBHConfig.m'))
    addpath(here);
end
% Installed wheel: this script is under matlab/samples/ih-temperature-bh.
matlabFolder=fileparts(fileparts(here));
if isfile(fullfile(matlabFolder,'radia_cuda_ih.mexw64')),addpath(matlabFolder);end
[ihConfig,ihAssembly]=makeTemperatureBHConfig(); %#ok<NASGU>
outputFolder=tempname;mkdir(outputFolder);
modelFile=radia.simulink.buildTemperatureBHResearchModel(ihConfig,...
    fullfile(outputFolder,'ih_temperature_bh_public.slx'),...
    CurrentA=1e4,AmbientK=300,Samples=5);
open_system(modelFile);
disp(modelFile);
