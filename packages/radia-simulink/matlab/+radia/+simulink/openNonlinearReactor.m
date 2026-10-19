function modelName = openNonlinearReactor()
%OPENNONLINEARREACTOR Open the tracked nonlinear HDiv-MMM reactor model.
matlabRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
modelPath = fullfile(matlabRoot,"radia_nonlinear_reactor.slx");
if ~isfile(modelPath)
    modelPath = radia.simulink.buildNonlinearReactorModel();
end
open_system(modelPath);
[~,modelName] = fileparts(modelPath);
end
