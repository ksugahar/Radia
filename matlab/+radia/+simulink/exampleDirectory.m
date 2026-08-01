function folder = exampleDirectory()
%EXAMPLEDIRECTORY Return the folder containing Radia Simulink samples.
folder = fileparts(fileparts(fileparts(mfilename("fullpath"))));
end
