function target = openIH(options)
%OPENIH Open the production native IH Simulink block.
%   radia.simulink.openIH() opens Applications/Induction Heating in the
%   installed Radia Simulink library.

arguments
    options.ModelName (1,1) string = "radia_ih"
    options.OutputDirectory (1,1) string = ""
    options.Save (1,1) logical = true
    options.Open (1,1) logical = true
end

libraryFile = string(which("radia_simulink_library.slx"));
if strlength(libraryFile) == 0
    if ispc
        outputDirectory = fullfile("C:\temp", "radia_simulink_library");
    else
        outputDirectory = fullfile(tempdir, "radia_simulink_library");
    end
    libraryFile = radia.simulink.buildLibrary( ...
        OutputDirectory=string(outputDirectory));
end
load_system(libraryFile);
target = "radia_simulink_library/Applications/Induction Heating";
if options.Open
    open_system(target);
end
end
