function target = openIH(options)
%OPENIH Open the standalone native IH Simulink preview model.

arguments
    options.ModelName (1,1) string = "radia_ih"
    options.OutputDirectory (1,1) string = ""
    options.Save (1,1) logical = true
    options.Open (1,1) logical = true
end

modelFile = "";
if options.ModelName == "radia_ih"
    modelFile = string(which("radia_ih.slx"));
end
if strlength(modelFile) == 0
    outputDirectory = options.OutputDirectory;
    if strlength(outputDirectory) == 0
        if ispc
            outputDirectory = fullfile("C:\temp", "radia_ih");
        else
            outputDirectory = fullfile(tempdir, "radia_ih");
        end
    end
    modelFile = radia.simulink.buildIHNativeModel( ...
        ModelName=options.ModelName, OutputDirectory=outputDirectory, ...
        Open=false);
end
load_system(modelFile);
target = options.ModelName;
if options.Open
    open_system(target);
end
end
