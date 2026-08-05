function target = openMagLev(options)
%OPENMAGLEV Open the standalone Radia MagLev Simulink model.

arguments
    options.ModelName (1,1) string = "radia_maglev"
    options.OutputDirectory (1,1) string = ""
    options.Open (1,1) logical = true
end

modelFile = "";
if options.ModelName == "radia_maglev"
    modelFile = string(which("radia_maglev.slx"));
end
if strlength(modelFile) == 0
    outputDirectory = options.OutputDirectory;
    if strlength(outputDirectory) == 0
        outputDirectory = fullfile("C:\temp", "radia_maglev");
    end
    modelFile = radia.simulink.buildMagLevModel( ...
        ModelName=options.ModelName, OutputDirectory=outputDirectory, ...
        Open=false);
end
load_system(modelFile);
target = options.ModelName;
if options.Open
    open_system(target);
end
end
